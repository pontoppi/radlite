@analyze_v2.pro
@analyze.pro
@make_abundance.pro
@xray_abundance.pro
@deriv
@int_tabulated
@uniq
@interpol

PRO excitation_setup, run_dir = run_dir

IF NOT KEYWORD_SET(run_dir) THEN run_dir = 'excitation_files'

@natconst.pro
@line_params.ini

maxpoints = 1000 ; Maximum number of theta (h) input points for the excitation code.
regN      = 100. ; number of grid points required for a regular theta grid (this defines the minimum dh.
dT_T      = 0.01 ; Required step in temperature in percent/100.

;
;Read the dust density and convert to H number density 
ddens   = read_dustdens()
density = ddens.rho / (mp * mu) * gtd * 2.0

;
;Read the dust temperature
dtemp = read_dusttemp()

;
;Read the mean intensity
;sfunc = read_source_function()
mint  = read_meanint()

;
;Mass density of the gas
rhogas  = ddens.rho[*,*,0]*gtd

;
;Number of grid points in RADMC setup
nr = n_elements(ddens.r)
nt = n_elements(ddens.theta)/2

;
;Write the gas density to a file 
openw,lund,'density.inp', /get_lun
printf,lund,nr,nt,1
for ir=0,nr-1 do begin
    for it=0,nt-1 do begin
        printf,lund,rhogas[ir,it]
    endfor
endfor
close,lund
free_lun, lund

;
;Now, there are a number of possibilities. The gas can be coupled to
;the dust or not. 
IF gas_decoup EQ 0 THEN BEGIN
   ;
   ;If the gas temp = dust temp, then we
   ;can just write the temperature right away.
   tgas   = dtemp.temp[*,*,0,0] 
   openw,lunt,'temperature.inp',/get_lun
   printf,lunt,nr,nt,1
   for ir=0,nr-1 do begin
      for it=0,nt-1 do begin
         printf,lunt,tgas[ir,it]
      endfor
   endfor
   close,lunt
   free_lun, lunt
ENDIF ELSE BEGIN
  ;In this case the gas
  ;temperature is calculated
  ;using Rowin's phenomenological
  ;model. We will make sure that there
  ;is no confusion by removing any
  ;existing gas temperature file.
 
   spawn, 'rm temperature.inp'   
   IF abun_str NE 5 THEN xray_abundance,abun=abun,tgas=tgas
   
ENDELSE

;
;There is a possibility that the gas temperature is calculated
;according to Rowin's model, but the user asks for a different
;(and inconsistent!) abundance structure. If so, calculate the new
;abundance structure here.
make_abundance,abun_str,PT_rel=PT_rel, abun=abun, tgas=tgas

;
;Now we can loop over each radius to produce all the excitation input
;directories. 
FOR ri=0,n_elements(ddens.r)-1 DO BEGIN
   
   writeu, -1, string(13b)
   writeu, -1, '---> '+strtrim(string(100d0*ri / (n_elements(ddens.r)-1),format='(f6.2)'),2) + '% complete'

   radius = ddens.r[ri]
   height = dblarr(maxpoints)

   ;
   ;cut out only the top part of the disk
   sel = where(0.5 * !pi - ddens.tt[ri,*] gt 0.0)

   gastemperature_dum = tgas[ri,sel]
   height_dum = (0.5*!pi-ddens.tt[ri,sel]) * radius

   ;
   ;What is the range of h covered?
   min_height = min(height_dum)
   max_height = max(height_dum)

   ;
   ;Not sure this really helps, but the
   ;idea is to smooth out the derivative
   ;a little.
   height_int = (max_height-min_height)*FINDGEN(10000)/9999+min_height
   gastemperature_int = INTERPOL(gastemperature_dum,height_dum,height_int)

   ;
   ;Calculate the derivative dh/dT
   dhdT = ABS(DERIV(gastemperature_int,height_int)) 

   ;
   ;Initialize the height at the minimum value in the grid.
   height[0] = min_height

   ;
   ;Calculate the minimum height step size.
   dhmin = max_height/regN

   ;
   ;Now loop until the maximum height is reached
   i = 0
   WHILE height[i] LT max_height DO BEGIN
      dh = INTERPOL(dhdT,height_int,height[i])*dT_T*INTERPOL(gastemperature_int,height_int,height[i])
      dh = MIN([dh,dhmin])
      height[i+1] = height[i] + dh
      i+=1
      IF i GT maxpoints THEN BEGIN
         PRINT, 'Not enough height grid points!!'
         STOP
      ENDIF
   ENDWHILE

   ;
   ;Truncate the output grid.
   ntot   = i
   height = height[0:ntot-1]
   
   ;
   ;Now interpolate all the required quantities on the new height grid
   height_radmc     = REFORM((ddens.tt[ri,*]-!pi/2.d0) * radius,ddens.ntheta)
   density_h        = INTERPOL(density[ri,*],height_radmc,height)
   gastemperature_h = INTERPOL(gastemperature_dum,height_dum,height)
   temperature_h    = INTERPOL(dtemp.temp[ri,*],height_radmc,height)
   abun_h           = INTERPOL(abun[ri,*],height_radmc,height)

   ;
   ;The mean intensity has to be interpolated for each dust frequency point
   meanintens_h = fltarr(ntot,mint.nnu)
   FOR inu=0,mint.nnu-1 DO BEGIN
      meanintens_h[*,inu]     = INTERPOL(mint.meanint[ri,*,inu],height_radmc,height)
   ENDFOR

   ;
   ;And write everything to file
   radius_dir = 'excit_dir' + STRTRIM(STRING(ri),2)
   spawn, 'mkdir '+ radius_dir, result,error
   openw, lun, radius_dir  + '/test_model.inp',/get_lun
   printf, lun, ntot
   openw, lunmi, radius_dir + '/test_model_mean_intens.inp',/get_lun  
   printf, lunmi, ntot

   mifstring = '('+STRTRIM(STRING(mint.nnu),2)+'(2x,e12.6))'
   ;
   ;Writing from high altitude to low
   FOR hi=ntot-1,0,-1 DO BEGIN   
      printf,lun, format='(8(e12.6, 2x))', height[hi], density_h[hi], gastemperature_h[hi], temperature_h[hi], $
             0.0, 0.5, 0.0, abun_h[hi]
      printf,lunmi, format=mifstring, height[hi], meanintens_h[hi,*]
   ENDFOR
   close,lun
   close,lunmi
   FREE_LUN, lun
   FREE_LUN, lunmi

endfor

;
;Copy into the excitation setup directory structure
spawn, 'mkdir ' + run_dir,result,error
spawn, 'cp -r excit_dir* ' + run_dir

;
;And clean up
spawn, 'rm -r excit_dir*'
end
