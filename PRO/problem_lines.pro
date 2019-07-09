PRO problem_lines, molfile, run_nlte=run_nlte
@line_params.ini
@natconst.pro

;====================================
; Messages
;====================================
IF N_ELEMENTS(VERBOSE) THEN print,' - Taking dust density x ', STRTRIM(STRING(gtd),2)

;====================================
; Make new radlite.inp
;====================================
write_radlite_line,niterline,/noisrf,cir_np=cir_np,b_per_r=b_per_r,$
                   b_extra=b_extra

;====================================
; Read the dust density and temperature
;====================================
ddens  = read_dustdens()
dtemp  = read_dusttemp()
nr     = n_elements(ddens.r)
nt     = n_elements(ddens.theta)/2

;====================================
;The gas temperature
;====================================
IF gas_decoup EQ 0 THEN BEGIN
   tgas   = dtemp.temp[*,*,0,0] ; gas temp = dust temp
   openw,lunt,'temperature.inp',/get_lun
   printf,lunt,nr,nt,1
   for ir=0,nr-1 do begin
      for it=0,nt-1 do begin
         printf,lunt,tgas[ir,it]
      endfor
   endfor
   close,lunt
   free_lun, lunt
ENDIF
IF gas_decoup EQ 1 THEN BEGIN
	   
  ;In this case the gas
  ;temperature is calculated
  ;using Rowin's phenomenological
  ;model. We will make sure that there
  ;is no confusion by removing any
  ;existing gas temperature file. 

  ;Note that xray_abundance is also
  ;called if make_abundance is called
  ;with abun_str = 5
   spawn, 'rm temperature.inp'
   xray_abundance, abun=abun,tgas=tgas
ENDIF
IF gas_decoup EQ 2 THEN BEGIN
    spawn, 'rm temperature.inp'
	parameterized_decoup_najita,tgas=tgas,mol_destruct=mol_destruct
ENDIF
IF gas_decoup EQ 3 THEN BEGIN
   spawn, 'rm temperature.inp'
   ; Right now, the only supported atomic species is oxygen. 
   IF isot EQ 341 THEN BEGIN
      mol_or_atom='atom'
   ENDIF ELSE BEGIN
      mol_or_atom='mol'
   ENDELSE

	parameterized_decoup_prodimo,tgas=tgas,speciesfrac=speciesfrac,mol_or_atom=mol_or_atom
ENDIF
;
;Some times the gas temperature comes
;out 0.0000 - have to check what goes wrong...
badspots = WHERE(tgas EQ 0.0000)
IF badspots[0] NE -1 THEN BEGIN
   tgas[badspots]=5d0           ;setting them to 5K
   print, 'WARNING: Temperature is 0 in some grid points - setting to 5K' 
ENDIF
;

;======================================
;Now create the gas velocity file.
;Some velocity fields will return an additional density component to
;satisfy the continuum equation.
;======================================
;
;Initialize add_dens
add_dens = 0.
make_velocity,ddens,vtype,add_dens=add_dens

;=======================================
; Create the gas density and temperature
;=======================================
rhogas = ddens.rho[*,*,0]*gtd;
gsubs = WHERE(add_dens NE 0)
IF gsubs[0] NE -1 THEN rhogas[gsubs] = add_dens[gsubs]

;======================================
; Write the gas density
;======================================

; Take only isize=1 and ispec=1
;
openw,lund,'density.inp', /get_lun
printf,lund,nr,nt,1
for ir=0,nr-1 do begin
    for it=0,nt-1 do begin
        printf,lund,rhogas[ir,it]
    endfor
endfor
close,lund
free_lun, lund


;======================================
;Make an abundance file
;======================================
make_abundance,abun_str,PT_rel=PT_rel,abun=abun, mol_destruc=mol_destruct,speciesfrac=speciesfrac

;====================================
;Test the velocity
;====================================
read_vel, 'velocity.inp',vel

;====================================
; Check passband width
;====================================
velmax = max(abs(vel.vphi))
IF N_ELEMENTS(VERBOSE) THEN PRINT,'Maximum velocity = ',velmax/1d5,' km/s'

;====================================
;Now that we have the molecular mass, 
;we can calculate the intrinsic
;line widths
;====================================

;====================================
;Check for user error
;====================================
IF turb_kep NE 0 and turb_sou NE 0 THEN BEGIN
   PRINT, 'Ambiguous turbulence description! Please set turb_kep or turb_sou to 0'
   STOP
ENDIF
IF turb_kep EQ 0 and turb_sou EQ 0 THEN BEGIN
   PRINT, 'No turbulence description set!'
   STOP
ENDIF
IF turb_kep NE 0 THEN BEGIN
   mol = READ_MOLECULE_LAMBDA(molfile)
   make_turbulence,ddens,2,alpha=0.0,kepler_frac=turb_kep,molmass=mol.mumol
ENDIF
IF turb_sou NE 0 THEN BEGIN
   mol = READ_MOLECULE_LAMBDA(molfile)
   make_turbulence,ddens,1,alpha=turb_sou,kepler_frac=0.0,molmass=mol.mumol
ENDIF

;====================================
;Compute level populations / extract previously computed populations
;====================================

make_levelpop, ddens=ddens, tgas=tgas, rhogas=rhogas, abun=abun, psum=psum, molfile=molfile, run_nlte=run_nlte

END
