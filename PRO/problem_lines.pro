PRO problem_lines, molfile
@line_params.ini
@natconst.pro

;
; Messages
;
print,' - Taking dust density x ', STRTRIM(STRING(gtd),2)
;
; Make new radlite.inp
;
write_radlite_line,niterline,/noisrf,cir_np=cir_np,b_per_r=b_per_r,$
                   b_extra=b_extra
;
; Read the dust density and temperature
;
ddens  = read_dustdens()
dtemp  = read_dusttemp()
nr     = n_elements(ddens.r)
nt     = n_elements(ddens.theta)/2

;
;The gas temperature
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
ENDIF ELSE BEGIN
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
   
ENDELSE
   
;
;Now create the gas velocity file.
;Some velocity fields will return an additional density component to
;satisfy the continuum equation.
;
add_dens = 0.
make_velocity,ddens,vtype,add_dens=add_dens

;
; Now create the gas density and temperature
;
rhogas = ddens.rho[*,*,0]*gtd;
gsubs = WHERE(add_dens NE 0)
IF gsubs[0] NE -1 THEN rhogas[gsubs] = add_dens[gsubs]

;
; Now create a gas density out of this

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


;
;And abundance
IF gas_decoup NE 0 THEN BEGIN
   IF isot NE 11 THEN BEGIN
      PRINT, 'WARNING: You are attempting to use Rowins phenomenological abundance + gas temperature model. This'
      PRINT, 'is intended for H2O only, but you are running some other molecule. '
;      stop
   ENDIF
   IF abun_str NE 5 THEN BEGIN
      PRINT, 'WARNING: We are using Rowins phenomenological gas temperature structure. This also includes'
      PRINT, 'an H2O abundance that is density dependent and is destroyed in the upper layers of the disk.'
      PRINT, 'This usually overrides any other abundance structure!'
      PRINT, 'Are you sure you want to use gas temp NE dust temp, but with some other abundance structure?'
;      stop
   ENDIF
ENDIF ELSE BEGIN
   make_abundance,abun_str,fr_temp=fr_temp,PT_rel=PT_rel,abun=abun
ENDELSE
;
;Test the velocity and plot for inspection
;
read_vel, 'velocity.inp',vel

;
; Check passband width
;
vmax = max(abs(vel.vphi))
print,'Maximum velocity = ',vmax/1d5,' km/s'
mol = read_molecule_lambda(molfile)

;====================================================================
;READING PARTITION SUMS
;====================================================================
dat = dblarr(3000,2)
openr, lunps, psumfile,/get_lun
str0 = '' & str1 = ''  & str2 = ''  & str3 = ''  & str4 = ''  & str5 = ''  & str6 = ''  & str7 = ''  & str8 = ''  
str9 = '' & str10 = ''  & str11 = ''  & str12 = ''  & str13 = ''  & str14 = ''  & str15 = ''  & str16 = ''  & str17 = ''  & str18 = ''
strdum = ''
i1 = '' & i2 = '' & i3 = '' & i4 = '' & i5 = '' & i6 = '' & i7 = '' & i8 = '' & i9 = ''
i10 = '' & i11 = '' & i12 = '' & i13 = '' & i14 = '' & i15 = '' & i16 = '' & i17 = '' & i18 = ''

readf,lunps,str0,str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17,str18,$ ;skip first line
  format='(a19,a28,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27)'
readf,lunps,str0,str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17,str18,$ 
  format='(a19,a24,a23,a24,a15,a25,a25,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27)'
psum_mols = STRTRIM([str0,str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17,str18],2)
readf,lunps,strdum,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,$ 
  format='(a19,a24,a23,a24,a15,a25,a25,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27)'
molmasses = STRTRIM([strdum,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18])

col_sub = WHERE(psum_mols eq mol.species)
molmass = float(molmasses[col_sub])

IF col_sub[0] EQ -1 THEN BEGIN
    print, 'No partition sum table found for the requested molecule: ', mol.species
    stop
ENDIF

i=0

WHILE NOT EOF(lunps) DO BEGIN
    readf, lunps, d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18
    psums = [d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18]
    dat[i,*] = [psums[0],psums[col_sub]]
    i=i+1
ENDWHILE
dat = dat[0:i-1,*]
close, lunps
free_lun, lunps
psum = dat

;
;Now that we have the molecular mass, we can calculate the intrinsic
;line widths

;
;Check for user error
IF turb_kep NE 0 and turb_sou NE 0 THEN BEGIN
   PRINT, 'Ambiguous turbulence description! Please set turb_kep or turb_sou to 0'
   STOP
ENDIF
IF turb_kep EQ 0 and turb_sou EQ 0 THEN BEGIN
   PRINT, 'No turbulence description set!'
   STOP
ENDIF
IF turb_kep NE 0 THEN BEGIN
   make_turbulence,ddens,2,alpha=0.0,kepler_frac=turb_kep,molmass=molmass
ENDIF
IF turb_sou NE 0 THEN BEGIN
   make_turbulence,ddens,1,alpha=turb_sou,kepler_frac=0.0,molmass=molmass
ENDIF

;====================================================================
;
; Put the populations to LTE, if requested
;
IF lte EQ 1 THEN BEGIN

   npop = dblarr(nr,nt,mol.nlevels)
   ntot = dblarr(nr,nt)
   
   ;
   ;Some times the gas temperature comes
   ;out 0.0000 - have to check what goes wrong...
   badspots = WHERE(tgas EQ 0.0000)
   IF badspots[0] NE -1 THEN BEGIN
      tgas[badspots]=5d0        ;setting them to 5K
      print, 'WARNING: Temperature is 0 in some grid points - setting to 5K' 
   ENDIF
   ;
   ntot = ntot + npop[*,*,0]
   FOR i=0,mol.nlevels-1 DO BEGIN
      npop[*,*,i] = mol.g[i] * $
                    exp(-(mol.energy_in_K[i]/(tgas[*,0:nt-1])))
   ENDFOR
   ;
   ;get partition sum:
   
   part = interpol(psum[*,1],psum[*,0],tgas[*,0:nt-1])
   FOR i=0,mol.nlevels-1 DO BEGIN
      npop[*,*,i] = npop[*,*,i] / part
   ENDFOR
   
   openw,lunl,'levelpop_'+molfile,/get_lun
   printf,lunl,nr,nt,mol.nlevels,1
   printf,lunl,mol.e[0:mol.nlevels-1]*hh*cc
   printf,lunl,mol.g[0:mol.nlevels-1]
   
   FOR ir=0,nr-1 DO BEGIN
      FOR it=0,nt-1 DO BEGIN
         printf,lunl,npop[ir,it,0:mol.nlevels-1]
      ENDFOR
   ENDFOR
   close,lunl
   free_lun, lunl
ENDIF ELSE BEGIN
   ;
   ;Check for existing non-lte level population file
   it_is_there = FILE_TEST('levelpop_nlte.fits')
   IF it_is_there THEN BEGIN
      PRINT, 'Existing level population file found - do you want to use it?'
      read, answer, prompt='[y/n]'
      
   ENDIF

   PRINT, 'You have selected non-LTE!'
   PRINT, '...starting detailed balance calculation...'
   
   nlte_main, tgas=tgas, rhogas=rhogas, abun=abun, ddens=ddens
ENDELSE

openw,lun,'levelpop.info',/get_lun
printf,lun,'-3'
printf,lun,'levelpop_'+molfile
printf,lun,0
close,lun
free_lun, lun

END
