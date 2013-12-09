;Procedure for extracting lines from a lamda file. Note that this may generate a lamda file with redundant
;levels (not connected by any of the included transitions). However,
;RADLite will not care about that, and it makes the bookkeeping much
;easier if we don't have to cull the population file. 
;
;

PRO lamda_extract_lines, isot=isot, lambdarange=lambdarange, max_energy=max_energy, $
                         lamda_path=lamda_path, molfile=molfile, freq=freq

@natconst.pro
@line_params.ini
;
;Set default keywords
;
IF ~KEYWORD_SET(lambdarange) THEN lambdarange=[4,5]
IF ~KEYWORD_SET(molfile) THEN molfile = 'moldata.dat'

CASE isot OF
   51: lamda_main='12CO_lamda.dat'
   52: lamda_main='13CO_Lamda.dat'
   141: lamda_main='HF_lamda.dat'
ENDCASE


molall = READ_MOLECULE_LAMBDA(lamda_path+lamda_main,/coll,/ghz,/kelvin)
mol    = LAMDA_EXTRACT_LEVELS(molall,vmax=vmax,jmax=jmax)

FREQRANGE=1d4/lambdarange

gsubs = WHERE(mol.freq GE FREQRANGE[1] AND mol.freq LE FREQRANGE[0] AND mol.eupper LT max_energy, nlines)

iup    = mol.iup[gsubs]
idown  = mol.idown[gsubs]
aud    = mol.aud[gsubs]
freq   = mol.freq[gsubs]
eupper = mol.eupper[gsubs]

OPENW,lun,molfile,/get_lun
printf,lun,'!MOLECULE'
printf,lun,mol.species
printf,lun,'!MOLECULAR WEIGHT'
printf,lun,mol.mumol,FORMAT='(f4.1)'
printf,lun,'!NUMBER OF ENERGY LEVELS'
printf,lun,mol.nlevels,FORMAT='(i6)'
printf,lun,'!LEVEL + ENERGIES(cm^-1) + WEIGHT + v + Q'
FOR i=0,mol.nlevels-1 DO BEGIN
    printf,lun,i+1,mol.e[i],mol.g[i],mol.lev_vib[i],mol.lev_rot[i],$
      FORMAT='(i5,d12.4,f7.1,a15,a15)'
ENDFOR
printf,lun,'!NUMBER OF RADIATIVE TRANSITIONS'
printf,lun,nlines,FORMAT='(i6)'
printf,lun,'!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(cm-1) + E_u(cm-1) + v_u + v_l + Q_p + Q_pp'

FOR i=0,nlines-1 DO BEGIN   
   printf, lun,i+1, iup[i],idown[i],aud[i],freq[i],eupper[i],mol.lev_vib[iup[i]-1],mol.lev_vib[idown[i]-1],mol.lev_rot[iup[i]-1],mol.lev_rot[idown[i]-1],$
           FORMAT='(i5,i5,i5,e12.3,f16.7,f12.5,a15,a15,a15,a15)'
ENDFOR
close,lun
free_lun, lun

END
