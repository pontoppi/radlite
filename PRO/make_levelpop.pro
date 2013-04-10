
;====================================
;Generate level population files for RADLite
;====================================

PRO make_levelpop, ddens=ddens, tgas=tgas, rhogas=rhogas, abun=abun, psum=psum, molfile=molfile, run_nlte=run_nlte
@line_params.ini
@natconst.pro

;====================================
; Set the populations to LTE, if requested
;====================================
IF lte EQ 1 THEN BEGIN
   ;====================================
   ;Read the molecular data
   ;====================================
   mol  = read_molecule_lambda(molfile)
   psum = read_psum(psumfile, mol)

   nlev   = mol.nlevels
   energy = mol.e[0:nlev-1]*hh*cc
   gunit  = mol.g[0:nlev-1]
   nr     = ddens.nr
   nt     = ddens.ntheta/2
   npop   = DBLARR(nr,nt,nlev)

   FOR i=0,nlev-1 DO BEGIN
      npop[*,*,i] = mol.g[i] * $
                    exp(-(mol.energy_in_K[i]/(tgas[*,0:nt-1])))
   ENDFOR
   ;
   ;get partition sum:
   part = interpol(psum.psum,psum.temp,tgas[*,0:nt-1])
   FOR i=0,nlev-1 DO BEGIN
      npop[*,*,i] = npop[*,*,i] / part
   ENDFOR

;===================================
;Or solve the detailed balance equations for NLTE
;===================================
ENDIF ELSE BEGIN
   mol    = read_molecule_lambda(molfile)
   CASE mol.species OF 
      '12CO': lamda_main='12CO_lamda.dat'
      '13CO': lamda_main='13CO_Lamda.dat'
      'HF': lamda_main='HF_lamda.dat'
   ENDCASE
   ;
   ;Read the full lamda file      
   molall = READ_MOLECULE_LAMBDA(main_path+'LAMDA/'+lamda_main,/coll,/ghz)
   ;
   ;use the Newton-Raphson global solver, if requested
   IF run_nlte EQ 1 THEN BEGIN
      PRINT, '...starting detailed balance calculation...'
      nlte_main, molall=molall, tgas=tgas, rhogas=rhogas, abun=abun, ddens=ddens
   ENDIF   
   it_is_there = FILE_TEST('levelpop_nlte.fits')
   IF it_is_there THEN BEGIN
      pop    = mrdfits('levelpop_nlte.fits',1,/SILENT)   
      mol    = mrdfits('levelpop_nlte.fits',2,/SILENT)
      nlev   = mol.nlevels
      energy = mol.energy_in_K * kk
      gunit  = mol.g

      nr = ddens.nr
      nt = ddens.ntheta/2

      npop   = DBLARR(nr,nt,nlev)
      FOR ir=0,nr-1 DO BEGIN
         FOR it=0,nt-1 DO BEGIN
            ;
            ;Convert to fractional level populations
            pop.npop_all[*,it,ir] = pop.npop_all[*,it,ir]/TOTAL(pop.npop_all[*,it,ir])
            FOR il=0,nlev-1 DO BEGIN 
               npop[ir,it,il] = pop.npop_all[il,it,ir]
            ENDFOR
         ENDFOR
      ENDFOR
   ENDIF ELSE BEGIN
      print, 'Error: Non-LTE level population file not found!'
   ENDELSE  
ENDELSE

lowsubs = WHERE(npop LT 1d-99)
npop[lowsubs] = 0.0

openw,lunl,'levelpop_'+molfile,/get_lun
printf,lunl,nr,nt,nlev,1
printf,lunl,energy
printf,lunl,gunit

FOR ir=0,nr-1 DO BEGIN
   FOR it=0,nt-1 DO BEGIN 
      printf,lunl,npop[ir,it,0:nlev-1]
   ENDFOR
ENDFOR
close,lunl
free_lun, lunl



openw,lun,'levelpop.info',/get_lun
printf,lun,'-3'
printf,lun,'levelpop_'+molfile
printf,lun,0
close,lun
free_lun, lun


END
