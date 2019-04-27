@analyze.pro
PRO parameterized_decoup_prodimo, tgas=tgas, mol_destruct=mol_destruct
@natconst.pro
@line_params.ini

PRINT, 'Calculating enhanced gas temperatures by scaling the typical ProDiMo gas temperature'
;from http://dianaproject.wp.st-andrews.ac.uk/data-results-downloads/an-example-disc-model/

;
;Read model
ddens = read_dustdens()
dtemp = read_dusttemp()

ndens  = N_ELEMENTS(dens)
ntemps = N_ELEMENTS(temps)

density  = ddens.rho / (mp * mu) * gtd
dusttemp = dtemp.temp(*,*,0,0)

column_density = DBLARR(ddens.nr,ddens.ntheta)

FOR i=1, ddens.nr-1 DO BEGIN
	FOR j=1,ddens.ntheta/2-1 DO BEGIN
    	  column_density[i,j] = INT_TABULATED(ddens.r[i]*ddens.theta[0:j], density[i,0:j], /sort, /double)
		  column_density[i,ddens.ntheta-1-j] = column_density[i,j]
   	ENDFOR
ENDFOR

;
;Read the phenomenological vertical structure
gd_xx = MRDFITS('gd_prodimo.fits',1, /silent)
gd_nh = MRDFITS('gd_prodimo.fits',3, /silent)
gd_td = MRDFITS('gd_prodimo.fits',4, /silent)
gd_tg = MRDFITS('gd_prodimo.fits',5, /silent)
nh2 = MRDFITS('gd_prodimo.fits',6, /silent)
gd_arr = gd_tg/gd_td

bsubs = WHERE(gd_td EQ 0)
gd_arr[bsubs] = 1

gd_r = gd_xx[*,0]

fac_r = FLTARR(N_ELEMENTS(gd_r))

gastemp = dusttemp
mol_destruct = INTARR(ddens.nr,ddens.ntheta)

FOR i=1,ddens.nr -1 DO BEGIN
   FOR j=1,ddens.ntheta-1 DO BEGIN

	  r_AU = ddens.r[i]/AU
	  NH = column_density[i,j]
	  FOR k=0,N_ELEMENTS(fac_r)-1 DO BEGIN
        IF NH LT MAX(gd_nh[k,*]) AND NH GT MIN(gd_nh[k,*]) THEN BEGIN
           fac_r[k] = INTERPOL(gd_arr[k,*],gd_nh[k,*],NH)
        ENDIF ELSE BEGIN
           fac_r[k] = 1
        ENDELSE
	  ENDFOR
	  fac = INTERPOL(fac_r,gd_r,ddens.r[i]/AU)
	  IF alog10(NH) LT 24. AND alog10(NH) GT 20 THEN BEGIN
		  gastemp[i,j] = gastemp[i,j] * fac
	  ENDIF
	  IF fac GT 4. THEN BEGIN
		  mol_destruct[i,j] = 1
	  ENDIF
   ENDFOR
ENDFOR

;
;Write the gas temperature.       
openw,lun,'temperature.inp',/get_lun
printf,lun, ddens.nr,ddens.ntheta/2,1

FOR ir=0,ddens.nr-1 DO BEGIN
   FOR it=0,ddens.ntheta/2-1 DO BEGIN
      printf, lun, gastemp[ir,it]
   ENDFOR
ENDFOR
close,lun
FREE_LUN, lun

tgas = gastemp

END
