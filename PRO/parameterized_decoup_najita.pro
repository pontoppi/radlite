@analyze.pro
PRO parameterized_decoup_najita, tgas=tgas, mol_destruct=mol_destruct
@natconst.pro
@line_params.ini

PRINT, 'Calculating enhanced gas temperatures by scaling the Najita et al. 2011 gas temperatures'

;
;Read model
ddens = read_dustdens()
dtemp = read_dusttemp()

ndens  = N_ELEMENTS(dens)
ntemps = N_ELEMENTS(temps)

density  = ddens.rho / (mp * mu) * gtd
dusttemp = dtemp.temp(*,*,0,0)

column_density = dblarr(ddens.nr,ddens.ntheta)
FOR ir=0,ddens.nr-1 DO BEGIN
    FOR it=1,ddens.ntheta/2-1 DO BEGIN
        column_density[ir,it] = column_density[ir,it-1] + $
          0.5 * ( density[ir,it] + density[ir,it-1] ) *$
          ( cos(ddens.theta[it-1]) - cos(ddens.theta[it]) ) * ddens.r[ir]
        ;
        ;Mirror
        column_density[ir,ddens.ntheta-it-1] = column_density[ir,it]
    ENDFOR
ENDFOR

;column_density = DBLARR(ddens.nr,ddens.ntheta)
;
;FOR i=1, ddens.nr-1 DO BEGIN
;	FOR j=1,ddens.ntheta/2-1 DO BEGIN
;    	  column_density[i,j] = INT_TABULATED(ddens.r[i]*ddens.theta[0:j], density[i,0:j], /sort, /double)
;		  column_density[i,ddens.ntheta-1-j] = column_density[i,j]
 ;  	ENDFOR
;ENDFOR

;
;Read the phenomenological vertical structure
readcol, 'gd.dat', gd_nh, gd_025, gd_05, gd_1, gd_4, gd_10, gd_20

gd_r = [0.01,0.25,0.5,1.,4.,10.,20.,2000.]
gd_arr = [[gd_025],[gd_025],[gd_05],[gd_1],[gd_4],[gd_10],[gd_20],[gd_20]]
fac_r = FLTARR(N_ELEMENTS(gd_r))

gastemp = dusttemp
mol_destruct = INTARR(ddens.nr,ddens.ntheta)

FOR i=1,ddens.nr -1 DO BEGIN
   FOR j=1,ddens.ntheta-1 DO BEGIN

	  r_AU = ddens.r[i]/AU
	  NH = column_density[i,j]
	  FOR k=0,N_ELEMENTS(fac_r)-1 DO BEGIN
		  fac_r[k] = INTERPOL(gd_arr[*,k],gd_nh,alog10(NH+1e-20))
	  ENDFOR
	  
	  fac = INTERPOL(fac_r,gd_r,ddens.r[i]/AU)
	  IF alog10(NH) LT 24. AND alog10(NH) GT 19 THEN BEGIN
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
