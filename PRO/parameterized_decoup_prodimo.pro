@analyze.pro
PRO parameterized_decoup_prodimo, tgas=tgas, mol_destruct=mol_destruct,speciesfrac=speciesfrac,mol_or_atom=mol_or_atom
@natconst.pro
@line_params.ini

PRINT, 'Calculating enhanced gas temperatures by scaling the typical ProDiMo gas temperature'
;from http://dianaproject.wp.st-andrews.ac.uk/data-results-downloads/an-example-disc-model/

IF ~KEYWORD_SET(mol_or_atom) THEN mol_or_atom = 'mol'

;
;Read model
ddens = read_dustdens()
dtemp = read_dusttemp()

nr = ddens.nr
nt = ddens.ntheta

density  = ddens.rho * gtd / (mp * mu) 
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
;FOR i=1, ddens.nr-1 DO BEGIN
;	FOR j=1,ddens.ntheta/2-1 DO BEGIN
;    	  column_density[i,j] = INT_TABULATED(ddens.r[i]*ddens.theta[0:j], density[i,0:j], /sort, /double)
;		  column_density[i,ddens.ntheta-1-j] = column_density[i,j]
;   	ENDFOR
;ENDFOR

;
;Read the phenomenological vertical structure
gd_xx = MRDFITS('gd_prodimo.fits',1, /silent)
gd_nh = MRDFITS('gd_prodimo.fits',3, /silent)
gd_td = MRDFITS('gd_prodimo.fits',4, /silent)
gd_tg = MRDFITS('gd_prodimo.fits',5, /silent)
nh2 = MRDFITS('gd_prodimo.fits',6, /silent)
nh = MRDFITS('gd_prodimo.fits',7, /silent)
nx_prodimo = (SIZE(gd_xx))[1]
ny_prodimo = (SIZE(gd_xx))[2]

gd_arr = gd_tg/gd_td
molfrac_arr = nh2/(nh+nh2)

bsubs = WHERE(gd_td EQ 0)
gd_arr[bsubs] = 1
gd_r = gd_xx[*,0]

gastemp = dusttemp
speciesfrac = FLTARR(ddens.nr,ddens.ntheta)

FOR i=1,ddens.nr-1 DO BEGIN
   r_AU = ddens.r[i]/AU
   NH = REFORM(column_density[i,*])

   fac_r = FLTARR(N_ELEMENTS(gd_r),nt)
   molfrac_r = FLTARR(N_ELEMENTS(gd_r),nt)

   FOR k=0,nx_prodimo-1 DO BEGIN
        LINTERP,REFORM(gd_nh[k,*]), REFORM(gd_arr[k,*]),NH, yout
        fac_r[k,*] = yout
        LINTERP,REFORM(gd_nh[k,*]),REFORM(molfrac_arr[k,*]),NH, yout
    	molfrac_r[k,*] = yout
   ENDFOR
   FOR j=1,ddens.ntheta-1 DO BEGIN
      LINTERP, gd_r, fac_r[*,j], ddens.r[i]/AU, fac
      LINTERP, gd_r, molfrac_r[*,j], ddens.r[i]/AU, molfrac_pt

	  gastemp[i,j] = gastemp[i,j] * MAX([fac,1])
     ; Is the species molecular or atomic?
     IF mol_or_atom EQ 'mol' THEN BEGIN
        speciesfrac[i,j] = molfrac_pt
     ENDIF ELSE BEGIN
        speciesfrac[i,j] = 1.0 - molfrac_pt
     ENDELSE
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
