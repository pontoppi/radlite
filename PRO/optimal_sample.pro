;===================================================
;This code smooths the mean intensity SED in the z=direction
;at each wavelength point using a Savitsky-Golay filter. 
;
;===================================================

PRO optimal_sample, z_col,tgas_col,rhogas_col,abun_col,JSED_col

nf = (SIZE(JSED_col))[2]

sg_filter = SAVGOL(3, 3, 0, 2)

FOR i=0,nf-1 DO BEGIN

   minJ    = MAX([1d-30,MIN(JSED_col[*,i])])
   lowsub = WHERE(JSED_col[*,i] LE minJ)
   JSED_col[lowsub,i] = minJ
   JSED_col[*,i] = ABS(CONVOL(MEDIAN(JSED_col[*,i],2),sg_filter,/EDGE_TRUNCATE))
ENDFOR

END
