PRO optimal_sample, z_col,tgas_col,rhogas_col,abun_col,JSED_col

ssubs      = SORT(z_col)
tgas_col   = tgas_col[ssubs]
rhogas_col = rhogas_col[ssubs]
abun_col   = abun_col[ssubs]
z_col      = z_col[ssubs]

nf = (SIZE(JSED_col))[2]

sg_filter = SAVGOL(3, 3, 0, 2)

FOR i=0,nf-1 DO BEGIN
   JSED_col[*,i] = JSED_col[ssubs,i]

   highsub = WHERE(JSED_col[*,i] GT 0)
   minJ    = MAX([1d-20,MIN(JSED_col[highsub,i])])

   lowsub = WHERE(JSED_col[*,i] LT minJ)
   JSED_col[lowsub,i] = minJ
   JSED_col[*,i] = ABS(CONVOL(MEDIAN(JSED_col[*,i],2),sg_filter,/EDGE_TRUNCATE))
ENDFOR

;tgas_col = ABS(CONVOL(tgas_col,sg_filter,/EDGE_TRUNCATE))


END
