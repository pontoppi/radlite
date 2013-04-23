PRO optimal_sample, z_col,tgas_col,rhogas_col,abun_col,JSED_col

ssubs = SORT(z_col)
tgas_col = tgas_col[ssubs]
rhogas_col = rhogas_col[ssubs]
abun_col = abun_col[ssubs]
nf = (SIZE(JSED_col))[2]
stop
FOR i=0,nf-1 DO BEGIN
   JSED_col[*,i] = JSED_col[ssubs,i]
ENDFOR

END
