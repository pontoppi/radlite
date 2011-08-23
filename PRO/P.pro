FUNCTION P, npop
COMMON coll,  Cul, TCul
COMMON const, dv
COMMON mol,  nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g
COMMON grid, np, z_col, tgas_col, abun_col, rhogas_col, J_col, JSED_col, nu_cont=nu_cont

@natconst

Cul_arr     = DBLARR(nlines,np)
Clu_arr     = DBLARR(nlines,np)

tau_arr     = DBLARR(nlines,np)
Pesc_arr    = DBLARR(nlines,np)
Ppump_arr   = DBLARR(nlines,np)
R_arr       = DBLARR(nlevels,nlevels,np)

P           = DBLARR(nlevels,np)

FOR i=0,nlines-1 DO BEGIN
   
   Cul_arr[i,*] = INTERPOL(Cul[i,*],TCul,tgas_col)
   Clu_arr[i,*] = Cul_arr[i,*] * gugl[i] * EXP(-hh*cc*freq[i]/(kk*tgas_col))
   
   int = INT_SIMPLE(z_col,npop[idown[i]-1,*]*gugl[i]-npop[iup[i]-1,*])
   tau_arr[i,*] = Aul[i]/(8d0*!pi*freq[i]^3d0*dv) * int
   ;
   ;Now calculate the escape probabilities
   FOR h=0,np-1 DO BEGIN
      Pesc_arr[i,h]  = Pesc(tau_arr[i,h])
      Ppump_arr[i,h] = Ppump(tau_arr[i,h])
   ENDFOR
   ;
   ;Now we can get the modified rate
   ;coefficients - note that the
   ;subscripts follow the
   ;IDL column/row matrix indexing
   R_arr[idown[i]-1,iup[i]-1,*] = Aul[i]*Pesc_arr[i,*] + Bul[i]*Ppump_arr[i,*]*J_col+Cul_arr[i,*]*rhogas_col  ;Rul
   R_arr[iup[i]-1,idown[i]-1,*] = Blu[i]*Ppump_arr[i,*]*J_col + Clu_arr[i,*]*rhogas_col                       ;Rlu
   
ENDFOR

FOR h=0,np-1 DO BEGIN
   ;
   ;1st term: rate into level j; 2nd term: rate out of level j
   ;Pj = Sj (Rji x nj) - ni x Sj(Rij)
   P[*,h] = TRANSPOSE(R_arr[*,*,h])##REFORM(npop[*,h]) - TOTAL(R_arr[*,*,h],1)*npop[*,h]                   
   ;
   ;Replacing the last equilibrium
   ;equation with mass conservation to
   ;form a linearly independent system
   P[nlevels-1,h]   = abun_col[h] * rhogas_col[h] - TOTAL(npop[*,h])
ENDFOR

RETURN, P

END
