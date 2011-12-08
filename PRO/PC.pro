FUNCTION PC, npop_in, col,m

@natconst

Cul_arr     = DBLARR(m.nctrans,m.np)
Clu_arr     = DBLARR(m.nctrans,m.np)

tau_arr     = DBLARR(m.nlines,m.np)
Pesc_arr    = DBLARR(m.nlines,m.np)
Ppump_arr   = DBLARR(m.nlines,m.np)
R_arr       = DBLARR(m.nlevels,m.nlevels,m.np)

P           = DBLARR(m.nlevels*m.np)

npop        = REFORM(npop_in,m.nlevels,m.np)

;======================================
;Add radiative terms to the rate matrix
;======================================
FOR i=0,m.nlines-1 DO BEGIN 
   int = INT_SIMPLE(col.z,npop[m.iup[i]-1,*] - npop[m.idown[i]-1,*]*m.gugl[i])  ;int(z->zmax)
   tau_arr[i,*] = m.Aul[i]/(8d0*!pi*m.freq[i]^3d0*m.dv) * int
   
   ;
   ;Now calculate the escape probabilities
   FOR h=0,m.np-1 DO BEGIN
      Pesc_arr[i,h]  = Pesc(tau_arr[i,h])
      Ppump_arr[i,h] = Ppump(tau_arr[i,h])
   ENDFOR
   ;
   ;Now we can get the modified rate
   ;coefficients - note that the
   ;subscripts follow the
   ;IDL column/row matrix indexing
   R_arr[m.idown[i]-1,m.iup[i]-1,*] = m.Aul[i]*Pesc_arr[i,*] + m.Bul[i]*Ppump_arr[i,*]*col.J[i,*]      ;Rul
   R_arr[m.iup[i]-1,m.idown[i]-1,*] = m.Blu[i]*Ppump_arr[i,*]*col.J[i,*]                             ;Rlu   
ENDFOR

;======================================
;Add collisional terms to the rate matrix
;======================================
FOR i=0,m.nctrans-1 DO BEGIN
   Cul_arr[i,*] = INTERPOL(REFORM(m.collrates[i,0:m.ntemps-1]),REFORM(m.coll_temps[0:m.ntemps-1]),col.tgas)
   Clu_arr[i,*] = Cul_arr[i,*] * m.g[m.coll_iup[i]-1]/m.g[m.coll_idown[i]-1] * $
                  EXP(-(m.energy_in_k[m.coll_iup[i]-1]-m.energy_in_k[m.coll_idown[i]-1])/col.tgas)
   R_arr[m.coll_idown[i]-1,m.coll_iup[i]-1,*] += Cul_arr[i,*]*col.rhogas  ;Rul
   R_arr[m.coll_iup[i]-1,m.coll_idown[i]-1,*] += Clu_arr[i,*]*col.rhogas  ;Rlu
ENDFOR

FOR h=0,m.np-1 DO BEGIN
   ;
   ;1st term: rate into level j; 2nd term: rate out of level j
   ;Pj = Sj (Rji x nj) - ni x Sj(Rij)
   P[h*m.nlevels:(h+1)*m.nlevels-1] = TRANSPOSE(R_arr[*,*,h])##REFORM(npop[*,h]) - TOTAL(R_arr[*,*,h],1)*npop[*,h]                   
   ;
   ;Replacing the last equilibrium
   ;equation with mass conservation to
   ;form a linearly independent system
   P[(h+1)*m.nlevels-1]   = col.abun[h] * col.rhogas[h] - TOTAL(npop[*,h])
ENDFOR

RETURN, P

END
