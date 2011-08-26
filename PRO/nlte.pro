PRO nlte, species=species, npop=npop
COMMON coll,  Cul, TCul
COMMON const, dv
COMMON mol,  nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g, collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans, partner
COMMON grid, np, z_col, tgas_col, abun_col, rhogas_col, J_col, JSED_col, nu_cont
@natconst


niter  = 20
dV     = 1d5  ;cm/s
frac    = 0.01

npop        = DBLARR(nlevels,np)
Jac         = DBLARR(nlevels,nlevels,np)
npop_new    = DBLARR(nlevels,np)
npop_iter   = DBLARR(nlevels,np,niter)
J_col       = DBLARR(nlines,np)

FOR h=0,np-1 DO BEGIN
   J_col[*,h] = INTERPOL(JSED_col[h,*],nu_cont/cc,freq)
ENDFOR

;
;Set initial level populations to LTE
FOR i=0,nlevels-1 DO BEGIN
   npop[i,*] = g[i] * exp(-(energy_in_K[i]/tgas_col))
ENDFOR
FOR h=0,np-1 DO BEGIN
   Ntot = TOTAL(npop[*,h])
   npop[*,h] = npop[*,h]/Ntot * abun_col[h] * rhogas_col[h] ;in cm^-3
ENDFOR

;
;Main iteration
FOR k=0,niter-1 DO BEGIN
   npop_iter[*,*,k] = npop
   ;
   ;Calculate the Jacobian
   dn = npop*frac
   ;
   ;i is the ith rate equation, jth level population
   FOR j=0,nlevels-1 DO BEGIN
      Pn          = P(npop)
      npopdn      = npop
      npopdn[j,*] = npopdn[j,*]+dn[j,*] 
      Pndn        = P(npopdn)
      FOR h=0,np-1 DO BEGIN
         Jac[j,*,h] = (Pndn[*,h]-Pn[*,h])/dn[j,h]
      ENDFOR
   ENDFOR
   
   FOR h=0,np-1 DO BEGIN
      npop_new[*,h] = npop[*,h] - LA_INVERT(REFORM(Jac[*,*,h]),/DOUBLE, STATUS=STATUS)##REFORM(Pn[*,h]) 
   ENDFOR
   bsubs = WHERE(FINITE(npop_new) NE 1)
   IF bsubs[0] NE -1 THEN BEGIN
      npop_new[bsubs] = npop[bsubs]*1.02
   ENDIF
;   print, (npop_new[*,1]-npop[*,1])/npop[*,1]
;   print, npop_new[*,1]
   npop = npop_new
   IF MEAN((npop_new[*,1]-npop[*,1])/npop[*,1]) LT 0.001 OR STATUS GT 0 THEN BREAK
ENDFOR 
 print, MEAN((npop_new[*,1]-npop[*,1])/npop[*,1])
   
END


