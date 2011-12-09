PRO nlte,  z_col, tgas_col, rhogas_col, abun_col, JSED_col, J_col, $
           dv, nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g, $
           collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans,$
           np, npop, ini_npop

@natconst


niter  = 8
dV     = 1d5  ;cm/s
frac    = 0.001

npop        = DBLARR(nlevels,np)
Jac         = DBLARR(nlevels,nlevels,np)
npop_new    = DBLARR(nlevels,np)
npop_iter   = DBLARR(nlevels,np,niter)

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
;Save the initial level pops
ini_npop = npop

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
      Pn          = P(npop,col,m)
      npopdn      = npop
      npopdn[j,*] = npopdn[j,*]+dn[j,*] 
      Pndn        = P(npopdn,col,m)
      FOR h=0,np-1 DO BEGIN
         Jac[j,*,h] = (Pndn[*,h]-Pn[*,h])/dn[j,h]
      ENDFOR
   ENDFOR
;   print, npop[*,1]
   FOR h=0,np-1 DO BEGIN
      npop_new[*,h] = npop[*,h] - LA_INVERT(REFORM(Jac[*,*,h]),/DOUBLE,STATUS=STATUS)##REFORM(Pn[*,h]) 
   ENDFOR
   bsubs = WHERE(FINITE(npop_new) NE 1)
   IF bsubs[0] NE -1 THEN BEGIN
      npop_new[bsubs] = npop[bsubs]*1.02
   ENDIF

   highsubs = WHERE(npop GT 100.)
   conv = ABS(MAX((npop_new[highsubs]-npop[highsubs])/npop[highsubs]))
   print, conv
   npop = npop_new
   IF conv LT 1d-12 THEN BEGIN
      PRINT, 'Converged in ' + STRTRIM(STRING(k+1),2) + ' iterations' 
      BREAK
   ENDIF
   IF STATUS GT 0 THEN BEGIN      
      PRINT, 'Warning: singular matrix detected!', STATUS
      BREAK
   ENDIF
  
ENDFOR 
   
END


