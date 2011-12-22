PRO nlteL,  z_col, tgas_col, rhogas_col, abun_col, JSED_col, J_col, $
           dv, nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g, $
           collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans,$
           np, npop, ini_npop

@natconst

niter  = 20
frac    = 0.00001

col = {z:z_col,tgas:tgas_col,rhogas:rhogas_col,abun:abun_col,JSED:JSED_col,J:J_col}
m   = {dv:dv,nlines:nlines,nlevels:nlevels,gugl:gugl,freq:freq,iup:iup,$
       idown:idown,Aul:Aul,Bul:Bul,Blu:Blu,energy_in_k:energy_in_k,g:g,collrates:collrates,$
       coll_iup:coll_iup,coll_idown:coll_idown,coll_temps:coll_temps,ntemps:ntemps,nctrans:nctrans,np:np}

npop        = DBLARR(nlevels,np)
Jac         = DBLARR(nlevels,nlevels,np)
lnpop_new   = DBLARR(nlevels,np)
lnpop_iter  = DBLARR(nlevels,np,niter)

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
lnpop = ALOG10(npop)
FOR k=0,niter-1 DO BEGIN
   lnpop_iter[*,*,k] = lnpop
   ;
   ;Calculate the Jacobian
   dln = lnpop*frac
   ;
   ;i is the ith rate equation, jth level population
   FOR j=0,nlevels-1 DO BEGIN
      Pln           = PL(lnpop,col,m)
      lnpopdln      = lnpop
      lnpopdln[j,*] = lnpopdln[j,*]+dln[j,*] 
      Plndln        = PL(lnpopdln,col,m)
      FOR h=0,np-1 DO BEGIN
         Jac[j,*,h] = (Plndln[*,h]-Pln[*,h])/dln[j,h]
      ENDFOR
   ENDFOR

   FOR h=0,np-1 DO BEGIN
      lnpop_new[*,h] = lnpop[*,h] - LA_INVERT(REFORM(Jac[*,*,h]),/DOUBLE,STATUS=STATUS)##REFORM(Pln[*,h]) 
   ENDFOR
   bsubs = WHERE(FINITE(lnpop_new) NE 1)
   IF bsubs[0] NE -1 THEN BEGIN
      lnpop_new[bsubs] = lnpop[bsubs]*1.02
   ENDIF

   highsubs = WHERE(lnpop GT 0.)
   conv = ABS(MAX((lnpop_new[highsubs]-lnpop[highsubs])/lnpop[highsubs]))
   print, conv
   lnpop = lnpop_new
   IF conv LT 1d-8 THEN BEGIN
      PRINT, 'Converged in ' + STRTRIM(STRING(k+1),2) + ' iterations' 
      BREAK
   ENDIF
   IF STATUS GT 0 THEN BEGIN      
      PRINT, 'Warning: singular matrix detected!', STATUS
      BREAK
   ENDIF
  
ENDFOR 
npop = 10^lnpop
   
END


