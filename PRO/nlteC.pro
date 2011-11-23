PRO nlteC, species=species, npop=npop, ini_npop=ini_npop
COMMON coll,  Cul, TCul
COMMON const, dv
COMMON mol,  nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g, collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans, partner
COMMON grid, np, z_col, tgas_col, abun_col, rhogas_col, J_col, JSED_col, nu_cont
@natconst


niter  = 8
dV     = 1d5  ;cm/s
frac    = 0.001

npop        = DBLARR(nlevels*np)
Jac         = DBLARR(nlevels*np,nlevels*np)
npop_new    = DBLARR(nlevels*np)
npop_iter   = DBLARR(nlevels*np,niter)

;
;Set initial level populations to LTE
FOR h=0,np-1 DO BEGIN
   npop[h*nlevels:(h+1)*nlevels-1] = g * exp(-(energy_in_K/tgas_col[h]))
   Ntot                            = TOTAL(npop[h*nlevels:(h+1)*nlevels-1])
   npop[h*nlevels:(h+1)*nlevels-1] = npop[h*nlevels:(h+1)*nlevels-1]/Ntot * abun_col[h] * rhogas_col[h] ;in cm^-3
ENDFOR
;
;Save the initial level pops
ini_npop = npop

;
;Main iteration
FOR k=0,niter-1 DO BEGIN
   npop_iter[*,k] = npop
   ;
   ;Calculate the Jacobian
   dn = npop*frac
   ;
   ;i is the ith rate equation, jth level population
;   FOR j=0,nlevels-1 DO BEGIN
;      Pn          = P(npop)
;      npopdn      = npop
;      npopdn[j,*] = npopdn[j,*]+dn[j,*] 
;      Pndn        = P(npopdn)
;      FOR h=0,np-1 DO BEGIN
;         Jac[j,*,h] = (Pndn[*,h]-Pn[*,h])/dn[j,h]
;      ENDFOR
;   ENDFOR
   FOR j=0,np*nlevels-1 DO BEGIN
      Pn         = PC(npop)
      npopdn     = npop
      npopdn[j] += dn[j]
      Pndn       = PC(npopdn)      
      Jac[j,*]   = (Pndn-Pn)/dn[j]
   ENDFOR
   
   npop_new = npop - LA_INVERT(REFORM(Jac),/DOUBLE,STATUS=STATUS)##REFORM(Pn) 
<<<<<<< local
   
=======
  
>>>>>>> other
   bsubs = WHERE(FINITE(npop_new) NE 1)
   IF bsubs[0] NE -1 THEN BEGIN
      npop_new[bsubs] = npop[bsubs]*1.02
   ENDIF

;   bsubs = WHERE(npop_new EQ 0)
;   npop_new[bsubs] = 1d-30
   highsubs = WHERE(npop GT 1.)
   conv = ABS(MAX((npop_new[highsubs]-npop[highsubs])/npop[highsubs]))
;   print, conv
   print, npop_new[highsubs]-npop[highsubs])/npop[highsubs]
   npop = npop_new
   IF conv LT 1d-12 THEN BEGIN
      PRINT, 'Converged in ' + STRTRIM(STRING(k+1),2) + ' iterations' 
      BREAK
   ENDIF
;   IF STATUS GT 0 THEN BEGIN      
;      PRINT, 'Warning: singular matrix detected!', STATUS
;      BREAK
;   ENDIF
  
ENDFOR 
;print, MAX((npop_new-npop)/npop[*,1])
   
END


