PRO nlteC, z_col, tgas_col, rhogas_col, abun_col, JSED_col, J_col, $
           dv, nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g, $
           collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans,$
           np, npop, ini_npop

@natconst

niter  = 16
frac   = 0.001

col = {z:z_col,tgas:tgas_col,rhogas:rhogas_col,abun:abun_col,JSED:JSED_col,J:J_col}
m   = {dv:dv,nlines:nlines,nlevels:nlevels,gugl:gugl,freq:freq,iup:iup,$
       idown:idown,Aul:Aul,Bul:Bul,Blu:Blu,energy_in_k:energy_in_k,g:g,collrates:collrates,$
       coll_iup:coll_iup,coll_idown:coll_idown,coll_temps:coll_temps,ntemps:ntemps,nctrans:nctrans,np:np}

npop        = DBLARR(nlevels*np)
Jac         = DBLARR(nlevels*np,nlevels*np)
npop_new    = DBLARR(nlevels*np)
npop_iter   = DBLARR(nlevels*np,niter)

;
;Set initial level populations to LTE
FOR h=0,np-1 DO BEGIN
   npop[h*nlevels:(h+1)*nlevels-1] = m.g * exp(-(m.energy_in_K/col.tgas[h]))
   Ntot                            = TOTAL(npop[h*nlevels:(h+1)*m.nlevels-1])
   npop[h*nlevels:(h+1)*nlevels-1] = npop[h*m.nlevels:(h+1)*m.nlevels-1]/Ntot * col.abun[h] * col.rhogas[h] ;in cm^-3
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

   FOR j=0,np*nlevels-1 DO BEGIN
      Pn         = PC(npop,col,m)
      npopdn     = npop
      npopdn[j] += dn[j]
      Pndn       = PC(npopdn,col,m)      
      Jac[j,*]   = (Pndn-Pn)/dn[j]
   ENDFOR
   
   npop_new = npop - LA_INVERT(REFORM(Jac),/DOUBLE,STATUS=STATUS)##REFORM(Pn) 
  
   bsubs = WHERE(FINITE(npop_new) NE 1)
   IF bsubs[0] NE -1 THEN BEGIN
      npop_new[bsubs] = npop[bsubs]*1.02
   ENDIF

   highsubs = WHERE(npop GT 1.)
   conv = ABS(MAX((npop_new[highsubs]-npop[highsubs])/npop[highsubs]))
   print, conv
   npop = npop_new
   IF conv LT 1d-8 THEN BEGIN
      PRINT, 'Converged in ' + STRTRIM(STRING(k+1),2) + ' iterations' 
      BREAK
   ENDIF
;   IF STATUS GT 0 THEN BEGIN      
;      PRINT, 'Warning: singular matrix detected!', STATUS
;      BREAK
;   ENDIF
  
ENDFOR 

END


