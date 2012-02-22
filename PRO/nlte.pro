PRO nlte,  z_col, tgas_col, rhogas_col, abun_col, JSED_col, J_col, $
           dv, nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g, $
           collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans,$
           np, npop, ini_npop, lte_npop

@natconst

niter  = 55
frac    = 0.001

col = {z:z_col,tgas:tgas_col,rhogas:rhogas_col,abun:abun_col,JSED:JSED_col,J:J_col}
m   = {dv:dv,nlines:nlines,nlevels:nlevels,gugl:gugl,freq:freq,iup:iup,$
       idown:idown,Aul:Aul,Bul:Bul,Blu:Blu,energy_in_k:energy_in_k,g:g,collrates:collrates,$
       coll_iup:coll_iup,coll_idown:coll_idown,coll_temps:coll_temps,ntemps:ntemps,nctrans:nctrans,np:np}

npop        = DBLARR(nlevels,np)
Jac         = DBLARR(nlevels,nlevels,np)
npop_new    = DBLARR(nlevels,np)
npop_two    = DBLARR(nlevels,np)
npop_iter   = DBLARR(nlevels,np,niter)

;
;Set initial level populations to LTE
FOR i=0,nlevels-1 DO BEGIN
   npop[i,*] = g[i] * exp(-(energy_in_K[i]/tgas_col))
ENDFOR

;
;Renormalize
FOR h=0,np-1 DO BEGIN
   Ntot = TOTAL(npop[*,h])
   npop[*,h] = npop[*,h]/Ntot * abun_col[h] * rhogas_col[h] ;in cm^-3
ENDFOR

lte_npop = npop

;
;Avoid starting with populations that are too low (0s will make a
;singular jacobian)
lsubs       = WHERE(npop LT 1d-30,lcount)
IF lcount GT 0 THEN npop[lsubs] = 1d-30

;
;Save the initial level pops
ini_npop = npop

;
;Main iteration
lam1 = FLTARR(np) + 1.
nf=0
FOR k=0,niter-1 DO BEGIN
   npop_iter[*,*,k] = npop
   ;
   ;Calculate the Jacobian
   dn = npop*frac
   zsub = WHERE(dn EQ 0,nzsub)
   IF nzsub GT 0 THEN dn[zsub] = 1d-5

   ;
   ;i is the ith rate equation, jth level population
   Pn = P(npop,col,m)
   FOR j=0,nlevels-1 DO BEGIN
      npopdn_neg      = npop
      npopdn_neg[j,*] = npopdn_neg[j,*]-dn[j,*] 
      Pndn_neg        = P(npopdn_neg,col,m)
      npopdn_pos      = npop
      npopdn_pos[j,*] = npopdn_pos[j,*]+dn[j,*] 
      Pndn_pos        = P(npopdn_pos,col,m)
      FOR h=0,np-1 DO BEGIN
         Jac[j,*,h] = (Pndn_pos[*,h]-Pndn_neg[*,h])/(2d0*dn[j,h])
         bsubs = WHERE(FINITE(Jac[j,*,h]) NE 1,bcount) 
         IF bcount GT 0 THEN stop
      ENDFOR
   ENDFOR

   FOR h=0,np-1 DO BEGIN
      newton = LA_INVERT(REFORM(Jac[*,*,h]),/DOUBLE,STATUS=STATUS)##REFORM(Pn[*,h]) 
      npop_new[*,h] = npop[*,h] - newton
   ENDFOR

   Pn_new = P(npop_new, col, m)
   ;
   ;Did the newton step fail to be closer to the solution?
   IF TOTAL(Pn_new^2) GT TOTAL(Pn^2) THEN BEGIN
      print,'newton failed', nf++
      lam  = FLTARR(np)
      hosubs = WHERE(lam1 NE 1,count)
      IF count EQ 0 THEN BEGIN
         FOR h=0,np-1 DO BEGIN
            gp_0 = TOTAL((REFORM(Jac[*,*,h])##REFORM(Pn[*,h])) * REFORM(newton))
            g_0  = 0.5d0*TOTAL(Pn[*,h]^2) 
            g_1  = 0.5d0*TOTAL(Pn_new[*,h]^2) 
            lam[h]  = -gp_0/(2d0*(g_1-g_0-gp_0))
         ENDFOR
      ENDIF ELSE BEGIN
         FOR h=0,np-1 DO BEGIN
            g_1  = 0.5d0*TOTAL(Pn_new[*,h]^2)
            a    = (1/(lam1[h]-lam2[h]))*((1/lam1[h]^2)*(g_1-gp_0*lam1[h]-g_0)-(1/lam2[h]^2)*(g_2-gp_0*lam2[h]-g_0))
            b    = (1/(lam1[h]-lam2[h]))*((-lam2[h]/(lam1[h]^2))*(g_1-gp_0*lam1[h]-g_0)+(lam1[h]/(lam2[h]^2))*(g_2-gp_0*lam2[h]-g_0))
            lam[h]  = (1/(3*a))*(-b+sqrt(b^2+3*a*gp_0))
         ENDFOR
      ENDELSE
      lsubs = WHERE(lam LT 0.1,lcount)
      hsubs = WHERE(lam GT 0.5,hcount)
      IF lcount GT 0 THEN    lam[lsubs] = 0.1*lam1[lsubs]
      IF hcount GT 0 THEN    lam[hsubs] = 0.5*lam1[hsubs]
      FOR h=0,np-1 DO BEGIN
         newton = LA_INVERT(REFORM(Jac[*,*,h]),/DOUBLE,STATUS=STATUS)##REFORM(Pn[*,h]) 
         npop_new[*,h] = npop[*,h] - lam[h]*newton
      ENDFOR
      lam2   = lam1
      lam1   = lam
      g_2     = g_1
   ENDIF
  FOR j=0,nlevels-1 DO BEGIN
      FOR h=0,np-1 DO BEGIN
         IF FINITE(npop_new[j,h]) NE 1 THEN BEGIN
            IF k GT 0 THEN npop_new[j,h] = npop_iter[j,h,k-1]
            IF k EQ 0 THEN npop_new[j,h] = lte_pop
           ; print,j,h,npop_new[j,h]
         ENDIF
      ENDFOR
   ENDFOR

;   bsubs = WHERE(FINITE(npop_new) NE 1)
;   IF bsubs[0] NE -1 THEN BEGIN
;      npop_new[bsubs] = npop[bsubs]*1.02
;   ENDIF

   highsubs = WHERE(npop GT 100.,highcount)
   IF highcount GT 0 THEN BEGIN
      conv = ABS(MAX((npop_new[highsubs]-npop[highsubs])/npop[highsubs]))
      print, conv
      npop = npop_new
      IF conv LT 1d-4 THEN BEGIN
         PRINT, 'Converged in ' + STRTRIM(STRING(k+1),2) + ' iterations' 
         BREAK
      ENDIF
      IF STATUS GT 0 THEN BEGIN      
         PRINT, 'Warning: singular matrix detected!', STATUS
         BREAK
      ENDIF
   ENDIF
   if k eq niter-1 then stop   
ENDFOR 

;return fractional populations
;FOR h=0,np-1 DO BEGIN
;   npop[*,h]     = npop[*,h]/(abun_col[h] * rhogas_col[h])
;   ini_npop[*,h] = ini_npop[*,h]/(abun_col[h] * rhogas_col[h])
;   lte_npop[*,h] = lte_npop[*,h]/(abun_col[h] * rhogas_col[h])
;ENDFOR

END


