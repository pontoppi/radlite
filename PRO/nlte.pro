PRO nlte,  z_col, tgas_col, rhogas_col, abun_col, JSED_col, J_col, $
           dv, nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g, $
           collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans,$
           np, npop, ini_npop, lte_npop

@natconst

niter   = 10
frac    = 0.0001
minlam  = 0.1
ALF     = 1d-4

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

;Avoid starting with populations that are too low (0s will make a
;singular jacobian)
lsubs       = WHERE(npop LT 1d-10,lcount)
IF lcount GT 0 THEN npop[lsubs] = 1d-10

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

      newton = -1*LA_INVERT(REFORM(Jac[*,*,h]),/DOUBLE,STATUS=STATUS)##REFORM(Pn[*,h]) 
      g_0    = 0.5d0*TOTAL(Pn[*,h]^2) 
      gp_0   = TOTAL((REFORM(Jac[*,*,h])##REFORM(Pn[*,h])) * REFORM(newton))
      
      alam   = 1.
      WHILE 1 DO BEGIN
         npop_new[*,h] = npop[*,h] + alam*newton
         Pn_new = P(npop_new, col, m)
         g_1    = 0.5d0*TOTAL(Pn_new[*,h]^2) 
         IF (alam LT minlam) THEN BEGIN
            npop_new[*,h] = npop[*,h]
            BREAK
         ENDIF 
         IF (g_1 LE g_0+ALF*alam*gp_0) THEN BEGIN
            BREAK
         ENDIF ELSE BEGIN
            IF alam EQ 1. THEN BEGIN ;First time.
               tmplam  = -gp_0/(2d0*(g_1-g_0-gp_0))
            ENDIF ELSE BEGIN
               rhs1 = g_1-g_0-alam*gp_0
               rhs2 = g_2-g_0-alam2*gp_0
               a    = (rhs1/alam^2.-rhs2/alam2^2.)/(alam-alam2)
               b    = (-alam2*rhs1/alam^2.+alam*rhs2/alam2^2.)/(alam-alam2)
               IF a EQ 0. THEN BEGIN
                  tmplam = -gp_0/(2.*b)
               ENDIF ELSE BEGIN
                  disc = b^2.-3.*a*gp_0
                  IF disc LT 0. THEN BEGIN
                     tmplam = 0.5*alam
                  ENDIF ELSE BEGIN
                     IF b LE 0. THEN BEGIN
                        tmplam =  (-b+sqrt(disc))/(3.0*a)
                     ENDIF ELSE BEGIN
                        tmplam = -gp_0/(b+SQRT(disc))
                     ENDELSE
                  ENDELSE
               ENDELSE
               IF tmplam GT 0.5*alam THEN tmplam = 0.5*alam
            ENDELSE
            if finite(tmplam) NE 1 then stop
            alam2 = alam
            g_2   = g_1
            alam  = MAX([tmplam, 0.1*alam])
         ENDELSE
      ENDWHILE
   ENDFOR
   
   highsubs = WHERE(npop_new GT 1.,highcount)
   
   IF highcount GT 0 THEN BEGIN
      conv = ABS(MAX((npop_new[highsubs]-npop[highsubs])/npop[highsubs]))
      print, conv
      npop = npop_new
      IF conv LT 1d-2 THEN BEGIN
         PRINT, 'Converged in ' + STRTRIM(STRING(k+1),2) + ' iterations' 
         BREAK
      ENDIF
      IF STATUS GT 0 THEN BEGIN
         PRINT, 'Warning: singular matrix detected!', STATUS
         frac = frac * 1.1
         
;         BREAK
      ENDIF
   ENDIF
ENDFOR 

END


