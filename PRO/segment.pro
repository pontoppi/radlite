PRO segment, theta, x0, y0, D, N=N, hex=hex, eta=eta, ksi=ksi, max_ksi=max_ksi

;
;Currently not using the shifts, since that can be applied after the
;heavy calculation.

;D = 800d0 ;cm

IF NOT KEYWORD_SET(N) THEN N=1500

;Define complex i
i = COMPLEX(0,1)
g = DBLARR(N,N)

IF NOT KEYWORD_SET(max_ksi) THEN max_ksi = 0.007
min_ksi = -max_ksi
ksi_p = (max_ksi-min_ksi)*findgen(N)/(N-1)+min_ksi
eta_p = ksi_p

;
;Rotate the pupil plane segment
ksi = ksi_p*cos(theta)+eta_p*sin(theta)
eta = -eta_p*sin(theta)+eta_p*cos(theta)

;
;Calculate analytical expression for trapezoid (Sabatke et al. 2005,
;Applied optics, 44, 1360).

FOR h=0,N-1 DO BEGIN
   FOR k=0,N-1 DO BEGIN
      g[h,k] = EXP(-i*!pi*D*(2d0*eta[k]/sqrt(3d0)+ksi[h])) / (4d0*!pi^2d0*(eta[k]^3d0-3d0*eta[k]*ksi[h]^2d0)) * $
         ((sqrt(3d0)*eta[k]-3d0*ksi[h]) * (EXP(i*!pi*D*sqrt(3d0)*eta[k])-EXP(i*!pi*D*((4d0/sqrt(3d0))*eta[k]+ksi[h]))) + $
         (sqrt(3d0)*eta[k]+3d0*ksi[h])*(EXP(i*!pi*D*eta[k]/sqrt(3d0))-EXP(i*!pi*D*ksi[h])))
   ENDFOR
ENDFOR

;
;Construct hexagonal
hex = g + ROTATE(g,5)



END
