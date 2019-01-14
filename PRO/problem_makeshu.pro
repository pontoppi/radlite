
FUNCTION DIFFERENTIAL, X, Y

   v     = Y[0]
   alpha = Y[1]
   DVDX = (alpha*(x-v)-2d0/x)*(x-v)/((x-v)^2d0-1d0)
   DALPHADX = (alpha-(2d0/x)*(x-v))*(x-v)*alpha/((x-v)^2d0-1)
 
   RETURN, [DVDX, DALPHADX]

END

FUNCTION make_shu,r,time,csenv=csenv,Aenv=Aenv

@natconst.pro
;Make initial density structure, rho(r,0):

Nx=10000.
Nr=10000d0
vs = fltarr(Nx)
alphas = fltarr(Nx)
Xs     = fltarr(Nx)
ms     = fltarr(Nx)
X0 = 1d0
H  = -X0/(Nx+1)

;
;Basic definitions

IF NOT KEYWORD_SET(Aenv) THEN BEGIN
   A    = 2.003
ENDIF ELSE BEGIN
   A    = Aenv
ENDELSE   

IF NOT KEYWORD_SET(csenv) THEN csenv = 0.2d5 ;cm/s

m_av = mu*mp
t0   = csenv                   ;s

ini_r   = findgen(Nr)/(Nr-1) * (MAX(r)-MIN(r)) + MIN(r)
ini_rho = csenv^2d0*A/(4d0*!pi*GG)*ini_r^(-2d0)
alpha0 = 4d0*!pi*GG*t0^2d0*ini_rho
v0     = fltarr(Nr)+csenv*1d-1

X = X0                            ;r[i]/(cs*t0)
Y = [-(A-2.)/X0,A/X0^2d0]         ;alpha0[i]]
FOR j=0L,Nx-1 DO BEGIN
    dydx = differential(X,Y)
    Y = RK4(Y,dydx,X,H,'differential',/DOUBLE)
    vs[j] = Y[0]
    alphas[j] = Y[1]
    Xs[j] = X
    ms[j] = X^2d0*Y[1]*(x-Y[0])
    X = X+H
ENDFOR

shu_r = Xs*csenv*time
shu_rho = alphas/(4d0*!pi*GG*time^2d0)
shu_v = vs*csenv


rho = INTERPOL(shu_rho,shu_r,r)
v   = INTERPOL(shu_v,shu_r,r)
outsubs = WHERE(r GT MAX(shu_r))
IF outsubs[0] NE -1 THEN BEGIN
   rho[outsubs] = INTERPOL(ini_rho,ini_r,r[outsubs])
   v[outsubs]   = 0.0
ENDIF

RETURN, {v:v,rho:rho,r:r}


END
