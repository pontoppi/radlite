;======================================================================
;      THESE ARE THE PAH OPACITY ROUTINES WRITTEN BY RUUD VISSER
;======================================================================


pro check_range,a,w,name=name,message=message,nostop=nostop

;a = available
;w = wanted

if not keyword_set(name) then name = 'Variable'
if not keyword_set(message) then message = 'Sorry ...'

if not ((min(w,/nan) ge min(a,/nan)) and (max(w,/nan) le max(a,/nan))) then begin
    print,'------------------------------------------------------------------------'
    print,'!! '+name+' requested between',min(w,/nan),' and',max(w,/nan)
    print,'   but only available between',min(a,/nan),' and',max(a,/nan),'  ....'
    print,'   '+message
    if not keyword_set(nostop) then stop else print,'  Continue anyway...'
    print,'------------------------------------------------------------------------'

endif

end


function get_index,xdrfile,wave_in
; purpose:
;  Get the refractive index from an generic xdr file
;  produced in the OPTICAL_DATA directory

nwav_in = n_elements(wave_in)


restore,xdrfile
; wav : wavelengths in microns
; m : complex optical index
; eps : complex dielectric index
; head : header containing informations about the optical constants

check_range,wav,wave_in,name='Optical Indexes'
re = INTERPOL(float(m),wav,wave_in) ; real part of m
im = INTERPOL(imaginary(m),wav,wave_in) ; imaginary part of m

;-- In case the interpolation produced negative optical indexes:
idre = where(re lt 0.)
if (size(idre))(0) gt 0 then re(idre) = 0.
idim = where(im lt 0.)
if (size(idim))(0) gt 0 then im(idim) = 0.

if nwav_in gt 1 then begin
  index = reform(complex(re,im))
  return,index
endif else return,complex(re,im)


end


pro mie,xgrain,mgrain,Qsca=Qsca,Qabs=Qabs,gmie=gmie,Qpr=Qpr,an=an,bn=bn, $
        phi=phi,phase=phase
; INPUTS :
;    xgrain = 2.*!Pi*radius / lambda : size parameter
;    mgrain : optical index
;
; OPTIONAL OUTPUTS:
;   Qsca
;   Qabs
;   gmie
;   Qpr
;   an
;   bn

;-- Constant values:
Qcste = double(2.d0/ (xgrain*xgrain))
mu = 1.d0 ; magnetic permitivity of the medium surrounding the grain
m = 1.d0 ; optical indice of the surrounding medium of the grain (vaccum)
mu_grain = 1.d0 ; magnetic permitivity of the grain
ro = xgrain * m
ro_grain = xgrain * mgrain
abcste = m*mu / mgrain / mu_grain
nstop = ceil(ro + 4.d0 * ro^(1.d0/3.d0)+4.d0)

;-- Computation  of Dn(ro_grain) put into D(n)
nmax = max([nstop,ceil(abs(ro_grain))]) + 15L
D = dcomplexarr(nmax+1,/nozero)
D(nmax) = dcomplex(0.d0,0.d0)
for i=nmax*1L,2,-1 do begin
    C = double(i)/ro_grain
    D(i-1) =  C - 1.d0/(D(i)+C)
endfor

;--- Computation of  an et bn:
;     xksi = xpsi - i*chi:
;     here: xpsi0=Psi(-1,ro),xpsi1=Psi(0,ro),chi0=Chi(-1,ro),chi1=Chi(0,ro):
xpsi0 = cos(ro)
xpsi1 = sin(ro)
chi0 = -xpsi1
chi1 = xpsi0
xksi0 = dcomplex(xpsi0,-chi0)
xksi1 = dcomplex(xpsi1,-chi1)

an = dcomplexarr(nstop+1,/nozero)
bn = dcomplexarr(nstop+1,/nozero)
for i=1L,nstop do begin
    x = double(i)
    xpsi = (2.d0*x-1.d0)/ro*xpsi1 - xpsi0
    xpsi0 = xpsi1
    xpsi1 = xpsi
    chi = (2.d0*x-1.d0)/ro*chi1 - chi0
    chi0 = chi1
    chi1 = chi
    xksi0 = dcomplex(xpsi0,-chi0)
    xksi1 = dcomplex(xpsi1,-chi1)

    ;      here: xpsi0=Psi(i-1,ro), xpsi1=Psi(i,ro), chi0=Chi(i-1,ro)
    ;          chi1=Chi(i,ro), an(i) = a(i,ro1,ro), bn(i) = b(i,ro1,ro):
    
    anum = (abcste*D(i)+x/ro)*xpsi1-xpsi0
    aden = (abcste*D(i)+x/ro)*xksi1-xksi0
    an(i) = anum / aden
    
    bnum = (D(i)/abcste+x/ro)*xpsi1-xpsi0
    bden = (D(i)/abcste+x/ro)*xksi1-xksi0
    bn(i) = bnum / bden
endfor

;-- Computation of  Qsca, Qabs, Qext, Qpr and gmie:
Qsca = 0.d0
Qext = 0.d0
gmie = 0.d0
nstop1 = nstop-1L
x = dindgen(nstop+1L)
for  i=1L,nstop1 do begin
    x = double(i)
    Qsca = Qsca+(2.*x+1.)*(an(i)*conj(an(i))+bn(i)*conj(bn(i)))
    Qext = Qext+(2.d0*x+1.d0)*(double(an(i)+bn(i)))
    g0 = an(i)*conj(an(i+1))+bn(i)*conj(bn(i+1))
    g1 = x*(x+2.d0)*double(g0)
    g2 = (2.d0*x+1.d0) * double(an(i)*conj(bn(i))) / x
    gmie = gmie + (g1 + g2) / (x+1.d0)
endfor
Qsca = Qsca*Qcste
Qext = Qext*Qcste
Qabs = Qext - Qsca
gmie = double(gmie * 4.d0/ro/ro / Qsca)
Qpr = Qext - gmie*Qsca

;-- Calculation of the phase function:
if keyword_set(phi) then begin
    nphi = n_elements(phi)
    phase = dblarr(nphi)
    for i=0L,nphi-1 do begin
        xcos = cos(phi(i))
        pi0 = 0.d0
        pi1 = 1.d0
        tau = xcos*pi1
        s1 = 3.d0/2.d0 * (an(1)*pi1 + bn(1)*tau)
        s2 = 3.d0/2.d0 * (an(1)*tau + bn(1)*pi1)
        for j=2,nstop1 do begin
            x = double(j)
            pi = ((2.d0*x-1.d0)*xcos*pi1 - x*pi0)/(x-1.d0)
            tau = x*xcos*pi - (x+1.d0)*pi1
            xn = (2.d0*x+1.d0)/x/(x+1.d0)
            s1 = s1 + xn*(an(j)*pi + bn(j)*tau)
            s2 = s2 + xn*(an(j)*tau + bn(j)*pi)
            pi0 = pi1
            pi1 = pi
        endfor
        phase(i) = (s1*conj(s1)+s2*conj(s2))/(xgrain*xgrain*2.d0*!Pi*Qsca)
    endfor
endif


end


;Name: CROSS2.PRO
;Author: Ruud Visser
;Date: July 1, 2004
;Description: This IDL procedure calculates the cross section per
;carbon atom for a set of reference energies according to Li & Draine
;(2001).
;To calculate the cross section for a single energy, call the function
;with size_ref = energy in erg. The procedure will recognize that
;size_ref < 1 and will subsequently treat size_ref as energy.
;References: LD - Li & Draine, 2001, ApJ 554:778-802
;Parameters: N_C - number of carbon atoms
;            N_H - number of hydrogen atoms
;            Z   - charge of the PAH
;            size_ref - size of the reference energy array OR the
;                       energy (see above)
;            E_ref - reference energies array
;            C_ref - reference cross sections

PRO cross2,N_C,N_H,Z,size_ref,E_ref=E_ref,C_ref=C_ref

;Ensure that N_C and N_H are DOUBLEs
N_C = DOUBLE(N_C)
N_H = DOUBLE(N_H)

;Global constants and other parameters
h = 6.626e-27                   ;erg s
c = 2.998e10                    ;cm s^-1
k = 1.381e-16                   ;erg K^-1
hxc = h * c                     ;erg cm (hxc = h x c = h times c)
xi_PAH = 0.99                   ;Weighting factor for including graphite properties in the PAH absorption cross section (LD:2,3)

;Hydrogenation fraction
H_C = N_H/N_C

;PAH radius in cm
radius = 1e-8 * (N_C/0.468)^(1./3)

;Number of benzenoid rings - LD:A3+
IF N_C LT 40 THEN Mbr = 0.3*N_C $
ELSE Mbr = 0.4*N_C

;Cutoff wavelength - LD:A3
IF Z LE 0 THEN lambda_c = 1 / (3.804 * Mbr^(-0.5) + 1.052) $ ;neutral
ELSE lambda_c = 1 / (2.282 * Mbr^(-0.5) + 0.889) ;ionized

;Initialisation of Drude profile parameters - LD:Table 1
gamma  = [0.195, 0.217, 0.012, 0.032, 0.091, 0.047, 0.018, 0.025, 0.024, 0.010, 0.036, 0.038, 0.046, 0.69]
lambda_central  = [0.0722, 0.2175, 3.3, 6.2, 7.7, 8.6, 11.3, 11.9, 12.7, 16.4, 18.3, 21.2, 23.1, 26.0]
IF Z LE 0 THEN sigma_int = [7.97e5, 1.23e5, 1.97*H_C, 0.196*3., 0.609*2., 0.347*2.*H_C, 4.27*H_C/3., 0.727*H_C/3., 1.67*H_C/3., 0.0552, 0.0604, 0.108, 0.0278, 0.152] $
ELSE sigma_int = [7.97e5, 1.23e5, 0.447*H_C, 1.57*3., 5.48*2., 2.42*2.*H_C, 4.00*H_C/3., 0.614*H_C/3., 1.49*H_C/3., 0.0552, 0.0604, 0.108, 0.0278, 0.152]
sigma_int = sigma_int * 1e-18
IF size_ref GE 1 THEN BEGIN
    gamma2 = REBIN(REFORM(gamma,1,14),size_ref,14) ;2D-array
    lambda_central2 = REBIN(REFORM(lambda_central,1,14),size_ref,14) ;2D-array
    sigma_int2 = REBIN(REFORM(sigma_int,1,14),size_ref,14) ;2D-array
ENDIF ELSE BEGIN
    gamma2 = gamma
    lambda_central2 = lambda_central
    sigma_int2 = sigma_int
ENDELSE

;Initialise wavelength array for which reference C_abs and u_E will be
;calculated
IF size_ref GE 1 THEN BEGIN
    lambda_ref = DINDGEN(size_ref) * 6./(size_ref-1) - 3
    lambda_ref = 10.^lambda_ref
    lambda_ref[0] = lambda_ref(0) * (1+1./size_ref)
    lambda_ref(size_ref-1) = lambda_ref(size_ref-1) / (1+1./size_ref)
    lambda2 = REBIN(REFORM(lambda_ref,size_ref,1),size_ref,15) ;2D-array
ENDIF ELSE BEGIN
    lambda_ref = 1e4 * h * c / size_ref
    lambda2 = lambda_ref
ENDELSE

;Get optical indices from a data file
m_par = get_index('Gra_Epar.xdr',lambda_ref)
m_per = get_index('Gra_Eper.xdr',lambda_ref)

;Convert wavelengths to energy and frequency
E_ref = 1e4*hxc/lambda_ref

;Cutoff function - LD:A2
y = lambda_c / lambda_ref
cutoff = (1/!Pi) * ATAN(1000.*(y-1)^3/y) + 0.5

;Drude profiles - LD:12
S = (2./!Pi) * gamma2 * 1e-4 * lambda_central2 * sigma_int2 / ((lambda2/lambda_central2 - lambda_central2/lambda2)^2 + gamma2^2) ;2D-array

;Conversion from wavelength to wavenumbers
x = 1./lambda_ref

;Graphite cross section (using Mie theory) - LD:Sec.2.2
IF size_ref GE 1 THEN BEGIN
    id = INDGEN(size_ref)
    idn = -1
    xg = 2. * !Pi * radius / (lambda_ref*1e-4)
    Q_par = DBLARR(size_ref)
    Q_per = DBLARR(size_ref)
    FOR i=0,(SIZE(id))[1]-1 DO BEGIN
        mie,xg[id[i]],m_par[i],Qabs=Qabs
        Q_par[id[i]] = FLOAT(Qabs)
        mie,xg[id[i]],m_per[i],Qabs=Qabs
        Q_per[id[i]] = FLOAT(Qabs)
    ENDFOR
    C_gra = DBLARR(size_ref,/nozero)
    IF (SIZE(id))[0] NE 0 THEN C_gra[id] = (1./N_C) * !Pi*radius^2 * (Q_par[id]/3. + 2.*Q_per[id]/3.)
    IF (SIZE(idn))[0] NE 0 THEN C_gra[idn] = (1./N_C) * !Pi*radius^2 * (Q_par[idn]/3. + 2.*Q_per[idn]/3.)
ENDIF ELSE BEGIN
    xg = 2. * !Pi * radius / (lambda_ref*1e-4)
    mie,xg,m_par,Qabs=Qabs
    Q_par = FLOAT(Qabs)
    mie,xg,m_per,Qabs=Qabs
    Q_per = FLOAT(Qabs)
    C_gra = (1./N_C) * !Pi*radius^2 * (Q_par/3. + 2.*Q_per/3.)
ENDELSE

;PAH crosss section per carbon atom - LD:5-11
id0 = WHERE(x GT 17.25)
id1 = WHERE(x GT 15. AND x LE 17.25)
id2 = WHERE(x GT 10. AND x LE 15.)
id3 = WHERE(x GT 7.7 AND x LE 10.)
id4 = WHERE(x GT 5.9 AND x LE 7.7)
id5 = WHERE(x GT 3.3 AND x LE 5.9)
id6 = WHERE(x GT 0.05 AND x LE 3.3)
id7 = WHERE(x LE 0.05)
IF size_ref GE 1 THEN BEGIN
    C_PAH = DBLARR(size_ref,/nozero)
    IF (SIZE(id0))[0] NE 0 THEN C_PAH[id0] = C_gra[id0]
    IF (SIZE(id1))[0] NE 0 THEN C_PAH[id1] = 1e-18 * (126.0 - 6.4943 * x[id1])
    IF (SIZE(id2))[0] NE 0 THEN C_PAH[id2] = 1e-18 * (-3.0 + 1.35 * x[id2]) + (REFORM(S[*,0]))[id2]
    IF (SIZE(id3))[0] NE 0 THEN C_PAH[id3] = 1e-18 * (66.302 - 24.367 * x[id3] + 2.950 * (x[id3])^2 - 0.1057 * (x[id3])^3)
    IF (SIZE(id4))[0] NE 0 THEN C_PAH[id4] = 1e-18 * (1.8687 + 0.1905 * x[id4] + 0.4175 * (x[id4] - 5.9)^2 + 0.04370 * (x[id4] - 5.9)^3) + (REFORM(S[*,1]))[id4]
    IF (SIZE(id5))[0] NE 0 THEN C_PAH[id5] = 1e-18 * (1.8687 + 0.1905 * x[id5]) + (REFORM(S[*,1]))[id5]
    IF (SIZE(id6))[0] NE 0 THEN C_PAH[id6] = 3.458 * 10.^(-17 - 3.431/x[id6]) * cutoff[id6] + (TOTAL(S[*,2:13],2))[id6]
    IF (SIZE(id7))[0] NE 0 THEN C_PAH[id7] = (TOTAL(S[*,2:13],2))[id7]
ENDIF ELSE BEGIN
    IF (SIZE(id0))[0] NE 0 THEN C_PAH = C_gra
    IF (SIZE(id1))[0] NE 0 THEN C_PAH = 1e-18 * (126.0 - 6.4943 * x)
    IF (SIZE(id2))[0] NE 0 THEN C_PAH = 1e-18 * (-3.0 + 1.35 * x) + S[0]
    IF (SIZE(id3))[0] NE 0 THEN C_PAH = 1e-18 * (66.302 - 24.367 * x + 2.950 * x^2 - 0.1057 * x^3)
    IF (SIZE(id4))[0] NE 0 THEN C_PAH = 1e-18 * (1.8687 + 0.1905 * x + 0.4175 * (x - 5.9)^2 + 0.04370 * (x - 5.9)^3) + S[1]
    IF (SIZE(id5))[0] NE 0 THEN C_PAH = 1e-18 * (1.8687 + 0.1905 * x) + S[1]
    IF (SIZE(id6))[0] NE 0 THEN C_PAH = 3.458 * 10.^(-17 - 3.431/x) * cutoff + TOTAL(S[2:13])
    IF (SIZE(id7))[0] NE 0 THEN C_PAH = TOTAL(S[2:13])
ENDELSE

;Reference cross sections
C_ref = xi_PAH * C_PAH + (1.-xi_PAH) * C_gra

;End of procedure
END


