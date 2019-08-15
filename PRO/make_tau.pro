;=======================================================================
;Calculate the optical depth along the line of sight to the central star
;=======================================================================

PRO make_tau, ref_freq, tau=tau, vtau=vtau, column=column
@natconst.pro

a  = read_dustdens()
s  = read_spectrum(file='starspectrum.inp')
ns = n_elements(a.rho[0,0,*]) 
nf = n_elements(s.freq)-1
nr = n_elements(a.r)
nt = n_elements(a.theta)
kappa = dblarr(nf,ns)
kappa_ref = dblarr(ns)
for is=0,ns-1 do begin
    file = strcompress(string(is+1,format='(I4)'),/remove_all)
    o = readopac(nr=file)
    kappa[*,is] = o.cabs[*] + o.csca[*]
    kappa_ref[is] = INTERPOL(kappa[*,is],s.freq[1:*],ref_freq)
endfor


;
; Computing the column from the star along the line of sight
;
column = dblarr(nr,nt,ns)
FOR it=0,nt-1 DO BEGIN
    FOR ir=1,nr-1 DO BEGIN
        column[ir,it,*] = column[ir-1,it,*] + $
          0.5 * ( a.rho[ir,it,*] + a.rho[ir-1,it,*] ) *$
          ( a.r[ir] - a.r[ir-1] )
    ENDFOR
ENDFOR

tau = dblarr(nr,nt)
FOR is=0,ns-1 DO BEGIN
    tau = tau + column[*,*,is]*kappa_ref[is]
ENDFOR

;
; Computing the vertical tau (really along constant r, not exactly right)
;

column = dblarr(nr,nt,ns)
FOR ir=0,nr-1 DO BEGIN
    FOR it=1,nt/2-1 DO BEGIN
        column[ir,it,*] = column[ir,it-1,*] + $
          0.5 * ( a.rho[ir,it,*] + a.rho[ir,it-1,*] ) *$
          ( cos(a.theta[it-1]) - cos(a.theta[it]) ) * a.r[ir]
        ;
        ;Mirror
        column[ir,nt-it-1] = column[ir,it,*]
    ENDFOR
ENDFOR


vtau = dblarr(nr,nt)
FOR is=0,ns-1 DO BEGIN
    vtau = vtau + column[*,*,is]*kappa_ref[is]
ENDFOR

END
