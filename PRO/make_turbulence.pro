;
;Procedure for setting up a turbulence grid.
;This is where additional turbulence mechanisms can be constructed.
;
;Currently included:
;
;1) Alpha-viscosity turbulence
;

PRO make_turbulence,a,ttype,alpha=alpha,tur=tur,kepler_frac=kepler_frac,molmass=molmass
@natconst.pro
@line_params.ini

iformat=1

gamma  = 1.4 ;Adiabatic constant for a diatomic gas
Mu     = 2.3 ;Mean molecular mass

IF NOT KEYWORD_SET(ttype) THEN BEGIN
    PRINT, 'No turbulence structure selected!' 
    stop
ENDIF

nr = n_elements(a.r)
nt = n_elements(a.theta)/2

Tstruct = read_temperature()
T = Tstruct.t

IF N_ELEMENTS(VERBOSE) THEN print, 'Using first dust component temperature to determine the turbulent velocities'

CASE ttype OF
    1: BEGIN                    ;Alpha viscosity
        IF NOT KEYWORD_SET(alpha) THEN BEGIN 
            alpha = 0.01
            IF N_ELEMENTS(VERBOSE) THEN print, 'No alpha set - assuming alpha = 0.01'
        ENDIF
        ;Calculating sound speed:
        cs = SQRT(gamma*kk*T/(Mu*mp))
        tur = alpha*cs
    END
    2: BEGIN                    ;fraction of kepler velocity
        IF NOT KEYWORD_SET(kepler_frac) THEN BEGIN 
            kepler_frac = 0.01
            IF N_ELEMENTS(VERBOSE) THEN print, 'No Kepler velocity fraction set, assuming 0.01'
        ENDIF
        ;Calculating sound speed:Kepler velocity:
        openr,1,'starinfo.inp'
        iformat=0
        readf,1,iformat
        rstar=0.d0
        mstar=0.d0
        readf,1,rstar
        readf,1,mstar
        close,1
        print, 'Used starinfo.inp for Keplerian velocity'
        vphi = sqrt(GG*mstar/a.rr)
        tur = kepler_frac*vphi
    END
ENDCASE

;
;Add thermal broadening (this sets a natural lower bound to the local
;line broadening):

tur = sqrt(tur^2.+2d0*kk*T/(molmass[0]*mp))

openw,1,'turbulence.inp'
printf,1,iformat
printf,1,nr,nt
for ir=0,nr-1 do begin
    for it=0,nt-1 do begin
        printf,1,tur[ir,it]/1d5  ;convert cm/s --> km/s
    endfor
endfor
close,1


END

PRO test_velocity


END

PRO read_tur, file, vel
;
;Reads a turbulence grid
;
openr,1,file
readf,1,nr,nt

tur = dblarr(nr,nt)

for ir=0,nr-1 do begin
    for it=0,nt-1 do begin
        readf,1,dum1
        tur[ir,it]     = dum1
    endfor
endfor
close,1


END
