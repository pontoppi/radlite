PRO xray_abundance,abun=h2o_abundance,tgas=tgas, ion_param = ion_param_dum
@natconst.pro
@line_params.ini

PRINT, 'Calculating enhanced gas temperatures using Rowins model'

;Definitions
;===================================================
; Integrated X-ray luminosity of the star 
Lx = 6.24150974d37 ; keV / s = 10^29 erg / s

; The energy needed to create one ion pair is 37 eV
ion_pair = 1e3 / 37.
Rx = 0.0571440

; Coefficients for the ionization cross sections from Morrison and
; McGammon 1983
c0 = [17.3, 34.6, 78.1, 71.4, 95.5, 308.9, 120.6, 141.3, 202.7, 342.7, 352.2, 433.9, 629.0, 701.2]
c1 = [608.1, 267.9, 18.8, 66.8, 145.8, -380.6, 169.3, 146.8, 104.7, 18.7, 18.7, -2.4, 30.9, 25.2]
c2 = [-2150., -476.1, 4.3, -51.4, -61.1, 294.0, -47.7, -31.5, -17.0, 0.0, 0.0, 0.75, 0.0, 0.0]

; Energy the intervals for the above mentioned coefficients
energies = [0.030, 0.100, 0.284, 0.400, 0.532, 0.707, 0.867, 1.303, 1.840, 2.471, 3.210, 4.038, $
            7.111, 8.331, 10.00]
;=====================================================

;
;Read model
ddens = read_dustdens()
dtemp = read_dusttemp()

; Energy grid between 0.03 and 10 keV
Emax = 10d0
Emin = 0.03
N_E  = 300.
EE   = (Emax-Emin)*findgen(N_E)/(N_E-1) + Emin

; Calculation of the cross sections
sigmaE = findgen(N_E)

FOR i=0,N_E-1 DO BEGIN
   k = 0
   FOR j=0,13 DO BEGIN
      IF ((EE[i] GE energies[j]) AND (EE[i] LE energies[j+1])) THEN BEGIN
         k = j
      ENDIF
   ENDFOR
   sigmaE[i] = EE[i]^(-3d0) * (c0[k] + c1[k] * EE[i] + c2[k] * EE[i]^2 ) * 1e-24
ENDFOR

restore, file=(freezefile)
ndens  = N_ELEMENTS(dens)
ntemps = N_ELEMENTS(temps)

i_logdens = alog10(dens)
j_temp    = temps

density  = ddens.rho / (mp * mu) * gtd
dusttemp = dtemp.temp(*,*,0,0)

column_density = DBLARR(ddens.nr,ddens.ntheta)

FOR th_i=0,ddens.ntheta-1 DO BEGIN
   FOR r_i=1, ddens.nr-1 DO BEGIN
      column_density[r_i,th_i] = int_tabulated(alog(ddens.rr[0:r_i,th_i]), $
                                           ddens.rr[0:r_i,th_i]*density[0:r_i,th_i], /double)
   ENDFOR
ENDFOR

; We take a thermal input spectrum
FE = exp(-EE / 1.0)

; Take an arbitrary cut off energy, which is due to absorption near
; star
sel = where(EE lt 0.1)

FE[sel] = 0d0

; Calculate the ionization rates from the ionization cross sections
; and the attenuating column densities calculated above.

; First normalize the integral

norm   = int_tabulated(EE, FE)
zeta_0 = DBLARR(ddens.nr, ddens.ntheta)

; And the flux
flux   = Lx / (4d0 * !pi * ddens.rr^2d0) ; 
zeta_0 = DBLARR(ddens.nr, ddens.ntheta)

FOR i=0,ddens.nr-1 DO BEGIN
   FOR j=0, ddens.ntheta-1 DO BEGIN
      NH = column_density[i,j]
      zeta_0[i,j] = INT_TABULATED(EE, sigmaE * FE * exp(-sigmaE * NH),/double,/sort) / $
                    norm * Flux[i,j] * ion_pair
   ENDFOR
ENDFOR

;
;Read the scale-free phenomenological vertical structure
fiducial_data = dblarr(3,100)
openr,1,'decoupled_gas_temp_H2O_abun.inp'
readf,1,fiducial_data
close,1

ion_param_dum  = fiducial_data[0,*]
h2o_param_dum  = fiducial_data[1,*]
temp_param_dum = fiducial_data[2,*]

Nion = 1000
maxion = -14d0
minion = -33d0
ion_param      = 10d0^((maxion-minion)*findgen(Nion) / (Nion-1) + minion)
temp_fiducial  = 10d0^INTERPOL(alog10(temp_param_dum), alog10(ion_param_dum), alog10(ion_param))
h2o_fiducial   = 10d0^INTERPOL(alog10(h2o_param_dum), alog10(ion_param_dum), alog10(ion_param))

gastemp = dusttemp
h2o_abundance = dblarr(ddens.nr, ddens.ntheta)
ion_param_dum = dblarr(ddens.nr, ddens.ntheta)

FOR i=0,ddens.nr -1 DO BEGIN
   FOR j=0,ddens.ntheta-1 DO BEGIN
      ion_param_dum[i,j] = zeta_0[i,j] / density[i,j]

      IF (ion_param_dum[i,j] GT ion_param[Nion-1]) THEN ion_param_dum[i,j] = ion_param[Nion-1]
      IF (ion_param_dum[i,j] LT ion_param[0]) THEN ion_param_dum[i,j] = ion_param[0]

      gastemp[i,j]       = gastemp[i,j] + 10d0^(INTERPOL(alog10(temp_fiducial), alog10(ion_param), alog10(ion_param_dum[i,j]), /spline))
      h2o_abundance[i,j] = 10d0^(INTERPOL(alog10(H2O_fiducial), alog10(ion_param), alog10(ion_param_dum[i,j]), /spline))

      IF (ion_param_dum[i,j] LT 1e-19) THEN h2o_abundance[i,j] = 1.0
   ENDFOR
ENDFOR

a_h2o_freeze = fltarr(ddens.nr,ddens.ntheta)
FOR i=0, ddens.nr - 1 DO BEGIN
   FOR j=0, ddens.ntheta - 1 DO BEGIN

      x=alog10(density[i,j])
      y=dusttemp[i,j]
      xnorm = (x - min(i_logdens)) / (max(i_logdens) - min(i_logdens)) * (n_elements(i_logdens) - 1.)
      ynorm = (y - min(j_temp)) / (max(j_temp) - min(j_temp)) * (n_elements(j_temp) - 1.)
      IF (y LT 273.) THEN BEGIN
         IF (xnorm lt 0.)  THEN xnorm = 0.
         IF (xnorm gt ndens-1) THEN xnorm = ndens-1
         IF (ynorm lt 0.)  THEN ynorm = 0.
         IF (ynorm gt ntemps-1) THEN ynorm = ntemps-1
         x_h2o_freeze = interpolate(abundance, xnorm, ynorm, cubic=-0.5)
      ENDIF ELSE BEGIN
         x_h2o_freeze = 0.0
      ENDELSE
      
      a_h2o_freeze[i,j]  = x_h2o_freeze
      h2o_abundance[i,j] = max_abun * (1.0 - x_h2o_freeze) * h2o_abundance[i,j]
      IF (h2o_abundance[i,j] LT h2o_fiducial[Nion-1]) THEN h2o_abundance[i,j] = max_abun * h2o_fiducial[Nion-1]
      IF (h2o_abundance[i,j] LT min_abun AND gastemp[i,j] LT 300d0) THEN h2o_abundance[i,j] = min_abun

   ENDFOR
ENDFOR

IF KEYWORD_SET(coldfinger) THEN BEGIN
   PRINT, 'Vertical cold finger effect activated!'
   mp_snowline = MAX(WHERE(a_h2o_freeze[*,ddens.ntheta/2] EQ 0))+1
   IF mp_snowline[0] EQ -1 THEN BEGIN
      PRINT, 'No midplane snow line found!'
      stop
   ENDIF
   mp_snowline_radius = ddens.r[mp_snowline]
   ssubs = WHERE(ddens.rr gt mp_snowline_radius)
   h2o_abundance[ssubs] = min_abun
ENDIF

tgas          = gastemp

;
;Write the gas temperature.       
openw,lun,'temperature.inp',/get_lun
printf,lun, ddens.nr,ddens.ntheta/2,1

FOR ir=0,ddens.nr-1 DO BEGIN
   FOR it=0,ddens.ntheta/2 -1 DO BEGIN
      printf, lun, tgas[ir,it]
   ENDFOR
ENDFOR
close,lun
FREE_LUN, lun

;
;Write the abundance.       

openw,lun,'abundance.inp',/get_lun
printf,lun,ddens.nr,ddens.ntheta/2
for ir=0,ddens.nr-1 do begin
    for it=0,ddens.ntheta/2-1 do begin
        printf,lun,h2o_abundance[ir,it], 0d0
    endfor
endfor
close,lun
FREE_LUN, lun



END
