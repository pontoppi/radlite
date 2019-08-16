;
;Procedure for setting up a turbulence grid.
;This is where additional turbulence mechanisms can be constructed.
;
;Currently included:
;
;1) Alpha-viscosity turbulence
;

PRO make_abundance,atype, PT_rel=PT_rel,abun=abun,tgas=tgas,mol_destruct=mol_destruct,speciesfrac=speciesfrac

@natconst.pro
@line_params.ini

iformat=1

gamma  = 1.4 ;Adiabatic constant for a diatomic gas

Tstruct = read_temperature()
Dstruct = read_density()
T       = Tstruct.t
moddens = Dstruct.rho
nr      = n_elements(Dstruct.r)
nt      = n_elements(Dstruct.theta)


CASE atype OF
    1: BEGIN                    ;Constant abundance
        IF N_ELEMENTS(VERBOSE) THEN PRINT, 'Setting constant abundance'
        abun = fltarr(nr,nt) + max_abun
        abun_collpartner = 0.d0
    END        
    2: BEGIN                    ;Freeze-out at some temperature
        IF N_ELEMENTS(VERBOSE) THEN PRINT, 'Setting abundance to ~0 below T='+STRTRIM(STRING(fr_temp),2)+' K'
        abun = fltarr(nr,nt)
        abun[*] = max_abun
        fsubs = WHERE(T lt fr_temp)
        IF fsubs[0] NE -1 THEN abun[fsubs] = min_abun
        abun_collpartner = 0.d0
    END
    3: BEGIN                    ;Step abundance at AV=3 (as seen from the star).
        IF N_ELEMENTS(VERBOSE) THEN PRINT, 'Setting abundance to 0 below AV=1...'
        abun = fltarr(nr,nt)
        make_tau, cc*1d4/0.12, tau=tau
        abun = (1-exp(-tau)) * max_abun * (1-(Dstruct.rr/Dstruct.r[0])^(-2))
        abun_collpartner = 0.d0
    END
    6: BEGIN                    ;Step abundance at col=5x10^18 cm-1
        IF N_ELEMENTS(VERBOSE) THEN PRINT, 'Setting abundance to 0 above a vertical H2 column of 5x10^18 cm-2'
        abun = fltarr(nr,nt)+max_abun
        make_tau, cc*1d4/0.12, tau=tau,column=column
		low_subs = WHERE(column*gtd/(2.4*1.66e-24) LT 5e18)
        abun[low_subs] = min_abun
        abun_collpartner = 0.d0
    END
    5: BEGIN                    ;infall model
        abun = fltarr(nr,nt)
        abun[*] = max_abun
        fsubs = WHERE(T lt fr_temp)
        IF fsubs[0] NE -1 THEN abun[fsubs] = 10.0^(alog10(Dstruct.rr[fsubs]/min(Dstruct.rr[fsubs])))*min_abun
        abun_collpartner = 0.d0
    END
    4: BEGIN                    ;Density-dependent freeze-out temperature
       IF N_ELEMENTS(VERBOSE) THEN PRINT, 'Using density-dependent freeze-out'
       IF NOT KEYWORD_SET(PT_rel) THEN BEGIN
          PT_rel_file = 'h2oabun_g2d_1280.sav'
       ENDIF ELSE BEGIN
          PT_rel_file = PT_rel
       ENDELSE
       IF N_ELEMENTS(VERBOSE) THEN BEGIN
          PRINT, 'Setting up density dependent abundance'
          PRINT, 'Using: ', PT_rel_file
       ENDIF

       restore, PT_rel_file
       logdens = alog10(dens)
       nd = N_ELEMENTS(temps)
       lsubs = WHERE(moddens EQ 0)
       IF lsubs[0] NE -1 THEN moddens[lsubs] = 1d-50
       h2dens  = alog10(moddens/(mu*mp))
       densnorm = (h2dens - min(logdens)) / (max(logdens) - min(logdens)) * (n_elements(logdens) - 1.)
       Tnorm    = (T - min(temps)) / (max(temps) - min(temps)) * (n_elements(temps) - 1.)

       subs     = WHERE(densnorm LT 0.)
       IF subs[0] NE -1 THEN densnorm[subs] = 0
       subs     = WHERE(densnorm GT nd)
       IF subs[0] NE -1 THEN densnorm[subs] = nd-1
       subs     = WHERE(Tnorm LT 0)
       IF subs[0] NE -1 THEN Tnorm[subs] = 0       
       subs     = WHERE(Tnorm GT nd)
       IF subs[0] NE -1 THEN Tnorm[subs] = nd-1
       subs     = WHERE(T GT 273)
       
       a_h2o_freeze = interpolate(abundance, densnorm, Tnorm,cubic=-0.5)
       a_h2o_freeze[subs] = 0d0

       abun = max_abun * (1.0 - a_h2o_freeze)
       subs = WHERE(abun LT min_abun)
       IF subs[0] NE -1 THEN abun[subs] = min_abun
       abun_collpartner = 0.d0
       IF KEYWORD_SET(coldfinger) THEN BEGIN
          IF N_ELEMENTS(VERBOSE) THEN PRINT, 'Vertical cold finger effect activated!'

          mp_snowline = MAX(WHERE(a_h2o_freeze[*,nt/2] EQ 0))+1
          IF mp_snowline[0] EQ -1 THEN BEGIN
             PRINT, 'No midplane snow line found!'
             stop
          ENDIF
          IF coldfinger EQ -1 THEN BEGIN
             mp_snowline_radius = dstruct.r[mp_snowline]
             ssubs = WHERE(dstruct.rr gt mp_snowline_radius)
          ENDIF ELSE BEGIN
             ssubs = WHERE(dstruct.rr gt coldfinger*AU)
          ENDELSE
          abun[ssubs] = min_abun
       ENDIF


    END
    7: BEGIN           ;Rowin's Xray temperature + destruction. In this case, the H2O abundance is calculated coupled with the gas temperature.
       xray_abundance,abun=abun,tgas=tgas, ion_param=ion_param
       abun_collpartner = 0d0
    END
ENDCASE

IF KEYWORD_SET(mol_destruct) THEN BEGIN
	dsubs = WHERE(mol_destruct EQ 1)
	abun[dsubs] = 1e-30
ENDIF

IF KEYWORD_SET(speciesfrac) THEN BEGIN
	abun = abun*speciesfrac
ENDIF

openw,1,'abundance.inp'
printf,1,nr,nt/2
for ir=0,nr-1 do begin
    for it=0,nt/2-1 do begin
        printf,1,abun[ir,it], abun_collpartner
    endfor
endfor
close,1


END

PRO read_abun, file, abun
;
;Reads a abundance grid
;
openr,1,file
readf,1,nr,nt

abun = dblarr(nr,nt)

for ir=0,nr-1 do begin
    for it=0,nt-1 do begin
        readf,1,dum1
        abun[ir,it]     = dum1
    endfor
endfor
close,1


END
