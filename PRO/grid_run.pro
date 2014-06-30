;=======================================================
;Example RADLite grid script. In this case, we vary 
;the freeze-out temperature in 6 steps and the gas-to-dust ratio,
;while keeping the dust mass constant. 
;
;KEYWORDS:
;clobber - overwrite previous run
;run_dirs - returns a named array of the run directories
;=======================================================

PRO grid_run, clobber=clobber, run_table=run_table,obsres=obsres

IF KEYWORD_SET(clobber) THEN spawn, 'rm -rf grid_*'
IF ~KEYWORD_SET(obsres) THEN obsres=300
IF ~KEYWORD_SET(run_table) THEN run_table='run_table.fits'

Ng2d = 6
g2d_max = alog10(10000.)
g2d_min = alog10(10.)
fr_temps = [1.,50.,100.,150.,200.,250.]
log_g2d      = (g2d_max-g2d_min)*FINDGEN(Ng2d)/(Ng2d-1)+g2d_min
g2d = 10^log_g2d

nx = N_ELEMENTS(fr_temps)
ny = N_ELEMENTS(g2d)
run_pars = []

MWRFITS, dummy, run_table, /CREATE

FOR i=0,N_ELEMENTS(fr_temps)-1 DO BEGIN
   FOR j=0,N_ELEMENTS(g2d)-1 DO BEGIN
      make_line_params, var_fr_temp=fr_temps[i],var_max_abun=3.6e-8*g2d[j]/100.,var_min_abun=1e-12*g2d[j]/100.
      ;We need to recompile all routines that inline line_params.ini
      RESOLVE_ROUTINE, ['line_run','problem_lines','lamda_extract_lines','make_abundance',$
                        'make_turbulence','make_velocity','nlte_main','xray_abundance']
      line_run, run_name='grid',rundir=rundir
      run_par = {dir:rundir,fr_temp:fr_temps[i],g2d:g2d[j]}
      run_pars = [run_pars,run_par]
      cd, rundir
      genspec, obsres=obsres
      cd, '..'
   ENDFOR
ENDFOR

MWRFITS,run_pars, run_table

END
