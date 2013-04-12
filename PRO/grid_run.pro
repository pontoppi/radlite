;=======================================================
;Example RADLite grid script. In this case, we vary 
;the freeze-out temperature in 6 steps and the gas-to-dust ratio,
;while keeping the dust mass constant. 
;
;KEYWORDS:
;clobber - overwrite previous run
;run_dirs - returns a named array of the run directories
;=======================================================

PRO grid_run, clobber=clobber, run_dirs=run_dirs

IF KEYWORD_SET(clobber) THEN spawn, 'rm -rf grid_*'

fr_temps = [1.,50.,100.,150.,200.,250.]
g2d      = [100.,1000.,10000.]

nx = N_ELEMENTS(fr_temps)
ny = N_ELEMENTS(g2d)
run_dirs = STRARR(nx,ny)

FOR i=0,N_ELEMENTS(fr_temps)-1 DO BEGIN
   FOR j=0,N_ELEMENTS(g2d)-1 DO BEGIN
      make_line_params, var_fr_temp=fr_temps[i],var_max_abun=3.6e-8*g2d[j]/100.,var_min_abun=1e-12*g2d[j]/100.
      RESOLVE_ROUTINE, ['line_run','problem_lines','lamda_extract_lines','make_abundance',$
                        'make_turbulence','make_velocity','nlte_main','xray_abundance']
      line_run, run_name='grid',wait_time=60.,rundir=rundir
      run_dirs[i,j] = rundir
   ENDFOR
ENDFOR

END
