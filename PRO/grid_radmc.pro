@make_line_params
@make_problem_params
;=======================================================
;Example RADMC grid script. In this case, we vary the disk mass and gas to dust ratio.
;
;KEYWORDS:
;clobber - overwrite previous run
;run_dirs - returns a named array of the run directories
;=======================================================

PRO grid_radmc, clobber=clobber, save_radmc=save_radmc, run_table=run_table, obsres=obsres

IF KEYWORD_SET(clobber) THEN spawn, 'rm -rf grid_*'
IF ~KEYWORD_SET(obsres) THEN obsres=300
IF ~KEYWORD_SET(run_table) THEN run_table='run_table.fits'

;fr_temps = [1.,50.,100.,150.,200.,250.]

mdisk=[1d-3];[1d-1,1d-2,1d-3,1d-4]
g2d=[1d3];[1d2,1d3,1d4]
hrgrid=[.4]
hrpl=[1./7.];[0.,1./7.,2./7.]
mstar=[1.0];[0.1,0.5,1.0]
rstar=[2.0]
tstar=[4275.0]
isot=[51]
maxabun=[1d-4]
minabun=[1d-4]

nx = N_ELEMENTS(mdisk)
ny = N_ELEMENTS(g2d)
nz = N_ELEMENTS(hrpl)
n4 = N_ELEMENTS(mstar)

MWRFITS, dummy, run_table, /CREATE
count=0

FOR i=0,N_ELEMENTS(mdisk)-1 DO BEGIN
   FOR j=0,N_ELEMENTS(g2d)-1 DO BEGIN
      FOR k=0,N_ELEMENTS(hrpl)-1 DO BEGIN
         FOR l=0,N_ELEMENTS(mstar)-1 DO BEGIN
            print,i,j,k,l
            make_problem_params, var_mdisk=mdisk[i],var_gtd=g2d[j],$
                                 var_plh=hrpl[k],var_mstar=mstar[l],$
                                 var_tstar=tstar[l],var_rstar=rstar[l],var_hrgrid=hrgrid
            problem_setup
            wait,2.0
            spawn,'radmc'
            FOR a=0,N_ELEMENTS(isot)-1 DO BEGIN
               count=count+1
               it_is_there = FILE_TEST('meanint_radmc.dat')
               IF it_is_there THEN BEGIN  
                                ;Only start RADLITE if RADMC was successful
                  make_line_params, var_gtd=g2d[j],var_isot=isot[a],var_max_abun=maxabun[a], var_min_abun=minabun[a]
                                ;We need to recompile all routines that inline line_params.ini
                  RESOLVE_ROUTINE, ['line_run','problem_lines','lamda_extract_lines','make_abundance',$
                                    'make_turbulence','make_velocity','nlte_main','xray_abundance']
                  line_run, run_name='grid',wait_time=60.,rundir=rundir
                  run_par = {dir:rundir,mdisk:mdisk[i],g2d:g2d[j],plh:hrpl[k],mstar:mstar[l]}
                  IF count EQ 1 THEN run_pars=run_par ELSE run_pars = [run_pars,run_par]
                  cd, rundir
                  genspec, obsres=obsres
                  cd, '..'
               ENDIF
            ENDFOR
            IF KEYWORD_SET(save_radmc) THEN BEGIN
               spawn, 'cp problem_params.pro '+run_dir
               spawn, 'cp meanint_radmc.dat '+run_dir
               spawn, 'cp scatsource.dat '+run_dir
               spawn, 'cp dusttemp_final.dat '+run_dir
               spawn, 'cp spectrum_all.dat '+run_dir
            ENDIF

         ENDFOR
      ENDFOR
   ENDFOR
ENDFOR

MWRFITS,run_pars, run_table
END
