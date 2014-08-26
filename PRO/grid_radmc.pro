@make_line_params
@make_problem_params
;=======================================================
;Example RADMC grid script. In this case, we vary the disk mass and gas to dust ratio.
;
;KEYWORDS:
;clobber - overwrite previous run
;run_dirs - returns a named array of the run directories
;=======================================================

PRO grid_radmc, clobber=clobber, save_radmc=save_radmc, run_table=run_table, obsres=obsres, start_index=start_index

IF KEYWORD_SET(clobber) THEN spawn, 'rm -rf grid_*'
IF ~KEYWORD_SET(obsres) THEN obsres=3.
IF ~KEYWORD_SET(run_table) THEN run_table='run_table.fits'
IF ~KEYWORD_SET(start_index) THEN start_index=[2,1,1,0]

;fr_temps = [1.,50.,100.,150.,200.,250.]

mdisk   = [1d-1,1d-2,1d-3,1d-4]
g2d     = [1d2,1d3,1d4]
hrgrid  = [.5,.43,.38,.32]
nphot   = [200000,200000,200000,200000]
hrpl    = [0,1./7.,2./7.]
mstar   = [0.5,1.0,2.0,3.5]
rstar   = [1.5,2.0,2.85,2.7]
tstar   = [3765.0,4275.0,4820.0,12045.0]

isot    = [51,52]
maxabun = [1d-4,1d-4/70.]
minabun = [1d-6,1d-6/70.]

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
			IF i GE start_index[0] AND j GE start_index[1] AND k GE start_index[2] AND l GE start_index[3] THEN BEGIN
				make_problem_params, var_mdisk=mdisk[i],var_gtd=g2d[j],$
                					 var_plh=hrpl[k],var_mstar=mstar[l],$
                                 	 var_tstar=tstar[l],var_rstar=rstar[l],$
                                     var_hrgrid=hrgrid[i],var_nphot=nphot[i]
            	RESOLVE_ROUTINE, ['problem_setup']
            	problem_setup
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
                      line_run, run_name='grid',rundir=rundir
                      run_par = {dir:rundir,mdisk:mdisk[i],g2d:g2d[j],plh:hrpl[k],mstar:mstar[l]}
                      IF count EQ 1 THEN run_pars=run_par ELSE run_pars = [run_pars,run_par]
                          cd, rundir
                          genspec, obsres=obsres
                          cd, '..'
                      ENDIF
                 ENDFOR
			
                spawn, 'cp problem_params.pro '+rundir
                IF KEYWORD_SET(save_radmc) THEN BEGIN
                   spawn, 'cp meanint_radmc.dat '+rundir
                   spawn, 'cp scatsource.dat '+rundir
                   spawn, 'cp dusttemp_final.dat '+rundir
                   spawn, 'cp spectrum_all.dat '+rundir
                ENDIF
			ENDIF
         ENDFOR
      ENDFOR
   ENDFOR
ENDFOR

MWRFITS,run_pars, run_table
END
