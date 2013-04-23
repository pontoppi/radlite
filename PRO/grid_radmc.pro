;=======================================================
;Example RADMC grid script. In this case, we vary the disk mass and gas to dust ratio.
;
;KEYWORDS:
;clobber - overwrite previous run
;run_dirs - returns a named array of the run directories
;=======================================================

PRO grid_radmc, clobber=clobber, run_dirs=run_dirs

;IF KEYWORD_SET(clobber) THEN spawn, 'rm -rf grid_*'

mdisk=[1d-1,1d-2,1d-3,1d-4]
g2d=[1d2,1d3,1d4,1d5]
hrgrid=[.5]

nx = N_ELEMENTS(mdisk)
ny = N_ELEMENTS(g2d)
run_dirs = STRARR(nx,ny)

FOR i=0,N_ELEMENTS(mdisk)-1 DO BEGIN
   FOR j=0,N_ELEMENTS(g2d)-1 DO BEGIN
      date = systime(0)
      run_dir='RUN_'+strmid(date,4,3)+'_'+STRTRIM(strmid(date,8,2),2)+'_'+strmid(date,11,8)
      run_dir_save=strcompress('DIRS_'+STRING(nx+ny)+run_dir+'.sav',/REMOVE_ALL)
      spawn, 'mkdir '+run_dir
      make_problem_params, var_mdisk=mdisk[i],var_gtd=g2d[j],var_hrgrid=hrgrid
      spawn, 'cp problem_params.pro '+run_dir
      spawn, 'cp silicate.Kappa '+run_dir
      spawn, 'cp carbon.Kappa '+run_dir
      spawn, 'cp forsterite.Kappa '+run_dir
     ; .r problem_setup.pro
     ; spawn,'radmc'
      run_dirs[i,j] = run_dir
      wait,1.5
   ENDFOR
ENDFOR
save, mdisk,g2d,run_dirs,filename=run_dir_save

END
