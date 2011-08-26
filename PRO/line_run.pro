@analyze.pro
@problem_lines
@hitran_extract.pro
@read_posvel.pro
@problem_files.pro
@problem_subroutines.pro
@problem_makeshu.pro
;@read_line.pro
@make_velocity.pro
@make_turbulence.pro
@make_tau.pro
@make_abundance.pro
@readopac.pro
@write_linespectrum_inp.pro
@Ppump.pro
@Pesc.pro
@P.pro
@nlte.pro
@nlte_main.pro

PRO line_run, run_name=run_name
@natconst.pro
@line_params.ini


IF image EQ 2 THEN BEGIN
   executable = 'time '+exe_path+'RADlite_imcir'
ENDIF ELSE BEGIN
   executable = 'time '+exe_path+'RADlite'
ENDELSE

IF NOT KEYWORD_SET(run_name) THEN run_name='run'

;First make a test run with hitran_extract to determine the total
;number of lines in the requested wavelength range. 
hitran_extract,cutoff=cutoff,lambdarange=[min_mu,max_mu],freq=freq,$
               isot=isot,molfile=molfile,H2O_OP=H2O_OP,$
               max_energy=max_energy, hitran_path=hit_path
;
;Sort in wavelength
wave = 1d4/freq
ssubs = SORT(wave)
wave  = wave[ssubs]
   
nlines = N_ELEMENTS(wave)
subN   = CEIL(nlines/FLOAT(ncores))
;
;Setting up the results directory
time = systime(0)
rundir = run_name+'_'+strmid(time,4,3)+'_'+STRTRIM(strmid(time,8,2),2)+'_'+strmid(systime(),11,5)
spawn, 'mkdir '+ rundir

PRINT, 'Splitting up in ',ncores, ' sub-runs'

FOR iii=0,ncores-1 DO BEGIN

    ;
    ;Some times the background process
    ;output files are not closed
    ;properly. I don't know why,
    ;but this check is intended to
    ;forcefully clear them.
    molfile = 'moldata_'+STRTRIM(STRING(iii),2)+'.dat'
    ;
    ;Find the boundaries of the subrun
    min_mu_run = min(wave[iii*subN:MIN([(iii+1)*subN-1,nlines-1])])
    max_mu_run = max(wave[iii*subN:MIN([(iii+1)*subN-1,nlines-1])])
    ;
    ;shift a little bit to make sure the
    ;lines on the edges are included.
    ;This may cause each run to not
    ;exactly have subN lines, but will
    ;make sure every line in interval is
    ;eventually included
    min_mu_run = min_mu_run*(1d0 - 1d-12) 
    max_mu_run = max_mu_run*(1d0 + 1d-12) 
    ;
    ;Now calculate moldata.dat for the subrun
    hitran_extract,cutoff=cutoff,lambdarange=[min_mu_run,max_mu_run],$
                   freq=freq,isot=isot,molfile=molfile,H2O_OP=H2O_OP,$
                   max_energy=max_energy, hitran_path=hit_path
    ;
    ;Setup the linespectrum.inp file
    write_linespectrum_inp,nlines=N_ELEMENTS(freq),incl=incl,dist=dist,image=image,$
                           vsampling=vsampling,ssampling=ssampling,passband=passband,$
                           imwidth=imwidth,molfile=molfile
    ;
    ;Run the setup script for the new moldata.dat
    problem_lines, molfile

    IF iii LT ncores-1 THEN BEGIN
       print, 'spawning background process for core: ', iii+1
       spawn, executable+'&'
       
    ENDIF ELSE BEGIN
       ;Were processes run in the background
       ;on some of the cores? if so, free the processes.
       print, 'spawning foreground process for core: ',iii+1
       spawn, executable,exit_status=exit3
    ENDELSE
    wait, 5.
ENDFOR
;
;Now save the result
;
;Copy the model setup parameters

IF ncores GT 1 THEN BEGIN
   PRINT, 'Make sure all processes are complete, then type .c to finish'
   stop
ENDIF

spawn, 'cp problem_params.pro '+rundir+'/.'
spawn, 'cp line_params.ini '+rundir+'/.'
FOR iii=0,ncores-1 DO BEGIN
   ;
   ;Save the molecular file to a unique name
   spawn, 'mv moldata_'+STRTRIM(STRING(iii),2)+'.dat '+rundir+'/.'
   spawn, 'mv levelpop_moldata_'+STRTRIM(STRING(iii),2)+'.dat '+rundir+'/.'
   ;
   ;And save the lines to a unique name
   IF image eq 0 or image EQ 2 THEN BEGIN
      spawn, 'mv linespectrum_moldata_'+STRTRIM(STRING(iii),2)+'.dat '+rundir+'/.'
   ENDIF
   IF image EQ 2 THEN BEGIN
      spawn, 'mv lineposvelcirc_moldata*.dat '+rundir+'/.'
   ENDIF
   IF image EQ 1 THEN BEGIN
      FOR lj=1,nlines DO BEGIN
         spawn, 'mv lineposvel_moldata_'+STRTRIM(STRING(iii),2)+'_'+STRTRIM(STRING(lj),2)+'.dat '+rundir+$
                '/lineposvel_moldata_'+STRTRIM(STRING(iii+lj),2)+'.dat'           
      ENDFOR
      IF KEYWORD_SET(dat2fits) THEN BEGIN
         cd, rundir
         FOR lj=1,nlines DO BEGIN
            read_posvel, 'lineposvel_moldata_'+STRTRIM(STRING(iii+lj),2)+'.dat', molecule=molecule
         ENDFOR
         cd, '..'
      ENDIF
   ENDIF
ENDFOR

END
