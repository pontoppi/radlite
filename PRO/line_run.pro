@analyze.pro
@problem_lines
@hitran_extract.pro
@read_posvel.pro
@problem_files.pro
@problem_subroutines.pro
@problem_makeshu.pro
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
@interpol
@read_molecule_lambda
@read_psum
@lamda_extract_lines
@lamda_extract_levels
@make_levelpop

PRO line_run, run_name=run_name, wait_time=wait_time, look_for_lpop_file=look_for_lpop_file, rundir=rundir, save_levelpop=save_levelpop, nodate=nodate
@natconst.pro
@line_params.ini

RESOLVE_ALL, /continue_on_error, /QUIET
spawn, 'clear'
;
;Welcome message
print, '------------------------------------------'
print, 'Welcome to RADLite Version 1.2 (2007-2011)'
print, '                                          '
print, 'Written by:                               '
print, 'Klaus Pontoppidan (pontoppi@stsci.edu)    '
print, 'Kees Dullemond                            '
print, 'Alex Lockwood                             '
print, 'Rowin Meijerink                           '
print, '------------------------------------------'
print, ' '
print, ' '

IF image EQ 2 THEN BEGIN
   executable = exe_path+'RADlite_imcir'
ENDIF ELSE BEGIN
   executable = exe_path+'RADlite'
ENDELSE

IF NOT KEYWORD_SET(run_name) THEN run_name='run'

;First make a test run with hitran_extract (for LTE) or
;lamda_extract_lines (for NLTE) to determine the total
;number of lines in the requested wavelength range. 
IF lte EQ 1 THEN BEGIN
   hitran_extract,cutoff=cutoff,lambdarange=[min_mu,max_mu],freq=freq,$
                  isot=isot,molfile=molfile,H2O_OP=H2O_OP,$
                  max_energy=max_energy, hitran_path=hit_path
ENDIF ELSE BEGIN  
   molfile = 'molfile.dat' 
   lamda_extract_lines, lambdarange=[min_mu,max_mu],$
                        isot=isot,molfile=molfile,max_energy=max_energy, $
                        lamda_path=lamda_path
   moldum = READ_MOLECULE_LAMBDA(molfile)
   freq   = moldum.freq
ENDELSE
;
;Sort in wavelength
wave = 1d4/freq
ssubs = SORT(wave)
wave  = wave[ssubs]
   
nlines = N_ELEMENTS(wave)
;
;Check that there are more lines than cores
IF nlines LT ncores THEN ncores = nlines
subN   = CEIL(nlines/FLOAT(ncores))
;
;Setting up the results directory
time = systime(0)

IF KEYWORD_SET(nodate) THEN BEGIN
	rundir = run_name+'_'+strmid(time,4,3)+'_'+STRTRIM(strmid(time,8,2),2)
ENDIF ELSE BEGIN
	rundir = run_name+'_'+strmid(time,4,3)+'_'+STRTRIM(strmid(time,8,2),2)+'_'+strmid(systime(),11,5)	
ENDELSE

spawn, 'mkdir '+ rundir

PRINT, 'Rendering lines on ',ncores, ' processor cores'

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
    min_mu_run = MIN(wave[MIN([iii*subN,nlines-1]):MIN([(iii+1)*subN-1,nlines-1])])
    max_mu_run = MAX(wave[MIN([iii*subN,nlines-1]):MIN([(iii+1)*subN-1,nlines-1])])
    ;
    ;shift a little bit to make sure the
    ;lines on the edges are included.
    ;This may cause each run to not
    ;exactly have subN lines, but will
    ;make sure every line in interval is
    ;eventually included
    min_mu_run = min_mu_run*(1d0 - 1d-4) 
    max_mu_run = max_mu_run*(1d0 + 1d-4) 
    ;
    ;Extract lines for LTE or nLTE
    IF lte EQ 1 THEN BEGIN
       ;
       ;Now calculate moldata.dat for the subrun
       hitran_extract,cutoff=cutoff,lambdarange=[min_mu_run,max_mu_run],$
                      freq=freq,isot=isot,molfile=molfile,H2O_OP=H2O_OP,$
                      max_energy=max_energy, hitran_path=hit_path
       
    ENDIF ELSE BEGIN
       lamda_extract_lines, lambdarange=[min_mu_run,max_mu_run],$
                            isot=isot,molfile=molfile,max_energy=max_energy, $
                            lamda_path=lamda_path, freq=freq
    ENDELSE
    ;
    ;Setup the linespectrum.inp file
    write_linespectrum_inp,nlines=N_ELEMENTS(freq),incl=incl,dist=dist,image=image,$
                           vsampling=vsampling,ssampling=ssampling,passband=passband,$
                           imwidth=imwidth,molfile=molfile
    ;
    ;Run the setup script for the new moldata.dat

    ;
    ;Check for existing non-lte level
    ;population file, and calculate it if
    ;this is the first core run (if
    ;it's a later core, it is
    ;assumed that the nlte module was
    ;already run.  
    run_nlte = -1
    IF iii EQ 0 THEN BEGIN
       it_is_there = FILE_TEST('levelpop_nlte.fits')
       IF it_is_there AND KEYWORD_SET(look_for_lpop_file) THEN BEGIN
          PRINT, 'Existing level population file found - do you want to use it?'
          answer=' '
          WHILE run_nlte EQ -1 DO BEGIN
             read, answer, prompt='[y/n]'
             CASE answer OF
                'y':  run_nlte = 0
                'n':  run_nlte = 1 
                ELSE: print, 'Please answer yes [y] or no [n]...'
             ENDCASE
          ENDWHILE
       ENDIF ELSE BEGIN
          run_nlte = 1   ;We didn't find a nlte file - calculate it!
       ENDELSE
    ENDIF
    
    problem_lines, molfile, run_nlte=run_nlte

	; Check that the setup has completed before proceeding (IO may not complete otherwise)
	WHILE FILE_TEST('levelpop_'+molfile) NE 1 DO BEGIN
		wait, 1.0
	ENDWHILE
	
    IF iii LT ncores-1 THEN BEGIN
       print, 'spawning background process for core: ', iii+1
       spawn, executable+' > RADLite_core'+STRTRIM(STRING(iii+1),2)+'.log&'
	   ; Check that we are done launching radlite before proceeding
   	   WHILE FILE_TEST('linespectrum'+'_'+molfile) NE 1 DO BEGIN
   	      wait, 1.0
   	   ENDWHILE
	   
    ENDIF ELSE BEGIN
       print, 'spawning foreground process for core: ',iii+1
       spawn, executable+' > RADLite_core'+STRTRIM(STRING(iii+1),2)+'.log',exit_status=exit3
    ENDELSE
ENDFOR

;
;Check for RADlite processes still running
IF ncores GT 1 and ~KEYWORD_SET(wait_time) THEN BEGIN
   WHILE 1 DO BEGIN
      spawn, 'ps cax | grep RADlite', radlite_running
      print, 'Waiting for all RADLite threads to finish', N_ELEMENTS(radlite_running), ' threads left'
	  writeu,-1,string(13b)
      IF radlite_running[0] EQ '' THEN BEGIN
		  wait, 5.0
		  BREAK
	  ENDIF
      wait, 5.0
   ENDWHILE
ENDIF ELSE BEGIN
   wait, 5.0
ENDELSE

;
;Now save the result
;
;Copy the model setup parameters

file_copy, 'problem_params.pro', rundir+'/.', /overwrite
file_copy, 'line_params.ini', rundir+'/.', /overwrite
file_copy, 'radius.inp', rundir+'/.', /overwrite
file_copy, 'theta.inp', rundir+'/.', /overwrite
file_copy, 'frequency.inp', rundir+'/.', /overwrite
file_copy, 'density.inp', rundir+'/.', /overwrite
file_copy, 'dustdens.inp', rundir+'/.', /overwrite
file_copy, 'dusttemp_final.dat', rundir+'/.', /overwrite
file_copy, 'dustopac.inp', rundir+'/.', /overwrite
file_copy, 'dustopac_*.inp', rundir+'/.', /overwrite
file_copy, 'abundance.inp', rundir+'/.', /overwrite
file_copy, 'temperature.inp', rundir+'/.', /overwrite

file_copy, 'RADLite_core*.log', rundir+'/.', /overwrite
spawn, 'rm RADLite_core*.log'

IF KEYWORD_SET(save_levelpop) THEN BEGIN
   file_copy, 'levelpop_nlte.fits', rundir+'/.', /overwrite
ENDIF


;
;Save the molecular file to a unique name
file_copy, 'moldata_*.dat', rundir+'/.', /overwrite
file_copy, 'levelpop_moldata_*.dat', rundir+'/.', /overwrite
spawn, 'rm moldata_*.dat'
spawn, 'rm levelpop_moldata_*.dat'
;
;And save the lines to a unique name
IF image eq 0 THEN BEGIN
   file_copy, 'linespectrum_moldata_*.dat', rundir+'/.', /overwrite
   spawn, 'rm linespectrum_moldata_*.dat'
ENDIF
IF image EQ 2 THEN BEGIN
   file_copy, 'lineposvelcirc_moldata_*.dat', rundir+'/.', /overwrite
   file_copy, 'linespectrum_moldata_*.dat', rundir+'/.', /overwrite
   spawn, 'rm lineposvelcirc_moldata_*.dat'
   spawn, 'rm linespectrum_moldata_*.dat'
ENDIF
   
FOR iii=0,ncores-1 DO BEGIN
   IF image EQ 1 THEN BEGIN
      FOR lj=1,nlines DO BEGIN
         file_move, 'lineposvel_moldata_'+STRTRIM(STRING(iii),2)+'_'+STRTRIM(STRING(lj),2)+'.dat', rundir+$
                '/lineposvel_moldata_'+STRTRIM(STRING(iii+lj),2)+'.dat', /overwrite   
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
