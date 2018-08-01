PRO rebuild_runtable, run_name=run_name

	IF ~KEYWORD_SET(run_table) THEN run_table='run_table.fits'

	run_pars = []
	line = ''
	
	spawn, 'ls | fgrep "'+run_name+'"',run_list
	FOR i=0,N_ELEMENTS(run_list)-1 DO BEGIN
		cd, run_list[i]
		openr, lun, 'line_params.ini', /get_lun
		WHILE (NOT EOF(lun)) DO BEGIN
			READF, lun, line
			words = STRTRIM(STRSPLIT(line, '=',/EXTRACT),2)
			
			IF words[0] EQ 'coldfinger' THEN coldfinger = FLOAT(words[1])
			IF words[0] EQ 'min_abun' THEN min_abun = FLOAT(words[1])
			IF words[0] EQ 'max_abun' THEN max_abun = FLOAT(words[1])
			IF words[0] EQ 'min_mu' THEN min_mu = FLOAT(words[1])
			IF words[0] EQ 'max_mu' THEN max_mu = FLOAT(words[1])			
		ENDWHILE
		close, lun
		free_lun, lun

		openr, lun, 'problem_params.pro', /get_lun
		WHILE (NOT EOF(lun)) DO BEGIN
			READF, lun, line
			words = STRTRIM(STRSPLIT(line, '=',/EXTRACT),2)
		
			IF words[0] EQ 'gastodust' THEN g2d = FLOAT(words[1])
			IF words[0] EQ 'plh' THEN plh = FLOAT(words[1])
			IF words[0] EQ 'mstar' THEN mstar = FLOAT(words[1])
			IF words[0] EQ 'mdmstr' THEN mdisk = mstar*FLOAT(words[1])
			IF words[0] EQ 'rstar' THEN rstar = FLOAT(words[1])
			IF words[0] EQ 'tstar' THEN tstar = FLOAT(words[1])
			
		ENDWHILE
		close, lun
		free_lun, lun
		
		run_par = {dir:run_list[i],coldfinger:coldfinger,min_abun:min_abun,max_abun:max_abun,linepos:min_mu+(max_mu-min_mu)/2.,$
				     mdisk:mdisk, g2d:g2d, plh:plh, mstar:mstar, rstar:rstar, tstar:tstar}
		run_pars = [run_pars,run_par]
		print, run_par
		cd, '..'
	ENDFOR
	
	MWRFITS, dummy, run_table, /CREATE
	MWRFITS, run_pars, run_table
	
END