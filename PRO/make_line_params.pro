PRO make_line_params,var_max_abun=var_max_abun,var_min_abun=var_min_abun,var_fr_temp=var_fr_temp,var_gtd=var_gtd,var_isot=var_isot,var_coldfinger=var_coldfinger,var_min_mu=var_min_mu,var_max_mu=var_max_mu,var_vsampling=var_vsampling,var_passband=var_passband,var_incl=var_incl

IF file_test('line_params_fixed.ini') THEN BEGIN
   @line_params_fixed.ini
ENDIF ELSE BEGIN
   @line_params.ini
ENDELSE

IF KEYWORD_SET(var_max_abun) THEN max_abun=var_max_abun
IF KEYWORD_SET(var_min_abun) THEN min_abun=var_min_abun
IF KEYWORD_SET(var_fr_temp) THEN fr_temp=var_fr_temp
IF KEYWORD_SET(var_gtd) THEN gtd=var_gtd
IF KEYWORD_SET(var_isot) THEN isot=var_isot
IF KEYWORD_SET(var_coldfinger) THEN coldfinger=var_coldfinger
IF KEYWORD_SET(var_min_mu) THEN min_mu=var_min_mu
IF KEYWORD_SET(var_max_mu) THEN max_mu=var_max_mu
IF KEYWORD_SET(var_vsampling) THEN vsampling=var_vsampling
IF KEYWORD_SET(var_passband) THEN passband=var_passband
IF KEYWORD_SET(var_incl) THEN incl=var_incl

openw,lun,'line_params.ini',/get_lun
printf, lun, 'main_path = "'+main_path+'"'
printf, lun, 'hit_path  = "'+hit_path+'"'
printf, lun, 'psumfile  = "'+psumfile+'"'
printf, lun, 'lamda_path= "'+lamda_path+'"'
printf, lun, 'exe_path  = "'+exe_path+'"'
printf, lun, 'freezefile= "'+freezefile+'"'
printf, lun, 'ncores    = '+STRING(ncores)
printf, lun, 'niterline = '+STRING(niterline)
printf, lun, 'image     = '+STRING(image)
printf, lun, 'gtd       = '+STRING(gtd)
printf, lun, 'lte       = '+STRING(lte)
printf, lun, 'gas_decoup= '+STRING(gas_decoup)
printf, lun, 'min_mu    = '+STRING(min_mu)
printf, lun, 'max_mu    = '+STRING(max_mu)
printf, lun, 'cutoff    = '+STRING(cutoff)
printf, lun, 'max_energy= '+STRING(max_energy)
printf, lun, 'turb_kep  = '+STRING(turb_kep)
printf, lun, 'turb_sou  = '+STRING(turb_sou)
printf, lun, 'fr_temp   = '+STRING(fr_temp)
printf, lun, 'H2O_OP    = '+STRING(H2O_OP)
printf, lun, 'isot      = '+STRING(isot)
printf, lun, 'max_abun  = '+STRING(max_abun)
printf, lun, 'min_abun  = '+STRING(min_abun)
printf, lun, 'abun_str  = '+STRING(abun_str)
printf, lun, 'coldfinger= '+STRING(coldfinger)
printf, lun, 'vtype     = '+STRING(vtype)
printf, lun, 'incl      = '+STRING(incl)
printf, lun, 'vsampling = '+STRING(vsampling)
printf, lun, 'passband  = '+STRING(passband)
printf, lun, 'vmax      = '+STRING(vmax)
printf, lun, 'jmax      = '+STRING(jmax)
printf, lun, 'parallel  = '+STRING(parallel)
printf, lun, 'cir_np    = '+STRING(cir_np)
printf, lun, 'b_per_r   = '+STRING(b_per_r)
printf, lun, 'b_extra   = '+STRING(b_extra)
close, lun
free_lun, lun

END

