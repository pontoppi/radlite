PRO make_problem_params, var_mdisk=var_mdisk, var_mstar=var_mstar, var_tstar=var_tstar, var_rstar=var_rstar,$
 var_gtd=var_gtd, var_plh=var_plh, var_nphot=var_nphot, var_hrgrid=var_hrgrid
;Mdisk is in units of Mstar
;Rstar is in units of Rsun
;Mstar is in units of Msun

@problem_natconst.pro
IF file_test('problem_params_fixed.pro') THEN BEGIN
   @problem_params_fixed.pro
ENDIF ELSE BEGIN
   @problem_params.pro
ENDELSE

IF KEYWORD_SET(var_mdisk) THEN mdmstr=var_mdisk
IF KEYWORD_SET(var_rstar) THEN rstar=var_rstar*RS
IF KEYWORD_SET(var_mstar) THEN mstar=var_mstar*MS
IF KEYWORD_SET(var_tstar) THEN tstar=var_tstar
IF KEYWORD_SET(var_gtd) THEN gastodust=var_gtd
IF KEYWORD_SET(var_plh) THEN plh=var_plh
IF KEYWORD_SET(var_nphot) THEN nphot=var_nphot
IF KEYWORD_SET(var_hrgrid) THEN hrgrid=var_hrgrid

openw,lun,'problem_params.pro',/get_lun
printf, lun, 'imakedisk = '+STRING(imakedisk)
printf, lun, ';-----------------------------------------'
printf, lun, "simser = '"+STRING(simser)+"'"
printf, lun, "simnr = '"+STRING(simnr)+"'"
printf, lun, "simweb = '"+STRING(simweb)+"'"
printf, lun, ';-----------------------------------------'
printf, lun, 'rstar = '+STRING(rstar)
printf, lun, 'mstar = '+STRING(mstar)
printf, lun, 'tstar = '+STRING(tstar)
printf, lun, 'kurucz ='+STRING(kurucz)
printf, lun, 'ifinstar = '+STRING(ifinstar)
printf, lun, ';-----------------------------------------'
printf, lun, 'abunc = '+STRING(abunc)
printf, lun, "mixnames = ['amorph_mix.Kappa']"  ;HARD CODED
printf, lun, "mixspecs = ['carbon.Kappa','silicate.Kappa']"  ;HARD CODED
printf, lun, 'mixabun  = [abunc,(1.0-abunc)]'  ;HARD CODED
printf, lun, 'fresmd = '+STRING(fresmd)
printf, lun, "infile =  ['amorph_mix.Kappa']"  ;HARD CODED
printf, lun, 'pll =[-1]'
printf, lun, 'scat =  '+STRING(scat)
printf, lun, ';-----------------------------------------'
printf, lun, 'gastodust = '+STRING(gastodust)
printf, lun, 'rhofloor = '+STRING(rhofloor)
printf, lun, 'run = '+STRING(run)
printf, lun, ';-----------------------------------------'
printf, lun, 'nr = '+STRING(nr)
printf, lun, 'nt = '+STRING(nt)
printf, lun, 'ntex = '+STRING(ntex)
printf, lun, 'rin = '+STRING(rin)
printf, lun, 'tin = '+STRING(tin)
printf, lun, 'rout = '+STRING(rout)
printf, lun, 'hrgrid = '+STRING(hrgrid)
printf, lun, 'hrgmax = '+STRING(hrgrid*1.4)
printf, lun, 'thintin = '+STRING(thintin)
printf, lun, ';-----------------------------------------'
printf, lun, 'rrefine = {nlevr:3, nspanr:3, nstepr:3}'  ;HARD CODED
printf, lun, 'drsm = '+STRING(drsm)
printf, lun, ';-----------------------------------------'
printf, lun, 'rdisk = '+STRING(rdisk)
printf, lun, 'sig0 = '+STRING(sig0)
printf, lun, 'mdmstr = '+STRING(mdmstr)
printf, lun, 'mdisk  = mdmstr * mstar'
printf, lun, 'plsig1 = '+STRING(plsig1)
printf, lun, 'plsig2 = '+STRING(plsig2)
printf, lun, 'bgdens = '+STRING(bgdens)
printf, lun, ';-----------------------------------------'
printf, lun, 'hrdisk = '+STRING(hrdisk)
printf, lun, 'hrmin = '+STRING(hrmin)
printf, lun, 'plh = '+STRING(plh)
printf, lun, ';-----------------------------------------'
printf, lun, 'rpfrin = '+STRING(rpfrin)
printf, lun, 'hrpuff = '+STRING(hrpuff)
printf, lun, ';-----------------------------------------'
printf, lun, 'nphot = '+STRING(nphot)
printf, lun, 'npdiff = '+STRING(npdiff)
printf, lun, 'errtol = '+STRING(errtol)
printf, lun, 'ifast = '+STRING(ifast)
printf, lun, ';-----------------------------------------'
printf, lun, 'nvstr = '+STRING(nvstr)
printf, lun, 'vserrt = '+STRING(vserrt)
printf, lun, 'ivstrt = '+STRING(ivstrt)
printf, lun, 'dostr  = [1, 1]'  ;HARD CODED
printf, lun, ';-----------------------------------------'
printf, lun, 'tauchop= 5.d3' 
printf, lun, 'lmbchop= 0.55d0'
printf, lun, 'idxchop= 1.d0'
printf, lun, ';-----------------------------------------'
printf, lun, 'xlevel   = 1'
printf, lun, 'cd,current=current'
printf, lun, "for i=1,xlevel do cd,'../'"
printf, lun, 'cd,current=maindir'
printf, lun, "cd,'~/'"
printf, lun, 'cd,current=homedir'
printf, lun, 'cd,current'
printf, lun, "bindir=maindir+'/bin/'"
printf, lun, ';-----------------------------------------'
printf, lun, "maindir = '/Users/pontoppi/WORK/RADLITE'"

printf, lun, "kuruczdir = maindir+'/KURUCZ/'"
printf, lun, "src_dust  = maindir+'/sources/makeopac/src/'"
printf, lun, "dustdata  = maindir+'/sources/makeopac/optconst/'"
printf, lun, "src_rayt  = maindir+'/sources/raytrace/src'"
printf, lun, "src_rdmc  = maindir+'/sources/radmc/src'"
printf, lun, "src_chop  = maindir+'/sources/chopdens/src'"
printf, lun, "src_pahc  = maindir+'/sources/pahcode/src'"
close, lun
free_lun, lun

END
