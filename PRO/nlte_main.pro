
;==========================================================
;Main procedure for the RADLite NLTE module. 
;
;Uses a Newton-Raphson global solver (e.g., Numerical Recipes,
;Chapter 9.7). The module currently has a framework for
;parallelization using the IDL_IDLBRIDGE class. 
;
;==========================================================

PRO nlte_main, molall=molall, tgas=tgas, rhogas=rhogas, abun=abun, species=species, ddens=ddens

@natconst
@line_params.ini

IF ~KEYWORD_SET(dV)       THEN dV     = 1d5  ;cm/s
IF ~KEYWORD_SET(parallel) THEN ncore  = 1
IF ~KEYWORD_SET(vmax)     THEN vmax   = 0
IF ~KEYWORD_SET(jmax)     THEN jmax   = 10
;
; Read the continuum mean intensity
;
mint   =  read_meanint()
np      = ddens.ntheta/2
nu_cont = ddens.nu

mol    = LAMDA_EXTRACT_LEVELS(molall,vmax=vmax,jmax=jmax)

nlines = N_ELEMENTS(mol.freq)

idown        = mol.idown
iup          = mol.iup
g            = mol.g
gugl         = g[iup-1]/g[idown-1]
energy_in_K  = mol.energy_in_K
freq         = mol.freq
Aul          = mol.Aud
Bul          = (1d0/(2d0*hh*cc*freq^3d0))*Aul
Blu          = gugl*Bul
nlevels      = mol.nlevels

collrates    = mol.collrates
coll_iup     = mol.coll_iup
coll_idown   = mol.coll_idown
coll_temps   = mol.temps
nctrans      = mol.nctrans
ntemps       = mol.ntemps

J_col        = DBLARR(nlines,np)

npop_all     = DBLARR(nlevels, np, ddens.nr)
npop_ini_all = DBLARR(nlevels, np, ddens.nr)
npop_lte_all = DBLARR(nlevels, np, ddens.nr)


npop           = DBLARR(nlevels,np)
ini_npop       = DBLARR(nlevels,np)
lte_npop       = DBLARR(nlevels,np)

IF KEYWORD_SET(parallel) THEN BEGIN
   bridges        = build_bridges(ncores)
   p_npop_all     = ptr_new(DBLARR(nlevels, np, ddens.nr), /no_copy)
   p_npop_ini_all = ptr_new(DBLARR(nlevels, np, ddens.nr), /no_copy)
   p_npop_lte_all = ptr_new(DBLARR(nlevels, np, ddens.nr), /no_copy)
   FOR i=0,ncores-1 DO BEGIN
      (bridges[i])->execute, '.compile Ppump.pro'
      (bridges[i])->execute, '.compile Pesc.pro'
      (bridges[i])->execute, '.compile P.pro'
      (bridges[i])->execute, '.compile nlte.pro'
      (bridges[i])->execute, '.compile callback.pro'
      (bridges[i])->setproperty, callback='callback'
   ENDFOR
ENDIF

FOR i=0,ddens.nr-1 DO BEGIN
  
   z_col      = REFORM((!pi/2-ddens.theta[0:np-1])*ddens.r[i])  ;exact - we are approximating photons to curve along circles in the polar system
   tgas_col   = REFORM(tgas[i,0:np-1])
   rhogas_col = REFORM(rhogas[i,0:np-1])/(mu*mp)
   abun_col   = REFORM(abun[i,0:np-1])
   JSED_col   = REFORM(mint.meanint[i,0:np-1,*])
   
   ;Change the input sampling to one that is optimal for the non-LTE calculation
   optimal_sample, z_col, tgas_col, rhogas_col, abun_col, JSED_col

   FOR h=0,np-1 DO BEGIN
      J_col[*,h] = INTERPOL(JSED_col[h,*],nu_cont/cc,freq)
   ENDFOR

   IF KEYWORD_SET(parallel) THEN BEGIN

      ud     = {i:i,p_npop_all:p_npop_all, p_npop_ini_all:p_npop_ini_all, p_npop_lte_all:p_npop_lte_all}
      
      bridge = get_idle_bridge(bridges)
      
      bridge->setproperty, userdata=ud
      bridge->setproperty
      bridge->setvar, 'rr', ddens.r[i]/AU
      bridge->setvar, 'z_col', z_col
      bridge->setvar, 'tgas_col', tgas_col
      bridge->setvar, 'rhogas_col', rhogas_col
      bridge->setvar, 'abun_col', abun_col
      bridge->setvar, 'JSED_col', JSED_col
      bridge->setvar, 'J_col', J_col
      bridge->setvar, 'dv',dv
      bridge->setvar, 'nlines',nlines
      bridge->setvar, 'nlevels',nlevels
      bridge->setvar, 'gugl',gugl
      bridge->setvar, 'freq', freq
      bridge->setvar, 'iup',iup
      bridge->setvar, 'idown',idown
      bridge->setvar, 'Aul',Aul
      bridge->setvar, 'Bul',Bul
      bridge->setvar, 'Blu',Blu
      bridge->setvar, 'energy_in_k',energy_in_k
      bridge->setvar, 'g',g
      bridge->setvar, 'collrates',collrates
      bridge->setvar, 'coll_iup',coll_iup
      bridge->setvar, 'coll_idown',coll_idown
      bridge->setvar, 'coll_temps',coll_temps
      bridge->setvar, 'ntemps',ntemps
      bridge->setvar, 'nctrans',nctrans
      bridge->setvar, 'np',np
      bridge->setvar, 'npop',npop
      bridge->setvar, 'ini_npop',ini_npop
      bridge->setvar, 'lte_npop',lte_npop

      bridge->execute, "PRINT, 'Radius = ', rr"
      bridge->execute, /nowait, 'nlte, z_col, tgas_col, rhogas_col, abun_col, JSED_col, J_col,'+$
                       'dv, nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g,'+$
                       'collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans, '+$
                       'np, npop, ini_npop, lte_npop' ;we can't pass an IDL structure - only arrays and scalars: hence this clunky arg list.

      writeu,-1,STRING(13b)
      writeu,-1,'Running statistical balance calculation in parallel... '+STRTRIM(STRING(100d0*i/(ddens.nr-1),format='(f6.2)'),2) + '% complete'


   ENDIF ELSE BEGIN 
      print, 'Radius: ', ddens.r[i]/AU, ' AU'
      nlte, z_col, tgas_col, rhogas_col, abun_col, JSED_col, J_col,$
             dv, nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g,$
             collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans, $
             np, npop, ini_npop, lte_npop
      
      npop_all[*,*, i]    = ABS(npop)
      npop_ini_all[*,*,i] = ini_npop
      npop_lte_all[*,*,i] = lte_npop
      
   ENDELSE

ENDFOR

IF KEYWORD_SET(parallel) THEN BEGIN
   barrier_bridges, bridges
   npop_all     = (*p_npop_all)
   npop_ini_all = (*p_npop_ini_all)
   npop_lte_all = (*p_npop_lte_all)
   ptr_free, p_npop_all
   ptr_free, p_npop_ini_all
   ptr_free, p_npop_lte_all
   burn_bridges, bridges
ENDIF

PRINT, ' '

MWRFITS, dum, 'levelpop_nlte.fits', /CREATE, /SILENT
MWRFITS, {npop_all:npop_all, npop_ini:npop_ini_all, npop_lte:npop_lte_all, theta:ddens.theta[0:np-1], radius:ddens.r}, $
         'levelpop_nlte.fits', /SILENT
MWRFITS, mol, 'levelpop_nlte.fits', /SILENT


END
