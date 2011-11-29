
;
;Main procedure for the RADLite NLTE module. 
;
;
;

PRO nlte_main, tgas=tgas, rhogas=rhogas, abun=abun, species=species, ddens=ddens, partner_name=partner_name

@natconst
@line_params.ini

IF ~KEYWORD_SET(dV)       THEN dV     = 1d5  ;cm/s
IF ~KEYWORD_SET(parallel) THEN ncore  = 1
;
; Read the continuum mean intensity
;
mint   =  read_meanint()
np      = ddens.ntheta/2
nu_cont = ddens.nu

CASE isot OF
   51: lamda_isotop='12CO_lamda.dat'
ENDCASE

CASE partner_name OF
   'H2': partner = 0
ENDCASE

molall = READ_MOLECULE_LAMBDA(main_path+'LAMDA/'+lamda_isotop,/coll,/ghz)
mol    = EXTRACT_LAMDA(molall,vmax=0,jmax=6)

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
nctrans      = mol.nctrans[partner]
ntemps       = mol.ntemps

J_col        = DBLARR(nlines,np)

npop_all     = DBLARR(nlevels, np, ddens.nr)
npop_ini_all = DBLARR(nlevels, np, ddens.nr)

bridges        = build_bridges(ncores)
p_npop_all     = ptr_new(DBLARR(nlevels, np, ddens.nr), /no_copy)
p_npop_ini_all = ptr_new(DBLARR(nlevels, np, ddens.nr), /no_copy)

npop           = DBLARR(nlevels*np)
ini_npop       = DBLARR(nlevels*np)

FOR i=0,ncores-1 DO BEGIN
   (bridges[i])->execute, '.compile nlteC.pro'
   (bridges[i])->execute, 'resolve_all'   
   (bridges[i])->execute, '.compile callback.pro'
  (bridges[i])->setproperty, callback='callback'
ENDFOR

FOR i=0,ddens.nr-1 DO BEGIN
   
   z_col      = REFORM((!pi/2-ddens.theta[0:np-1])*ddens.r[i])  ;exact - we are approximating photons to curve along circles in the polar system
   tgas_col   = MEDIAN(REFORM(tgas[i,0:np-1]),5)
   rhogas_col = REFORM(rhogas[i,0:np-1])/(mu*mp)
   minsub     = WHERE(rhogas_col LT 1.) ;Densities cannot be too close to 0.
   rhogas_col[minsub] = 1.
   abun_col   = REFORM(abun[i,0:np-1])
   JSED_col   = REFORM(mint.meanint[i,0:np-1,*])

   FOR h=0,np-1 DO BEGIN
      J_col[*,h] = INTERPOL(JSED_col[h,*],nu_cont/cc,freq)
   ENDFOR


   IF KEYWORD_SET(parallel) THEN BEGIN
      ud     = {i:i,p_npop_all:p_npop_all, p_npop_ini_all:p_npop_ini_all}
      
      bridge = get_idle_bridge(bridges)
      
      bridge->setproperty, userdata=ud
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
      bridge->setvar, 'partner',partner
      bridge->setvar, 'np',np
      bridge->setvar, 'npop',npop
      bridge->setvar, 'ini_npop',ini_npop
      
      bridge->execute, nowait=0, 'nlteC, z_col, tgas_col, rhogas_col, abun_col, JSED_col, J_col,'+$
                       'dv, nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g,'+$
                       'collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans, partner,'+$
                       'np, npop, ini_npop' ;we can't pass an IDL structure - only arrays and scalars

   ENDIF ELSE BEGIN 
      nlteC, z_col, tgas_col, rhogas_col, abun_col, JSED_col, J_col,$
             dv, nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g,$
             collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans, partner,$
             np, npop, ini_npop
      
      npop_all[*,*, i]    = npop
      npop_ini_all[*,*,i] = ini_npop
   ENDELSE

ENDFOR

IF KEYWORD_SET(parallel) THEN BEGIN
   barrier_bridges, bridges
   npop_all     = (*p_npop_all)
   npop_ini_all = (*p_npop_ini_all)
ENDIF

MWRFITS, dum, 'levelpop_nlte.fits', /CREATE
MWRFITS, {npop_all:npop_all, npop_ini:npop_ini_all, idown:idown, iup:iup, g:g, gugl:gugl, energy_in_K:energy_in_K, $
          freq:freq, Aul:Aul, Bul:Bul, Blu:Blu, nlevels:nlevels, theta:ddens.theta[0:np-1], radius:ddens.r}, $
         'levelpop_nlte.fits'



END
