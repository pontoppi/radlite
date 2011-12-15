
;
;Main procedure for the RADLite NLTE module. 
;
;
;


@analyze_v2
PRO nlte_main_old, tgas=tgas, rhogas=rhogas, abun=abun, species=species, ddens=ddens, partner_name=partner_name
COMMON coll,  Cul, TCul
COMMON mol,  nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g, collrates, coll_iup, coll_idown, coll_temps, ntemps, nctrans, partner
COMMON grid, np, z_col, tgas_col, abun_col, rhogas_col, J_col, JSED_col, nu_cont
@natconst
@line_params.ini
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
mol    = EXTRACT_LAMDA(molall,vmax=2,jmax=10)

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

npop_all     = DBLARR(nlevels, np, ddens.nr)
npop_ini_all = DBLARR(nlevels, np, ddens.nr)

FOR i=0,ddens.nr-1 DO BEGIN
   
   z_col      = REFORM((!pi/2-ddens.theta[0:np-1])*ddens.r[i])  ;exact - we are approximating photons to curve along circles in the polar system
   tgas_col   = MEDIAN(REFORM(tgas[i,0:np-1]),5)
   rhogas_col = REFORM(rhogas[i,0:np-1])/(mu*mp)
   minsub     = WHERE(rhogas_col LT 100.) ;Densities cannot be too close to 0.
   rhogas_col[minsub] = 100.
   abun_col   = REFORM(abun[i,0:np-1])
   JSED_col   = MEDIAN(REFORM(mint.meanint[i,0:np-1,*]),5)
   
   nlte, species=species, npop=npop, ini_npop=ini_npop
   npop_all[*,*, i]    = npop
   npop_ini_all[*,*,i] = ini_npop
   print, i
ENDFOR

MWRFITS, dum, 'levelpop_nlte.fits', /CREATE
MWRFITS, {npop_all:npop_all, npop_ini:npop_ini_all, idown:idown, iup:iup, g:g, gugl:gugl, energy_in_K:energy_in_K, $
          freq:freq, Aul:Aul, Bul:Bul, Blu:Blu, nlevels:nlevels, theta:ddens.theta[0:np-1], radius:ddens.r}, $
         'levelpop_nlte.fits'

stop

END
