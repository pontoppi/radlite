@analyze_v2
PRO nlte_main, tgas=tgas, rhogas=rhogas, abun=abun, species=species, ddens=ddens
COMMON coll,  Cul, TCul
COMMON mol,  nlines, nlevels, gugl, freq, iup, idown, Aul, Bul, Blu, energy_in_k, g
COMMON grid, np, z_col, tgas_col, abun_col, rhogas_col, J_col, JSED_col, nu_cont
@natconst
;
; Read the continuum mean intensity
;
mint   =  read_meanint()
np      = ddens.ntheta/2
nu_cont = ddens.nu

;HITRAN_EXTRACT, isotop = 51, lambdarange=[800,4000.74], max_energy = 400., vmax=1
mol    = READ_MOLECULE_LAMBDA('moldata.dat',/coll)
nlines = N_ELEMENTS(mol.freq)
Cul    = FLTARR(nlines,10)
TCul   = [2.,5.,10.,50.,200,400,800,1500,3000.,10000.]

PureRot = WHERE(mol.freq LT 1000.)
RoVib   = WHERE(mol.freq GE 1000.)
Cul[PureRot,*] = 3d-11
Cul[RoVib,*]   = 1d-15

idown       = mol.idown
iup         = mol.iup
g           = mol.g
gugl        = g[iup-1]/g[idown-1]
energy_in_K = mol.energy_in_K
freq        = mol.freq
Aul         = mol.Aud
Bul         = (1d0/(2d0*hh*cc*freq^3d0))*Aul
Blu         = gugl*Bul
nlevels     = mol.nlevels


npop_all = DBLARR(nlevels, np, ddens.nr)

FOR i=0,ddens.nr-1 DO BEGIN
   
   z_col      = REFORM((!pi/2-ddens.theta[0:np-1])*ddens.r[i])  ;exact - we are approximating photons to curve along circles in the polar system
   tgas_col   = MEDIAN(REFORM(tgas[i,0:np-1]),5)
   rhogas_col = REFORM(rhogas[i,0:np-1])/(mu*mp)
   minsub = WHERE(rhogas_col LT 100.) ;Densities cannot be too close to 0.
   rhogas_col[minsub] = 100.
   abun_col   = REFORM(abun[i,0:np-1])
   JSED_col   = MEDIAN(REFORM(mint.meanint[i,0:np-1,*]),5)
   
   nlte, species=species, npop=npop
   npop_all[*,*, i] = npop
   print, i
ENDFOR

MWRFITS, dum, 'levelpop_nlte.fits', /CREATE
MWRFITS, {npop_all:npop_all, idown:idown, iup:iup, g:g, gugl:gugl, energy_in_K:energy_in_K, $
          freq:freq, Aul:Aul, Bul:Bul, Blu:Blu, nlevels:nlevels, theta:ddens.theta[0:np-1], radius:ddens.r}, $
         'levelpop_nlte.fits'

stop

END
