;=======================================================================
;                     COLLECTION OF MODEL-SETUPS
;
; This file contains ready-to-use model-setups, including everything.
; These routines call various other subroutines in problem_disk.pro,
; problem_shuulrich.pro, problem_files.pro etc, which contain standard
; model components. Here these model components are glued together into
; a full things, and the files are written out, and optionally the
; RADMC code is actually called. 
;
; THESE ARE THE ROUTINES THAT THE USER SHOULD USE!
; 
; The following models are available:
;
;   diskenv_radmc():         Disk and envelope model, in which the sigma(R)
;                            is given, for instance as output from another 
;                            code such as diskevol or so. The envelope is
;                            simply parameterized according to the Shu and
;                            Ulrich models.
;
;   makedisk_vertstruct():   A routine that is used in the initial setups 
;                            of dust coagulation. It is actually nearly
;                            the same as diskenv_radmc(), but with an 
;                            interface back to the dustcoag routines.
;
;   simpledisk_vertstruct(): This is the D&D05 setup in which the basic
;                            disk structure is parameterized. But the
;                            vertical structure can be computed
;
;=======================================================================

@problem_subroutines.pro
@problem_makeopac.pro
@problem_kurucz.pro
@problem_makeshu.pro
@problem_grid.pro
@problem_radmc_files.pro
@problem_disk.pro

pro simplediskenv,rstar=rstar,tstar=tstar,mstar=mstar,$
         ifinstar=ifinstar,fresmd=fresmd,$
         scat=scat,ntex=ntex,nrr=nrr,ntt=ntt,rin=rin,tin=tin,$
         rout=rout,hrgrid=hrgrid,hrgmax=hrgmax,rrefine=rrefine,$
         drsm=drsm,rdisk=rdisk,sigdust0=sigdust0,mdisk=mdisk,$
         plsig1=plsig1,plsig2=plsig2,$
         opacnames=opacnames,pllongs=pllongs,$
         schmidt=schmidt,ab_r0=ab_r0,ab_ab0=ab_ab0,ab_pl=ab_pl,$
         gastodust=gastodust,ab_min=ab_min,hrdisk=hrdisk,hrmin=hrmin,$
         plh=plh,rpfrin=rpfrin,hrpuff=hrpuff,nvstr=nvstr,$
         ivstrt=ivstrt,vserrtol=vserrtol,nphot=nphot,$
         npdiff=npdiff,errtol=errtol,tauchop=tauchop,lmbchop=lmbchop,$
         idxchop=idxchop,rhofloor=rhofloor,pnc=pnc,pnh=pnh,pz=pz,$
         imakedisk=imakedisk,run=run,hrstore=hrstore,$
         thintin=thintin,tt=tt,radius=r,theta=theta,rhodust=rhodust,$
         sigdust=sigdust,kurucz=kurucz,kurdir=kurdir,bindir=bindir,$
         ifast=ifast,dostr=dostr,csenv=csenv,time=time,env=env
@problem_natconst.pro

; Copy nr and nt to local variables
;
; NOTE: Internally the nt is always the full theta array size,
;       hence the nt = ntt+ntex
;
nr = nrr
nt = ntt+ntex
;
; ...Grid refinement near inner edge
;
if not keyword_set(rrefine) then begin
    rrefine={nlevr:3,nspanr:2,nstepr:2}
endif
;
; Find the outer radius of the disk
;
rin  = min(rdisk)
rout = max(rdisk)

;
;               SET UP THE 2-D RADIATIVE TRANSFER STUFF
;
; Make the R-grid
;
r      = make_rgrid(rin,rout,nr,rrefine=rrefine)
;
; Make the Theta-grid
;
theta  = make_tgrid(hrgrid,nt,hrgmax=hrgmax,ntex=ntex,$
                    hrlg=hrlg,zrefine=zrefine)
;
; Check the grid sizes
;
if nr ne n_elements(r) then stop
if nt ne n_elements(theta) then stop
;
; Now create the disk model using the routines of problem_disk.pro
; NOTE: Single dust species here, and we convert back to gas density
;       because this is what we are interested in here.
;
if keyword_set(sigdust0) then begin
   ;;
   ;; Directly from sigdust0
   ;;
   rhodusttot = disk_model_1(r,theta,rdisk,sigdust0,plsig1,plsig2,$
                             hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                             hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdust,gap=gap,gdrop=gdrop)
endif else begin
   ;;
   ;; Compute sigdust0 from mdisk
   ;;
   sigdust00  = 1.d0
   rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                             hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                             hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot,gap=gap,gdrop=gdrop)
   mddum      = integrate(r,2*pi*r*sigdusttot)*gastodust
   sigdust00  = mdisk / mddum
   rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                             hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                             hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot,gap=gap,gdrop=gdrop)
endelse

rhodisk = rhodisk * gastodust
;
; Compute the local qplus from the vertically integrated qpldisk
;
if keyword_set(qpldisk) then begin
   qloc   = dblarr(nr,nt)
   for ir=0,irdisk do begin
      ;;
      ;; Get the z coordinate
      ;;
      z          = (!dpi/2-theta)*r[ir]
      ;;
      ;; Basically take the qloc to follow the shape of the 
      ;; density profile
      ;;
      qloc[ir,*] = rhodisk[ir,*]
      if keyword_set(qsig) then begin
         ;;
         ;; If qsig is set, then require all the energy release to be
         ;; in a surface layer of the disk down to sigma=qsig. This 
         ;; can be useful to speed up the MC code.
         ;;
         sig= integrate(z,rhodisk[ir,*],/prim)
         ii = where(sig gt qsig and max(sig)-sig gt qsig)
         if ii[0] ne -1 then qloc[ir,ii] = 0.d0
      endif
      ;;
      ;; Find the qplus here
      ;;
      qpld = interpol(qpldisk,rdisk,r[ir])
      ;;
      ;; Normalize the qloc such that the total released 
      ;; energy equals qpldisk
      ;;
      qq         = 2*integrate(z,transpose(qloc[ir,*]))
      qloc[ir,*] = qloc[ir,*] * ( qpld / qq )
   endfor
endif
;
; Add the floor density
;
rhodisk = rhodisk + rhofloor

;
; If requested, include the Shu infall model
;
stop
IF KEYWORD_SET(env) THEN BEGIN
   IF n_elements(time) NE 1 OR n_elements(csenv) NE 1 OR $
      n_elements(csenv) NE 1 then begin
      print,'ERROR: If you include the envelope, you must also define'
      print,'   time and csenv...'
      stop
   endif
   q    = make_shu(r,time,csenv=csenv)
   
   ;; Mutually exclude the densities
   ;;
   ii = WHERE( rhodisk GT q.rho ) 
   IF ii[0] GE 0 THEN q.rho[ii] = 0.d0
   ii = WHERE( rhodisk LE q.rho ) 
   IF ii[0] GE 0 THEN rhodisk[ii] = 0.d0
   PRINT, 'Included a Shu envelope!'
ENDIF
;
;---------------------------------------------------------------------
;             NOW ASSEMBLE EVERYTHING AND WRITE IT OUT
;---------------------------------------------------------------------
;
; Now create the full density profile with all species
; (when envelope is switched off, this is same as rhodisk/gastodust)
;
nspec    = 1
rhodust  = dblarr(nr,nt,nspec)
IF KEYWORD_SET(env) THEN BEGIN
   rhodust[*,*,0] = (rhodisk+q.rho) / gastodust
ENDIF ELSE BEGIN
   rhodust[*,*,0] = rhodisk / gastodust
ENDELSE
;
; Write the grid and the density profiles for RADMC
;
write_radmc_density,r,theta,rhodust
;
; Write the star spectrum and info for RADMC
;
write_star,mstar,rstar,tstar,kurucz=kurucz,kurdir=kurdir
;
; If chopping is active, then write the chop file, else delete
;
if keyword_set(tauchop) and keyword_set(lmbchop) $
  and keyword_set(idxchop) then begin
    write_chopdens,tauchop,lmbchop,idxchop
    print,'*** CHOPPING ACTIVE! ***'
    spawn,bindir+'/chopdens'
endif else begin
    spawn,'rm -f chopdens.inp'
endelse
;
; Write the qplus.inp file
;
if keyword_set(qpldisk) then begin
   openw,1,'qplus.inp'
   printf,1,nr,nt
   for ir=0,nr-1 do begin
      for it=0,nt-1 do begin
         printf,1,qloc[ir,it]+1d-99
      endfor
   endfor
   close,1
endif else begin
   spawn,'rm -f qplus.inp'
endelse
;
; Write the raytrace.inp for the ray-tracing (spectra and images)
;
write_raytrace,nextra=nextra,nrref=nrref
;
; Write the radmc.inp file
;
; ADDED 09.01.07: ifinstar
;
write_radmc,tdec=tdec,nphot=nphot,ifinstar=ifinstar,npdiff=npdiff,$
            nvstr=nvstr,ivstrt=ivstrt,vserrtol=vserrtol,ifast=ifast
;
; If requested, install an opacity
;
if keyword_set(opacname) then begin
    if not keyword_set(pllongs) then pllongs=[-2]
    useopac,[opacname],pllongs,fresmd=fresmd,scat=scat
endif
;
; Return some information concerning the envelope
;
if n_elements(q) ne 0 then info = q
;
; If we want to run this, then do it
;
if keyword_set(run) then begin
    spawn,'radmc'
    spawn,'raytrace spectrum'
endif
;
end



;------------------------------------------------------------------------
;      MAKE 2-D VERTICAL STRUCTURE MODEL OF 1-D RADIAL DISK MODEL
;
;               ***** COMPATIBILITY ROUTINE ******
;
; Given a 1-D radial disk model, this routine produces a 2-D disk setup,
; runs RADICAL and VERTSTRUCT to produce the vertical structure of the
; disk, and then re-converts (upon request) back to 1+1D structure.
;
; (This routine must at some point be simplified or integrated with
; the diskenvsed or so; But for now it is needed for the coupling
; to dustcoag.F)
;
; ARGUMENTS:
;  mstar             Mass of the star [g]
;  rstar             Radius of the star [cm]
;  tstar             Effective temperature of the star [K]
;  rdisk             Array of radial grid for the input disk model [cm]
;  sigmagasdisk      Array of Sigma_gas(R) of the disk [g/cm^2]
;  qpldisk           The viscous dissipation (vert.integr; both sides) 
;                    [erg/cm^2/s]
;  qsig              If set, then release all Qplus energy in surface
;                    layer of DeltaSigma=qsig
;  rrim              Radius within which the dust is assumed evaporated [cm]
;  nnr               Nr of radial grid points for 2-D/3-D RT problem
;  nnt               Basic nr of theta grid points for 2-D/3-D RT problem
;  hrgrid            Height (pi/2-theta) of main Theta grid 
;  hrgmax            Height (pi/2-theta) of total Theta grid 
;  ntex              Extra grid points in between hrgrid and hrgmax
;  zrref             (For refinement near the midplane): refinement parameter
;  nrzef             (For refinement near the midplane): nr of extra z points
;  gastodust         Gas-to-dust mass ratio. Default 100.
;  rhofloor          The density floor value
;  nib               [Technical] Number of rays inward of Rin for SED.
;                    Normally put this to 4 or so. But sometimes this
;                    causes problems when the star radius is not too much
;                    smaller than Rin. Then put it to 0.
;  nphot             Nr of photon packages
;  flang             If set: flaring angle to take for initial guess
;  mugas             Mean molecular weight of gas
;  save              If set, then save all intermediate VET results
;                    (all vertical structure density results are always
;                    saved)
;  ifinstar          =1 -> Use finite-size star
;  zrefine           A structure for making a very refined grid in 
;                    theta near the midplane. See problem_grid.pro the
;                    routine make_tgrid(). Note that if this option
;                    is used, the nr of theta grid points will exceed
;                    the value ntt+ntex.
; (For returned array:)
;  zdisk             The full z-coordinates (for each rdisk) of the disk
;  nzdisk            If not defined zdisk, then define nr of z points
;  rrefine           R-grid refinement structure
;  maindir           Where can we find the codes? Default ../
;  run               =1 --> Run the RADMC code
;                    =0 --> Write a dorun script to do this
;  ifast             [Careful!!] When this is 1 or 2 then the calculations
;                    of the dust temperature are not done upon each
;                    interaction of the photon with a grid cell, but only if
;                    the energy in the cell is appreciable larger than the
;                    last update. NOTE: This is not related to the new
;                    improved RADMC version, which is anyway much faster
;                    than the old one. Perhaps the ifast is even no longer
;                    necessary, but that remains to be seen.
;
; RETURNS:
;  struct            The 1+1D structure of the disk in r,z, but with
;                    interpolation of z along theta.
;
; NOTE:
;  - In this routine, for backward compatilibity reasons, the ntt is
;    the basic nr of theta points. The total nr is ntt+ntex.
;
;------------------------------------------------------------------------
pro makedisk_vertstruct,mstar,rstar,tstar,rdisk,sigmagasdisk,$
               rrim=rrim,nrr=nrr,ntt=ntt,ntex=ntex,$
               gastodust=gastodust,zref=zrref,nref=nzref,$
               rhofloor=rhofloor,qpldisk=qpldisk,zrefine=zrefine,$
               hrgrid=hrgrid,hrgmax=hrgmax,info=info,$
               nib=nib,nphot=nphot,qsig=qsig,flang=flang,mugas=mugas,$
               nvstr=nvstr,save=save,struct=struct,$
               nzdisk=nzdisk,zdisk=zdisk,rrefine=rrefine,maindir=maindir,$
               ifinstar=ifinstar,ifast=ifast
@problem_natconst.pro
;
; Default
;
if not keyword_set(flang) then flang=0.05
if not keyword_set(mugas) then mugas=2.3
if not keyword_set(nvstr) then nvstr=3
if not keyword_set(gastodust) then gastodust=100.
if not keyword_set(maindir) then maindir='../'
;
; First make an estimate of the temperature
;
tgas0  = flang^0.25 * tstar * (rstar/rdisk)^0.5 / sqrt(2.0)
hpdisk = sqrt(kk*tgas0*rdisk^3/(mugas*mp*GG*mstar))
;
; Set up the initial 2-D setup
;
diskenv_radmc,mstar,rstar,tstar,rdisk,sigmagasdisk,hpdisk,$
               rrim=rrim,nrr=nrr,ntt=ntt,ntex=ntex,$
               gastodust=gastodust,zref=zrref,nref=nzref,$
               rhofloor=rhofloor,qpldisk=qpldisk,$
               hrgrid=hrgrid,hrgmax=hrgmax,info=info,$
               nib=nib,nphot=nphot,qsig=qsig,rrefine=rrefine,$
               zrefine=zrefine,npdiff=npdiff,$
               ifinstar=ifinstar,nvstr=nvstr,/run
;
; Now read the final structure
;
a=read_dustdens()
b=read_dusttemp()
;
; Retrieve nr and nt
;
nr = n_elements(a.r)
nt = n_elements(a.theta)/2
;
; Make z-grid, if not given
;
nrdisk=n_elements(rdisk)
if not keyword_set(zdisk) then begin
    if not keyword_set(hrgrid) then hrgrid = 0.6
    if not keyword_set(nzdisk) then nzdisk = 60
    zdisk = dblarr(nzdisk,nrdisk)
    for irdisk=0,nrdisk-1 do begin
        zdisk[*,irdisk]=hrgrid*rdisk[irdisk]*dindgen(nzdisk)/(nzdisk-1.d0)
    endfor
endif
nzdisk=n_elements(zdisk[*,0])
;
; Interpolate onto 1+1D disk grid, using standard routine
;
from_2d_to_11d,a.r,a.theta,a.rho,rdisk,transpose(zdisk),rhogas
from_2d_to_11d,b.r,b.theta,b.temp,rdisk,transpose(zdisk),tgas
rhogas = rhogas * gastodust
rhogas = transpose(rhogas)
tgas   = transpose(tgas)
;
; Return...
;
struct={rdisk:rdisk,zdisk:zdisk,rhogas:rhogas,tgas:tgas}
;
end




;------------------------------------------------------------------------
;               BASIC DISK STRUCTURE: THE D&D05 SETUP
;
; If one looks for a simple routine to compute the D&D05 models, of 
; course with the newest techniques, then this is the routine!
; Note that the original routines were different. But is should be
; reasonably compatible.
;
; ARGUMENTS:
;  imakedisk         = 1   (This is normally what you want):
;                          Make an ad-hoc disk structure, which can be
;                          iterated on (set nvstr > 0), so that a vertically
;                          hydrostatic model comes out. If nvstr=0 then
;                          only the ad-hoc model is produced. Note that
;                          nvstr is written to the radmc.inp file which 
;                          tells whether or not to iterate on the structure.
;                    =-1   Read dust density (1 species) from old run
;                          to serve as the basis upon which to put
;                          different dust species. You may want to do 
;                          this if you do not want to recompute the 
;                          vertical structure all the time, but instead
;                          take a precalculated structure from the local
;                          directory. This is an expert option.
;  mstar             Mass of the star [g]
;  rstar             Radius of the star [cm]
;  tstar             Effective temperature of the star [K]
;  kurucz            If set, then use Kurucz model
;  kurdir            The directory where the Kurucz models are located
;  rdisk             Array of radial grid for the input disk model [cm]
;  nrr               Nr of radial grid points for 2-D/3-D RT problem
;  ntt               Basic nr of theta grid points for 2-D/3-D RT problem
;                      NOTE: Here it is the full nr, so
;                            we shall verify that n_elements(theta).eq.ntt
;  rin               Radius of inner wall
;  tin               The inner wall temperature (---> rin computed from this)
;  thintin           If set --> Use real opacities to compute rin(T=tin)
;  rout              Outer radius of the grid
;  hrgrid            Height (pi/2-theta) of main Theta grid 
;  hrgmax            Height (pi/2-theta) of total Theta grid 
;  ntex              Extra grid points in between hrgrid and hrgmax
;  zrref             (For refinement near the midplane): refinement parameter
;  nrzef             (For refinement near the midplane): nr of extra z points
;  gastodust         Gas-to-dust mass ratio. Default 100.
;  rhofloor          The density floor value
;  tdec              Enable thermal decoupling of temperatures?
;  nib               [Technical] Number of rays inward of Rin for SED.
;                    Normally put this to 4 or so. But sometimes this
;                    causes problems when the star radius is not too much
;                    smaller than Rin. Then put it to 0.
;  nphot             Nr of photon packages
;  npdiff            Nr photons per cell below wich the diffusion is used
;  tauchop           The max allowed vert optical depth before this is 
;                    going to be capped by the chopping routine.
;                    (the chopping routine simply removes mass from the 
;                    equator until tau<=chop. this is necessary in order 
;                    not to stall the Monte Carlo simulation)
;  lmbchop           Wavelength at which this tau is computed
;  ifinstar          =1 -> uses finite size star (iso point source)
;                    (ADDED 09.01.07)
;  rrefine           R-grid refinement structure
;  hrdisk            Value of H_p(R=Rdisk)/Rdisk for parameterized setup
;                      (NOTE: This value gets lost when iterating on struct)
;  plh               The powerlaw for H_p(R)/R
;                      (NOTE: This value gets lost when iterating on struct)
;  hrmin             Floor value for H_p/R
;                      (NOTE: This value gets lost when iterating on struct)
;  hrpuff            [optional] The artificially puffed-up rim: Hrim/Rrim
;                      (NOTE: This value gets lost when iterating on struct)
;  rpfrin            [optional] The artificially puffed-up rim: Rrim=rpfrin*Rin
;                      (NOTE: This value gets lost when iterating on struct)
;  ab_r0             Radius of pivot for simple powerlaw abundance profile
;  ab_ab0            Abundances (of species 2,3,4...) at R=ab_r0
;  ab_pl             Powerlaw index of abundance profile for species 2,3,4...
;  ab_min            Floor value for these abundances
;  opacnames         If set, then call useopac with these names. Else use
;                    the opacity that is currently in the files 
;                    frequency.inp and dustopac_1.inp and dustopac.inp
;  pllongs           [only if opacnames; see routine useopac()]
;  fresmd            [only if opacnames; see routine useopac()]
;  scat              [only if opacnames; see routine useopac()]
;  nvstr             >0 --> Iterate on the vertical structure nvstr times
;  ivstrt            Dust species that is representative for T_gas (default=1)
;  dostr             [optional] This allows you to switch on/off the
;                    participation of each dust species to the vertical 
;                    structure calculation. Must be array of size nspec,
;                    containing 0 or 1 for each species.
;  vserrtol          The convergence criterion for the vertical structure
;                    iteration. If 0.d0, then it will always iterate 
;                    precisely nvstr times.
;  run               =1 --> Run the RADMC code
;                    =0 --> Write a dorun script to do this
;  bindir            If the executables are not located in ../bin, then 
;                    they are located in bindir.
;  ifast             [Careful!!] When this is 1 or 2 then the calculations
;                    of the dust temperature are not done upon each
;                    interaction of the photon with a grid cell, but only if
;                    the energy in the cell is appreciable larger than the
;                    last update. NOTE: This is not related to the new
;                    improved RADMC version, which is anyway much faster
;                    than the old one. Perhaps the ifast is even no longer
;                    necessary, but that remains to be seen.
;
; RETURNS:
;  hrstore           The Hp(R)/R array from the simple setup
;  tt                The tau array, for plotting if required.
;  r                 The radial grid (NOTE: Is already saved to disk!)
;  theta             The theta grid (NOTE: Is already saved to disk!)
;  rhodust           The dust density (NOTE: Is already saved to disk!)
;  sigdust           The dust surface density  (NOTE: Is already saved to disk!)
;
; NOTE:
;  - In this routine the 'ntt' variable is the TOTAL nr of theta points 
;    in the upper quadrant (old style: n_elements(theta) = ntt + ntex).
;
;------------------------------------------------------------------------
pro simpledisk_vertstruct,rstar=rstar,tstar=tstar,mstar=mstar,$
         ifinstar=ifinstar,fresmd=fresmd,$
         scat=scat,ntex=ntex,nrr=nrr,ntt=ntt,rin=rin,tin=tin,$
         rout=rout,hrgrid=hrgrid,hrgmax=hrgmax,rrefine=rrefine,$
         drsm=drsm,rdisk=rdisk,sigdust0=sigdust0,mdisk=mdisk,$
         plsig1=plsig1,plsig2=plsig2,$
         opacnames=opacnames,pllongs=pllongs,$
         schmidt=schmidt,ab_r0=ab_r0,ab_ab0=ab_ab0,ab_pl=ab_pl,$
         gastodust=gastodust,ab_min=ab_min,hrdisk=hrdisk,hrmin=hrmin,$
         plh=plh,rpfrin=rpfrin,hrpuff=hrpuff,nvstr=nvstr,$
         ivstrt=ivstrt,vserrtol=vserrtol,nphot=nphot,$
         npdiff=npdiff,errtol=errtol,tauchop=tauchop,lmbchop=lmbchop,$
         idxchop=idxchop,rhofloor=rhofloor,pnc=pnc,pnh=pnh,pz=pz,$
         imakedisk=imakedisk,run=run,hrstore=hrstore,$
         thintin=thintin,tt=tt,radius=r,theta=theta,rhodust=rhodust,$
         sigdust=sigdust,kurucz=kurucz,kurdir=kurdir,bindir=bindir,$
         ifast=ifast,dostr=dostr,csenv=csenv,time=time,env=env,gap=gap,gdrop=gdrop
;
@problem_natconst.pro
;
; Defaults
;
if not keyword_set(rstar) then rstar=RS
if not keyword_set(tstar) then tstar=TS
if not keyword_set(mstar) then mstar=MS
if n_elements(ifinstar) eq 0 then ifinstar=1
if n_elements(fresmd) eq 0 then fresmd=3
if n_elements(scat) eq 0 then scat=1
if not keyword_set(nrr) then nrr=130
if not keyword_set(ntt) then ntt=60
nr = nrr
nt = ntt
if not keyword_set(ntex) then ntex=10
if keyword_set(rin) and keyword_set(tin) then begin
    print,'ERROR: Cannot set both tin and rin simultaneously'
    stop
endif
if not keyword_set(rin) and not keyword_set(tin) then tin=1500.
if not keyword_set(rout) then rout = 100*AU
if not keyword_set(hrgrid) then hrgrid = 0.6
if not keyword_set(hrgmax) then hrgmax = 0.9*!pi/2.d0
if not keyword_set(drsm) then drsm=0.d0
if not keyword_set(rdisk) then rdisk=100*AU
if keyword_set(sigdust0) and keyword_set(mdisk) then begin
    print,'ERROR: Either specify sigdust0 or mdisk. Not both.'
    stop
endif
if not keyword_set(sigdust0) and not keyword_set(mdisk) then mdisk=0.01*mstar
if not keyword_set(plsig1) then plsig1=-1.d0
if not keyword_set(plsig2) then plsig2=-12.d0
if not keyword_set(schmidt) then schmidt=1.d0
if not keyword_set(ab_r0) then ab_r0=0
if not keyword_set(ab_ab0) then ab_ab0=0
if not keyword_set(ab_pl) then ab_pl=0
if not keyword_set(ab_min) then ab_min=0
if keyword_set(ab_r0) then begin
    if n_elements(ab_r0) ne n_elements(ab_ab0) then begin
        print,'ERROR: Nr of elements of ab_* not mutually consistent'
        stop
    endif
    if n_elements(ab_r0) ne n_elements(ab_pl) then begin
        print,'ERROR: Nr of elements of ab_* not mutually consistent'
        stop
    endif
    if n_elements(ab_r0) ne n_elements(ab_min) then begin
        print,'ERROR: Nr of elements of ab_* not mutually consistent'
        stop
    endif
endif
if not keyword_set(hrdisk) then begin
    print,'WARNING: hrdisk is not specified, so I put it arbitrarily at 0.1...'
    hrdisk = 0.1
endif
if not keyword_set(hrmin) then hrmin=0.02
if not keyword_set(plh) then plh=0.1
if not keyword_set(hrpuff) then hrpuff=0.d0
if not keyword_set(nvstr) then nvstr=0
if not keyword_set(nphot) then nphot=100000L
if n_elements(npdiff) eq 0 then npdiff=30
if not keyword_set(errtol) then errtol=1d-10
if not keyword_set(tauchop) then tauchop=0.d0
if not keyword_set(lmbchop) then lmbchop=0.55
if not keyword_set(idxchop) then tdxchop=1.d0
if n_elements(imakedisk) eq 0 then begin
    print,'WARNING: imakedisk is not specified. As default I will do'
    print,'         the full vertical structure and RT (imakedisk=1).'
    imakedisk = 1
endif
if not keyword_set(opacnames) then begin
    print,'ERROR: No opacities specified! Set opacnames[].'
    stop
endif
if keyword_set(ab_r0) then begin
    if n_elements(opacnames) ne n_elements(ab_ab0)+1 then begin
        print,'ERROR: Nr of opacities specified does not match nr of abundances'
        print,'       of additional species + 1.'
        stop
    endif
endif
if n_elements(opacnames) ne n_elements(pllongs) and $ 
     keyword_set(pllongs) then begin
    print,'ERROR: nr of elements in opacnames must equal that of pllongs'
    stop
endif
if not keyword_set(pllongs) then pllongs=-2.d0+dblarr(n_elements(opacnames))
if n_elements(gastodust) eq 0 then gastodust = 100.
if n_elements(rhofloor) eq 0 then rhofloor = 1d-26
if not keyword_set(bindir) then bindir='../bin/'
if nvstr gt 0 and tauchop gt 0.d0 then begin
   print,'ERROR: The vertical structure iteration module is'
   print,'       at present not compatible with vertical structure'
   print,'       iteration. Either put nvstr or tauchop to 0.'
   stop
endif
;;
;;--------------------------------------------------------------------
;;
;; ...Grid refinement near inner edge
;;
if not keyword_set(rrefine) then begin
    rrefine={nlevr:3,nspanr:2,nstepr:2}
endif
;;
;; Install the opacities
;;
useopac,opacnames,pllongs,fresmd=fresmd,scat=scat,$
        nf=nf,pnc=pnc,pnh=pnh,pz=pz,$
        nspec=nspec,npah=npah,ntherm=ntherm
;
; Write the star spectrum and info
;
write_star,mstar,rstar,tstar,kurucz=kurucz,kurdir=kurdir,$
         inupk=inupk,lampk=lampk
;
; Find the opacity at the peak of the stellar spectrum
;
kappa = dblarr(nspec)
for i=0,nspec-1 do begin
    o = readopac(nr=strcompress(string(i+1),/remove_all))
    kappa[i] = o.cabs(inupk) + o.csca(inupk)
endfor
print,"lambda_peak = ",lampk," micron,  kappa = ",kappa
;
; Estimate rin from tin using a simple blackbody wall thing
;
; NOTE: In contrast to the original D&D models we now do not iterate
;       on the self-irradiation.
;
; NOTE: We now have a new option: pre-compute the optically thin
;       temp, and from there try to find the rin
;
if not keyword_set(thintin) then begin
    ;;
    ;; Use the D&D04 formula, but this time without self-irradiation
    ;;
    if keyword_set(tin) then begin
        rin  = rstar*(tstar/tin)^2
    endif
endif else begin
    ;;
    ;; Use the optically thin dust temperature to estimate the 
    ;; radius rin belonging to tin.
    ;;
    print,'PROBLEM: The optically thin real temperature mode for tin'
    print,'         is still under construction....'
    stop
endelse
;;
;; Compute the puffing-up radius
;;
if keyword_set(rpfrin) then begin
    if(rpfrin lt 1.d0) then begin
        print,'ERROR: rpfrin should not be smaller than 1 '
        print,'       (unless you put it to 0, i.e. deactivate puffing)'
        stop
    endif
    rpuff = rpfrin * rin
endif else begin
    rpuff = 0.d0
endelse
;
; NOW DECIDE WHETHER TO IMPORT A BASIC SINGLE-DUST-SPECIES STRUCTURE AS
; A STARTING POINT OR MAKE AN AD-HOC MULTI-DUST-SPECIES STRUCTURE WHICH
; CAN OPTIONALLY BE ITERATED IN VERTICAL STRUCTURE.
; 
case imakedisk of 
    -1: begin
        ;;
        ;;---------------------------------------------------------------
        ;;     USE PREVIOUSLY GENERATED DISK STRUCTURE
        ;;  (dustdens_struct.inp, radius_struct.inp, theta_struct.inp) 
        ;;---------------------------------------------------------------
        ;;
        ;; Read the density structure
        ;;
        openr,1,'radius_struct.inp'
        nr=0
        readf,1,nr
        r=dblarr(nr)
        readf,1,r
        close,1
        openr,1,'theta_struct.inp'
        nt=0
        idum=0
        readf,1,nt,idum
        theta=dblarr(nt)
        readf,1,theta
        close,1
        openr,1,'dustdens_struct.inp'
        nnr=0
        nnt=0
        nsp=0
        idum=0
        readf,1,nsp,nnr,nnt,idum
        if nsp ne 1 then begin
            print,'ERROR: Input density structure must have only 1 dust species'
            stop
        endif
        if nnr ne nr or nnt ne nt then begin
            print,'ERROR: The radial and/or theta grid size are not equal'
            print,'       to that of the dustdens_radical.inp'
            stop
        endif
        rhodusttot=dblarr(nnt,nnr)
        readf,1,rhodusttot
        rhodusttot = transpose(rhodusttot)
        close,1
        ;;
        ;; Now split this dust density up in the various abundances
        ;;
        ;;
        ;; [[ HERE ONE CAN DO THING DIFFERENTLY IF ONE WISHES DIFFERENT 
        ;;    ABUNDANCE PROFILES ]]
        ;; 
        
        if nspec gt 1 then begin
           abun = dblarr(nnr,nnt,nspec) + 1.d0
           for ispec=1,nspec-1 do begin
              q  = dblarr(nnr) + ab_ab0[ispec-1]
              ii = where(r gt ab_r0[ispec-1])
              q[ii] = ab_ab0[ispec-1]*(r[ii]/ab_r0[ispec-1])^ab_pl[ispec-1]
              q = q > ab_min[ispec-1]
              abun[*,*,ispec] = rebin(q,nnr,nnt)
              abun[*,*,0] = abun[*,*,0] - abun[*,*,ispec]
           endfor
           if min(abun[*,*,0]) lt 0.d0 then begin
              print,'ERROR: Total abundances dont add up!'
              stop
           endif
        endif
        ;;
        ;; Now convert this into dustdens
        ;;
        rhodust = dblarr(nnr,nnt,nspec)
        for ispec=0,nspec-1 do begin
            rhodust[*,*,ispec] = rhodusttot * abun[*,*,ispec]
        endfor
        ;;
        ;; Write the stuff to the input files
        ;;
        write_radmc_density,r,theta,rhodust
        ;;
    end
    1: begin
        ;;
        ;;---------------------------------------------------------------
        ;;     MAKE AN AD-HOC DISK STRUCTURE FROM THE PARAMETERS
        ;;       WITH OR WITHOUT VERTICAL STRUCTURE ITERATION
        ;;---------------------------------------------------------------
        ;;
        ;; Make the R-grid
        ;;
        r      = make_rgrid(rin,rout,nr,rrefine=rrefine)
        ;;
        ;; Make the Theta-grid
        ;;
        theta  = make_tgrid(hrgrid,nt,hrgmax=hrgmax,ntex=ntex,$
                            hrlg=hrlg,zrefine=zrefine)
        ;;
        ;; Just for safety, get nnr and nnt
        ;;
        nnr    = n_elements(r)
        nnt    = n_elements(theta)
        ;;
        ;; Now make a disk model
        ;;
        if keyword_set(sigdust0) then begin
            ;;
            ;; Directly from sigdust0
            ;;
            rhodusttot = disk_model_1(r,theta,rdisk,sigdust0,plsig1,plsig2,$
                          hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                          hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdust,gap=gap,gdrop=gdrop)
        endif else begin
            ;;
            ;; Compute sigdust0 from mdisk
            ;;
            sigdust00  = 1.d0
            rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                            hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                            hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot,gap=gap,gdrop=gdrop)
            mddum      = integrate(r,2*pi*r*sigdusttot)*gastodust
            sigdust00  = mdisk / mddum
            rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                         hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                         hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot,gap=gap,gdrop=gdrop)
        endelse



        ;;
        ;; Compute the mass
        ;;
        md = mass(rhodusttot,r,theta)
        print,'Mdisk = ',md/MS,' Msun (using gas-to-dust=100)'

        ;; Now add an envelope

        IF KEYWORD_SET(env) THEN BEGIN
           IF n_elements(time) NE 1 OR n_elements(csenv) NE 1 OR $
              n_elements(csenv) NE 1 then begin
              print,'ERROR: If you include the envelope, you must also define'
              print,'   time and csenv...'
              stop
           endif
           q    = make_shu(r,time*3600.*24.*365.,csenv=csenv)

           FOR ith = 0,nt-1 DO BEGIN
               ;; Mutually exclude the densities
              ;;
              gsubs = WHERE(r LT rdisk) 
              IF gsubs[0] NE -1 THEN q.rho[gsubs] = 0.d0
              gsubs = WHERE( rhodusttot[*,ith] LT q.rho/gastodust) 
              rhodusttot[gsubs,ith] = q.rho[gsubs] / gastodust
            ENDFOR
           PRINT, 'Included a Shu envelope!'
        ENDIF
;  


        ;;
        ;; Now split this dust density up in the various abundances
        ;;
        if nspec gt 1 then begin
           freeze_out=1
           IF freeze_out THEN BEGIN
              abun = dblarr(nnr,nnt,nspec) + 1.d0
              abun[WHERE(r GT 2.*AU),*,0] = 0.
              abun[WHERE(r LE 2.*AU),*,1] = 0.
           ENDIF ELSE BEGIN
              
              abun = dblarr(nnr,nnt,nspec) + 1.d0
              for ispec=1,nspec-1 do begin
                 q  = dblarr(nnr) + ab_ab0[ispec-1]
                 ii = where(r gt ab_r0[ispec-1])
                 q[ii] = ab_ab0[ispec-1]*(r[ii]/ab_r0[ispec-1])^ab_pl[ispec-1]
                 q = q > ab_min[ispec-1]
                 abun[*,*,ispec] = rebin(q,nnr,nnt)
                 abun[*,*,0] = abun[*,*,0] - abun[*,*,ispec]
              endfor
              if min(abun[*,*,0]) lt 0.d0 then begin
                 print,'ERROR: Total abundances dont add up!'
                 stop
              endif
           ENDELSE
        endif
        ;;
        ;; Now convert this into dustdens
        ;;
        if keyword_set(abun) then begin
            rhodust = dblarr(nr,nt,nspec)
            for ispec=0,nspec-1 do begin
                rhodust[*,*,ispec] = rhodusttot * abun[*,*,ispec]
            endfor
        endif else begin
            rhodust = rhodusttot
        endelse
        ;;
        ;; Now smooth the inner rim if necessary
        ;;
        if keyword_set(drsm) then begin
            smooth_rim,r,theta,rhodust,kappa,sigdust=sigdust,$
                    drsm=drsm,tautol=tautol
        endif
        ;;
        ;; Write the dust density to the input files of RADMC
        ;;
        write_radmc_density,r,theta,rhodust
        ;;
    end
    else: begin
        print,'ERROR: Must set imakedisk to -1 or 1'
        stop
    end
endcase
;;
;;
;; Make optical depth array
;;
b  = {r:r,theta:theta,rho:rhodust}
tt =maketau(b,kappa)
;;
;; Write a message
;;
ddr=r(nr-1)/r(nr-2)-1.d0
print,'TAUR  = ',tt.taur(nr-1,nt-1)
print,'DR/R  = ',ddr
;;
;; Check the optical depth at the inner edge
;;
tr0 = tt.taur(1,nt-1)
if tr0 gt 1.d0 then begin
    print,'Optical depth error! Inner grid cell at equator optically thick'
    print,' Is this okay? (type 1)'
    read,i
    if i ne 1 then stop
endif
;;
;; Remove radical.inp  (obsolete)
;;
spawn,'\mv -f radical.inp radical.inp_obsolete >& /dev/null'
;;
;; Write the raytrace.inp for the ray-tracing (spectra and images)
;;
write_raytrace,nextra=nextra,nrref=nrref
;;
;; Write the radmc.inp file
;;
;if keyword_set(ifinstar) then begin
;    print,'Note: New mode active: Finite size of the star'
;endif
write_radmc,nphot=nphot,ifinstar=ifinstar,npdiff=npdiff,$
            nvstr=nvstr,ivstrt=ivstrt,vserrtol=vserrtol,ifast=ifast
;;
;; Make a chopdens file and immediately do the chopping!!
;;
if tauchop ne 0 then begin
    if lmbchop eq 0 or idxchop eq 0 then begin
        print,'ERROR: When tauchop, also define lmbchop and idxchop'
        stop
    endif
    write_chopdens,tauchop,lmbchop,idxchop
    spawn,'nice '+bindir+'/chopdens'
endif
;;
;; Make the vstruct.inp file
;;
if keyword_set(nvstr) and n_elements(dostr) eq nspec then begin
   write_vstruct,dostr
endif else begin
   spawn,'rm -f vstruct.inp'
endelse
;
; Verify
;
if nr ne n_elements(r) then stop
if nt ne n_elements(theta) then stop
;
; If we want to run this, then do it
;
if keyword_set(run) then begin
   if not keyword_set(bindir) then begin
      print,'ERROR: When /run, then also specify bindir'
      stop
   endif
   spawn,'nice '+bindir+'/radmc'
   spawn,'nice raytrace spectrum'
endif
;
end

