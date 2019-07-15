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
;*NOENV*
;@problem_shuulrich.pro
@problem_makeshu.pro
@problem_grid.pro
@problem_radmc_files.pro
@problem_disk.pro

PRO problem_models
   
END
;------------------------------------------------------------------------
;           DO THE RADIATIVE TRANSFER FOR DISK + ENVELOPE
;
; This subroutine sets up the 2-D/3-D density distribution for the disk +
; envelope around a young star, and then subsequently calls RADMC to do
; the Monte Carlo simulation. It can also do the vertical structure
; iteration using the updated RADMC with diffusion and vertical structure.
;
; The envelope is computed internally using the Shu
; model coupled to the Ulrich rotating infall model. The disk model 
; has to be computed externally (for instance from a viscous accretion
; disk model) and the resulting Sigma(R) inserted here in this routine.
; Of course, also a simple powerlaw disk model can be inserted, if
; desired.
;
; ARGUMENTS:
;  mstar             Mass of the star [g]
;  rstar             Radius of the star [cm]
;  tstar             Effective temperature of the star [K]
;  kurucz            If set, then use Kurucz model
;  kurdir           The directory where the Kurucz models are located
;  rdisk             Array of radial grid for the input disk model [cm]
;  sigmagasdisk      Array of Sigma_gas(R) of the disk [g/cm^2]
;  hpdisk            Array of pressure scale heigh H_p(R) of disk [cm]
;
; KEYWORDS:
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
;  chop              Upper limit to vertical optical depth at 0.55 micron
;                    (the chopping routine simply removes mass from the 
;                    equator until tau<=chop. this is necessary in order 
;                    not to stall the Monte Carlo simulation)
;  gastodust         Gas-to-dust mass ratio. Default 100.
;  rhofloor          The density floor value
;  tdec              Enable thermal decoupling of temperatures?
;  nib               [Technical] Number of rays inward of Rin for SED.
;                    Normally put this to 4 or so. But sometimes this
;                    causes problems when the star radius is not too much
;                    smaller than Rin. Then put it to 0.
;  nphot             Nr of photon packages
;  ifinstar          =1 -> uses finite size star (iso point source)
;                    (ADDED 09.01.07)
;  zrefine           A structure for making a very refined grid in 
;                    theta near the midplane. See problem_grid.pro the
;                    routine make_tgrid(). Note that if this option
;                    is used, the nr of theta grid points will exceed
;                    the value ntt+ntex.
;  opacname          If set, then call useopac with this name. Else use
;                    the opacity that is currently in the files 
;                    frequency.inp and dustopac_1.inp and dustopac.inp
;  pllongs           [only if opacname; see routine useopac()]
;  fresmd            [only if opacname; see routine useopac()]
;  scat              [only if opacname; see routine useopac()]
;  nvstr             >0 --> Iterate on the vertical structure nvstr times
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
; ARGUMENTS FOR THE SHU/ULRICH INFALL MODEL:
;  time              Time after onset of collapse of the cloud core [s]
;  csenv             Sound speed in the cloud core (see Shu model!)
;
; RETURNS:
;
; NOTE:
;  - In this routine, for backward compatilibity reasons, the ntt is
;    the basic nr of theta points. The total nr is nt=ntt+ntex.
;
;------------------------------------------------------------------------
; BUGFIX 08.12.06: I accidently RETURN the new values of nr and nt,
;   so that they become different upon return than they were when the
;   diskenv_radmc() was called. This can be extremely dangerous when 
;   this routine is called multiple times!!!!
;   This is now fixed!!
;------------------------------------------------------------------------
pro diskenv_radmc, rstar=rstar,tstar=tstar,mstar=mstar,ifinstar=ifinstar,$
			       fresmd=fresmd,scat=scat,ntex=ntex,nrr=nrr,ntt=ntt,rin=rin,tin=tin,$
			       out=rout,hrgrid=hrgrid,hrgmax=hrgmax,rrefine=rrefine,drsm=drsm,$
			       rdisk=rdisk,sigdust0=sig0,mdisk=mdisk,$
	               plsig1=plsig1,plsig2=plsig2,kurucz=kurucz,kurdir=kurdir,$
        	       opacnames=opacnames,pllongs=pllongs,$
			       gastodust=gastodust,hrdisk=hrdisk,hrmin=hrmin,$
			       plh=plh,rpfrin=rpfrin,hrpuff=hrpuff,nvstr=nvstr,$
	               nphot=nphot,npdiff=npdiff,errtol=errtol,tauchop=tauchop,lmbchop=lmbchop,$
			       idxchop=idxchop,rhofloor=rhofloor,pnc=pnc,pnh=pnh,pz=pz,ref2=ref2,bindir=bindir,$
				   tt=tt,radius=r,theta=theta,csenv=csenv,time=time,env=env,cav=cav,opening=opening, Aenv=Aenv
				   
				   
				          
@problem_natconst.pro
;
; Set default parameters
;
if n_elements(gastodust) eq 0 then gastodust = 100.
if n_elements(ntt) eq 0 then ntt=50
if n_elements(nrr) eq 0 then nrr=90
if n_elements(ntex) eq 0 then ntex=10
if n_elements(rhofloor) eq 0 then rhofloor = 1d-90
if n_elements(hrgrid) eq 0 then hrgrid = 0.6  
if n_elements(hrgmax) eq 0 then hrgmax = 0.9*!pi/2.d0
if n_elements(nib) eq 0 then nib = 4

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
; Install the opacities
;
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
; Find the outer radius of the disk
;
if keyword_set(tin) then begin
    rin  = rstar*(tstar/tin)^2
endif

;
; Make the R-grid
;

r      = make_rgrid(rin,rout,nr,rrefine=rrefine)

;
; If requested, then make another radius of extra R-refinement.
; This is useful if you have a "second wall" somewhere, or 
; equivalently the outer part of an annular gap or something like
; this. 
;
if keyword_set(ref2) then begin
   refine,r,ref2.irstart,nlevr=ref2.nlevr,nspanr=ref2.nspanr,nstepr=ref2.nstepr
endif
;
; Make the Theta-grid
;
theta  = make_tgrid(hrgrid,nt,hrgmax=hrgmax,ntex=ntex,$
                    hrlg=hrlg,zrefine=zrefine)


if keyword_set(sigdust0) then begin
     ;;
     ;; Directly from sigdust0
	 ;;
     rhodusttot = disk_model_1(r,theta,rdisk,sigdust0,plsig1,plsig2,$
                   hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                   hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdust)
endif else begin
     ;;
     ;; Compute sigdust0 from mdisk
     ;;
     sigdust00  = 1.d0
     rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                     hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                     hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot)
     mddum      = integrate(r,2*pi*r*sigdusttot)*gastodust
     sigdust00  = mdisk / mddum
     rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                  hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                  hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot)
endelse

;
;Make the Shu infalling envelope
;
if keyword_set(env) then begin
	if ~keyword_set(cav_dens) then cav_dens=0.
   envelope = make_shu(r,time*3600.*24.*365.,csenv=csenv, Aenv=Aenv)
	env_subs = WHERE(r GT rdisk)
	for i=0,nt-1 do begin
		envelope_at_theta = envelope.rho
		;
		;Carve out a cavity, if requested
		;
		;In polar coordinates:
		;r = 2a/(1+cos(theta))
		;
		if keyword_set(cav) then begin
			cav_subs = WHERE(r LT 2.*opening/(1.+cos(!pi-theta[i])))
			envelope_at_theta[cav_subs] = cav_dens
		endif	

		rhodusttot[env_subs,i] += envelope_at_theta[env_subs]/gastodust
	endfor
endif

;
; Compute the mass
;
md = mass(rhodusttot,r,theta)
print,'Mdisk = ',md/MS*(gastodust/100),' Msun'
;
; Now convert this into dustdens
;
if keyword_set(abun) then begin
    rhodust = dblarr(nnr,nt,nspec)
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

;; Add the density background
rhodust += rhofloor

;
; Write the dust density to the input files of RADMC
;
write_radmc_density,r,theta,rhodust

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
;
; Make optical depth array
;
b  = {r:r,theta:theta,rho:rhodust}
tt =maketau(b,kappa)
;
; Check the optical depth at the inner edge
;
tr0 = tt.taur(1,nt-1)
if tr0 gt 1.d0 then begin
    print,'Optical depth error! Inner grid cell at equator optically thick'
    print,' Is this okay? (type 1)'
    read,i
    if i ne 1 then stop
endif

;
; Write the raytrace.inp for the ray-tracing (spectra and images)
;
write_raytrace,nextra=nextra,nrref=nrref
;
; Write the radmc.inp file
;
write_radmc,nphot=nphot,ifinstar=ifinstar,npdiff=npdiff,$
            nvstr=nvstr,ivstrt=ivstrt,vserrtol=vserrtol,ifast=ifast

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
         ifast=ifast,dostr=dostr,ref2=ref2,snowline=snowline
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
if n_elements(rhofloor) eq 0 then rhofloor = 1d-90
if not keyword_set(snowline) then snowline=0
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

        rhodust += rhofloor
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
        ;; If requested, then make another radius of extra R-refinement.
        ;; This is useful if you have a "second wall" somewhere, or 
        ;; equivalently the outer part of an annular gap or something like
        ;; this. 
        ;;
        if keyword_set(ref2) then begin
           refine,r,ref2.irstart,nlevr=ref2.nlevr,nspanr=ref2.nspanr,nstepr=ref2.nstepr
        endif
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
                          hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdust)
        endif else begin
            ;;
            ;; Compute sigdust0 from mdisk
            ;;
            sigdust00  = 1.d0
            rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                            hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                            hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot)
            mddum      = integrate(r,2*pi*r*sigdusttot)*gastodust
            sigdust00  = mdisk / mddum
            rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                         hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                         hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot)
			if keyword_set(ref2) then begin	
				rhodusttot[0:ref2.irstart,*] = rhodusttot[0:ref2.irstart,*]/ref2.deplete_inner
			endif
        endelse
        ;;
        ;; Compute the mass
        ;;
        md = mass(rhodusttot,r,theta)
        print,'Mdisk = ',md/MS,' Msun (using gas-to-dust=100)'
        ;;
        ;; Now split this dust density up in the various abundances
        ;;
        if (nspec gt 1) and not snowline then begin
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
        
        if (nspec eq 2) and snowline then begin
            abun = dblarr(nnr,nnt,nspec) + 1.d0
            ice = r gt snowline*AU
            plane = rebin(ice,nnr,nnt)
            abun[*,*,0] = plane
            noice = r le snowline*AU
            plane = rebin(noice,nnr,nnt)
            abun[*,*,1] = plane
        endif
        ;;
        ;; Now convert this into dustdens
        ;;
        if keyword_set(abun) then begin
            rhodust = dblarr(nnr,nt,nspec)
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

        ;; Add the density background
        rhodust += rhofloor

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
ddr=r(nnr-1)/r(nnr-2)-1.d0
print,'TAUR  = ',tt.taur(nnr-1,nt-1)
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
;if nr ne n_elements(r) then stop
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

