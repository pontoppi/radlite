;==========================================================================
;               MAIN SETUP ROUTINE FOR THE 2-D DISK MODELS
;          EXECUTE THIS SCRIPT IN ORDER TO SET UP A SIMULATION
;                 (PARAMETERS ARE IN PROBLEM_PARAMS.PRO)
;==========================================================================
@problem_models.pro
@problem_mixopacities.pro
;
PRO problem_setup, show=show
@natconst.pro
;
; Checks
;
if n_elements(itstr) ne 0 then begin
   print,'ERROR: In this version itstr ---> nvstr (i.e. equal to'
   print,'       input variable of RADMC).'
   stop
endif
if n_elements(dnphot) ne 0 then begin
   print,'ERROR: In this version dnphot ---> npdiff (i.e. similar to'
   print,'       input variable of RADMC: nphotdiff).'
   stop
endif
;
; Defaults
;
mugas  = 2.3
close,/all
if n_elements(show) eq 0 then begin
   show   = 1           ; If 1, then show plots of setup
endif
;
; Get date, hostname and current directory
;
spawn,'echo `hostname`',hostname
hostname=strcompress(hostname,/remove_all)
hostname=hostname[0]
cd,current=directory
date=systime()
;
; Open and clear the log file
;
openw,1,'diskmodel.log'
printf,1,'Setup of 2-D Disk Model Run:'
printf,1,'    Date      : '+date
printf,1,'    Hostname  : '+hostname
printf,1,'    Directory : '+directory
printf,1,'-------------------------------------------------'
close,1
;
; Read the parameters
;
@problem_params.pro
;
; Checks
;
if rdisk gt rout then begin
   print,'PROBLEM: Outer disk radius is larger than the outer grid radius!'
   stop
endif
;
;
; kuruczdir conflicts with the kurucz keyword because of IDL disambiguation
if KEYWORD_SET(kuruczdir) then begin
	kurdir = kuruczdir
endif
;
;
; Mix dust opacities
;
if keyword_set(mixspecs) then begin
   sz = size(mixspecs)
   if sz[0] eq 1 then nmix=1 else nmix=sz[2]
   for imix=0,nmix-1 do mixopacities,mixspecs[*,imix],$
                   mixnames[imix],mixabun[*,imix]
endif
;
; Call the model setup
;
IF keyword_set(csenv) THEN BEGIN
	diskenv_radmc,rstar=rstar,tstar=tstar,mstar=mstar,ifinstar=ifinstar,$
       fresmd=fresmd,scat=scat,ntex=ntex,nrr=nr,ntt=nt,rin=rin,tin=tin,$
       out=rout,hrgrid=hrgrid,hrgmax=hrgmax,rrefine=rrefine,drsm=drsm,$
       rdisk=rdisk,sigdust0=sig0,mdisk=mdisk,$
       plsig1=plsig1,plsig2=plsig2,kurucz=kurucz,kurdir=kurdir,$
       opacnames=infile,pllongs=pll,$
       gastodust=gastodust,hrdisk=hrdisk,hrmin=hrmin,$
       plh=plh,rpfrin=rpfrin,hrpuff=hrpuff,nvstr=nvstr,$
       nphot=nphot,npdiff=npdiff,errtol=errtol,tauchop=tauchop,lmbchop=lmbchop,$
       idxchop=idxchop,rhofloor=rhofloor,bindir=bindir,tt=tt,radius=r,$
	   theta=theta,csenv=csenv,time=time,env=env,cav=cav,opening=opening, Aenv=Aenv
ENDIF ELSE BEGIN
	simpledisk_vertstruct,rstar=rstar,tstar=tstar,mstar=mstar,ifinstar=ifinstar,$
       fresmd=fresmd,scat=scat,ntex=ntex,nrr=nr,ntt=nt,rin=rin,tin=tin,$
       rout=rout,hrgrid=hrgrid,hrgmax=hrgmax,rrefine=rrefine,drsm=drsm,$
       rdisk=rdisk,sigdust0=sig0,mdisk=mdisk,$
       plsig1=plsig1,plsig2=plsig2,kurucz=kurucz,kurdir=kurdir,$
       opacnames=infile,pllongs=pll,$
       schmidt=schmidt,ab_r0=ab_r0,ab_ab0=ab_ab0,ab_pl=ab_pl,$
       gastodust=gastodust,ab_min=ab_min,hrdisk=hrdisk,hrmin=hrmin,$
       plh=plh,rpfrin=rpfrin,hrpuff=hrpuff,nvstr=nvstr,$
       ivstrt=ivstrt,vserrtol=vserrt,nphot=nphot,$
       npdiff=npdiff,errtol=errtol,tauchop=tauchop,lmbchop=lmbchop,$
       idxchop=idxchop,rhofloor=rhofloor,pnc=pnc,pnh=pnh,pz=pz,$
       imakedisk=imakedisk,run=run,hrstore=hrstore,$
       thintin=thintin,tt=tt,radius=r,theta=theta,bindir=bindir,$
       dostr=dostr,ref2=ref2,snowline=snowline
ENDELSE
;
;Added the accretion description of Nikoletta Sipos. If there is time
;as some point, this should be cleaned up and implemented as a proper subroutine
;
if N_ELEMENTS(do_acc) gt 0 then begin
   if do_acc then begin
      @problem_accretion.pro
   endif
endif
;
; Write a data file with the names simser,simnr,simweb
;
openw,1,'diskmodel.names'
printf,1,simser
printf,1,simnr
printf,1,simweb
close,1

;
; Make plots and write info
;
if show eq 1 then begin
   pos = [0.1,0.55,0.95,0.95]
   xtitle = 'R [AU]'
   ytitle = '!4p!X/2-!4H!X = H/R'
   nnt = n_elements(theta)
   window,1,xsize=800,ysize=700
   !p.multi(2)=2
;   surface,tt.taur gt 1,r,theta,/xl,ax=50,/xs,/ys
   surface,tt.taur gt 1,r/au,!pi/2-theta,/xl,ax=90,az=0,/xs,/ys,position=pos,xtitle=xtitle,ytitle=ytitle
   contour,tt.taur gt 1,r/au,!pi/2-theta,/xl,/xs,/ys,position=pos,/noerase,/fill,$
          title='Grid and the !4s!X!DV!N>1 region'
;   surface,tt.taur gt 1,r/au,!pi/2-theta,xtitle=xtitle,ytitle=ytitle,/xl,ax=90,az=0,/xs,/ys,position=pos,/noerase
   plot,tt.taur(0:20,nnt-1),yrange=[0,6],psym=6,yminor=1,yticks=6,$
            xtitle='Radial grid points from inner edge',$ 
            ytitle='!4s!X!DV!N at midplane'
   oplot,tt.taur(0:20,nnt-1)
   !p.multi(2)=1
endif

END
