;
;Procedure for setting up a velocity grid.
;This is where additional velocity structures can be constructed.
;
;Currently included:
;
;1) Keplerian velocity
;
PRO make_velocity,a,vtype,vr=vr,vtheta=vtheta,vphi=vphi,add_dens=add_dens
@natconst.pro
@problem_params.pro
@line_params.ini

IF NOT KEYWORD_SET(vtype) THEN BEGIN
    PRINT, 'No velocity structure selected!' 
    stop
ENDIF

nr = n_elements(a.r)
nt = n_elements(a.theta)/2

CASE vtype OF
    1: BEGIN    ;Keplerian velocity structure
        openr,1,'starinfo.inp'
        iformat=0
        readf,1,iformat
        rstar=0.d0
        mstar=0.d0
        readf,1,rstar
        readf,1,mstar
        close,1
        IF N_ELEMENTS(VERBOSE) THEN print, 'Used starinfo.inp for Keplerian velocity'
        vphi = sqrt(GG*mstar/a.rr)
        vr   = fltarr(nr,nt)
        vth  = fltarr(nr,nt)
        IF KEYWORD_SET(lefthand) THEN vphi = -vphi
    END
    2: BEGIN    ;Keplerian velocity structure for the disk and infall for the 
                ;envelope
        openr,1,'starinfo.inp'
        iformat=0
        readf,1,iformat
        rstar=0.d0
        mstar=0.d0
        readf,1,rstar
        readf,1,mstar
        close,1
        IF N_ELEMENTS(VERBOSE) THEN print, 'Used starinfo.inp for Keplerian velocity'
        vphi = sqrt(GG*mstar/a.rr);*(a.rr le 50*1.5d13)
        IF KEYWORD_SET(env) THEN BEGIN
           IF n_elements(time) NE 1 OR n_elements(csenv) NE 1 OR $
              n_elements(csenv) NE 1 then begin
              print,'ERROR: If you include the envelope, you must also define'
              print,'   time and csenv...'
              stop
           endif
           q = make_shu(a.rr,time*3600.*24.*365.,csenv=csenv, Aenv=Aenv)
        ENDIF ELSE BEGIN
           PRINT, 'ERROR: If you want to have an infalling envelope, you must define an envelope in problem_params.pro!'
           stop
        ENDELSE
		
        vr   = q.v
        dsubs = WHERE(a.rr LE rdisk)
		  esubs = WHERE(a.rr GT rdisk)

		;No infall within the disk.
        vr[dsubs] = 0.d0
		;Uncomment this to have no rotation within the envelope.
	     ;vphi[esubs] = 0.d0
		
        vth  = fltarr(nr,nt) 
		  
        IF KEYWORD_SET(lefthand) THEN vphi = -vphi
    END
    3: BEGIN   ;Magnetospheric accretion (Hartmann, Hewett and Calvet, 1994)
       openr,1,'starinfo.inp'
       iformat=0
       readf,1,iformat
       rstar=0.d0
       mstar=0.d0
       readf,1,rstar
       readf,1,mstar
       close,1
       IF N_ELEMENTS(VERBOSE) THEN print, 'Used starinfo.inp for Funnel flow/Keplerian velocity'
 
       th = a.tt
       
       r  = a.rr
       rm = r/sin(th)^2d0
       y  = sin(th)^2.

       Rfl = rm*sin(th)*sin(th)^2.
       
       Rvec = 3.*y^0.5*(1.-y)^(0.5)/(4.-3.*y)^(0.5)
       zvec = (2.-3.*y)/(4.-3.*y)^(0.5)
       
       vp = (2.*GG*mstar/rstar*(rstar/Rfl - rstar/rm))^(0.5);/1d5

       ;The velocity vectors in cartesian coordinates       
       velx = -vp*Rvec
       vely = -vp*zvec

       ;make sure the field is properly mirrored:
       vely[*,a.ntheta/2:a.ntheta-1] = -vely[*,a.ntheta/2:a.ntheta-1]

       ;Now we must transform the cartesian
       ;velocity vectors to the polar basis.
       vth = cos(!pi/2.-th)*velx-sin(!pi/2.-th)*vely
       vr  = sin(!pi/2.-th)*velx+cos(!pi/2.-th)*vely
       vphi = SQRT(GG*mstar/a.rr)*0.
       IF KEYWORD_SET(lefthand) THEN vphi = -vphi
;Debug plot
;       gsubs = where(finite(velx) and finite(vely) AND Rfl GT rstar)
;       partvelvec,velx[gsubs],vely[gsubs], Rfl[gsubs]/AU, zfl[gsubs]/AU,/isotropic,xrange=[0,10],yrange=[-10,10],length=0.01

       
;
;Calculate the density in the flow, given a mass accretion rate
       Mdot = 1d-9*Mstar/(3600.*24.*365.)
       rmi = min(r)
       rmo = max(r)
       add_dens = Mdot * rstar/(4.*!pi*(Rstar/rmi - Rstar/rmo))*r^(-5./2.)/(2.*GG*mstar)^0.5*(4.-3.*y)^0.5/(1.-y)^0.5

    END
    4: BEGIN   ;Disk wind (Kurosawa et al. 2006)
       openr,1,'starinfo.inp'
       iformat=0
       readf,1,iformat
       rstar=0.d0
       mstar=0.d0
       readf,1,rstar
       readf,1,mstar
       close,1
       IF N_ELEMENTS(VERBOSE) THEN print, 'Used starinfo.inp for Keplerian velocity'
       vphi = sqrt(GG*mstar/a.rr)
       vth  = fltarr(nr,nt)
       vr   = fltarr(nr,nt)
       vstream = fltarr(nr,nt)
       l    = fltarr(nr,nt)
       rhowind  = fltarr(nr,nt)
       mdots    = fltarr(nr,nt)
       add_dens = fltarr(nr,nt*2)
;       tau_surf = 0.75
;       d = 5d0 * Rstar
;       v_esc = 20.d0*1d5 ;km/s
;       f     = 1.d0
;       p     = -7./2.
;       beta  = 2.0
;       scale_fac = 150.
;       mdot_sc   = 2.d-9        ;g/cm^2/s
       d = d_rstar * Rstar

       ;
       ;Define disk surface (the launching point of the wind - below this, the
       ;disk is Keplerian)
       make_tau, cc*1d4/0.55, vtau=tau
       surf_ith = fltarr(nr)
       
       FOR ir=0,nr-1 DO BEGIN
          surf_ith[ir] = MIN(WHERE(tau[ir,0:nt-1] GT tau_surf))
          IF surf_ith[ir] EQ -1 THEN surf_ith[ir] = nt-1
       ENDFOR
       surf_th = a.theta[surf_ith]
       surf_z  = a.r*sin(!pi/2.-surf_th)
       surf_w  = a.r*cos(!pi/2.-surf_th)
       ;
       ;Now loop over each grid point
       Tstruct = read_temperature()
       T = Tstruct.t
       gamma  = 1.4 ;Adiabatic constant for a diatomic gas
       Mu     = 2.3 ;Mean molecular mass

       cs = SQRT(gamma*kk*T/(Mu*mp))                   
       FOR ir=0,nr-1 DO BEGIN
          FOR ith=0,surf_ith[ir]-1 DO BEGIN
             w0 = a.r[ir]*cos(!pi/2.-a.theta[ith])
             z0 = a.r[ir]*sin(!pi/2.-a.theta[ith])
             diff = ((z0+d)/w0 * surf_w - d) - surf_z
             h=0
             REPEAT BEGIN
                h = h+1
             ENDREP UNTIL diff[h] GT 0
             cross_ir = h-1
             mindiff = diff[cross_ir]
             
;             mindiff  = MIN(diff)
;             absmindiff = MIN(ABS(diff),cross_ir)
 ;               plot, surf_w/AU, ((z0+d)/w0 * surf_w - d)/AU,xrange=[0,2],yrange=[0,2],/isotropic
 ;               oplot, [w0,w0]/AU,[z0,z0]/AU,psym=1
 ;               wait,0.05
 ;               print, mindiff
                ;
                ;The point where the stream line
                ;crosses the disk surface           

             IF mindiff LE 0 THEN BEGIN ;Check that the stream line crosses the disk surface (or comes close)
                w_cross = surf_w[cross_ir]
                z_cross = surf_z[cross_ir]

                ;
                ;Calculate l (the distance from the
                ;current grid point to the launch surface)
                l[ir,ith] = SQRT((w0-w_cross)^2. + (z0-z_cross)^2.)
                ;
                ;Then we can get the velocity along
                ;the stream line at the current grid point
                Rscale = scale_fac*w_cross[0]
                vstream[ir,ith] = cs[ir,ith] + (f*v_esc-cs[ir,ith])*(1.-Rscale/(l[ir,ith]+Rscale))^beta
                ;
                ;Angular momentum is conserved along the streamline (kinda)
                vphi[ir,ith] = vphi[cross_ir,surf_ith[cross_ir]] * w_cross/w0
                ;
                ;Now transform the velocity vector
                ;back on the polar basis
                ;
                ;Distance between streamline apex and grid point:
                psi = !pi-a.theta[ith]
                S   = SQRT(d^2.+a.r[ir]^2d0 - 2d0*d*a.r[ir]*cos(psi))
                ;
                ;Angle between theta and streamline is ksi:
                ksi = asin(sin(psi)*d/S)
                ;
                ;The radial velocity vector is
                ;projected from the streamline to the
                ;radial unit vector.
                vr[ir,ith]  = vstream[ir,ith]*cos(ksi)
                ;
                ;The theta unit vector is
                ;perpendicular to the radial vector, so:
                vth[ir,ith]  = vstream[ir,ith]*sin(ksi)
                ;
                ;The wind velocity field requires a
                ;certain density profile. delta is the
                ;angle between the streamline and the
                ;disk normal.
                delta = a.theta[ith] - ksi
                mdot = mdot_sc * (w_cross/a.r[0])^p
                mdots[ir,ith] = mdot
                rhowind[ir,ith] = mdot/(vstream[ir,ith]*ABS(cos(delta)))*(d/(S*cos(delta)))^2d0
             ENDIF ELSE BEGIN ;If the streamline does not cross the disk surface, the velocity is 0 (and the wind density should be 0 too!)
                vstream[ir,ith] = 0d0
                l[ir,ith] = alog(-1)
                ;
                ;And vphi remains Keplerian, the other
                ;velocity components are 0. 
             ENDELSE
          ENDFOR
       ENDFOR
       ;multiplied by 2 because the disk has two surfaces
       IF N_ELEMENTS(VERBOSE) THEN $
          print, 2d0*INT_TABULATED(a.r, 2.*!pi*a.r  *  mdot_sc*(a.r/a.r[0])^p,/sort,/double) * 3600.*24.*365./ms, ' Msol/yr, mass loss' 
       add_dens[*,0:nt-1] = rhowind
       add_dens[*,nt:*]   = ROTATE(rhowind,7)
       IF KEYWORD_SET(lefthand) THEN vphi = -vphi
   END
ENDCASE

openw,1,'velocity.inp'
printf,1,nr,nt
for ir=0,nr-1 do begin
    for it=0,nt-1 do begin
        printf,1,vr[ir,it],vth[ir,it],vphi[ir,it]
    endfor
endfor
close,1


END


PRO read_vel, file, vel
;
;Reads a velocity grid
;
openr,1,file
readf,1,nr,nt

vr = dblarr(nr,nt)
vtheta = dblarr(nr,nt)
vphi = dblarr(nr,nt)

for ir=0,nr-1 do begin
    for it=0,nt-1 do begin
        readf,1,dum1,dum2,dum3
        vr[ir,it]     = dum1
        vtheta[ir,it] = dum2
        vphi[ir,it]   = dum3
    endfor
endfor
close,1

vel = {vr:vr,vtheta:vtheta,vphi:vphi}

END
