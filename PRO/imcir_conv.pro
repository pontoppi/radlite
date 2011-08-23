@imdisp.pro

PRO imcir_conv, filename,xs=xs,ys=ys,plgrid=plgrid,ifr=ifr,nx=nx,ny=ny,$
                logsc = logsc,maxi=maxi,scale=scale,addstar=addstar,stopi=stopi,$
                saveit=saveit,plotit=plotit,zlog=zlog,wnum=wnum,sqrt=sqrt,$
                specast=specast,posang=posang,specres=specres,pvit=pvit

@natconst.pro
IF NOT KEYWORD_SET(xs)  THEN xs  = 10.0 ;AU
IF NOT KEYWORD_SET(ys)  THEN ys  = 10.0 ;AU
IF NOT KEYWORD_SET(ifr) THEN ifr = 0
IF NOT KEYWORD_SET(nx)  THEN nx = 500
IF NOT KEYWORD_SET(ny)  THEN ny = 500
IF NOT KEYWORD_SET(wnum) THEN wnum = 3
IF NOT KEYWORD_SET(posang) THEN posang = 0d0
IF NOT KEYWORD_SET(addstar) THEN addstar = 1
IF NOT KEYWORD_SET(saveit) THEN saveit='test.fits'
IF NOT KEYWORD_SET(spacast) THEN specast='specast.fits'

posang_rad = !pi*posang/180.d0   ;
xs_au = xs * AU
ys_au = ys * AU

IF KEYWORD_SET(addstar) THEN BEGIN
   nr_star = 30.
ENDIF ELSE BEGIN
   nr_star = 0.
ENDELSE

IF KEYWORD_SET(stopi) THEN ifr = 0

PRINT, 'Reading data'
imcir = read_imcir(filename,stopi=stopi)  ;We can stop the reader at some frequency (since this is a slow step).
rstar = imcir.ri[0]*0.99                  ;exactly 1 confuses the interpolation

IF ifr LT 0 THEN BEGIN
   first_fr = 0
   last_fr  = imcir.nfr-1
ENDIF ELSE BEGIN
   IF ifr GT imcir.nfr-1 THEN BEGIN
      PRINT, 'Requested frequency not in data!'
      stop
   ENDIF
   first_fr = ifr
   last_fr  = ifr
ENDELSE
IF ifr EQ 0 THEN BEGIN
   first_fr = 0
   last_fr = imcir.nfr-1
ENDIF

IF KEYWORD_SET(specres) THEN BEGIN
   fwhm   = specres/(imcir.velo[1]-imcir.velo[0])
   KERNEL = exp(-(findgen(imcir.nfr)-(imcir.nfr-1)/2.)^2.*2.77/fwhm^2)
   KERNEL = KERNEL/TOTAL(KERNEL)
   PRINT, 'In IMCIR_CONV.PRO: Convolving to requested spectral resolution: ', specres, ' km/s'
   FOR i=0,imcir.nphi-1 DO BEGIN
      FOR j=0,imcir.nr-1 DO BEGIN 
         imcir.imcir[*,i,j] = CONVOL(REFORM(imcir.imcir[*,i,j],imcir.nfr),KERNEL,TOTAL(KERNEL),/edge_truncate)
      ENDFOR
   ENDFOR
ENDIF

PRINT, 'Calculating spectrum from circular image'
fl_cir = fltarr(imcir.nfr)
FOR iifr=first_fr,last_fr DO BEGIN
   surf = !pi*rstar^2d0
   fl_cir[iifr] = fl_cir[iifr] + surf*imcir.radray[iifr]
   FOR ir=0,imcir.nr-1 DO BEGIN
      surf = (imcir.ri[ir+1]^2.d0-imcir.ri[ir]^2.d0)*!pi/imcir.nphi
      FOR iphi=0,imcir.nphi-1 DO BEGIN
         fl_cir[iifr] = fl_cir[iifr]+surf*imcir.imcir[iifr,iphi,ir]
      ENDFOR
   ENDFOR
ENDFOR
fl_cir = fl_cir/pc^2d0  ; move to a distance of 1 pc

imcir_im = fltarr(imcir.nphi,imcir.nr+nr_star)
phi = 2d0*!pi*findgen(imcir.nphi)/(imcir.nphi)
x        = fltarr(imcir.nphi,imcir.nr+nr_star)
y        = fltarr(imcir.nphi,imcir.nr+nr_star)
imrect   = fltarr(nx,ny,imcir.nfr)

   
FOR ip=0,imcir.nphi-1 DO BEGIN
   FOR ir=1,imcir.nr-1 DO BEGIN
      x[ip,ir]   = imcir.r[ir]*cos(phi[ip]-posang_rad)
      y[ip,ir]   = imcir.r[ir]*sin(phi[ip]-posang_rad)      
   ENDFOR
ENDFOR


PRINT, 'Interpolating on rectangular grid'
FOR iifr=first_fr,last_fr DO BEGIN
   imcir_im[0:imcir.nphi-1,0:imcir.nr-1] = imcir.imcir[iifr,*,*]
   
;
;Now add the star, if requested
 
   IF KEYWORD_SET(addstar) THEN BEGIN
      ;the last two points are outside the star to make sure the interpolation is grounded
      star_rad = rstar*(findgen(nr_star))/(nr_star-2-1) 
      
      FOR ip=0,imcir.nphi-1 DO BEGIN
         FOR ir=0,nr_star-1 DO BEGIN
            x[ip,ir+imcir.nr]     = star_rad[ir]*cos(phi[ip])
            y[ip,ir+imcir.nr]     = star_rad[ir]*sin(phi[ip])      
            P = asin(star_rad[ir]/rstar)
            IF star_rad[ir] LE rstar THEN BEGIN
               imcir_im[ip,ir+imcir.nr] = imcir.radray[iifr]*(1.d0-0.3*(1d0-cos(P)))
            ENDIF ELSE BEGIN
               imcir_im[ip,ir+imcir.nr] = imcir_im[ip,1]
            ENDELSE
         ENDFOR
      ENDFOR
   ENDIF
   ;
                                ;The rectangular image is a xs by ys
                                ;AU box with nx by ny pixels. For
                                ;example, if the user sets xs = ys =
                                ;1.0 AU, then the image will span -0.5
                                ;to 0.5 AU. 
   triangulate, x,y,triangles
   imrect[*,*,iifr] = trigrid(x,y,imcir_im[*,*],triangles,xgrid=xgrid,ygrid=ygrid,$
                    [xs_au/(nx-1),ys_au/(ny-1)], [-xs_au/2.,-ys_au/2.,xs_au/2.,ys_au/2.])
ENDFOR

spec2d = fltarr(nx,imcir.nfr)
IF KEYWORD_SET(pvit) THEN BEGIN
   PRINT, 'Calculating simulated 2D spectrum'
   FOR i=0,imcir.nfr-1 DO BEGIN
      spec2d[*,i] = TOTAL(imrect[*,*,i],2)
   ENDFOR
ENDIF

sa = fltarr(imcir.nfr)
fl = fltarr(imcir.nfr)
IF KEYWORD_SET(specast) THEN BEGIN
   PRINT, 'Calculating the spectro-astrometric signal at PA: ', posang, ' degrees'
   FOR ifr=first_fr,last_fr DO BEGIN
      im_coll = TOTAL(imrect[*,*,ifr],2,/DOUBLE)
      sa[ifr] = TOTAL(im_coll*xgrid,/DOUBLE)/TOTAL(im_coll,/DOUBLE)     
      fl[ifr] = TOTAL(imrect[*,*,ifr],/DOUBLE)*(xgrid[1]-xgrid[0])*(ygrid[1]-ygrid[0])
   ENDFOR
   fl = fl/pc^2d0
   
   PRINT, 'Saving spectro-astrometry to file: ', specast
   mwrfits, {xgrid:xgrid,ygrid:ygrid,velo:imcir.velo,sa:sa,fl:fl,fl_cir:fl_cir,$
             addstar:addstar,spec2d:spec2d,nu0:imcir.nu0,posang:posang,nx:nx,ny:ny,xs:xs,ys:ys},specast


   loadct,0
   @plot_setup.h
   set_plot, 'ps'
   device, filename='imcir_sa.eps',/encapsulated,color=0
   plot, imcir.velo, (sa-sa[0])/AU, psym=10,xtitle='!6Velocity [km s!u-1!N]',ytitle='Offset [AU / (1+C/L)]',$
         /xs,/ys
   device, /close
   set_plot, 'x'
   @plot_clear.h

ENDIF

IF KEYWORD_SET(saveit) THEN BEGIN
   imcir_file = saveit
   PRINT, 'Saving to .fits file: ', imcir_file
   mwrfits, imrect, imcir_file, /create
   mwrfits, {xgrid:xgrid,ygrid:ygrid,velo:imcir.velo,sa:sa,fl:fl,fl_cir:fl_cir,$
             addstar:addstar,spec2d:spec2d,nu0:imcir.nu0},imcir_file
ENDIF


IF KEYWORD_SET(plotit) THEN BEGIN
   PRINT, 'Plotting to eps, as requested.'
   @plot_setup.h
   set_plot, 'ps'
   device, filename='imcir.ps',/color,bits_per_pixel=24,encapsulated=0;,xsize=4,ysize=4
   loadct, 3
   IF KEYWORD_SET(zlog) THEN BEGIN
      imrect = alog10(imrect)
   ENDIF 
   IF KEYWORD_SET(sqrt) THEN BEGIN
      imrect = sqrt(imrect)
   ENDIF 

   FOR ifr=first_fr,last_fr DO BEGIN
      imrect_cg = imrect[*,*,ifr]
      
      
      mini = MIN(imrect[*,*,0])
      IF NOT KEYWORD_SET(maxi) THEN maxi = MAX(imrect[*,*,0])
      IF NOT KEYWORD_SET(scale) THEN scale = 1d0
      imrect_cg = ((imrect_cg-mini)/(maxi-mini))* 255. * scale
      msubs = WHERE(imrect_cg GT 255)
      IF msubs[0] NE -1 THEN imrect_cg[msubs] = 255
      msubs = WHERE(imrect_cg LT 0)
      IF msubs[0] NE -1 THEN imrect_cg[msubs] = 0
      
      minx = min(xgrid)/AU
      maxx = max(xgrid)/AU
      miny = min(ygrid)/AU
      maxy = max(ygrid)/AU

      imdisp, imrect_cg, /axis, xrange=[minx,maxx],$
              yrange=[miny,maxy],/xs,/ys,/erase,$
              xtitle='!6AU',ytitle='!6AU',/negative,/noscale,color=140
      xyouts, (maxx-minx)*0.1+minx,(maxy-miny)*0.9+miny,$
              strtrim(string(imcir.velo[ifr],format='(f6.2)'),2)+' km s!u-1!N',color=5000,/data
      
      IF KEYWORD_SET(plgrid) THEN BEGIN
;         oplot, x/AU,y/AU,psym=1,symsize=0.01,color=500
         cmasksubs = WHERE(imcir.cmask[ifr,*,*] EQ 1)
         IF cmasksubs[0] NE -1 THEN $
            oplot, x[cmasksubs]/AU,y[cmasksubs]/AU,psym=4,color=5600,symsize=0.01
      ENDIF
   ENDFOR
   device,/close
   set_plot, 'x'
   @plot_clear.h
ENDIF

END
