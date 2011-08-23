@imdisp.pro

PRO imcir_conv, filename,xs=xs,ys=ys,plgrid=plgrid,ifr=ifr,nx=nx,ny=ny,$
                logsc = logsc,maxi=maxi,scale=scale,addstar=addstar,stopi=stopi,$
                saveit=saveit,plotit=plotit,zlog=zlog,wnum=wnum

@natconst.pro
IF NOT KEYWORD_SET(xs)  THEN xs  = 1.0 ;AU
IF NOT KEYWORD_SET(ys)  THEN ys  = 1.0 ;AU
IF NOT KEYWORD_SET(ifr) THEN ifr = 0
IF NOT KEYWORD_SET(nx)  THEN nx = 500
IF NOT KEYWORD_SET(ny)  THEN ny = 500
IF NOT KEYWORD_SET(wnum) THEN wnum = 3

xs = xs * AU
ys = ys * AU

IF KEYWORD_SET(addstar) THEN BEGIN
   nr_star = 30.
ENDIF ELSE BEGIN
   nr_star = 0.
ENDELSE

IF KEYWORD_SET(stopi) THEN ifr = 0

imcir = read_imcir(filename,stopi=stopi)  ;We can stop the reader at some frequency (since this is a slow step).
rstar = imcir.imcir_ri[1]*0.99           ;exactly 1 confuses the interpolation


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

imcir_im = fltarr(imcir.nphi,imcir.nr+nr_star)
phi = 2d0*!pi*findgen(imcir.nphi)/(imcir.nphi)
x        = fltarr(imcir.nphi,imcir.nr+nr_star)
y        = fltarr(imcir.nphi,imcir.nr+nr_star)
imrect   = fltarr(nx,ny,imcir.nfr)

   
FOR ip=0,imcir.nphi-1 DO BEGIN
   FOR ir=1,imcir.nr-1 DO BEGIN
      x[ip,ir]   = imcir.imcir_ri[ir]*cos(phi[ip])
      y[ip,ir]   = imcir.imcir_ri[ir]*sin(phi[ip])      
   ENDFOR
ENDFOR
FOR iifr=first_fr,last_fr DO BEGIN
   imcir_im[0:imcir.nphi-1,0:imcir.nr-1] = imcir.imcir[iifr,*,*]
   
;
;Now add the star, if requested
   IF KEYWORD_SET(addstar) THEN BEGIN
      star_rad = rstar*(findgen(nr_star))/(nr_star-1)
      
      FOR ip=0,imcir.nphi-1 DO BEGIN
         FOR ir=0,nr_star-1 DO BEGIN
            x[ip,ir+imcir.nr]     = star_rad[ir]*cos(phi[ip])
            y[ip,ir+imcir.nr]     = star_rad[ir]*sin(phi[ip])      
            P = asin(star_rad[ir]/rstar)
            imcir_im[ip,ir+imcir.nr] = imcir.radray*(1.d0-0.3*(1d0-cos(P)))
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
                    [xs/(nx-1),ys/(ny-1)], [-xs/2.,-ys/2.,xs/2.,ys/2.])
   
ENDFOR

IF KEYWORD_SET(plotit) THEN BEGIN
   @plot_setup.h
   set_plot, 'ps'
   device, filename='imcir.ps',/color,xsize=4,ysize=4
;   device,decomposed=0
   loadct, 3
;   window,wnum,xsize=1000,ysize=1000
   FOR ifr=first_fr,last_fr DO BEGIN
      contour, imrect[*,*,ifr],xgrid/AU,ygrid/AU,/fill,/isotropic,$
               xrange=[-xs,xs]/AU/2d0,yrange=[-ys,ys]/AU/2d0,nlevels=2,/xs,/ys,$
               charsize=0.1,xtitle='!6AU', ytitle='!6AU',/nodata
      PX2=!X.WINDOW*!D.X_VSIZE
      PY2=!Y.WINDOW*!D.Y_VSIZE
      SX2 = PX2[1]-PX2[0]+1
      SY2 = PY2[1]-PY2[0]+1
      imrect_cg = CONGRID(imrect[*,*,ifr],sx2,sy2)
      
      IF KEYWORD_SET(zlog) THEN BEGIN
         imrect_cg = alog10(imrect_cg)
      ENDIF 
      
      mini = MIN(imrect_cg)
      IF NOT KEYWORD_SET(maxi) THEN maxi = MAX(imrect_cg)
      IF NOT KEYWORD_SET(scale) THEN scale = 1d0
      imrect_cg = ((imrect_cg-mini)/(maxi-mini))* 255. * scale
      msubs = WHERE(imrect_cg GT 255)
      IF msubs[0] NE -1 THEN imrect_cg[msubs] = 255
;   imrect_cg = smooth(imrect_cg,5)
      
      tv,imrect_cg,PX2[0],PY2[0],xsize=SX2,ysize=SY2
      
      contour, imrect[*,*,ifr],xgrid/AU,ygrid/AU,/fill,/isotropic,$
               xrange=[-xs,xs]/AU/2d0,yrange=[-ys,ys]/AU/2d0,nlevels=2,/xs,/ys,$
               charsize=0.1,xtitle='!6AU', ytitle='!6AU',/nodata,/noerase
      
      IF KEYWORD_SET(plgrid) THEN BEGIN
         oplot, x/AU,y/AU,psym=1,symsize=0.05
         cmasksubs = WHERE(imcir.cmask[ifr,*,*] EQ 1)
         IF cmasksubs[0] NE -1 THEN $
            oplot, x[cmasksubs]/AU,y[cmasksubs]/AU,psym=4,color=5600,symsize=0.1
      ENDIF
   endfor
   device,/close
   set_plot, 'x'
   @plot_clear.h
ENDIF

IF KEYWORD_SET(saveit) THEN BEGIN
   writefits, imrect,filename=saveas
ENDIF



stop
END
