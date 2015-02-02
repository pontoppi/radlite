;==============================================================
;Procedure to simulate a telescopic line image, including noise.
;The input is an image data cube from RADLite. The original 
;purpose was to simulate line imaging with ELT-METIS.
;
;Klaus Pontoppidan, December 2008
;
;===============================================================


PRO telsim, lineposvel,lineposvel_small=lineposvel_small,dist=dist,$
            telescope=telescope,outname=outname,psf=psf,vsamp=vsamp,$
            slitwidth=slitwidth,specres=specres,noiseless=noiseless,$
            add_model=add_model,pixel_aspect=pixel_aspect
@natconst.pro

IF ~KEYWORD_SET(slitwidth) THEN slitwidth = 0.2 ;arcsec

IF ~KEYWORD_SET(pixel_aspect) THEN pixel_aspect = 1.

IF ~KEYWORD_SET(dist) THEN BEGIN
   dist = 125d0   ;pc
ENDIF

IF ~KEYWORD_SET(telescope) THEN BEGIN
   telescope = 'EELT'
ENDIF
IF ~KEYWORD_SET(outname) THEN BEGIN
   outname = 'telescope_image'
ENDIF

fidelity = 500.

CASE telescope OF
   'VLT_K' : BEGIN
      eltpix  = 0.020           ;arcsec
      eltres  = 0.056           ;arcsec
      sigma   = 0.6e-6          ;Jy/beam
   END
   'VLT' : BEGIN
      eltpix  = 0.043           ;arcsec
      eltres  = 0.120           ;arcsec
      sigma   = 0.002           ;Jy/beam
   END
   'EELT_CO' : BEGIN
      eltpix  = 0.009          ;arcsec
      eltres  = 0.034           ;arcsec
      sigma   = 0.00023          ;Jy/beam
   END
   'EELT_N' : BEGIN
      eltpix  = 0.0172          ;arcsec
      eltres  = 0.100           ;arcsec
      sigma   = 0.4*1d-3        ;Jy/beam
   END
   'EELT_LW' : BEGIN
      eltpix  = 0.020           ;arcsec
      eltres  = 0.060           ;arcsec
      sigma   = 3.70*1d-4        ;Jy/beam
   END
    'MIRI' : BEGIN
      eltpix  = 0.15             ;arcsec
      eltres  = 0.5              ;arcsec
      sigma   = 7.0*1d-5        ;Jy/beam (6d-21 W/m^2 at 22.5)
   END
    'NIRSPEC' : BEGIN
      eltpix  = 0.08             ;arcsec
      eltres  = 0.20             ;arcsec
      sigma   = 7.0*1d-6        ;Jy/beam (6d-21 W/m^2 at 22.5)
   END
   'ATLAST' : BEGIN
     eltpix  = 0.005            ;arcsec
     eltres  = 0.010             ;arcsec
     sigma   = 7.*1d-9        ;Jy/beam 
  	END   
    'SUPER_OWL' : BEGIN
      eltpix  = 0.0005
      eltres  = 0.0010
      sigma   = 1d-8
    END
  ELSE : BEGIN
      PRINT, 'Unknown telescope!'
      stop
   ENDELSE
ENDCASE

a    = readfits(lineposvel)
info = mrdfits(lineposvel,1)
; 
;Get the size of a pixel (in cm) and central frequency (in Hz).

xpixsize = info.xgrid[1]-info.xgrid[0]
ypixsize = info.ygrid[1]-info.ygrid[0]
freq_hz  = info.nu0

IF xpixsize NE ypixsize THEN BEGIN
   PRINT, 'Non-square pixels not yet supported!'
   stop
ENDIF
pixsize_cm = xpixsize
nx = (size(a))[1]
nf = (size(a))[3]
line = fltarr(nf)
findex = findgen(nf)

;
;Allow for babuschka images
IF KEYWORD_SET(lineposvel_small) THEN BEGIN
   a_s=readfits(lineposvel_small,hdr)
   info_s = mrdfits(lineposvel_small,1)
   ; 
   ;Get the size of a pixel (in cm) and central frequency (in Hz).
   xpixsize_s = info_s.xgrid[1]-info_s.xgrid[0]
   ypixsize_s = info_s.ygrid[1]-info_s.ygrid[0]
   freq_hz_s  = info_s.nu0
   
   IF xpixsize NE ypixsize THEN BEGIN
      PRINT, 'Non-square pixels not yet supported!'
      stop
   ENDIF
   pixsize_s_cm = xpixsize_s
   nx_s = (size(a_s))[1]
   nf_s = (size(a_s))[3]
   IF nf NE nf_s THEN BEGIN
      PRINT, 'Something is wrong, babuschka images do not match'
      stop
   ENDIF
   factor = pixsize_s_cm/pixsize_cm
   IF abs(nx*factor - ROUND(nx*factor))/(nx*factor) GT 1d-2 THEN BEGIN
      PRINT, 'Image size must be divisible by the babuschka size!'
      stop
   ENDIF
   n_int = 2*(ROUND(nx_s*factor)/2)+1
   FOR i=0,nf-1 DO BEGIN
      tmp = FREBIN(a_s[*,*,i],n_int,n_int)
      a[nx/2-n_int/2:nx/2+n_int/2,nx/2-n_int/2:nx/2+n_int/2,i] = tmp
   ENDFOR
ENDIF

;
;determine angular pixel size
pixsize = pixsize_cm/AU/dist ;arcsec/pixel

;
;convert image to Jy/arcsec^2
;Old conversion from the rectangular RADLite imager
;a = a * freq_hz^2./3.25465503368d36
;a = a * 2.35d-11 *1d23

a = a * 2.35d-11 * 1d23 ;erg/s/cm^2/Hz/ster --> Jy/arcsec^2

;
;Read telescope convolution kernel

kernel    = readfits(psf,hdr)
k_pixsize = sxpar(HDR, 'xpix')
nk        = (size(kernel))[2]

;
;Rebin the kernel on the image grid
factor = k_pixsize/pixsize

kernel_int = FREBIN(kernel,ROUND(nk*factor),ROUND(nk*factor))
error = abs(1-ROUND(nk*factor)/(nk*factor))
print, 'Kernel rebin error: ', error*100, '%'
kernel_int = kernel_int/total(kernel_int)

size_im     = (SIZE(a[*,*,0]))[1]
size_kernel = (SIZE(kernel_int))[1]
IF size_im LT size_kernel THEN BEGIN
   kernel_int = kernel_int[(size_kernel-size_im)/2:(size_kernel+size_im)/2,(size_kernel-size_im)/2-1:(size_kernel+size_im)/2-1]
ENDIF
 
;
;Rebin to ELT pixels
n_elt_pix = floor(nx*pixsize/eltpix)
im_elt = fltarr(n_elt_pix, n_elt_pix,nf)

;
;Convolve images
FOR i=0,nf-1 DO BEGIN
   a[*,*,i] = CONVOLVE(a[*,*,i],kernel_int[*,*])         
   im_elt[*,*,i] = FREBIN(a[*,*,i],n_elt_pix,n_elt_pix)
   WRITEU,-1,string(13b)
   WRITEU,-1, 'convolving frequency '+STRTRIM(STRING(i+1),2), ' of ', STRTRIM(STRING(nf),2)
ENDFOR
WRITEU,-1, ' '

IF KEYWORD_SET(specres) THEN BEGIN
   fwhm   = specres/(info.velo[1]-info.velo[0])
   sp_KERNEL = exp(-(findgen(nf)-(nf-1)/2.)^2.*2.77/fwhm^2)
   sp_KERNEL = sp_KERNEL/TOTAL(sp_KERNEL)
   PRINT, 'In TELSIM.PRO: Convolving to requested spectral resolution: ', specres, ' km/s'
   FOR ix=0,n_elt_pix-1 DO BEGIN
      FOR iy=0,n_elt_pix-1 DO BEGIN
          im_elt[ix,iy,*] = CONVOL(REFORM(im_elt[ix,iy,*],nf),sp_KERNEL,TOTAL(sp_KERNEL),/edge_truncate)        
      ENDFOR
   ENDFOR
ENDIF

;
;Resample
IF KEYWORD_SET(vsamp) THEN BEGIN
   new_nf = ROUND(nf*(info.velo[1]-info.velo[0])/vsamp)
   im_resamp = fltarr(n_elt_pix,n_elt_pix,new_nf)
   FOR i=0,n_elt_pix-1 DO BEGIN
      im_resamp[i,*,*] = $
         FREBIN(REFORM(im_elt[i,*,*],n_elt_pix,nf),n_elt_pix,new_nf)
   ENDFOR
   velo = FREBIN(info.velo,new_nf)
   im_elt = im_resamp
   nf = new_nf
ENDIF ELSE BEGIN
   velo = info.velo
ENDELSE

;There is an option to add a noiseless model cube at this point. The
;original purpose was to add a planet calculation. It is the users
;responsibility that the added model flux cube actually makes sense.
IF KEYWORD_SET(add_model) THEN BEGIN
   xoff = 13
   yoff = -3
   fsh  = -2
   add_data = MRDFITS(add_model)
   add_data = SHIFT(add_data,0,0,fsh)
   nf_add = (SIZE(add_data))[3]
   nx_add  = (SIZE(add_data))[1]
   ny_add  = (SIZE(add_data))[2]
   IF nf_add NE nf THEN BEGIN
      print, 'Shape of added model cube does not match!'
      stop
   ENDIF
   im_elt[(n_elt_pix-nx_add)/2+xoff:(n_elt_pix+nx_add)/2+xoff-1,(n_elt_pix-ny_add)/2+yoff:(n_elt_pix+ny_add)/2+yoff-1,*] $
      += add_data

ENDIF

FOR i=0,nf-1 DO BEGIN
   ;
   ;Add noise
   IF ~KEYWORD_SET(noiseless) THEN BEGIN
      ns = randomn(seed,n_elt_pix,n_elt_pix) * $
           (sigma/(!pi*(eltres/2.)^2)+im_elt[*,*,i]/fidelity)
      im_elt[*,*,i] = im_elt[*,*,i] + ns
   ENDIF

   ;
   ;and calculate line
   line[i] = total(im_elt[*,*,i])*eltpix^2.
ENDFOR

;
;We can change the pixel aspect ratio at this point. Note that noise
;per pixel will not be preserved. 
IF pixel_aspect NE 1 THEN BEGIN
   ny_new = ROUND(n_elt_pix/pixel_aspect)
   im_elt_rebin = fltarr(n_elt_pix,ny_new,nf)
   FOR i=0,nf-1 DO BEGIN
      im_elt_rebin[*,*,i] = FREBIN(im_elt[*,*,i],n_elt_pix,ny_new)
   ENDFOR
   im_elt_original = im_elt
   im_elt = im_elt_rebin
ENDIF

mwrfits, im_elt, outname+'.fits', /CREATE
line_im = fltarr(n_elt_pix,n_elt_pix)
FOR i=0,nf-1 DO BEGIN
   line_im = line_im + (im_elt[*,*,i] - (im_elt[*,*,0]+im_elt[*,*,nf-1])/2d0)
ENDFOR
mwrfits, line_im,outname+'.fits'
mwrfits, {totspec:total(total(im_elt,1),1)*eltpix^2,velo:velo, pixsize:eltpix, pixunit:'Arcsec/pixel'},outname+'.fits'

IF pixel_aspect NE 1 THEN BEGIN
   im_elt = im_elt_original
ENDIF
;
;Create slit image
slit_im = total(im_elt,2)
mwrfits, slit_im,outname+'_slit.fits',/CREATE
;
;Create astrometric spectrum
slitradius_pix = ROUND(slitwidth/eltpix/2.0)
sa = fltarr(nf)
sa_index = findgen(slitradius_pix*2+1)-slitradius_pix
fl = fltarr(nf)
FOR i=0,nf-1 DO BEGIN
   sa[i] = TOTAL(sa_index*slit_im[n_elt_pix/2-slitradius_pix:n_elt_pix/2+slitradius_pix,i]) / $
           TOTAL(slit_im[n_elt_pix/2-slitradius_pix:n_elt_pix/2+slitradius_pix,i])
   fl[i] = TOTAL(slit_im[n_elt_pix/2-slitradius_pix:n_elt_pix/2+slitradius_pix,i])
ENDFOR
sa    = sa * eltpix
mwrfits, {velo:velo, fl:fl, sa:sa},outname+'_slit.fits'

;
;Create continuum subtracted slit image
im_elt_csub = fltarr(n_elt_pix,n_elt_pix,nf)
FOR i=0,nf-1 DO BEGIN
   im_elt_csub[*,*,i] = im_elt[*,*,i] - im_elt[*,*,0]
ENDFOR

IF KEYWORD_SET(slitwidth) THEN BEGIN
   slitim = total(im_elt_csub[*,n_elt_pix/2-slitradius_pix:n_elt_pix/2+slitradius_pix,*],2)
ENDIF ELSE BEGIN
   slitim = total(im_elt_csub,2)
ENDELSE

mwrfits, slitim, outname+'_slit_csub.fits',/CREATE



PRINT, ' '
PRINT, 'Success!'

END
