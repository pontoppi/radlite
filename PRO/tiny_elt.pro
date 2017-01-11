;======================================================
;Procedure to calculate the diffracted PSF of a hexagonal 
;segmented telescope. The procedure uses the analytical
;description of (Sabatke et al. 2005,
;Applied optics, 44, 1360)
;
;Klaus Pontoppidan, December 2008
;
;======================================================

@segment.pro

PRO tiny_elt,fitsfile=fitsfile,lambda=lambda,tel=tel, N=N,max_ksi=max_ksi, single_photon=single_photon

IF NOT KEYWORD_SET(lambda) THEN lambda = 3.0 ;micron
IF NOT KEYWORD_SET(tel) THEN tel = 'EELT'
IF NOT KEYWORD_SET(fitsfile) THEN fitsfile = tel
IF NOT KEYWORD_SET(N) THEN N = 1000
IF NOT KEYWORD_SET(normalize) THEN normalize = 0
lambda_cm = lambda*1d-4   ;convert to cgs units (cm)

fitsfile = fitsfile+'_psf_'+STRTRIM(STRING(lambda),2)+'.fits'
fitsfile_pupil = fitsfile+'_pupil.fits'
;
;Define segment centers for the ELT:


IF tel eq 'EELT' THEN BEGIN
   dd      = 0d3                ;cm
   Rmirror = 1940d0             ;cm
   Rinner  = 555d0              ;cm
   dseg    = 140.*sqrt(3.)/2d0   ;125.574d0 ;cm
   nside   = 39.
ENDIF
IF tel eq 'Keck' THEN BEGIN
   dd      = 0.3                ;cm
   Rmirror = 600d0              ;cm
   Rinner  = 90d0              ;cm
   dseg    = 180.*sqrt(3.)/2d0   ;125.574d0 ;cm
   nside   = 7.
ENDIF
IF tel eq 'VLT' THEN BEGIN
   dd      = 0.03              ;cm
   Rmirror = 370d0             ;cm
   Rinner  = 90d0              ;cm
   dseg    = 60.*sqrt(3.)/2d0   ;125.574d0 ;cm
   nside   = 21.
ENDIF
IF tel eq 'JWST' THEN BEGIN
   dd      = 0.03              ;cm
   Rmirror = 325d0            ;cm
   Rinner  = 50d0              ;cm
   dseg    = 130.*sqrt(3.)/2d0   ;
   nside   = 5.
ENDIF
IF tel eq 'SUPEROWL' THEN BEGIN
   dd      = 0.03              ;cm
   Rmirror = 15000d0            ;cm
   Rinner  = 1500d0              ;cm
   dseg    = 2000.*sqrt(3.)/2d0   ;
   nside   = 40.
ENDIF
IF tel eq 'FIRS' THEN BEGIN
	dd      = 0.3
	Rmirror = 650.
	Rinner  = 100.
	dseg    = 200.*sqrt(3.)/2d0
	nside   = 21.
ENDIF


nsegs = nside^2

xcen = fltarr(nside,nside)
ycen = fltarr(nside,nside)


FOR i=-(nside-1)/2,(nside-1)/2 DO BEGIN
   FOR j=-(nside-1)/2,(nside-1)/2 DO BEGIN
      xcen[i+(nside-1)/2,j+(nside-1)/2] = i*(dseg*sqrt(3.)/2.+dd)
      ycen[i+(nside-1)/2,j+(nside-1)/2] = (j+0.5*(i mod 2))*(dseg+dd)
   ENDFOR
ENDFOR
radii = SQRT(xcen^2+ycen^2)
gsubs = WHERE((radii lt Rmirror-dseg/2.) AND (radii gt Rinner))

xcen = xcen[gsubs]
ycen = ycen[gsubs]

PRINT, 'Constructed primary with ', n_elements(gsubs), ' elements'

segment, 0.,0.,0.,dseg,hex=hex,N=N,eta=eta,ksi=ksi,max_ksi=max_ksi

FF = DBLARR(N,N)

i = COMPLEX(0,1)
FOR h=0,N-1 DO BEGIN
   FOR k=0,N-1 DO BEGIN
      FF[h,k] = TOTAL(EXP(i*2d0*!pi*(xcen*eta[k]+ycen*ksi[h]))*hex[h,k])
   ENDFOR
   writeu,-1,string(13b)
   writeu,-1,string(100d0*h/(N-1),format='(f5.1)')+'%'
ENDFOR

int = abs(FF)^2d0

pupil = ABS((FFT(int, /center, /inverse))^2d0)

IF single_photon NE 0 THEN BEGIN
	int = int/TOTAL(int)
ENDIF

;
;Create fits file and write it
mkhdr, HDR,int
xpixsize = max_ksi*2d0/N * lambda_cm * (180.d0/!PI) * 3600d0
ypixsize = xpixsize
sxaddpar, HDR, 'WAVEL',lambda, 'Wavelength [micron]'
sxaddpar, HDR, 'XPIX',xpixsize, 'Size of a pixel in the x direction [arcsec]'
sxaddpar, HDR, 'YPIX',ypixsize, 'Size of a pixel in the y direction [arcsec]'
writefits, fitsfile,int,HDR
writefits, fitsfile_pupil,pupil,HDR

plot, xcen, ycen, psym=4,/isotropic

END
