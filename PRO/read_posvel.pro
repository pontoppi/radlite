PRO read_posvel,file,molecule=molecule

IF NOT KEYWORD_SET(molecule) then molecule = 'unknown'

dum = 0.
mol = ''
molfile = ''

openr,lun,file,/get_lun

readf,lun,dum
readf,lun,mol
readf,lun,molfile
readf,lun,dist,radvelo,incl
readf,lun,lev_up,lev_down
readf,lun,cent_freq
readf,lun,nfreq
readf,lun,imr_nx,imr_ny,imr_spx,imr_spy,imr_phioff,imr_xoff,imr_yoff

vel = fltarr(nfreq)
readf,lun,vel

Ints = fltarr(imr_nx,imr_ny,nfreq)
Taus = fltarr(imr_nx,imr_ny,nfreq)
FOR i=0,nfreq-1 DO BEGIN
    writeu,-1,string(13b)
    writeu,-1,string(100*i/(nfreq-1))+'%'
    FOR j=0,imr_nx-1 DO BEGIN
        FOR h=0,imr_ny-1 DO BEGIN
            readf,lun,dum1,dum2
            Ints[j,h,i] = dum1
            Taus[j,h,i] = dum2
        ENDFOR
    ENDFOR
ENDFOR
close,lun

;
;Create fits header
mkhdr, HDR,Ints
sxaddpar, HDR, 'XPIX',imr_spx, 'Size of a pixel in the x direction [cm]'
sxaddpar, HDR, 'YPIX',imr_spy, 'Size of a pixel in the y direction [cm]'
sxaddpar, HDR, 'INCL',incl, 'Inclination [degrees]'
sxaddpar, HDR, 'FREQ',cent_freq, 'Central frequency [Hz]'
sxaddpar, HDR, 'MOLEC',molecule, 'Molecular species'

;
;Create file name
fitsfile = STRMID(file,0,STRLEN(file)-4)+'.fits'

writefits, fitsfile,Ints,HDR
;
;Add the velocity grid as an extension
mwrfits, {velocity:vel}, fitsfile

;
;Create tau image
fitsfile = STRMID(file,0,STRLEN(file)-4)+'_tau.fits'
writefits, fitsfile,Taus,HDR
mwrfits, {velocity:vel}, fitsfile

END
