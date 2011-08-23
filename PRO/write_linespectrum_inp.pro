;===========================================
;
;Procedure to write the linespectrum.inp input file 
;to RADLite.
;
;===========================================


PRO write_linespectrum_inp,nlines=nlines,incl=incl,dist=dist,image=image,$
                           vsampling=vsampling,ssampling=ssampling,passband=passband,$
                           imwidth=imwidth,molfile=molfile

IF NOT KEYWORD_SET(incl) THEN incl = 45.0
IF NOT KEYWORD_SET(vlsr) THEN vlsr = 0.0
IF NOT KEYWORD_SET(dist) THEN dist = 1.0
IF NOT KEYWORD_SET(vsampling) THEN vsampling = 1.2 ;km/s
IF NOT KEYWORD_SET(ssampling) THEN ssampling = 0.1 ;AU/pixel
IF NOT KEYWORD_SET(passband) THEN passband = 30.0
IF NOT KEYWORD_SET(imwidth) THEN imwidth = 10. ;AU
npix = CEIL(imwidth/ssampling)
IF image EQ 0 OR image EQ 2 THEN BEGIN
    radmode = 0 
    telesc_command = 4
ENDIF ELSE BEGIN
    radmode = 2
    telesc_command = 6
ENDELSE

AU = 1.496e+13   ;cm
imwidth_cm = imwidth*AU

openw,lun,'linespectrum.inp',/get_lun

printf,lun, '1       Format number'
printf,lun, '1       Spectrum output style'
printf,lun, '---------------------------------------------------------------'
printf,lun, '2       Format number'
printf,lun, '---------------------------------------------------------------'
printf,lun, STRTRIM(STRING(passband,format='(f6.2)'),2),'Width of line passband [km/s]',          format='(a-8,a-50)'
printf,lun, STRTRIM(STRING(vsampling,format='(f6.2)'),2),'Velocity sampling [km/s]',              format='(a-8,a-50)'
printf,lun, '-----------------------------------------------------------------'
printf,lun, '2       Format number'
printf,lun, STRTRIM(STRING(molfile),2),'Molecular data file',                                                        format='(a-30,a-50)'
printf,lun, STRTRIM(STRING(radmode,format='(i1)'),2),'Command (0=spectrum,2=image[3-D P/V cube])',format='(a-8,a-50)'
printf,lun, STRTRIM(STRING(dist,format='(f7.1)'),2),'Distance in [pc]',                           format='(a-8,a-50)'
printf,lun, STRTRIM(STRING(incl,format='(f6.1)'),2),'Inclination [deg]',                          format='(a-8,a-50)'
printf,lun, STRTRIM(STRING(vlsr,format='(f6.1)'),2),'Radial velocity, rel. to local standard of rest [km/s]',format='(a-8,a-60)'
printf,lun, STRTRIM(STRING(nlines,format='(i4)'),2),'Nr of lines to make spectrum/image',         format='(a-8,a-50)'
printf,lun, '1       Starting line to make spectrum/image'
IF radmode eq 2 THEN BEGIN
    printf,lun, STRTRIM(STRING(npix,format='(i4)'),2),'Nr of x pixels',                           format='(a-8,a-50)'
    printf,lun, STRTRIM(STRING(npix,format='(i4)'),2),'Nr of y pixels',                           format='(a-8,a-50)'
    printf,lun, '1       image size in cm?'
    printf,lun, STRTRIM(STRING(imwidth_cm,format='(e8.1E2)'),2),'size x direction',                    format='(a-8,a-50)'
    printf,lun, STRTRIM(STRING(imwidth_cm,format='(e8.1E2)'),2),'size y direction',                    format='(a-8,a-50)'
    printf,lun, '0       Phi offset?'
    printf,lun, '0       x offset?'
    printf,lun, '0       y offset?'
    printf,lun, '1       add star?'
ENDIF
close,lun
free_lun, lun

END
