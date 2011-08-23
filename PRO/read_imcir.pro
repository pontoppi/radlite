function read_imcir, filename,stopi=stopi

nrmax = 1000.
nphimax = 1000.
str = ' '
dum = 0d0

openr, lun, filename, /get_lun
readf,lun,nfr
velo = fltarr(nfr)
IF KEYWORD_SET(stopi) THEN BEGIN
   imcir = fltarr(1,nphimax,nrmax)
   imcir_cmask = bytarr(1,nphimax,nrmax)
   radray = fltarr(1)
ENDIF ELSE BEGIN
   imcir = fltarr(nfr,nphimax,nrmax)
   imcir_cmask = bytarr(nfr,nphimax,nrmax)
   radray = fltarr(nrmax)
ENDELSE

readf,lun,nu0
readf,lun,nphi,nr
ri = fltarr(nr+1)
r = fltarr(nr+1)
FOR ir=0,nr DO BEGIN
   readf,lun,dum
   ri[ir] = dum
ENDFOR
FOR ir=0,nr DO BEGIN
   readf,lun,dum
   r[ir] = dum
ENDFOR

;In some cases, one may not want to read all the frequencies, 
;e.g., if RADLite is still running. 
IF KEYWORD_SET(stopi) THEN BEGIN
   nfr = stopi
ENDIF

FOR ifr=0,nfr-1 DO BEGIN
   readf,lun,str
   readf,lun,dumvelo
   velo[ifr] = dumvelo
   readf,lun,str
   readf,lun,radray_dum
   radray[ifr] = radray_dum
   FOR iphi=0,nphi-1 DO BEGIN
      FOR ir=0,nr-1 DO BEGIN
         readf,lun,dum,flag        
         imcir[ifr,iphi,ir] = dum
         imcir_cmask[ifr,iphi,ir] = flag
      ENDFOR
   ENDFOR
ENDFOR
close,lun
free_lun, lun

imcir = imcir[*,0:nphi-1,0:nr-1]
imcir_cmask = imcir_cmask[*,0:nphi-1,0:nr-1]
radray      = radray[0:nfr-1]
return, {imcir:imcir,cmask:imcir_cmask,$
         ri:ri,r:r,nphi:nphi,nr:nr,nfr:nfr,$
         radray:radray,velo:velo,nu0:nu0}

END
