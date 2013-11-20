FUNCTION read_abundance,abunfile=abunfile

IF NOT KEYWORD_SET(abunfile) THEN abunfile = 'abundance.inp'

openr, lun, abunfile, /get_lun
readf,lun,nr,nt
abun = dblarr(nr,nt,2)

FOR ir = 0,nr-1 DO BEGIN
    FOR it = 0,nt-1 DO BEGIN
        readf,lun,dum1,dum2
        abun[ir,it,*] = [dum1,dum2]
    ENDFOR
ENDFOR

close,lun
free_lun,lun
;
;Also read the radius and theta

openr, lun, 'radius.inp', /get_lun
readf,lun,nr
rr = dblarr(nr)
FOR ir=0,nr-1 DO BEGIN
    readf,lun,dum
    rr[ir] = dum
ENDFOR
close,lun
free_lun,lun


openr, lun, 'theta.inp', /get_lun
readf,lun,nt
tt = dblarr(nt)
FOR it=0,nt-1 DO BEGIN
    readf,lun,dum
    tt[it] = dum
ENDFOR
close,lun
free_lun,lun


return, {abun:abun,rr:rr,tt:tt,nr:nr,nt:nt}

END

