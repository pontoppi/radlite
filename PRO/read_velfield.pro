FUNCTION read_velfield,velfile=velfile

IF NOT KEYWORD_SET(velfile) THEN velfile = 'velocity.inp'

openr, lun, velfile, /get_lun
readf,lun,nr,nt
vel=dblarr(nr,nt,3)

FOR ir = 0,nr-1 DO BEGIN
    FOR it = 0,nt-1 DO BEGIN
        readf,lun,dum1,dum2,dum3
        vel[ir,it,*] = [dum1,dum2,dum3]
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


return, {vel:vel,rr:rr,tt:tt,nr:nr,nt:nt}

END

