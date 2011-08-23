FUNCTION read_meanint, meanint_file=meanint_file

IF NOT KEYWORD_SET(meanint_file) THEN meanint_file = 'meanint_radmc.dat'

openr,lun,meanint_file, /get_lun
readf,lun,nf,nr,nt,dum
MI = dblarr(nr,nt,nf)
FOR fi=0,nf-1 DO BEGIN
    FOR ti=0,nt-1 DO BEGIN
        FOR ri=0,nr-1 DO BEGIN
            readf,lun,dum
            MI[ri,ti,fi]=dum
        ENDFOR
    ENDFOR
ENDFOR
close,lun
free_lun,lun

close
stop
return, MI


END
