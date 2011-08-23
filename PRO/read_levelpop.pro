FUNCTION read_levelpop,lpop_file=lpop_file

IF NOT KEYWORD_SET(lpop_file) THEN lpop_file='levelpop_moldata.dat'

openr,lun,lpop_file,/get_lun

;
;First read array sizes
readf,lun,nr,nt,nl,dum
energy = fltarr(nl)
weights = fltarr(nl)
readf,lun,energy
readf,lun,weights

lpop = fltarr(nr,nt,nl)

FOR ir=0,nr-1 DO BEGIN
    FOR it=0,nt-1 DO BEGIN
        FOR il=0,nl-1 DO BEGIN
            readf,lun,dum
            lpop[ir,it,il]=dum
        ENDFOR
    ENDFOR
ENDFOR

close,lun
free_lun,lun

return,{nr:nr,nt:nt,nl:nl,lpop:lpop}

END
