FUNCTION read_line, file


str=''

openr, lun, file,/get_lun
FOR i=1,4 DO readf,lun,str
readf,lun,dum,nfreq
nlines = dum

flux = dblarr(nfreq,nlines)
velo = dblarr(nfreq,nlines)
cfreqs = dblarr(nlines)
lcount = 0L

FOR i=1,1 DO readf,lun,str

j=0
WHILE NOT EOF(lun) DO BEGIN
;FOR j=0,nlines-1 DO BEGIN
   FOR i=1,2 DO readf,lun,str
   readf,lun,f0
   cfreqs[j] = f0
   FOR i=1,2 DO readf,lun,str
   FOR k=0,nfreq-1 DO BEGIN
      readf,lun,dum1,dum2
      velo[k,j]=dum1
      flux[k,j]=dum2
   ENDFOR
;ENDFOR
   j=j+1
ENDWHILE
close, lun
free_lun, lun

return, {velo:velo,flux:flux,cfreqs:cfreqs}


END
