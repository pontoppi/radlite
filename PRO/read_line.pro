FUNCTION read_line, file

fileformat1 = 0
fileformat2 = 0
molname     = ''
molfile     = ''
nlines      = 0
nfreq       = 0
levup       = 0
levdown     = 0
f0          = 0.
beamsize    = 0.
nfreq       = 0
velo_dum    = 0.
flux_dum    = 0.
str         = ' '

openr, lun, file,/get_lun

readf, lun, fileformat1, format='(i2)'
readf, lun, fileformat2, format='(i2)'
readf, lun, molname,  format='(a80)'
readf, lun, molfile,  format='(a80)'
readf, lun, nlines, format='(i10)'
readf, lun, maxfreq, format='(i10)'

readf, lun, distance, radvelo, inclination, format='(e12.4,e12.4,f7.3)'

flux = dblarr(maxfreq,nlines)
velo = dblarr(maxfreq,nlines)
cfreqs = dblarr(nlines)
lcount = 0L

FOR j=0,nlines-1 DO BEGIN
   
   readf, lun, str,format='(a1)' ;empty line
   readf, lun, levup, levdown, format='(i5,i5)'
   readf, lun, f0, format='(f14.9)'
   readf, lun, beamsize, format='(f10.5)'
   readf, lun, nfreq, format='(i5)'
   readf, lun, str,format='(a1)' ;empty line            

   cfreqs[j] = f0


   FOR k=0,nfreq-1 DO BEGIN
      readf,lun,velo_dum,flux_dum, format='(E13.6,1X,E13.6)'
      velo[k,j] = velo_dum
      flux[k,j] = flux_dum
   ENDFOR
ENDFOR
close, lun
free_lun, lun

velo = velo[0:nfreq-1,*]
flux = flux[0:nfreq-1,*]


return, {velo:velo,flux:flux,cfreqs:cfreqs,nlines:nlines,nfreq:nfreq,inclination:inclination,distance:distance,radvelo:radvelo,maxfreq:maxfreq}


END
