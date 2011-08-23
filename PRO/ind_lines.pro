!p.multi=[0,1,1]

set_plot, 'ps'

device, filename='individual_lines.ps'

close,1

spawn, 'ls linespectrum_moldata_*dat', linefiles
nfiles = N_ELEMENTS(linefiles)

IF KEYWORD_SET(maxnfile) THEN nfiles = maxnfile

nlines= fltarr(nfiles)
str= ''

openr,1, 'linespectrum_moldata_0.dat'
FOR i=1,4 DO BEGIN
   readf,1,str
   print, str
ENDFOR

readf,1,dum,nfreq
print, dum,nfreq
close,1

lines = dblarr(2000, nfreq, 2)
cfreqs = dblarr(2000)
lcount=0L

FOR ff=0,nfiles-1 DO BEGIN

   close,1
   openr,1,'linespectrum_moldata_'+STRTRIM(STRING(ff),2)+'.dat'

   FOR i=1,4 DO BEGIN
      readf,1,str
      print, str
   ENDFOR
   readf,1,dum,nfreq
   print, dum, nfreq

   nlines[ff] = dum
   FOR i=1,1 DO readf,1,str

   FOR j=0,nlines[ff] - 1 DO BEGIN 
      FOR i=1,2 DO readf,1,str
      readf,1,f0
      cfreqs[lcount] = f0
      FOR i=1,2 DO readf,1,str
      FOR k=0,nfreq-1 DO BEGIN
         readf,1,dum1, dum2
         lines[lcount,k,*]=[dum1,dum2]
      ENDFOR
      plot, lines[lcount,*,0], smooth(lines[lcount,*,1],2) * 1e23, xtitle='Velocity [km/s]', ytitle='Flux [Jy]', $
            title=string(cfreqs[lcount]) / 2.998e10, $
            yrange=[0.999 * min(lines[lcount,*,1] * 1e23), 1.001 *  max(lines[lcount,*,1] * 1e23)], $
            xrange=[-80,80], xstyle=1
      FOR k=0,nfreq-1 DO print, lines[lcount,k,0], lines[lcount,k,1]
      lcount = lcount + 1
   ENDFOR
   close,1
ENDFOR


device, /close

set_plot, 'x'

end
