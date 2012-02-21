PRO genspec,obsres=obsres,maxnfile=maxnfile,xrange=xrange,$
              psym=psym,yrange=yrange,xlog=xlog,ylog=ylog,$
              scale_spec=scale_spec,sedcomp=sedcomp,sampling=sampling,$
              multiruns=multiruns,noise=noise,dist=dist

LMAX = 5000
n_x_max = 1d7
c=2.99792458d14

IF NOT KEYWORD_SET(dist) THEN dist = 150. ;pc

IF NOT KEYWORD_SET(obsres) THEN obsres=3.
IF NOT KEYWORD_SET(scale_spec) THEN scale_spec = 1d0

IF NOT KEYWORD_SET(sampling) THEN BEGIN
   reswidth = 1.0              ;km/s
ENDIF ELSE BEGIN
   reswidth = sampling
ENDELSE

obsres   = obsres/reswidth
max_vel  = 500. ;km/s : extend the box width to this

IF KEYWORD_SET(sedcomp) THEN BEGIN
;For testing purposes only
   readcol, '../spectrum.dat',x,y
   x = c/x
   y = y*1d23
ENDIF

IF KEYWORD_SET(multiruns) THEN BEGIN
   FOR i=0,N_ELEMENTS(multiruns)-1 DO BEGIN
      spawn, 'ls '+multiruns[i]+'/linespectrum*.dat',linefiles_subrun
      IF linefiles_subrun[0] EQ '' THEN BEGIN
         PRINT, 'No files in specified subdirectory: ', multiruns[i]
         STOP
      ENDIF
;      FOR j=0,N_ELEMENTS(linefiles_subrun)-1 DO BEGIN
;         stop
;         linefiles_subrun[j] = multiruns[i]+'/'+linefiles_subrun[j]
;      ENDFOR
      IF i NE 0 THEN BEGIN
         linefiles = [linefiles,linefiles_subrun]
      ENDIF ELSE BEGIN
         linefiles = linefiles_subrun
      ENDELSE
   ENDFOR
   nfiles = N_ELEMENTS(linefiles)
ENDIF ELSE BEGIN
;
;Determine number of line spectrum files
   spawn, 'ls -tr linespectrum_moldata_*.dat', linefiles
   nfiles = N_ELEMENTS(linefiles)
ENDELSE

IF KEYWORD_SET(maxnfile) THEN nfiles = maxnfile

nlines = fltarr(nfiles)
str=''

openr, 1, linefiles[0]
FOR i=1,4 DO readf,1,str
readf,1,dum,nfreq
close,1

lines  = dblarr(LMAX,nfreq)
velos  = dblarr(LMAX,nfreq)
cfreqs = dblarr(LMAX)
lcount = 0L

FOR ff=0,nfiles-1 DO BEGIN
   dum_str = read_line(linefiles[ff])
   dum_nl  = N_ELEMENTS(dum_str.cfreqs)
   cfreqs[lcount:lcount+dum_nl-1]   = dum_str.cfreqs
   FOR i=0,dum_nl-1 DO BEGIN
      velos[lcount+i,*]  = dum_str.velo[*,i]
      lines[lcount+i,*]  = dum_str.flux[*,i]
   ENDFOR
   lcount = lcount + dum_nl
ENDFOR

;
;remove duplicate lines if present
cfreqs = cfreqs[0:lcount-1]

uniqsubs = UNIQ(cfreqs)
cfreqs = cfreqs[uniqsubs]
velos = velos[uniqsubs,*]
lines = lines[uniqsubs,*]
lcount = N_ELEMENTS(uniqsubs)

lines = lines[0:lcount-1,*]
velos = velos[0:lcount-1,*]
vel   = velos[0,*]
res_el = vel[1]-vel[0]
N_vel = LONG(2*Max_vel/res_el+1L)
lines_int = fltarr(lcount, N_vel,2)
new_vel = 2*Max_vel*findgen(N_vel)/(N_vel-1)-Max_vel

c_lines = fltarr(lcount)

specx = dblarr(lcount*N_vel)
specy = dblarr(lcount*N_vel)
 
FOR i=0L,lcount-1 DO BEGIN
    line = fltarr(nfreq)
    line[*] = lines[i,*] 
    vel  = fltarr(nfreq)
    vel[*]  = velos[i,*]

    ;
    ;Calculate the line continuum (assuming a linear function)
    cont  = interpol([line[0],line[nfreq-1]],[vel[0],vel[nfreq-1]],vel)
    ;And save the 0 velocity continuum for later
    c_lines[i] = interpol([line[0],line[nfreq-1]],[vel[0],vel[nfreq-1]],0.) 
    ;And then subtract the continuum
    line  = line - cont

    line_int = interpol([line[0],line,line[nfreq-1]],$
                        [-max_vel,vel,max_vel],new_vel)
    lines_int[i,*,1] = line_int
ENDFOR

;
;Determine minimum and maximum wavelength
gsubs = where(cfreqs NE 0)
max_mu = MAX(c/cfreqs[gsubs])
max_mu = max_mu + max_vel*1d9*max_mu/c
min_mu = MIN(c/cfreqs[gsubs])
min_mu = min_mu - max_vel*1d9*min_mu/c
;
;Define new wavelength grid
;n_x_all = FLOOR((max_mu-min_mu)/ ((min_mu+(max_mu-min_mu)/2d0) * reswidth*1d9/c))
x_all   = dblarr(n_x_max)

x_all[0] = min_mu
n_x_all = 0d0
i       = 0L
WHILE x_all[i] LT max_mu DO BEGIN
   x_all[i+1] = x_all[i] * (1d0 + reswidth / 2.99792458d5)
   n_x_all = n_x_all + 1
   IF n_x_all GT n_x_max THEN BEGIN
      PRINT, 'Too many wavelength points!!!'
      stop
   ENDIF
   i = i+1
ENDWHILE
x_all   = x_all[0:n_x_all-1]
y_all   = fltarr(n_x_all)

FOR i=0,lcount-1 DO BEGIN
   
   x_mu = new_vel*1d9/cfreqs[i]+c/cfreqs[i]
   print,new_vel
   y_Jy = lines_int[i,*,1]
   gsubs = where(x_all ge min(x_mu) and x_all le max(x_mu))   
   y_all[gsubs] = y_all[gsubs] + INTERPOL(y_jy,x_mu,x_all[gsubs])
   
ENDFOR

gauss=exp(-(findgen(3001)-1500.)^2./obsres^2.*2./alog(2.))
ssubs = sort(c/cfreqs)

c_all = interpol(c_lines[ssubs],c/cfreqs[ssubs],x_all)
l_only = y_all * 1d23/dist^2 
;
;Remember to add the continuum back
y_all = (y_all+c_all) * 1d23/dist^2.
;
;and convolve to requested resolving power
y_all = CONVOL(y_all,gauss,total(gauss),/edge_truncate)


IF NOT KEYWORD_SET(xrange) THEN xrange = [MIN(x_all),MAX(x_all)]
rsubs = where(x_all gt xrange[0] and x_all lt xrange[1])

IF NOT KEYWORD_SET(yrange) THEN yrange=[min(y_all[rsubs]),max(y_all[rsubs])]

IF KEYWORD_SET(noise) THEN BEGIN
   ns = randomn(seed,N_ELEMENTS(x_all))*noise
   y_all = y_all+ns
ENDIF

plot, x_all,y_all, xtitle='Wavelength [micron]', ytitle='Flux [Jy]',$
  yrange=yrange,xrange=xrange,psym=psym,/xs,ylog=ylog,xlog=xlog

IF KEYWORD_SET(sedcomp) THEN $
   oplot, x,y/dist^2d0*scale_spec,color=56300

@plot_setup.h
set_plot, 'ps'
device, filename='model.eps',xsize=20,ysize=16,/encapsulated
plot, x_all,y_all, xtitle='Wavelength [micron]', ytitle='Flux [Jy]',$
  yrange=yrange,xrange=xrange,/xs,ylog=ylog,xlog=xlog
device,/close
set_plot, 'x'
@plot_clear.h

mwrfits, dum, 'model.fits',/create
mwrfits, {wave:x_all,spec:y_all,lines:l_only},'model.fits'

END
