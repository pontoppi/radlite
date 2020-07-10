PRO genspec,obsres=obsres,maxnfile=maxnfile,xrange=xrange,$
              psym=psym,yrange=yrange,xlog=xlog,ylog=ylog,$
              scale_spec=scale_spec,sedcomp=sedcomp,sampling=sampling,$
              multiruns=multiruns,noise=noise,dist=dist,plotit=plotit,outfile=outfile

LMAX = 500000
n_x_max = 1d8
c = 2.99792458d14
obsres = obsres[0] ;make sure this is a scalar

IF NOT KEYWORD_SET(dist) THEN dist = 150. ;pc
IF NOT KEYWORD_SET(obsres) THEN obsres=3.
IF NOT KEYWORD_SET(scale_spec) THEN scale_spec = 1d0
IF NOT KEYWORD_SET(outfile) THEN outfile = 'model'

IF NOT KEYWORD_SET(sampling) THEN BEGIN
   sampling = obsres/2.              ;km/s
ENDIF 

print, 'Spectral resolving power set to: ', STRTRIM(STRING(obsres),2), ' km/s (change with keyword *obsres*)'
print, 'Spectral sampling set to: ', STRTRIM(STRING(sampling),2), ' km/s (change with keyword *sampling*)'

IF KEYWORD_SET(sedcomp) THEN BEGIN
;For testing purposes only
   readcol, '../spectrum.dat',x,y
   x = c/x
   y = y*1d23
ENDIF

IF KEYWORD_SET(multiruns) THEN BEGIN
   FOR i=0,N_ELEMENTS(multiruns)-1 DO BEGIN
      spawn, 'ls '+multiruns[i]+'/linespectrum*.dat',linefiles_subrun
      spawn, 'ls '+multiruns[i]+'/moldata_*.dat',molfiles_subrun

      IF linefiles_subrun[0] EQ '' THEN BEGIN
         PRINT, 'No files in specified subdirectory: ', multiruns[i]
         STOP
      ENDIF
      IF i NE 0 THEN BEGIN
         linefiles = [linefiles,linefiles_subrun]
         molfiles  = [molfiles,molfiles_subrun]
      ENDIF ELSE BEGIN
         linefiles = linefiles_subrun
         molfiles  = molfiles_subrun
      ENDELSE
   ENDFOR
   nfiles = N_ELEMENTS(linefiles)
ENDIF ELSE BEGIN
;
;Determine number of line spectrum files
   spawn, 'ls linespectrum_moldata_*.dat', linefiles
   spawn, 'ls moldata_*.dat', molfiles
   nfiles = N_ELEMENTS(linefiles)
ENDELSE

IF nfiles NE N_ELEMENTS(molfiles) THEN BEGIN
   PRINT, 'number of line files does not equal number of molecular data files!'
   STOP
ENDIF

IF KEYWORD_SET(maxnfile) THEN nfiles = maxnfile

nlines = fltarr(nfiles)
nfreqs = FLTARR(nfiles)
max_vels = FLTARR(nfiles) 


FOR ff=0,nfiles-1 DO BEGIN
	test_read = read_line(linefiles[ff])
	nfreqs[ff]  = test_read.nfreq
	max_vels[ff]  = MAX([3.*ABS(test_read.velo[0,0]),obsres*3.]) ;km/s : extend the box width to this
ENDFOR

max_vel = MAX(max_vels)
highest_res = (WHERE(nfreqs EQ max(nfreqs)))[0]
high_read = read_line(linefiles[highest_res])
nfreq  = high_read.nfreq
velo_master = high_read.velo[*,0]

lines   = dblarr(LMAX,nfreq)
velos   = dblarr(LMAX,nfreq)
cfreqs  = dblarr(LMAX)
species = strarr(LMAX)
trans   = strarr(LMAX)
eupper  = dblarr(LMAX)
aud     = dblarr(LMAX)
gupper  = dblarr(LMAX)
glower  = dblarr(LMAX)


lcount = 0L

FOR ff=0,nfiles-1 DO BEGIN
   dum_mol = read_molecule_lambda(molfiles[ff])
   dum_str = read_line(linefiles[ff])
   dum_nl  = dum_str.nlines
   nfreq_file = dum_str.nfreq

   ;RADLite calculates the line frequency as a the difference between energy levels, but this is not very accurate in HITRAN, so can lead to very imprecise line centers!!!
   ;   cfreqs[lcount:lcount+dum_nl-1] = dum_str.cfreqs   
   cfreqs[lcount:lcount+dum_nl-1] = c/(1e4/dum_mol.freq)

   trans[lcount:lcount+dum_nl-1]  = 'v='+STRCOMPRESS(dum_mol.lin_vib)+$
                                    ' J='+STRCOMPRESS(dum_mol.lin_rot)
   eupper[lcount:lcount+dum_nl-1] = dum_mol.eupper
   aud[lcount:lcount+dum_nl-1]    = dum_mol.aud
   gupper[lcount:lcount+dum_nl-1] = dum_mol.g[dum_mol.iup-1]
   glower[lcount:lcount+dum_nl-1] = dum_mol.g[dum_mol.idown-1]

   FOR i=0,dum_nl-1 DO BEGIN
      velos[lcount+i,0:nfreq-1]  = velo_master
      lines[lcount+i,0:nfreq-1]  = INTERPOL(dum_str.flux[*,i],dum_str.velo[*,i],velo_master)
      species[lcount+i]          = dum_mol.species
   ENDFOR
   lcount = lcount + dum_nl
   
ENDFOR

;
;remove duplicate lines if present
cfreqs  = cfreqs[0:lcount-1]
trans   = trans[0:lcount-1]
eupper  = eupper[0:lcount-1]
aud     = aud[0:lcount-1]
species = species[0:lcount-1]
gupper  = gupper[0:lcount-1]
glower  = glower[0:lcount-1]

total_name = species+trans+string(gupper)+string(glower)+string(eupper)
;uniqsubs  = UNIQ(cfreqs, SORT(cfreqs))
uniqsubs  = UNIQ(total_name, SORT(total_name))
cfreqs    = cfreqs[uniqsubs]
trans   = trans[uniqsubs]
eupper  = eupper[uniqsubs]
aud     = aud[uniqsubs]
species = species[uniqsubs]
gupper  = gupper[uniqsubs]
glower  = glower[uniqsubs]

velos     = velos[uniqsubs,*]
lines     = lines[uniqsubs,*]
print, 'removed', lcount- N_ELEMENTS(uniqsubs), ' duplicate lines'
lcount    = N_ELEMENTS(uniqsubs)

vel       = velos[0,*]
res_el    = vel[1]-vel[0]
N_vel     = LONG(2*Max_vel/res_el+1L)
lines_int = fltarr(lcount, N_vel,2)
new_vel   = 2*Max_vel*findgen(N_vel)/(N_vel-1)-Max_vel

c_lines     = fltarr(lcount)
line_fluxes = fltarr(lcount)
line_widths = fltarr(lcount)
l2cs        = fltarr(lcount)

specx = dblarr(lcount*N_vel)
specy = dblarr(lcount*N_vel)

FOR i=0L,lcount-1 DO BEGIN
    line = fltarr(nfreq)
    line[*] = lines[i,*] 
    vel  = fltarr(nfreq)
    vel[*]  = velos[i,*]

    ;
    ;Calculate the line continuum (assuming a linear function)
    cont       = interpol([line[0],line[nfreq-1]],[vel[0],vel[nfreq-1]],vel)
    ;And save the 0 velocity continuum for later
    c_lines[i] = interpol([line[0],line[nfreq-1]],[vel[0],vel[nfreq-1]],0.) 
    ;And then subtract the continuum
    line       = line - cont

    line_int    = INTERPOL([line[0],line,line[nfreq-1]],$
                           [-max_vel,vel,max_vel],new_vel)
    lines_int[i,*,1] = line_int

    fnu       = cfreqs[i]
    fnus      = (1.+vel*1d9/c)*fnu
    
    line_flux = INT_TABULATED(fnus,line,/DOUBLE,/SORT)/dist^2.
    line_fluxes[i] = line_flux
	
	
	 gsubs = WHERE(abs(line) gt MAX(abs(line)*0.2)) ;Avoid extended line wings
	 line_widths[i] = 2.35482*SQRT(TOTAL(line[gsubs]*vel[gsubs]^2)/TOTAL(line[gsubs]))  ;Moment width -> FWHM
	
	 l2cs[i] = MAX(line)/MEAN(cont)
    
ENDFOR

;
;Determine minimum and maximum wavelength
gsubs = where(cfreqs NE 0)
max_mu = MAX(c/cfreqs[gsubs])
max_mu = max_mu + max_vel*1d9*max_mu/c
min_mu = MIN(c/cfreqs[gsubs])
min_mu = min_mu - max_vel*1d9*min_mu/c
;
;Define wavelength grid with constant velocity resolution
x_all    = dblarr(n_x_max)
x_all[0] = min_mu
n_x_all  = 0d0
i        = 0L
WHILE x_all[i] LT max_mu DO BEGIN
   x_all[i+1] = x_all[i] * (1d0 + res_el / 2.99792458d5)
   n_x_all = n_x_all + 1
   IF n_x_all GT n_x_max THEN BEGIN
      PRINT, 'Too many wavelength points!!!'
      stop
   ENDIF
   i = i+1
ENDWHILE
x_all   = x_all[0:n_x_all-1]
y_all   = fltarr(n_x_all)

;
;Define wavelength grid on requested output velocity sampling
x_out    = dblarr(n_x_max)
x_out[0] = min_mu
n_x_out  = 0d0
i        = 0L
WHILE x_out[i] LT max_mu DO BEGIN
   x_out[i+1] = x_out[i] * (1d0 + sampling / 2.99792458d5)
   n_x_out = n_x_out + 1
   IF n_x_out GT n_x_max THEN BEGIN
      PRINT, 'Too many wavelength points!!!'
      stop
   ENDIF
   i = i+1
ENDWHILE
x_out   = x_out[0:n_x_out-1]

FOR i=0,lcount-1 DO BEGIN
   x_mu = new_vel*1d9/cfreqs[i]+c/cfreqs[i]
   y_Jy = lines_int[i,*,1]
   gsubs = where(x_all ge min(x_mu) and x_all le max(x_mu))   
   y_all[gsubs] = y_all[gsubs] + INTERPOL(y_jy,x_mu,x_all[gsubs])
ENDFOR

ngauss = (CEIL(3.*obsres/res_el))[0]  ;At some IDL version CEIL began returning an array, rather than a scalar. Extremely dangerous.
obsres_sampling = obsres/res_el

gauss = exp(-(FINDGEN(ngauss)-(ngauss-1)/2.)^2./obsres_sampling^2.*2./alog(2.))
ssubs = sort(c/cfreqs)

c_all = INTERPOL(c_lines[ssubs],c/cfreqs[ssubs],x_all)

l_only = y_all * 1d23/dist^2 
;
;Remember to add the continuum back
y_all = (y_all+c_all) * 1d23/dist^2.

;
;and convolve to requested resolving power
y_all = CONVOL(y_all,gauss,total(gauss),/edge_truncate)
l_only = CONVOL(l_only,gauss,total(gauss),/edge_truncate)

;
;Resample to requested output grid
y_out = INTERPOL(y_all,x_all,x_out) ;Total spectrum
l_only = INTERPOL(l_only,x_all,x_out) ;Lines only
c_all = INTERPOL(c_all,x_all,x_out) * 1d23/dist^2. ;Continuum only

IF NOT KEYWORD_SET(xrange) THEN xrange = [MIN(x_out),MAX(x_out)]
rsubs = where(x_out gt xrange[0] and x_out lt xrange[1])

IF NOT KEYWORD_SET(yrange) THEN yrange=[min(y_out[rsubs]),max(y_out[rsubs])]

IF KEYWORD_SET(noise) THEN BEGIN
   ns = randomn(seed,N_ELEMENTS(x_out))*noise
   y_out = y_out+ns
ENDIF

IF KEYWORD_SET(plotit) THEN BEGIN	
	@plot_setup.h
	set_plot, 'ps'
	device, filename=outfile+'.eps',xsize=20,ysize=16,/encapsulated
	plot, x_out,y_out, xtitle='Wavelength [micron]', ytitle='Flux [Jy]',$
  	  yrange=yrange,xrange=xrange,/xs,ylog=ylog,xlog=xlog
	device,/close
	set_plot, 'x'
	@plot_clear.h
	
	plot, x_out,y_out, xtitle='Wavelength [micron]', ytitle='Flux [Jy]',$
  	  yrange=yrange,xrange=xrange,psym=psym,/xs,ylog=ylog,xlog=xlog
	
	IF KEYWORD_SET(sedcomp) THEN $
   	 oplot, x,y/dist^2d0*scale_spec,color=56300
	
ENDIF

mkhdr, header, 2, 0, /extend
sxaddpar, header, 'WAVEUNIT', 'micron'
sxaddpar, header, 'FLUXUNIT', 'Jy'

mwrfits, dum, outfile+'.fits',header,/create
mwrfits, {wave:x_out,spec:y_out,lines:l_only,continuum:c_all},outfile+'.fits'
table = REPLICATE({fluxes:line_fluxes[0],trans:trans[0],species:species[0],eupper:eupper[0],$
    	  		   aud:aud[0],gupper:gupper[0],glower:glower[0],wavelength:c/cfreqs[0],freq:cfreqs[0],width:line_widths[0],l2c:l2cs[0]},$
				   N_ELEMENTS(line_fluxes))
FOR i=1,N_ELEMENTS(line_fluxes)-1 DO BEGIN
	table[i] = {fluxes:line_fluxes[i],trans:trans[i],species:species[i],eupper:eupper[i],$
    	        aud:aud[i],gupper:gupper[i],glower:glower[i],wavelength:c/cfreqs[i],freq:cfreqs[i],width:line_widths[i],l2c:l2cs[i]}
ENDFOR

mwrfits, table,outfile+'.fits'


END
