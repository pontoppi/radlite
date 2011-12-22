PRO plot_rotdia, popfile=popfile, xrange=xr, yrange=yr, radius=radius

IF ~KEYWORD_SET(popfile) THEN popfile = 'levelpop_nlte.fits'
IF ~KEYWORD_SET(xr) THEN xr=[0,5000]
IF ~KEYWORD_SET(yr) THEN yr=[-10,40]
IF ~KEYWORD_SET(radius) THEN radius=0

pop = MRDFITS(popfile,1)
mol = MRDFITS(popfile,2)
nr = N_ELEMENTS(pop.radius)
nt = N_ELEMENTS(pop.theta)

!p.multi=[0,2,1]
cgloadct, 10
cgplot, [0], xr=xr, yr=yr
FOR j=0,nt-1 DO BEGIN
   cgplot, mol.energy_in_k, ALOG(pop.npop_all[*,j,radius]/mol.g), psym=16, /overplot,color=FLOOR(255.*j/nt)
   cgplot, mol.energy_in_k, ALOG(pop.npop_all[*,j,radius]/mol.g), psym=0, /overplot,color=FLOOR(255.*j/nt)   
ENDFOR

cgplot, [0], xr=xr, yr=yr
FOR j=0,nt-1 DO BEGIN
   cgplot, mol.energy_in_k, ALOG(pop.npop_ini[*,j,radius]/mol.g), psym=16, /overplot,color=FLOOR(255.*j/nt)
   cgplot, mol.energy_in_k, ALOG(pop.npop_ini[*,j,radius]/mol.g), psym=0, /overplot,color=FLOOR(255.*j/nt)   
ENDFOR

stop




END
