PRO plot_nltepop, popfile=popfile, level=level, xrange=xrange, yrange=yrange

@natconst

IF ~KEYWORD_SET(popfile) THEN popfile = 'levelpop_nlte.fits'
IF ~KEYWORD_SET(level)   THEN level   = 0

pop = MRDFITS(popfile,1)

lte_ratio = REFORM(pop.npop_all[level,*,*]/pop.npop_ini[level,*,*])

nr = N_ELEMENTS(pop.radius)
nt = N_ELEMENTS(pop.theta)
x  = FLTARR(nr,nt)
y  = FLTARR(nr,nt)

FOR i=0,nr-1 DO BEGIN
   FOR j=0,nt-1 DO BEGIN
      x[i,j] = sin(pop.theta[j])*pop.radius[i]
      y[i,j] = cos(pop.theta[j])*pop.radius[i]
   ENDFOR
ENDFOR

log_ratio = ALOG10(ROTATE(lte_ratio,1))

cgloadct,4
cgcontour, log_ratio, x/AU, y/AU, /fill, nlevels=255,xrange=xrange,yrange=yrange,/xs,/ys, $
           title='!6', xtitle='!6Radius [AU]', ytitle='Height [AU]',/xl,/yl
<<<<<<< local
cgcolorbar, range=[Min(log_ratio,/nan), MAX(log_ratio, /nan)],format='(f5.2)', /VERTICAL, $
=======
cgcolorbar, range=[MIN(log_ratio,/nan), MAX(log_ratio,/nan)],format='(f5.2)', /VERTICAL, $
>>>>>>> other
            POSITION=[0.95, 0.10, 0.98, 0.90]
stop
END
