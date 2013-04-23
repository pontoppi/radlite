PRO plot_nltepop, popfile=popfile, level=level, xrange=xrange, yrange=yrange, zrange=zrange, log=log

@natconst

IF ~KEYWORD_SET(popfile) THEN popfile = 'levelpop_nlte.fits'
IF ~KEYWORD_SET(level)   THEN level   = 0

pop = MRDFITS(popfile,1)

lsubs = WHERE(pop.npop_all LT 1.)
pop.npop_all[lsubs] = 0.

lte_ratio = REFORM(pop.npop_all[level,*,*]/pop.npop_lte[level,*,*])

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

IF KEYWORD_SET(log) THEN BEGIN
   ratio = ALOG10(ROTATE(lte_ratio,4))
ENDIF ELSE BEGIN
   ratio = ROTATE(lte_ratio,4)
ENDELSE


cgloadct,4
cgwindow
!P.Multi = [0, 2, 1]
levels = findgen(100)/20.-1

cgcontour, ratio, x/AU, y/AU, /fill, levels=levels,xrange=xrange,yrange=yrange,/xs,/ys, $
           xtitle='!6Radius [AU]', ytitle='Height [AU]',/xl,/yl,title='!6Departure from LTE',/addcmd
cgcolorbar, range=[MIN(levels), MAX(levels)],format='(f6.2)', /VERTICAL, $
            POSITION=[0.95/2., 0.10, 0.98/2., 0.90],/addcmd

levels = findgen(100)/10.
cgcontour, ALOG10(ROTATE(REFORM(pop.npop_all[level,*,*]),4)), x/AU, y/AU, /fill, levels=levels,xrange=xrange,yrange=yrange,/xs,/ys, $
           title='!6', xtitle='!6Radius [AU]', ytitle='Height [AU]',/xl,/yl,/addcmd
cgcolorbar, range=[MIN(levels), MAX(levels)],format='(e9.2)', /VERTICAL, $
            POSITION=[0.95, 0.10, 0.98, 0.90],/addcmd

END
