PRO plot_abundance, abunfile=abunfile,xrange=xrange,yrange=yrange
@natconst

IF ~KEYWORD_SET(xrange) THEN xrange=[0,100]
IF ~KEYWORD_SET(yrange) THEN yrange=[0,100]

abunstr = read_abundance(abunfile=abunfile)

abunstr.abun[*,*,0] = alog10(abunstr.abun[*,*,0])

min_lev = -16
max_lev = 0 
levels = (max_lev-min_lev)*findgen(150)/149+min_lev

x = fltarr(abunstr.nr,abunstr.nt)
y = fltarr(abunstr.nr,abunstr.nt)

FOR i=0,abunstr.nr-1 DO BEGIN
   FOR j=0,abunstr.nt-1 DO BEGIN
      x[i,j] = abunstr.rr[i]*cos(!pi/2.-abunstr.tt[j])
      y[i,j] = abunstr.rr[i]*sin(!pi/2.-abunstr.tt[j])
   ENDFOR
ENDFOR

set_plot, 'ps'
@plot_setup.h
loadct,4
device, filename='abundance_plot.eps', /encapsulated, /color

contour, REFORM(abunstr.abun[*,*,0]), x/AU,y/AU, /fill,levels=levels,yrange=yrange,xrange=xrange,$
         /xs,/ys

;colorbar, range=[min_lev,max_lev]/1d5,position = [0.8, 0.1, 0.85, 0.4],/vertical,/right,title='!6km/s';,ncolors=50
device, /close
@plot_clear.h
!p.multi=[0,1,1]
set_plot, 'x'
END
