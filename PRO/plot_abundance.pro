PRO plot_abundance, abunfile=abunfile,xrange=xrange,yrange=yrange,comps=comps
@natconst
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
!p.multi=[0,1,TOTAL(comps)]
@plot_setup.h
loadct,4
device, filename='abundance_plot.eps', /encapsulated, /color,xsize=10,ysize=14 

IF comps[0] THEN BEGIN
   contour, REFORM(abunstr.abun[*,*,0]), x/AU,y/AU, /fill,levels=levels,yrange=yrange,xrange=xrange,$
            /xs,/ys,ymargin=[0,1],xmargin=[6,9],xtickname=[' ',' ',' ', ' ',' ',' ',' ',' ']         
   xyouts, (xrange[1]-xrange[0])/2.+xrange[0],  (yrange[1]-yrange[0])*0.9+yrange[0],'!6v!lr!N', color=255,charsize=2 ;, /xl, /yl
ENDIF
IF comps[1] THEN BEGIN
   contour, REFORM(abunstr.abun[*,*,1]), x/AU,y/AU, /fill,levels=levels,yrange=yrange,xrange=xrange,$
            /xs,/ys,ymargin=[0,0],xmargin=[6,9],xtickname=[' ',' ',' ', ' ',' ',' ',' ',' ']                      ;, /xl, /yl
   xyouts, (xrange[1]-xrange[0])/2.+xrange[0],  (yrange[1]-yrange[0])*0.9+yrange[0],'!6v!7!lh!N', color=255,charsize=2 ;, /xl, /yl
ENDIF
IF comps[2] THEN BEGIN
   contour, REFORM(abunstr.abun[*,*,2]), x/AU,y/AU, /fill,levels=levels,yrange=yrange,xrange=xrange,$
            /xs,/ys,ymargin=[3,0],xmargin=[6,9], xtitle='!6R [AU]', ytitle='!6z [AU]'                             ;, /xl, /yl
   xyouts, (xrange[1]-xrange[0])/2.+xrange[0],  (yrange[1]-yrange[0])*0.9+yrange[0],'!6v!7!lu!N', color=255,charsize=2 ;, /xl, /yl
ENDIF
colorbar, range=[min_lev,max_lev]/1d5,position = [0.8, 0.1, 0.85, 0.4],/vertical,/right,title='!6km/s';,ncolors=50
device, /close
@plot_clear.h
!p.multi=[0,1,1]
set_plot, 'x'
END
