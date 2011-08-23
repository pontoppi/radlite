PRO plot_vfield, velfile=velfile,xrange=xrange,yrange=yrange,comps=comps
@natconst
velstr = read_velfield(velfile=velfile)

min_lev = 0. * 1d5
max_lev = 30. * 1d5
levels = (max_lev-min_lev)*findgen(150)/149+min_lev

x=fltarr(velstr.nr,velstr.nt)
y=fltarr(velstr.nr,velstr.nt)

FOR i=0,velstr.nr-1 DO BEGIN
   FOR j=0,velstr.nt-1 DO BEGIN
      x[i,j] = velstr.rr[i]*cos(!pi/2.-velstr.tt[j])
      y[i,j] = velstr.rr[i]*sin(!pi/2.-velstr.tt[j])
   ENDFOR
ENDFOR

set_plot, 'ps'
!p.multi=[0,1,TOTAL(comps)]
@plot_setup.h
loadct,4
device, filename='velfield_plot.eps', /encapsulated, /color,xsize=10,ysize=14 
;contour, REFORM(velstr.vel[*,*,0]), velstr.rr/AU, 90-180*velstr.tt/!pi, /fill, /xl, /yl,nlevels=50
;contour, REFORM(velstr.vel[*,*,1]), velstr.rr/AU, 90-180*velstr.tt/!pi, /fill, /xl, /yl,nlevels=50
;contour, REFORM(velstr.vel[*,*,2]), velstr.rr/AU, 90-180*velstr.tt/!pi, /fill, /xl, /yl,nlevels=50

IF comps[0] THEN BEGIN
   contour, REFORM(velstr.vel[*,*,0]), x/AU,y/AU, /fill,levels=levels,yrange=yrange,xrange=xrange,$
            /xs,/ys,ymargin=[0,1],xmargin=[6,9],xtickname=[' ',' ',' ', ' ',' ',' ',' ',' ']         
   xyouts, (xrange[1]-xrange[0])/2.+xrange[0],  (yrange[1]-yrange[0])*0.9+yrange[0],'!6v!lr!N', color=255,charsize=2 ;, /xl, /yl
ENDIF
IF comps[1] THEN BEGIN
   contour, REFORM(velstr.vel[*,*,1]), x/AU,y/AU, /fill,levels=levels,yrange=yrange,xrange=xrange,$
            /xs,/ys,ymargin=[0,0],xmargin=[6,9],xtickname=[' ',' ',' ', ' ',' ',' ',' ',' ']                      ;, /xl, /yl
   xyouts, (xrange[1]-xrange[0])/2.+xrange[0],  (yrange[1]-yrange[0])*0.9+yrange[0],'!6v!7!lh!N', color=255,charsize=2 ;, /xl, /yl
ENDIF
IF comps[2] THEN BEGIN
   contour, REFORM(velstr.vel[*,*,2]), x/AU,y/AU, /fill,levels=levels,yrange=yrange,xrange=xrange,$
            /xs,/ys,ymargin=[3,0],xmargin=[6,9], xtitle='!6R [AU]', ytitle='!6z [AU]'                             ;, /xl, /yl
   xyouts, (xrange[1]-xrange[0])/2.+xrange[0],  (yrange[1]-yrange[0])*0.9+yrange[0],'!6v!7!lu!N', color=255,charsize=2 ;, /xl, /yl
ENDIF
colorbar, range=[min_lev,max_lev]/1d5,position = [0.8, 0.1, 0.85, 0.4],/vertical,/right,title='!6km/s';,ncolors=50
device, /close
@plot_clear.h
!p.multi=[0,1,1]
set_plot, 'x'
END
