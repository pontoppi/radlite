@analyze.pro	

PRO plot_density,xrange=xrange,yrange=yrange,min_dens=min_dens,max_dens=max_dens,nomirror=nomirror
@natconst

IF ~KEYWORD_SET(xrange) THEN xrange=[0,100]
IF ~KEYWORD_SET(yrange) THEN yrange=[-100,100]
IF ~KEYWORD_SET(min_dens) THEN min_dens = -3
IF ~KEYWORD_SET(max_dens) THEN max_dens = 16

densstr = read_dustdens()
number_dens = alog10(densstr.rho/(2.3*mp))

levels = (max_dens-min_dens)*findgen(150)/149+min_dens

x = fltarr(densstr.nr,densstr.ntheta)
y = fltarr(densstr.nr,densstr.ntheta)

FOR i=0,densstr.nr-1 DO BEGIN
   FOR j=0,densstr.ntheta-1 DO BEGIN
      x[i,j] = densstr.r[i]*cos(!pi/2.-densstr.theta[j])
      y[i,j] = densstr.r[i]*sin(!pi/2.-densstr.theta[j])
   ENDFOR
ENDFOR

IF ~KEYWORD_SET(nomirror) THEN BEGIN
	sh = size(number_dens)
	
	number_dens_mirror = fltarr(sh[1]*2,sh[2])
	number_dens_mirror[0:sh[1]-1,*] = number_dens
	number_dens_mirror[sh[1]:2*sh[1]-1,*] = number_dens

	x_mirror = fltarr(sh[1]*2,sh[2])
	x_mirror[0:sh[1]-1,*] = x
	x_mirror[sh[1]:2*sh[1]-1,*] = -x

	y_mirror = fltarr(sh[1]*2,sh[2])
	y_mirror[0:sh[1]-1,*] = y
	y_mirror[sh[1]:2*sh[1]-1,*] = y
ENDIF ELSE BEGIN
	number_dens_mirror = number_dens
	x_mirror = x
	y_mirror = y
ENDELSE


set_plot, 'ps'
@plot_setup.h
loadct,4
device, filename='density_plot.eps', /encapsulated, /color

contour, number_dens_mirror, x_mirror/AU,y_mirror/AU,/fill,/irregular,levels=levels,yrange=yrange,xrange=xrange,$
		  /xs,/ys,/isotropic

;colorbar, range=[min_lev,max_lev]/1d5,position = [0.8, 0.1, 0.85, 0.4],/vertical,/right,title='!6km/s';,ncolors=50
device, /close
@plot_clear.h
!p.multi=[0,1,1]
set_plot, 'x'
END
