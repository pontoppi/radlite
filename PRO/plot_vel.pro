;=====================================================
;Procedure to plot the velocity field of a RADLite model
;
;Klaus Pontoppidan 20/8/07
;
;======================================================

PRO plot_vel

AU = 1.495d13
vel=read_velfield()

loadct,5
!p.multi=[0,2,1]
@plot_setup.h
set_plot, 'ps'
device, filename='velocities.ps',/color,/landscape,xsize=24,ysize=12

;top-down midplane view
xpos = REPLICATE_ARRAY(vel.rr,vel.nr)
ypos = TRANSPOSE(REPLICATE_ARRAY(vel.rr,vel.nr))
rad  = SQRT(xpos^2.+ypos^2.)
absvel  = INTERPOL(vel.vel[*,vel.nt-1,2],vel.rr,rad)
angle = atan(ypos/xpos)
velx = -sin(angle);*absvel
vely = cos(angle);*absvel

partvelvec, velx,vely,xpos/AU,ypos/AU,/isotropic,length=0.03,$
  fraction=1.,veccolors=255*sqrt(absvel)/max(sqrt(absvel)),$
  xtitle='!6AU',ytitle='!6AU'

;r, theta view
;xpos = REPLICATE_ARRAY(vel.rr,vel.nr)
;ypos = fltarr(nr,nt)
;FOR i=0,nr-1 DO BEGIN
;    ypos[i,*] = vel.rr[i]*sin(vel.tt)
;ENDFOR



device,/close
set_plot, 'x'
@plot_clear.h

END
