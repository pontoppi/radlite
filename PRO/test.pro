
t = (findgen(100)+1)/100.*!pi/2.
rm = 2.
r  = rm*sin(t)^2
th = t

velx = -r*sin(th)+2.*sin(th)*cos(th)^2.*rm
vely = r*cos(th)+2.*sin(th)^2.*cos(th)*rm
norm = -sqrt(velx^2+vely^2)


partvelvec, vely/norm,velx/norm, r*sin(th), r*cos(th),/nan, length=0.01,xrange=[0,0.1]
stop
END
