PRO make_rotdia

lpop = MRDFITS('levelpop_nlte.fits',1)

nradius = N_ELEMENTS(lpop.radius)
ntheta  = N_ELEMENTS(lpop.theta)

loadct,34
set_plot, 'ps'
device, filename='rotdia.eps'
FOR j=0,nradius-1 DO BEGIN
   cgplot, lpop.energy_in_k, ALOG(lpop.npop_all[*,0,j]/lpop.g),psym=16,yrange=[-50,50]
   FOR i=1,ntheta-1 DO BEGIN
      cgplot, lpop.energy_in_k, ALOG(lpop.npop_all[*,i,j]/lpop.g),psym=16,/overplot,color=i*4
   ENDFOR
ENDFOR
device, /close
set_plot, 'x'

END
