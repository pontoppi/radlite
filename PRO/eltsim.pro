@imcir_conv
@telsim

PRO eltsim

radlite_circim =  'H2O_12.85_May_17_15:14/lineposvelcirc_moldata_0_1.dat'
lsize = 100d0 ;AU
ssize =   2d0 ;AU
specres = 6d0 ;km/s
dist  = 125.d0 ;pc
telescope = 'EELT_LW'
psf = '~/WORK/PAPERS/ELTSIM/elt_psf12.0000.fits'

imcir_conv, radlite_circim,saveit='imcir.fits',$
            /addstar,specres=specres,/specast,xs=lsize,ys=lsize,nx=1000,ny=1000
imcir_conv, radlite_circim,saveit='imcir_s.fits',$
            /addstar,specres=specres,/specast,xs=ssize,ys=ssize,nx=1000,ny=1000
telsim, 'imcir.fits', lineposvel_small='imcir_s.fits', telescope=telescope, psf=psf,dist=dist

END
