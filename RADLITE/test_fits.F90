program test_fits

implicit none
real image(1:1024,1:1024),r
integer i, j, status, unit, blocksize,bitpix,naxis,naxes(2)
integer group,fpixel,nelements
logical simple,extend
character filename*80

DO i=1,1024
   DO j=1,1024
      r = SQRT((i-512)**2d0+(j-512)**2d0)
      image(i,j) = dexp(-r**2/30d0**2)
   ENDDO
ENDDO


status = 0
filename = 'test.fits'

call ftgiou(unit,status)

blocksize = 1
call ftinit(unit,filename,blocksize,status)

simple=.true.
bitpix=-32
naxis=2
naxes(1)=1024
naxes(2)=1024
extend=.true.

call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

group  = 1
fpixel = 1
nelements = naxes(1)*naxes(2)
call ftppre(unit,group,fpixel,nelements,image,status)

call ftpkyj(unit,'EXPOSURE',1500,'total exp time',status)

call ftclos(unit,status)
call ftfiou(unit,status)

end program
