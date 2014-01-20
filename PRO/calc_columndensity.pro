@analyze.pro
PRO calc_columndensity,outname=outname
	;
	;This procedure calculates and saves the total column density of the 
	;appropriate molecular species a function of radius.
	;
@natconst.pro
@problem_params.pro
@line_params.ini

if ~KEYWORD_SET(outname) THEN outname = 'columndensity.fits'
	
abunstr = read_abundance(abunfile=abunfile)
densstr = read_dustdens()

nt = abunstr.nt
nr = abunstr.nr
rr = abunstr.rr
tt = abunstr.tt
moldens = abunstr.abun[*,*,0]*densstr.rho[*,nt:nt*2-1]*gastodust/mp/mu

columns = fltarr(nr)

FOR i=0,nr-1 DO BEGIN
	;times 2 to get both sides of the disk
	columns[i] = 2*INT_TABULATED(tt*rr[i],moldens[i,*],/DOUBLE, /SORT)
ENDFOR

outstr = {columndensity:columns,rr:rr,isot:isot}
mwrfits, dummy, outname 
mwrfits, outstr, outname

END