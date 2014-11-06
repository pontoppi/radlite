FUNCTION read_psum, psumfile, mol

;====================================================================
;READING PARTITION SUMS
;====================================================================
psum_temp = DBLARR(3000)
psum = DBLARR(3000)
OPENR, lunps, psumfile,/get_lun
row = ' ' ;A string
READF, lunps, row ;skip the first line
READF, lunps, row
row = STRCOMPRESS(row)
pos = 0
nspecies = 0
species = STRARR(100)
posits  = INTARR(100)
WHILE pos NE -1 DO BEGIN
   pos = STRPOS(row, ' ', pos+1)
   posits[nspecies+1] = pos
   species[nspecies]  = STRMID(row,posits[nspecies],posits[nspecies+1]-posits[nspecies])
   nspecies++   
ENDWHILE
species = STRCOMPRESS(species, /remove_all)
molmass = FLTARR(nspecies-1)
READF, lunps, row
row = STRCOMPRESS(row)
READS, row, molmass

i=0
floatrow = FLTARR(nspecies)
gsub = WHERE(mol.species EQ species)

IF gsub[0] EQ -1 THEN BEGIN
   PRINT, 'No partition sum table found for the requested molecule!'
   STOP
ENDIF

WHILE NOT EOF(lunps) DO BEGIN
   READF, lunps, row
   row = STRCOMPRESS(row)
   READS, row, floatrow
   psum_temp[i] = floatrow[0]
   psum[i]      = floatrow[gsub+1]
   i++
ENDWHILE

psum_temp = psum_temp[0:i-1] 
psum = psum[0:i-1] 

CLOSE, lunps
FREE_LUN, lunps

RETURN, {temp:psum_temp,psum:psum,molmass:molmass[gsub]}

END

