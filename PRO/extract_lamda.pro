;Function to extract a subset of lines and molecules from a lamda IDL structure
;
;INPUT: A lamda structure
;OUTPUT: Also a lamda structure, edited to match the constraints

FUNCTION extract_lamda, mol, vmax=vmax, jmax=jmax, emax=emax, lambdarange=lambdarange

glines            = WHERE(mol.eupper LT emax,nlines) ;emax in Kelvin
giup              = mol.iup[glines]
gidown            = mol.idown[glines]
all_levels        = [giup,gidown]
all_levels_sorted = all_levels(SORT(all_levels))

lind_nonconsec = all_levels_sorted[UNIQ(all_levels_sorted)] ;This is no longer necessarily a consecutive index
nlevels = N_ELEMENTS(lind_nonconsec)

energy      = FLTARR(nlevels)
energy_in_k = FLTARR(nlevels)
e           = FLTARR(nlevels)
g           = FLTARR(nlevels)
lev_vib     = STRARR(nlevels)
lev_rot     = STRARR(nlevels)

;Extract level information
FOR i=0,nlevels-1 DO BEGIN
   gsub = WHERE(lind_nonconsec[i] EQ mol.lind)
   energy[i]      = mol.energy[gsub]
   energy_in_k[i] = mol.energy_in_k[gsub]
   e[i]           = mol.e[gsub]
   g[i]           = mol.g[gsub]
   lev_vib[i]     = mol.lev_vib[gsub]
   lev_rot[i]     = mol.lev_rot[gsub]
ENDFOR

aud         = mol.aud[glines]
freq        = mol.freq[glines]
lin_vib     = mol.lin_vib[glines]
lin_rot     = mol.lin_rot[glines]
eupper      = mol.eupper[glines]
nctrans     = mol.nctrans

;Now we need to calculate a new set of iup and idown to match a
;consecutive level index, 1,2,3,4,...
lind  = INDGEN(nlevels)+1
iup   = INTARR(nlines)
idown = INTARR(nlines)
FOR i=0,nlines-1 DO BEGIN
   iup[i]   = WHERE(mol.iup[glines[i]] EQ lind_nonconsec)+1 ;+1 to match the 1-indexed lamda table
   idown[i] = WHERE(mol.idown[glines[i]] EQ lind_nonconsec)+1
ENDFOR

collrates   = FLTARR(MAX(mol.nctrans),MAX(mol.ntemps),mol.nr_coll_partners)
coll_iup    = FLTARR(MAX(mol.nctrans),mol.nr_coll_partners)
coll_idown  = FLTARR(MAX(mol.nctrans),mol.nr_coll_partners)

FOR k=0,mol.nr_coll_partners-1 DO BEGIN
   gcount = 0
   FOR h=0,nctrans[k]-1 DO BEGIN
      iupdum   = WHERE(mol.coll_iup[h,k] EQ lind_nonconsec,nup)
      idowndum = WHERE(mol.coll_idown[h,k] EQ lind_nonconsec,ndown)
      IF nup GT 0 AND ndown GT 0 THEN BEGIN
         coll_iup[gcount,k]    = iupdum+1   ;+1 to match the 1-indexed lamda table
         coll_idown[gcount,k]  = idowndum+1
         collrates[gcount,0:mol.ntemps[k]-1,k] = mol.collrates[h,0:mol.ntemps[k]-1,k]
         gcount++
      ENDIF
   ENDFOR
   nctrans[k] = gcount ;Update to new number of collisional transitions
ENDFOR

;Free unused rate entries
true_maxcoll = MAX(nctrans)
collrates    = collrates[0:MAX(nctrans)-1,*,*]
coll_iup     = coll_iup[0:MAX(nctrans)-1,*]
coll_idown   = coll_idown[0:MAX(nctrans)-1,*]

return,{lind:lind,energy:energy,energy_in_K:energy_in_K,e:e,g:g,iup:iup,idown:idown,$
        aud:aud,freq:freq,species:mol.species,mumol:mol.mumol,nlevels:nlevels,lev_vib:lev_vib,lev_rot:lev_rot,$
        lin_vib:lin_vib,lin_rot:lin_rot,eupper:eupper,ntemps:mol.ntemps,temps:mol.temps,collrates:collrates,$
        coll_iup:coll_iup,coll_idown:coll_idown,nr_coll_partners:mol.nr_coll_partners,$
        partner_name:mol.partner_name,nctrans:nctrans}



END
