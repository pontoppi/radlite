;Function to extract a subset of lines and molecules from a lamda IDL structure
;
;INPUT: A lamda structure
;OUTPUT: Also a lamda structure, edited to match the constraints

FUNCTION extract_lamda, mol, vmax=vmax, emax=emax, lambdarange=lambdarange

glines  = WHERE(mol.eupper*1.4387973 LT emax) ;emax in Kelvin
giup    = mol.iup[glines]
gidown  = mol.idown[glines]
all_levels = [giup,gidown]
all_levels_sorted = all_levels(SORT(all_levels))
unique_lev = all_levels_sorted[UNIQ(all_levels_sorted)]

nlevels = N_ELEMENTS(unique_lev)

lind        = mol.lind[unique_lev]
energy      = mol.energy[unique_lev]
energy_in_k = mol.energy_in_k[unique_lev]
e           = mol.e[unique_lev]
g           = mol.g[unique_lev]
iup         = mol.iup[glines]
idown       = mol.idown[glines]
aud         = mol.aud[glines]
freq        = mol.freq[glines]
lev_vib     = mol.lev_vib[unique_lev]
lev_rot     = mol.lev_rot[unique_lev]
lin_vib     = mol.lin_vib[glines]
lin_rot     = mol.lin_rot[glines]
eupper      = mol.eupper[glines]
nctrans     = mol.nctrans

collrates   = FLTARR(MAX(mol.nctrans),MAX(mol.ntemps),mol.nr_coll_partners)
coll_iup    = FLTARR(MAX(mol.nctrans),mol.nr_coll_partners)
coll_idown  = FLTARR(MAX(mol.nctrans),mol.nr_coll_partners)

FOR k=0,mol.nr_coll_partners-1 DO BEGIN
   gcoll = INTARR(nctrans[k])
   FOR h=0,nctrans[k]-1 DO BEGIN
      dum = WHERE(mol.coll_iup[h] EQ lind,nup)
      dum = WHERE(mol.coll_idown[h] EQ lind,ndown)
      IF nup GT 0 AND ndown GT 0 THEN gcoll[h] = 1
   ENDFOR
   gsubs = WHERE(gcoll,gcount)
   nctrans[k] = gcount
   coll_iup[0:nctrans[k]-1,k]    = mol.coll_iup[gsubs,k]
   coll_idown[0:nctrans[k]-1,k]  = mol.coll_idown[gsubs,k]
   collrates[0:nctrans[k]-1,0:mol.ntemps[k]-1,k] = mol.collrates[gsubs,0:mol.ntemps[k]-1,k]
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
