;Function to extract a subset of lines and molecules from a lamda IDL structure
;
;INPUT: A lamda structure
;OUTPUT: Also a lamda structure, edited to match the constraints

FUNCTION extract_lamda, mol, vmax=vmax, jmax=jmax, lambdarange=lambdarange

OPR     = 2.0 ;H2 ortho-to-para ratio

nlines  = N_ELEMENTS(mol.iup)
glevels = WHERE(mol.lev_vib LE vmax AND mol.lev_rot LE jmax, nlevels_new)

energy      = FLTARR(nlevels_new)
energy_in_k = FLTARR(nlevels_new)
e           = FLTARR(nlevels_new)
g           = FLTARR(nlevels_new)
lev_vib     = INTARR(nlevels_new)
lev_rot     = INTARR(nlevels_new)

;
;Extract level information
FOR i=0,nlevels_new-1 DO BEGIN
   energy[i]      = mol.energy[glevels[i]]
   energy_in_k[i] = mol.energy_in_k[glevels[i]]
   e[i]           = mol.e[glevels[i]]
   g[i]           = mol.g[glevels[i]]
   lev_vib[i]     = mol.lev_vib[glevels[i]]
   lev_rot[i]     = mol.lev_rot[glevels[i]]
ENDFOR

;
;Which lines are still ok?
keep = INTARR(nlines)
FOR i=0,nlines-1 DO BEGIN
   IF (WHERE(glevels+1 EQ mol.iup[i]))[0] NE -1 AND (WHERE(mol.idown[i] EQ glevels))[0] NE -1 THEN BEGIN
      keep[i] = 1
   ENDIF
ENDFOR
glines  = WHERE(keep,nlines_new)
freq    = mol.freq[glines]
aud     = mol.aud[glines]
eupper  = mol.eupper[glines]
lin_vib = mol.lin_vib[glines]
lin_rot = mol.lin_rot[glines]

;Now we need to calculate a new set of iup and idown to match a
;consecutive level index, 1,2,3,4,...
lind  = INDGEN(nlevels_new)+1
iup   = INTARR(nlines_new)
idown = INTARR(nlines_new)
FOR i=0,nlines_new-1 DO BEGIN
   iup[i]   = WHERE(mol.iup[glines[i]] EQ mol.lind[glevels])+1 ;+1 to match the 1-indexed lamda table
   idown[i] = WHERE(mol.idown[glines[i]] EQ mol.lind[glevels])+1
ENDFOR

;
;Do the same thing for the collisional rates
collrates   = FLTARR(MAX(mol.nctrans),MAX(mol.ntemps),mol.nr_coll_partners)
coll_iup    = FLTARR(MAX(mol.nctrans),mol.nr_coll_partners)
coll_idown  = FLTARR(MAX(mol.nctrans),mol.nr_coll_partners)

nctrans_new     = mol.nctrans
FOR k=0,mol.nr_coll_partners-1 DO BEGIN
   gcount = 0
   FOR h=0,mol.nctrans[k]-1 DO BEGIN
      iupdum   = WHERE(mol.coll_iup[h,k] EQ mol.lind[glevels],nup)
      idowndum = WHERE(mol.coll_idown[h,k] EQ mol.lind[glevels],ndown)
      IF nup GT 0 AND ndown GT 0 THEN BEGIN
         coll_iup[gcount,k]    = iupdum+1   ;+1 to match the 1-indexed lamda table
         coll_idown[gcount,k]  = idowndum+1
         collrates[gcount,0:mol.ntemps[k]-1,k] = mol.collrates[h,0:mol.ntemps[k]-1,k]
         gcount++
      ENDIF
   ENDFOR
   nctrans_new[k] = gcount ;Update to new number of collisional transitions
ENDFOR

;Free unused rate entries
true_maxcoll = MAX(nctrans_new)
collrates    = collrates[0:true_maxcoll-1,*,*]
coll_iup     = coll_iup[0:true_maxcoll-1,*]
coll_idown   = coll_idown[0:true_maxcoll-1,*]
master_collrates  = FLTARR(TOTAL(nctrans_new),MAX(mol.ntemps))
master_coll_iup   = FLTARR(TOTAL(nctrans_new))
master_coll_idown = FLTARR(TOTAL(nctrans_new))

count = 0

FOR p=0,mol.nr_coll_partners-1 DO BEGIN
;Identify collisional rates by partner
;for now we consider H2,H, and He

   IF (STRPOS(mol.partner_name[p],'H') NE -1) THEN BEGIN
      IF (STRPOS(mol.partner_name[p],'H2') NE -1) THEN BEGIN
         IF (STRPOS(mol.partner_name[p],'p-H2') NE -1) THEN BEGIN
            collrates[*,*,p] *= (1d0/(OPR+1d0)) ;weighted para-H2 rates
         ENDIF ELSE BEGIN
            IF (STRPOS(mol.partner_name[p],'o-H2') NE -1) THEN BEGIN
               collrates[*,*,p] *= (OPR/(OPR+1d0)) ;weighted ortho-H2 rates
            ENDIF
         ENDELSE
      ENDIF ELSE BEGIN
         IF (STRPOS(mol.partner_name[p],'He') NE -1) THEN BEGIN
            collrates[*,*,p] *= 0.2 ;cosmic abundance relative to number of H2 molecules
         ENDIF ELSE BEGIN
            collrates[*,*,p] *= 2d0 ;collisions with H, assuming density is number of H2 molecules
         ENDELSE
      ENDELSE
   ENDIF
;Make master list of collisional transitions
master_collrates[count:count+nctrans_new[p]-1,*] = collrates[0:nctrans_new[p]-1,*,p]
master_coll_iup[count:count+nctrans_new[p]-1]    = coll_iup[0:nctrans_new[p]-1]
master_coll_idown[count:count+nctrans_new[p]-1]  = coll_idown[0:nctrans_new[p]-1]
count += nctrans_new[p]
ENDFOR
 
collrates  = master_collrates
coll_iup   = master_coll_iup
coll_idown = master_coll_idown
stop
return,{lind:lind,energy:energy,energy_in_K:energy_in_K,e:e,g:g,iup:iup,idown:idown,$
        aud:aud,freq:freq,species:mol.species,mumol:mol.mumol,nlevels:nlevels,lev_vib:lev_vib,lev_rot:lev_rot,$
        lin_vib:lin_vib,lin_rot:lin_rot,eupper:eupper,ntemps:mol.ntemps,temps:mol.temps,collrates:collrates,$
        coll_iup:coll_iup,coll_idown:coll_idown,nr_coll_partners:mol.nr_coll_partners,$
        partner_name:mol.partner_name,nctrans:nctrans_new}



END
