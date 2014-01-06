;Function to extract a subset of levels from a lamda IDL structure
;
;INPUT: A lamda structure
;OUTPUT: Also a lamda structure, edited to match the constraints

FUNCTION lamda_extract_levels, mol, vmax=vmax, jmax=jmax, OPR=OPR, temp_type=temp_type

IF ~KEYWORD_SET(OPR) THEN OPR             = 2.0
IF ~KEYWORD_SET(temp_type) THEN temp_type   = 1

CASE temp_type OF
   1: BEGIN
      temp_grid = mol.temps[UNIQ(mol.temps, SORT(mol.temps))] ;Combined set of all temperatures 
      temp_grid = temp_grid[WHERE(temp_grid GT 0)]
      ntemps    = N_ELEMENTS(temp_grid)
   END
ENDCASE

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
   IF (WHERE(mol.iup[i] EQ glevels+1))[0] NE -1 AND (WHERE(mol.idown[i] EQ glevels+1))[0] NE -1 THEN BEGIN
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
   iup[i]   = WHERE(mol.iup[glines[i]]-1 EQ glevels)+1 ;+1 to match the 1-indexed lamda table
   idown[i] = WHERE(mol.idown[glines[i]]-1 EQ glevels)+1
ENDFOR

;
;Do the same thing for the collisional rates
collrates   = FLTARR(MAX(mol.nctrans),MAX(mol.ntemps),mol.nr_coll_partners)
coll_iup    = FLTARR(MAX(mol.nctrans),mol.nr_coll_partners)
coll_idown  = FLTARR(MAX(mol.nctrans),mol.nr_coll_partners)

nctrans_culled     = mol.nctrans
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
   nctrans_culled[k] = gcount ;Update to new number of collisional transitions
ENDFOR

;Free unused rate entries
true_maxcoll = MAX(nctrans_culled)
collrates    = collrates[0:true_maxcoll-1,*,*]
coll_iup     = coll_iup[0:true_maxcoll-1,*]
coll_idown   = coll_idown[0:true_maxcoll-1,*]

FOR p=0,mol.nr_coll_partners-1 DO BEGIN
   ;
   ;get collision partner string according to the LAMDA standard
   pos_left  = STRPOS(mol.partner_name[p], mol.species)+STRLEN(mol.species)+2 
   pos_right = STRPOS(mol.partner_name[p], '!')
   partner_string = STRCOMPRESS(STRMID(mol.partner_name[p], pos_left, pos_right-pos_left),/REMOVE_ALL)
   ;
   ;Scale rates appropriately
   CASE partner_string OF 
      'p-H2' : collrates[*,*,p] *= (1./(OPR+1.))  ;weighted para-H2 rates
      'o-H2' : collrates[*,*,p] *= (OPR/(OPR+1.)) ;weighted ortho-H2 rates
      'H2'   : collrates[*,*,p] *= 1d0            ;do nothing
      'He'   : collrates[*,*,p] *= 0.2            ;cosmic abundance relative to number of H2 molecules
      'H'    : collrates[*,*,p] *= 1d-4           ;collisions with H
      'e-'   : collrates[*,*,p] *= 1d-6           ;collisions with electrons
      ELSE   : BEGIN
         PRINT, 'Unknown collision partner: ', partner_string
         STOP
      END                       ;no other species supported!
   ENDCASE
ENDFOR

;
;Interpolate all rates onto same temperature grid
collrates_int = FLTARR(true_maxcoll,ntemps,mol.nr_coll_partners)
FOR p=0,mol.nr_coll_partners-1 DO BEGIN
   FOR i=0,nctrans_culled[p]-1 DO BEGIN
      collrates_int[i,*,p] = INTERPOL(REFORM(collrates[i,0:mol.ntemps[p]-1,p]),REFORM(mol.temps[0:mol.ntemps[p]-1,p]),temp_grid)  
      ;
      ;Make sure we don't make any crazy extrapolations
      mini = MIN(mol.temps[0:mol.ntemps[p]-1,p])
      maxi = MAX(mol.temps[0:mol.ntemps[p]-1,p])
      lsubs = WHERE(temp_grid LT mini)
      usubs = WHERE(temp_grid GT maxi)
      IF lsubs[0] NE -1 THEN collrates_int[i,lsubs,p] = 0
      IF usubs[0] NE -1 THEN collrates_int[i,usubs,p] = 0
   ENDFOR
ENDFOR


master_collrates  = FLTARR(TOTAL(nctrans_culled),ntemps)
master_coll_iup   = FLTARR(TOTAL(nctrans_culled))
master_coll_idown = FLTARR(TOTAL(nctrans_culled))
count             = 0L

idown_reform = REFORM(coll_idown)
iup_reform   = REFORM(coll_iup)

FOR i=0,nlevels_new-1 DO BEGIN
   FOR j=0,nlevels_new-1 DO BEGIN
      gsubs = WHERE(coll_idown EQ i+1 AND coll_iup EQ j+1,nterms)
      IF nterms GT 0 THEN BEGIN
         index         = ARRAY_INDICES(coll_idown, gsubs)
         index_trans   = REFORM(index[0,*])
         index_partner = REFORM(index[1,*])
         
         FOR k=0,nterms-1 DO master_collrates[count,*] += collrates_int[index[0,k],*,index[1,k]]

         master_coll_idown[count] = i+1
         master_coll_iup[count]   = j+1

         count++
      ENDIF
   ENDFOR
ENDFOR

;
;clean up
collrates_new  = master_collrates[0:count-1,*]
coll_iup_new   = master_coll_iup[0:count-1]
coll_idown_new = master_coll_idown[0:count-1]
nctrans_new    = count

return,{lind:lind,energy:energy,energy_in_K:energy_in_K,e:e,g:g,iup:iup,idown:idown,$
        aud:aud,freq:freq,species:mol.species,mumol:mol.mumol,nlevels:nlevels_new,lev_vib:lev_vib,lev_rot:lev_rot,$
        lin_vib:lin_vib,lin_rot:lin_rot,eupper:eupper,temps:temp_grid,ntemps:ntemps,collrates:collrates_new,$
        coll_iup:coll_iup_new,coll_idown:coll_idown_new,nr_coll_partners:1,$
        partner_name:mol.partner_name,nctrans:nctrans_new}


END
