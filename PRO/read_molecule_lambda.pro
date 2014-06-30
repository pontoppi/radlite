;-----------------------------------------------------------------
;              READ THE LAMBDA MOLECULE DATABASE FORMAT
;-----------------------------------------------------------------
;ghz    - frequencies in data file in GHz (standard lamda)
;kelvin - Line Eupper in Kelvin (standard lamda
;

FUNCTION read_molecule_lambda,file, coll=coll,ghz=ghz,kelvin=kelvin
;
;Fixed parameters
MAXTEMPS = 30
MAXCOLL  = 10000 ;Maximum number of collisional transitions
;
;Definitions
cc     = 2.99792458d10
hh     = 6.62620755d-27
kk     = 1.380658d-16   

;
;Declarations

str     = ''
mumol   = 0.d0
nlev    = 0
nline   = 0

dum   = 0
edum  = 0d0
gdum  = 0d0
vibdum = ' '
rotdum = ' '
llinddum = 0

iupdum   = 0
idowndum = 0
auddum   = 0d0
freqdum  = 0d0

;print,file
openr,lunm,file, /get_lun
readf,lunm,str
readf,lunm,str
species = STRCOMPRESS(str) ;Name of the molecule
readf,lunm,str
readf,lunm,mumol     ;Molecular weight
readf,lunm,str
readf,lunm,nlev,format='(i6)'
readf,lunm,str

e     = dblarr(nlev)
g     = dblarr(nlev)
lind  = intarr(nlev)
lev_vib = strarr(nlev)   ;The vibrational quantum number(s) of the level
lev_rot = strarr(nlev)   ;The rotational quantum number(s)
FOR i=0,nlev-1 DO BEGIN
    readf,lunm,linddum,edum,gdum,vibdum,rotdum,format='(i5,f12.4,f7.1,a15,a15)'
 ;   print,linddum,edum,gdum,vibdum
    e[i]       = edum
    g[i]       = gdum
    lind[i]    = linddum
    lev_vib[i] = vibdum
    lev_rot[i] = rotdum
ENDFOR

readf,lunm,str
readf,lunm,nline
readf,lunm,str

iup   = intarr(nline)
idown = intarr(nline)
aud   = dblarr(nline)
freq  = dblarr(nline)
eupper = dblarr(nline)
lin_vib = strarr(nline)   ;The vibrational quantum number(s) of the line (upper to lower). Can be used as a redundant sanity check with the level numbers
lin_rot = strarr(nline)   ;The rotational quantum number(s) of the line (upper to lower). Can be used as a redundant sanity check with the level numbers

FOR i=0,nline-1 DO BEGIN
   readf,lunm,dum,iupdum,idowndum,auddum,freqdum,eupperdum,vibdum,rotdum,format='(i5,i5,i5,e12.3,f16.7,f12.5,a30,a30)'
 ;  print,lunm,dum,iupdum,idowndum,auddum,freqdum,eupperdum
   iup[i]   = iupdum
   idown[i] = idowndum
   aud[i]   = auddum
   freq[i]  = freqdum
   eupper[i] = eupperdum
   lin_vib[i] = vibdum   
   lin_rot[i] = rotdum
ENDFOR

IF KEYWORD_SET(ghz) THEN freq = freq*1d9/cc ;->cm-1 if originally in GHz (the standard radio LAMDA format).
IF KEYWORD_SET(kelvin) THEN eupper = eupper * kk/(hh*cc) ;->cm-1 if originally in Kelvin (the standard radio LAMDA format).

;
;
IF KEYWORD_SET(coll) THEN BEGIN ;Read collisional rates?
   readf, lunm, str,format='(a100)'
   readf, lunm, nr_coll_partners
   partner_name = STRARR(nr_coll_partners)
   nctrans      = INTARR(nr_coll_partners)
   ntemps       = INTARR(nr_coll_partners)
   temps        = FLTARR(MAXTEMPS,nr_coll_partners)
   collrates    = FLTARR(MAXCOLL,MAXTEMPS,nr_coll_partners)
   coll_iup     = FLTARR(MAXCOLL,nr_coll_partners)
   coll_idown   = FLTARR(MAXCOLL,nr_coll_partners)
   namedum = ' '
   ntrandum = 0
   ntempdum = 0
   FOR k=0,nr_coll_partners-1 DO BEGIN
      readf, lunm, str, format='(a100)'
      readf, lunm, namedum,format='(a100)'
      partner_name[k] = namedum
      readf, lunm, str,format='(a100)'
      readf, lunm, ntrandum,format='(i4)'
      nctrans[k] = ntrandum
      readf, lunm, str,format='(a100)'
      readf, lunm, ntempdum,format='(i4)'
      ntemps[k] = ntempdum
      readf, lunm, str,format='(a100)'
      dumtemps = fltarr(ntemps[k])
      readf, lunm, dumtemps
      temps[0:ntemps[k]-1,k]=dumtemps
      readf, lunm, str,format='(a100)'
      ratedum = FLTARR(ntemps[k])

      FOR h=0,nctrans[k]-1 DO BEGIN
         readf, lunm, indexdum, iupdum, idowndum, ratedum
         collrates[h,0:ntemps[k]-1,k] = ratedum 
         coll_iup[h,k]   = iupdum
         coll_idown[h,k] = idowndum
      ENDFOR
   ENDFOR
   ;Free unused rate entries
   true_maxcoll = MAX(nctrans)
   collrates    = collrates[0:true_maxcoll-1,*,*]
   coll_iup     = coll_iup[0:true_maxcoll-1,*]
   coll_idown   = coll_idown[0:true_maxcoll-1,*]
ENDIF ELSE BEGIN
   ;
   ;These definitions for compatibility when coll=0
   nr_coll_partners = 1
   partner_name     = STRARR(nr_coll_partners)
   nctrans          = INTARR(nr_coll_partners)
   ntemps           = INTARR(nr_coll_partners)
   temps            = FLTARR(nr_coll_partners,MAXTEMPS)
   collrates        = 0
   coll_iup         = 0
   coll_idown       = 0
ENDELSE

close,lunm
free_lun, lunm
energy = e * hh * cc ;energy in K * k (formerly used by RADLite, but now changed to energy_in_K, as below)
energy_in_K = e * hh * cc / kk ;energy in Kelvin

return,{lind:lind,energy:energy,energy_in_K:energy_in_K,e:e,g:g,iup:iup,idown:idown,$
        aud:aud,freq:freq,species:species,mumol:mumol,nlevels:nlev,lev_vib:lev_vib,lev_rot:lev_rot,$
        lin_vib:lin_vib,lin_rot:lin_rot,eupper:eupper,ntemps:ntemps,temps:temps,collrates:collrates,$
        coll_iup:coll_iup,coll_idown:coll_idown,nr_coll_partners:nr_coll_partners,$
        partner_name:partner_name,nctrans:nctrans}
end
