FUNCTION read_psum, psumfile, mol

;====================================================================
;READING PARTITION SUMS
;====================================================================
psum_temp = DBLARR(3000)
psum_psum = DBLARR(3000)
openr, lunps, psumfile,/get_lun
str0 = '' & str1 = ''  & str2 = ''  & str3 = ''  & str4 = ''  & str5 = ''  & str6 = ''  & str7 = ''  & str8 = ''  
str9 = '' & str10 = ''  & str11 = ''  & str12 = ''  & str13 = ''  & str14 = ''  & str15 = ''  & str16 = ''  & str17 = ''  & str18 = ''
strdum = ''
i1 = '' & i2 = '' & i3 = '' & i4 = '' & i5 = '' & i6 = '' & i7 = '' & i8 = '' & i9 = ''
i10 = '' & i11 = '' & i12 = '' & i13 = '' & i14 = '' & i15 = '' & i16 = '' & i17 = '' & i18 = ''

readf,lunps,str0,str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17,str18,$ ;skip first line
  format='(a19,a28,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27)'
readf,lunps,str0,str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17,str18,$ 
  format='(a19,a24,a23,a24,a15,a25,a25,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27)'
psum_mols = STRTRIM([str0,str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17,str18],2)
readf,lunps,strdum,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,$ 
  format='(a19,a24,a23,a24,a15,a25,a25,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27,a27)'
molmasses = STRTRIM([strdum,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18])

col_sub = WHERE(psum_mols eq mol.species)
molmass = FLOAT(molmasses[col_sub])

IF col_sub[0] EQ -1 THEN BEGIN
    print, 'No partition sum table found for the requested molecule: ', mol.species
    stop
ENDIF

i=0

WHILE NOT EOF(lunps) DO BEGIN
   readf, lunps, d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18
   psums = [d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18]
   psum_temp[i] = psums[0]
   psum_psum[i] = psums[col_sub]
   i=i+1
ENDWHILE

psum_temp = psum_temp[0:i-1]
psum_psum = psum_psum[0:i-1]

CLOSE, lunps
FREE_LUN, lunps

RETURN, {temp:psum_temp,psum:psum_psum,molmass:molmass}

END
