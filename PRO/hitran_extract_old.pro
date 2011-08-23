PRO hitran_extract, isotopologue=isotopologue,vrange=vrange,Jrange=Jrange,$
                    Branches=Branches,hitran_file=hitran_file


MAX_LINES = 10000
ISOT   = INTARR(MAX_LINES)
FREQ   = DBLARR(MAX_LINES)
XX1    = DBLARR(MAX_LINES)
A      = DBLARR(MAX_LINES)
XX2    = DBLARR(MAX_LINES)
XX3    = DBLARR(MAX_LINES)
Elower = DBLARR(MAX_LINES)
Eupper = DBLARR(MAX_LINES)
XX4    = DBLARR(MAX_LINES)
XX5    = DBLARR(MAX_LINES)
vu     = DBLARR(MAX_LINES)
vl     = DBLARR(MAX_LINES)
Branch = STRARR(MAX_LINES)
Jl     = DBLARR(MAX_LINES)
Ju     = DBLARR(MAX_LINES)
XX6    = DBLARR(MAX_LINES)
XX7    = DBLARR(MAX_LINES)
XX8    = DBLARR(MAX_LINES)
XX9    = DBLARR(MAX_LINES)
XX10   = DBLARR(MAX_LINES)
XX11   = DBLARR(MAX_LINES)
XX12   = DBLARR(MAX_LINES)
gl     = DBLARR(MAX_LINES)
gu     = DBLARR(MAX_LINES)
;
;Set default keywords
;
IF NOT KEYWORD_SET(isotopologue) THEN isotopologue=51
IF NOT KEYWORD_SET(vrange) THEN vrange=[0,0]
IF NOT KEYWORD_SET(Jrange) THEN Jrange=[0,40]
IF NOT KEYWORD_SET(Branches) THEN Branches=['R','P']
IF NOT KEYWORD_SET(hitran_file) THEN hitran_file='05_hit04.par'
;
;
;
CASE isotopologue OF
    51:molweight=28.
    52:molweight=29.
    53:molweight=30.
    54:molweight=29.
    55:molweight=31.
    56:molweight=30.
ENDCASE
;
;HITRAN04 file format
;
FORMAT='(i3,d12.6,e10.3,e10.3,d5.4,d5.4,d11.5,d3.2,d8.6,i15,i15,a21,i3,i12,i2,i2,i2,i2,i2,i2,d8.1,d7.1)'
;FORMAT='(i3,d12.6,e10.3,e10.3,d5.4,d5.4,d10.4,d4.2,d8.6,i15,i15,a21,i3,i12,i2,i2,i2,i2,i2,i2,d8.1,d7.1)'

d1  = 0
d2  = 0.d0
d3  = 0.d0 
d4  = 0.d0
d5  = 0.d0
d6  = 0.d0
d7  = 0.d0
d8  = 0.d0
d9  = 0.d0
d10  = 0
d11 = 0
d12 = ''
d13 = 0
d14 = 0.0
d15 = 0
d16 = 0
d17 = 0
d18 = 0
d19 = 0
d20 = 0
d21 = 0.d0
d22 = 0.d0

i=0
openr,lun, hitran_file,/GET_LUN
WHILE NOT EOF(lun) DO BEGIN
    readf,lun,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11, $
      d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,FORMAT=FORMAT
    ISOT[i]   = d1
    FREQ[i]   = d2
    XX1[i]    = d3
    A[i]      = d4
    XX2[i]    = d5
    XX3[i]    = d6
    Elower[i] = d7
    Eupper[i] = Elower[i]+FREQ[i]
    XX4[i]    = d8
    XX5[i]    = d9
    vu[i]     = d10
    vl[i]     = d11
    Branch[i] = STRTRIM(d12,2)
    Jl[i]     = d13
    XX5[i]    = d14
    XX6[i]    = d15
    XX7[i]    = d16
    XX8[i]    = d17
    XX9[i]    = d18
    XX10[i]   = d19
    XX11[i]   = d20
    gu[i]     = d21
    gl[i]     = d22
    IF Branch[i] EQ 'R' THEN Ju[i] = Jl[i]+1
    IF Branch[i] EQ 'P' THEN Ju[i] = Jl[i]-1
    IF Branch[i] EQ 'Q' THEN Ju[i] = Jl[i]

    i = i+1
ENDWHILE
close,1
free_lun,lun
Nlines = i
;
;Determine selection - chosen parameters + ro-vibrational lines only,
;                      i.e. no purely rotational lines
;
outsubs = WHERE(ISOT EQ isotopologue AND vu GE vrange[0] AND vu LE vrange[1] $
                AND Jl GE Jrange[0] AND Jl LE Jrange[1] AND $
                (Branch EQ Branches[0] OR Branch EQ Branches[1]) AND vl-vu NE 0)
ISOT   = ISOT[outsubs]
FREQ   = FREQ[outsubs]
A      = A[outsubs]
Elower = Elower[outsubs]
Eupper = Eupper[outsubs]
vu     = vu[outsubs]
vl     = vl[outsubs]
Branch = Branch[outsubs]
Jl     = Jl[outsubs]
Ju     = Ju[outsubs]
gl     = gl[outsubs]
gu     = gu[outsubs]
;
;
;
ssubs  = SORT(Elower)
ISOT   = ISOT[ssubs]
FREQ   = FREQ[ssubs]
A      = A[ssubs]
Elower = Elower[ssubs]
Eupper = Eupper[ssubs]
vu     = vu[ssubs]
vl     = vl[ssubs]
Branch = Branch[ssubs]
Jl     = Jl[ssubs]
Ju     = Ju[ssubs]
gl     = gl[ssubs]
gu     = gu[ssubs]
;
;Find unique energy levels
;
E_all = [Elower,Eupper]
v_all = [vl,vu]
J_all = [Jl,Ju]
g_all = [gl,gu]
lsubs    = SORT(E_all)
E_all = E_all[lsubs]
v_all = v_all[lsubs]
J_all = J_all[lsubs]
g_all = g_all[lsubs]
E_levels = [[E_all],[v_all],[J_all],[g_all]]

E_duplic = INTARR(N_ELEMENTS(E_levels[*,0]))
FOR i=1,N_ELEMENTS(E_levels[*,0])-1 DO BEGIN
    IF E_levels[i-1,0] EQ E_levels[i,0] OR (E_levels[i-1,1] EQ E_levels[i,1] AND E_levels[i-1,2] EQ E_levels[i,2]) THEN E_duplic[i]=-1
ENDFOR 
gsubs = WHERE(E_duplic NE -1)
E_levels = E_levels[gsubs,*]

OPENW,lun,'moldata.dat',/get_lun
printf,lun,'!MOLECULE'
printf,lun,'CO'
printf,lun,'!MOLECULAR WEIGHT'
printf,lun,molweight,FORMAT='(f4.1)'
printf,lun,'!NUMBER OF ENERGY LEVELS'
printf,lun,N_ELEMENTS(E_levels[*,0]),FORMAT='(i6)'
printf,lun,'!LEVEL + ENERGIES(cm^-1) + WEIGHT + v + J'
FOR i=0,N_ELEMENTS(E_levels[*,0])-1 DO BEGIN
    printf,lun,i+1,E_levels[i,0],E_levels[i,3],E_levels[i,1],E_levels[i,2],$
      FORMAT='(i5,d11.4,f5.1,i6,i6)'
ENDFOR
printf,lun,'!NUMBER OF RADIATIVE TRANSITIONS'
printf,lun,N_ELEMENTS(FREQ)-1,FORMAT='(i6)'
printf,lun,'!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(cm-1) + E_u(cm-1) + v_u + v_l + J_u + J_l'
FOR i=0,N_ELEMENTS(FREQ)-2 DO BEGIN
    levelu = WHERE(vu[i] EQ E_levels[*,1] AND Ju[i] EQ E_levels[*,2])
    levell = WHERE(vl[i] EQ E_levels[*,1] AND Jl[i] EQ E_levels[*,2])
    printf, lun,i+1,levelu[0]+1,levell[0]+1,A[i],FREQ[i],Eupper[i],vu[i],vl[i],Ju[i],Jl[i],Branch[i],$
      FORMAT='(i5,i5,i5,e12.3,f16.7,f10.2,i3,i3,i3,i3,a3)'
ENDFOR

close,lun

free_lun, lun

END
