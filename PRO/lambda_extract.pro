;only set for CO now - does not have the ortho/para
;functionality of hitran_extract.pro

PRO lambda_extract, max_energy=max_energy, molfile=molfile,$
                    lambdarange=lambdarange,$
                    vmax=vmax,jmax=jmax, isot=isotop

@natconst.pro
@line_params.ini

IF NOT KEYWORD_SET(molfile) THEN molfile = 'moldata.dat'
IF NOT KEYWORD_SET(max_energy) THEN max_energy = 1d33

FREQRANGE=1d4/lambdarange

;
;Set default keywords
;
IF NOT KEYWORD_SET(isotop) THEN isotop=51
IF NOT KEYWORD_SET(cutoff) THEN cutoff = 0
IF NOT KEYWORD_SET(lambdarange) THEN lambdarange=[12,13]
IF NOT KEYWORD_SET(vmax) THEN vmax = 4
IF NOT KEYWORD_SET(jmax) THEN jmax = 40

;
CASE isotop OF
   51:BEGIN
      molweight=28. 
      molecule='CO'
      lamda_file='12CO_lamda.dat'
   END
ENDCASE
;
;HITRAN08 file format
;
FORMAT='(i3,d12.6,d10.3,d10.3,d5.4,d5.3,d10.4,d4.2,d8.6,a15,a15,a15,a15,6I1,6I2,A1,d7.1,d8.1)'

mol_all = READ_MOLECULE_LAMBDA(main_path+'LAMDA/'+lamda_file,/coll,/ghz)
mol     = EXTRACT_LAMDA(mol_all,vmax=vmax,jmax=jmax)
ntrans  = size(mol.aud,/N_ELEMENTS)
print,mol.iup
ISOT        = INTARR(ntrans)
ISOT        = ISOT + 51
FREQ        = mol.freq
S           = DBLARR(ntrans)
A           = mol.aud
lambda_air  = DBLARR(ntrans)
lambda_self = DBLARR(ntrans)
Elower      = mol.e[mol.iup]
Eupper      = mol.e[mol.idown]
n_air       = DBLARR(ntrans) 
delta_air   = DBLARR(ntrans)
vu          = STRTRIM(STRING(mol.lev_vib[mol.iup]),2)
vl          = STRTRIM(STRING(mol.lev_vib[mol.idown]),2)
Qp          = STRTRIM(STRING(mol.lev_rot[mol.iup]),2)
Qpp         = STRTRIM(STRING(mol.lev_rot[mol.idown]),2)
Ierr        = 0
Iref        = 0
star        = STRARR(ntrans)
gu          = mol.g[mol.iup]
gl          = mol.g[mol.idown]
print,gu
print,freq
print,freqrange
;
;Determine selection based on constraints from line_params.ini
outsubs = WHERE(FREQ LT FREQRANGE[0] AND FREQ GT FREQRANGE[1] AND Eupper LT max_energy)

ISOT   = ISOT[outsubs]
FREQ   = FREQ[outsubs]
A      = A[outsubs]
Elower = Elower[outsubs]
Eupper = Eupper[outsubs]
vu     = vu[outsubs]
vl     = vl[outsubs]
Qp     = Qp[outsubs]
Qpp    = Qpp[outsubs]
gu     = gu[outsubs]
gl     = gl[outsubs]
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
Qp     = Qp[ssubs]
Qpp    = Qpp[ssubs]
gu     = gu[ssubs]
gl     = gl[ssubs]
;
;Combine all energy levels
;
E_all = [Elower,Eupper]
v_all = [vl,vu]
Q_all = [Qpp,Qp]
g_all = [gl,gu]
lsubs    = SORT(E_all)
E_all = E_all[lsubs]
v_all = v_all[lsubs]
Q_all = Q_all[lsubs]
g_all = g_all[lsubs]

;
;Remove duplicate levels - will also take care of degenerate energy levels
E_duplic = INTARR(N_ELEMENTS(E_all))
FOR i=1,N_ELEMENTS(E_all)-1 DO BEGIN
    IF ABS(E_all[i-1]-E_all[i])/(E_all[i]+0.1) LT 0.0001 AND $
       (v_all[i-1] EQ v_all[i]) AND (g_all[i-1] EQ g_all[i]) THEN E_duplic[i]=-1
ENDFOR 

usubs = WHERE(E_duplic NE -1)
E_all = E_all[usubs]
v_all = v_all[usubs]
Q_all = Q_all[usubs]
g_all = g_all[usubs]


OPENW,lun,molfile,/get_lun
printf,lun,'!MOLECULE'
printf,lun,molecule
printf,lun,'!MOLECULAR WEIGHT'
printf,lun,molweight,FORMAT='(f4.1)'
printf,lun,'!NUMBER OF ENERGY LEVELS'
printf,lun,N_ELEMENTS(E_all),FORMAT='(i6)'
printf,lun,'!LEVEL + ENERGIES(cm^-1) + WEIGHT + v + Q'
FOR i=0,N_ELEMENTS(E_all)-1 DO BEGIN
    printf,lun,i+1,E_all[i],g_all[i],v_all[i],Q_all[i],$
      FORMAT='(i5,d11.4,f6.1,a15,a15)'
ENDFOR
printf,lun,'!NUMBER OF RADIATIVE TRANSITIONS'
printf,lun,N_ELEMENTS(FREQ),FORMAT='(i6)'
printf,lun,'!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(cm-1) + E_u(cm-1) + v_u + v_l + Q_p + Q_pp'

FOR i=0,N_ELEMENTS(FREQ)-1 DO BEGIN
    levelu = WHERE(ABS(E_all-Eupper[i])/Eupper[i] LT 0.0001)
    levell = WHERE(ABS(E_all-Elower[i])/(Elower[i]+0.1) LT 0.0001) ;+0.1 in case Elower is the E=0 ground level
   
    printf, lun,i+1,levelu[0]+1,levell[0]+1,A[i],FREQ[i],Eupper[i],vu[i],vl[i],Qp[i],Qpp[i],$
      FORMAT='(i5,i5,i5,e12.3,f16.7,f12.5,a15,a15,a15,a15)'
ENDFOR
close,lun
free_lun, lun



END
