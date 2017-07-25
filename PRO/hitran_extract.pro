;
;
;Change history:
;
;KMP - August 2011: Added vmax keyword
;

PRO hitran_extract, isotop=isotop,lambdarange=lambdarange,$
                    cutoff=cutoff,FREQ=FREQ,molfile=molfile,$
                    hitran_path=hitran_path, H2O_OP=H2O_OP, max_energy=max_energy,$
                    vmax=vmax,verbose=verbose

IF NOT KEYWORD_SET(molfile) THEN molfile = 'moldata.dat'
IF NOT KEYWORD_SET(hitran_path) THEN hitran_path = './'
IF NOT KEYWORD_SET(max_energy) THEN max_energy = 1d33

c = 2.9979246d14 ;micron/s
MAX_LINES = 3000000L

freqrange = 1d4/lambdarange

ISOT        = INTARR(MAX_LINES)
FREQ        = DBLARR(MAX_LINES)
S           = DBLARR(MAX_LINES)
A           = DBLARR(MAX_LINES)
lambda_air  = DBLARR(MAX_LINES)
lambda_self = DBLARR(MAX_LINES)
Elower      = DBLARR(MAX_LINES)
Eupper      = DBLARR(MAX_LINES)
n_air       = DBLARR(MAX_LINES)
delta_air   = DBLARR(MAX_LINES)
vu          = STRARR(MAX_LINES)
vl          = STRARR(MAX_LINES)
Qp          = STRARR(MAX_LINES)
Qpp         = STRARR(MAX_LINES)
star        = STRARR(MAX_LINES)
gl          = DBLARR(MAX_LINES)
gu          = DBLARR(MAX_LINES)

IF KEYWORD_SET(H2O_OP) THEN BEGIN
   v1up = INTARR(MAX_LINES)
   v2up = INTARR(MAX_LINES)
   v3up = INTARR(MAX_LINES)

   v1do = INTARR(MAX_LINES)
   v2do = INTARR(MAX_LINES)
   v3do = INTARR(MAX_LINES)

   Jup  = INTARR(MAX_LINES)
   KAup = INTARR(MAX_LINES)
   KCup = INTARR(MAX_LINES)

   Jdo  = INTARR(MAX_LINES)
   KAdo = INTARR(MAX_LINES)
   KCdo = INTARR(MAX_LINES)
ENDIF

;
;Set default keywords
;
IF NOT KEYWORD_SET(isotop) THEN isotop=11
IF NOT KEYWORD_SET(cutoff) THEN cutoff = 0
IF NOT KEYWORD_SET(lambdarange) THEN lambdarange=[12,13]
;
CASE isotop OF
   51:BEGIN
      molweight=28. 
      molecule='12CO'
      hitran_file=hitran_path+'05_hit12.par'
   END
   52:BEGIN
      molweight=29. 
      molecule='13CO'
      hitran_file=hitran_path+'05_hit12.par'
   END
   53:BEGIN
      molweight=30. 
      molecule='C18O'
      hitran_file=hitran_path+'05_hit12.par'
   END
   54:BEGIN
      molweight=29. 
      molecule='C17O'
      hitran_file=hitran_path+'05_hit12.par'
   END
   11:BEGIN
      molweight=18. 
      molecule='H2O'
      hitran_file=hitran_path+'01_hit12.par'
   END
   12:BEGIN
      molweight=20. 
      molecule='H218O'
      hitran_file=hitran_path+'01_hit12.par'
   END
   13:BEGIN
      molweight=19. 
      molecule='H217O'
      hitran_file=hitran_path+'01_hit12.par'
   END
   14:BEGIN
      molweight=19. 
      molecule='HDO'
      hitran_file=hitran_path+'01_hit12.par'
   END
   131:BEGIN
      molweight=17.
      molecule='OH'
      hitran_file=hitran_path+'13_HITEMP2010.par'
   END
   231:BEGIN
      molweight=29.
      molecule='HCN'
      hitran_file=hitran_path+'23_hit12.par'
   END
   21:BEGIN
      molweight=48.
      molecule='CO2'
      hitran_file=hitran_path+'02_hit12.par'
   END
   61:BEGIN
      molweight=16.
      molecule='CH4'
      hitran_file=hitran_path+'06_hit12.par'
   END
   111:BEGIN
      molweight=17.
      molecule='NH3'
      hitran_file=hitran_path+'11_hit12.par'
   END
   141:BEGIN
      molweight=15.
      molecule='HF'
      hitran_file=hitran_path+'14_hit12.par'
   END
   201:BEGIN
      molweight=30.
      molecule='H2CO'
      hitran_file=hitran_path+'20_hit12.par'
   END
   221:BEGIN
      molweight=28.
      molecule='N2'
      hitran_file=hitran_path+'22_hit12.par'
   END
   261:BEGIN
      molweight=26.
      molecule='C2H2'
      hitran_file=hitran_path+'26_hit12.par'
   END
   271:BEGIN
      molweight=30.
      molecule='C2H6'
      hitran_file=hitran_path+'27_hit12.par'
   END
   341:BEGIN
      molweight=16.
      molecule='O'
      hitran_file=hitran_path+'34_hit12.par'
   END
   391:BEGIN
      molweight=32.
      molecule='CH3OH'
      hitran_file=hitran_path+'39_hit12.par'
   END
   321:BEGIN
      molweight=46.
      molecule='HCOOH'
      hitran_file=hitran_path+'32_hit12.par'
   END
   381:BEGIN
      molweight=28.
      molecule='C2H4'
      hitran_file=hitran_path+'38_hit12.par'
   END
   422:BEGIN
      molweight=3.
      molecule='HD'
      hitran_file=hitran_path+'HD.par'
   END
   151:BEGIN
      molweight=36.
      molecule='HCL'
      hitran_file=hitran_path+'15_hit12.par'
   END
ENDCASE
;
;HITRAN08 file format
;
FORMAT='(i3,d12.6,d10.3,d10.3,d5.4,d5.3,d10.4,d4.2,d8.6,a15,a15,a15,a15,6I1,6I2,A1,d7.1,d8.1)'


d1  = 0
d2  = 0.d0
d3  = 0.d0 
d4  = 0.d0
d5  = 0.d0
d6  = 0.d0
d7  = 0.d0
d8  = 0.d0
d9  = 0.d0
d10  =''
d11 = ''
d12 = ''
d13 = ''
d14 = intarr(6)
d15 = intarr(6)
d16 = ''
d17 = 1.d0
d18 = 1.d0

i=0L
openr,lun, hitran_file,/GET_LUN
WHILE ~EOF(lun) DO BEGIN
    readf,lun,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11, $
      d12,d13,d14,d15,d16,d17,d18,FORMAT=FORMAT
    ISOT[i]        = d1
    FREQ[i]        = d2
    S[i]           = d3
    A[i]           = d4
    lambda_air[i]  = d5
    lambda_self[i] = d6
    Elower[i]      = d7
    Eupper[i]      = Elower[i]+FREQ[i]
    n_air[i]       = d8
    delta_air[i]   = d9
    vu[i]          = STRTRIM(d10,2)
    vl[i]          = STRTRIM(d11,2)
    Qp[i]          = STRTRIM(d12,2)
    Qpp[i]         = STRTRIM(d13,2)
    Ierr           = d14
    Iref           = d15
    star[i]        = d16
    gu[i]          = d17
    gl[i]          = d18

    ;
    ;If the molecule is water and the user
    ;wishes to extract only ortho or para lines
    IF KEYWORD_SET(H2O_OP) AND isotop EQ 11 THEN BEGIN
       reads,vu[i],d1,d2,d3
       v1up[i] = d1
       v2up[i] = d2
       v3up[i] = d3
       
       reads,vl[i],d1,d2,d3
       v1do[i] = d1
       v2do[i] = d2
       v3do[i] = d3
       
       reads,Qp[i],d1,d2,d3
       Jup[i]  = d1
       KAup[i] = d2
       KCup[i] = d3
       
       reads,Qpp[i],d1,d2,d3
       Jdo[i]  = d1
       KAdo[i] = d2
       KCdo[i] = d3
    ENDIF

    IF d17 EQ 0 THEN BEGIN
		IF KEYWORD_SET(verbose) THEN PRINT, 'Warning: incomplete molecular data in the HITRAN table at frequency: ', freq[i]
    ENDIF ELSE BEGIN
	   i = i+1
    ENDELSE
	   
ENDWHILE
close,1
free_lun,lun
Nlines = i

;
;Determine selection
;

IF KEYWORD_SET(H2O_OP) AND isotop EQ 11 THEN BEGIN
   
   CASE H2O_OP OF 
      'ortho': BEGIN
          outsubs = WHERE(ISOT eq isotop AND FREQ lt FREQRANGE[0] AND FREQ gt FREQRANGE[1] AND $
                         S gt cutoff AND ((KAup+KCup+v3up) mod 2) EQ 1 AND Eupper LT max_energy)
      END
      'para': BEGIN
         outsubs = WHERE(ISOT eq isotop AND FREQ lt FREQRANGE[0] AND FREQ gt FREQRANGE[1] AND $
                         S gt cutoff AND ((KAup+KCup+v3up) mod 2) EQ 0 AND Eupper LT max_energy)
      END
      ELSE: BEGIN
         PRINT, 'Value of H2O_OP must be ortho or para!!'
         STOP
      END
   ENDCASE
ENDIF ELSE BEGIN
   IF NOT KEYWORD_SET(vmax) THEN BEGIN
      outsubs = WHERE(ISOT eq isotop AND FREQ lt FREQRANGE[0] AND FREQ gt FREQRANGE[1] AND S gt cutoff AND Eupper LT max_energy)      
   ENDIF ELSE BEGIN
      outsubs = WHERE(ISOT eq isotop AND FREQ lt FREQRANGE[0] AND FREQ gt FREQRANGE[1] AND S gt cutoff AND Eupper LT max_energy $
                      AND vu LE vmax)
   ENDELSE
ENDELSE

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
;    IF (E_all[i-1] EQ E_all[i]) AND (v_all[i-1] EQ v_all[i]) AND (g_all[i-1] EQ g_all[i]) THEN E_duplic[i]=-1
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
      FORMAT='(i5,d12.4,f7.1,a15,a15)'
ENDFOR
printf,lun,'!NUMBER OF RADIATIVE TRANSITIONS'
printf,lun,N_ELEMENTS(FREQ),FORMAT='(i6)'
printf,lun,'!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(cm-1) + E_u(cm-1) + v_u + v_l + Q_p + Q_pp'

FOR i=0,N_ELEMENTS(FREQ)-1 DO BEGIN
;15/2/2009 BUG FIX. Identifying levels by their
;hitran spectroscopic notation was not unique in some cases.
;Changed to identification by energy.
;
;    levelu = WHERE(vu[i] EQ v_all AND Qp[i] EQ Q_all)
;    levell = WHERE(vl[i] EQ v_all AND Qpp[i] EQ Q_all)
;    levelu = WHERE(Eupper[i] EQ E_all AND Qp[i] EQ Q_all)
;    levell = WHERE(Elower[i] EQ E_all AND Qpp[i] EQ Q_all)
    levelu = WHERE(ABS(E_all-Eupper[i])/Eupper[i] LT 1e-4)
    IF Elower[i] NE 0 THEN BEGIN
       levell = WHERE(ABS(E_all-Elower[i])/Elower[i] LT 1e-4)
    ENDIF ELSE BEGIN
       levell = WHERE(E_all EQ 0.) 
    ENDELSE
    
;    IF (N_ELEMENTS(levelu) GT 1) OR (N_ELEMENTS(levell) GT 1) THEN STOP
    printf, lun,i+1,levelu[0]+1,levell[0]+1,A[i],FREQ[i],Eupper[i],vu[i],vl[i],Qp[i],Qpp[i],$
      FORMAT='(i5,i5,i5,e12.3,f16.7,f12.5,a15,a15,a15,a15)'
ENDFOR

close,lun

free_lun, lun

END
