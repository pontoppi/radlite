pro write_conf_radical,nr,nt,nf,nref,nspec,onlyray=onlyray,$
        nximage=nximage,ntempdb=ntempdb,ndiff=ndiff
safemode = 1
savemeanint = 0
openw,1,'configure.h'
if keyword_set(onlyray) then begin
   printf,1,'#define ONLY_RAY_TRACING'
endif
if not keyword_set(nximage) then nximage=200
if keyword_set(onlyray) then printf,1,'#define SAVE_MEANINT'
printf,1,'#define RADGRID_TWODIM'
printf,1,'#define INCLUDE_DUST'
printf,1,'#define COORD_SPHERICAL'
printf,1,'#define ESC_MINIMAL_EXTENSION '
if keyword_set(ntempdb) then begin
    printf,1,'#define DB_NTEMP_MAX ',ntempdb
endif
if keyword_set(ndiff) then begin
    printf,1,'#define FRSIZE_NDIFF ',ndiff
endif
printf,1,'#define FRSIZE_X ',nr
printf,1,'#define FRSIZE_Y ',2*(nt)
printf,1,'#define FRSIZE_FREQ ',nf
printf,1,'#define FRSIZE_MU   36'
printf,1,'#define FRSIZE_PHI  16'
printf,1,'#define FRSIZE_CHAR 24'
printf,1,'#define FRSIZE_DRR  ',1+nref
printf,1,'#define FRSIZE_INF_PHI  40'
printf,1,'#define FRSIZE_INF_R    (2*FRSIZE_X+100)'
printf,1,'#define FRSIZE_INF_X ',nximage
printf,1,'#define FRSIZE_INF_Y ',nximage
printf,1,'#define DUST_SPECIES_MAX ',nspec
printf,1,'#define DUST_SIZE_MAX    1'
printf,1,'#define DUST_TRANGE_MAX  1'
printf,1,'#define DUST_OPAC_TEMPDEP'
printf,1,'#define SMALL_MEMORY'
if keyword_set(safemode) then printf,1,'#define CHECK_NUMBERS'
printf,1,'#define MIRROR_THETA'
printf,1,'#define MIRROR_PHI'
printf,1,'#define INTERPOL_PHI_3'
printf,1,'#define INTERPOL_THETA_3'
printf,1,'#define INTERPOL_SRC_1'
printf,1,'#define INTERPOL_FREQ_1'
printf,1,'#define SRCQDR_PNT_2 '
printf,1,'#define ALI_PNT_2'
;printf,1,'#define WRITE_CONVERGENCE'
printf,1,'#define SAFETY_CHECKS_ACTIVE'
printf,1,'#define SAFETY_CHECK_TAUMAX 1.d5'
printf,1,'#define SHORTCHAR_EPS  1.d-10'
printf,1,'#define TELESC_EPS     1.d-10'
printf,1,'#define SHORTCHAR_NEWSTYLE'
printf,1,'#define FRSIZE_SRC_BK  4'
printf,1,'#define MAXRAYNR       256000'
printf,1,'#define FRSIZE_CAM_MU  100'
printf,1,'#define FRSIZE_CAM_PHI 128'
printf,1,'#define CENTRAL_SOURCE'
;printf,1,'#define WRITEINT'
;printf,1,'#define ALI_FILE'
printf,1,'#define DEBUG_FILLINTENS'
printf,1,'#define LINUXVERSION'
printf,1,'#define INCLUDE_EMISQUANT'
printf,1,'#define INCLUDE_QUANTSOURCENEW'
close,1
end



;pro write_conf_vertstruct,nr,nt
;openw,1,'configure.h'
;printf,1,'#define FRSIZE_X ',nr
;printf,1,'#define FRSIZE_Y ',2*(nt)
;printf,1,'#define DUST_SPECIES_MAX 3'
;close,1
;end


;pro write_conf_diffusion,nr,nt
;openw,1,'configure.h'
;printf,1,'#define DIFF_SIZE_X ',nr+10
;printf,1,'#define DIFF_SIZE_Y ',nt+10
;printf,1,'#define DIFF_SIZE_Z 1'
;printf,1,'#define LOW_TAU      (1.d-4)'
;close,1
;end


pro write_conf_dust,nf,nsize,dustdata1,dustdata2
openw,1,'configure.h'
printf,1,'#define DUST_DATABASE_PATH "'+dustdata1+'"'
printf,1,'#define DUST_DATABASE_PATHE "'+dustdata2+'"'
printf,1,'#define FRSIZE_FREQ ',nf+10
printf,1,'#define DUST_SPECIES_MAX 10'
printf,1,'#define DUST_SIZE_MAX ',nsize
close,1
end


pro write_conf_chop,nr,nt,nf,nspec
openw,1,'configure.h'
printf,1,'#define FRSIZE_R ',nr
printf,1,'#define FRSIZE_THETA ',2*(nt)
printf,1,'#define FRSIZE_FREQ ',nf
printf,1,'#define DUST_SPECIES_MAX ',nspec
printf,1,'#define DUST_SIZE_MAX    1'
printf,1,'#define DUST_TRANGE_MAX  1'
close,1
end

pro write_conf_pah,nr,nt,nf,nspec
openw,1,'configure.h'
printf,1,'#define INCLUDE_DUST '
printf,1,'#define QUANTUM_NTEMP 2000'
printf,1,'#define FRSIZE_X ',nr
printf,1,'#define FRSIZE_Y ',2*(nt)
printf,1,'#define FRSIZE_FREQ ',nf
printf,1,'#define DUST_SPECIES_MAX ',nspec
printf,1,'#define DUST_SIZE_MAX    1'
printf,1,'#define DUST_TRANGE_MAX  1'
printf,1,'#define DUST_OPAC_TEMPDEP'
close,1
end
