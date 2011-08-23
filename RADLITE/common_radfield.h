c     --------------------------------------------------------------
c              THE DATA SPECIFYING THE RADIATION FIELD
c     --------------------------------------------------------------
#ifndef NO_INTENSITY_STORAGE
      common/intensity/intens
      doubleprecision intens  
      dimension intens(0:FRSIZE_PHI_SMALL,
     %    -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_X)
#endif 
c
#ifndef JSRC_NO_ARRAYS
      common/moment0/intmom_0
      doubleprecision intmom_0
      dimension intmom_0(FREQ_START:FRSIZE_FREQ,
     %                   0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      common/imoment0/intmom_iterscat
      integer intmom_iterscat
      common/cummom/cumul_intmom
      doubleprecision cumul_intmom
      dimension cumul_intmom(FREQ_START:FRSIZE_FREQ,
     %                   0:FRSIZE_Y_SMALL,0:FRSIZE_X)
#endif
c
      common/eddtens/feddij
      doubleprecision feddij
      dimension feddij(1:3,1:3,0:FRSIZE_Y_SMALL,0:FRSIZE_X)
c
      common/transerr/transerr_overshoot,transerr_starresolution
      integer transerr_overshoot,transerr_starresolution
c
#ifndef NO_INTENSITY_STORAGE
#ifdef INTENSITY_ARRAY
#ifdef INTARRAY_FREQ
      common/intensarr/intensity_array
      doubleprecision intensity_array(1:INTARRAY_FREQ,
     %         1:FRSIZE_PHI_SMALL,-FRSIZE_MU_HALF:FRSIZE_MU_HALF,
     %         1:FRSIZE_Y_SMALL,1:FRSIZE_X)
#endif
#endif
#endif
c
      common/radtrwarn/warn_bug_centrbeam,warn_mubin_theta
      integer warn_bug_centrbeam,warn_mubin_theta
c
      common/theluminosity/tlum_lumnu,tlum_lumtot
      doubleprecision tlum_lumnu(FREQ_START:FRSIZE_FREQ,1:FRSIZE_X)
      doubleprecision tlum_lumtot(1:FRSIZE_X)
c     
      common/theflux/fl_hr,fl_ht,fl_divh,fl_divhc,fl_corrfact
      doubleprecision fl_hr(FRSIZE_Y_SMALL,FRSIZE_X)
      doubleprecision fl_ht(FRSIZE_Y_SMALL,FRSIZE_X)
      doubleprecision fl_divh(FRSIZE_Y_SMALL,FRSIZE_X)
      doubleprecision fl_divhc(FRSIZE_Y_SMALL,FRSIZE_X)
      doubleprecision fl_corrfact(FRSIZE_Y_SMALL,FRSIZE_X)
c
