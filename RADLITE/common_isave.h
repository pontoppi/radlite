c     --------------------------------------------------------------
c            THE INFORMATION ABOUT THE ITERATION SAVE LOOP
c     --------------------------------------------------------------
      common/isave/isave_short_char,isave_long_char,
     %    isave_alioper,isave_intens_inu,isave_source,
     %    isave_physvar,isave_scatnonlte,isave_meanint,
     %    isave_dusttemp,isave_levpop,isave_electemp,
     %    isave_fluxcons
      integer isave_short_char,isave_meanint,isave_levpop
      integer isave_long_char,isave_intens_inu,isave_source
      integer isave_physvar,isave_scatnonlte,isave_dusttemp
      integer isave_electemp,isave_alioper,isave_fluxcons
c
      common/icount/icount_dusttemp,icount_meanint,icount_levpop,
     %              icount_electemp,icount_source
      integer icount_dusttemp,icount_meanint,icount_levpop,
     %              icount_electemp,icount_source
c
      common/lastsave/lastsave_dusttemp,lastsave_meanint,
     %     lastsave_levpop,lastsave_electemp,lastsave_source
      integer lastsave_dusttemp,lastsave_meanint,lastsave_levpop,
     %     lastsave_electemp,lastsave_source
c
