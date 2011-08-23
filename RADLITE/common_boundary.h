c     --------------------------------------------------------------
c               THE VARIABLES SPECIFYING THE BOUNDARIES
c     --------------------------------------------------------------
      common/radbound/radbnd_rstar,radbnd_starspec,radbnd_startemp,
     %                radbnd_lstar,radbnd_interstellfield,
     %                radbnd_starspec_line
      doubleprecision radbnd_rstar,radbnd_starspec,radbnd_startemp,
     %                radbnd_lstar,radbnd_interstellfield,
     %                radbnd_starspec_line
      dimension radbnd_starspec(FREQ_START:FRSIZE_FREQ)
      dimension radbnd_starspec_line(FREQ_START:FRSIZE_FREQ)
      dimension radbnd_interstellfield(FREQ_START:FRSIZE_FREQ)
      common/iradbound/iradbnd_in_itype,iradbnd_out_itype,
     %                 iradbnd_in_spectype,iradbnd_read_rstar,
     %                 iradbnd_read_starspec
      integer iradbnd_in_itype,iradbnd_out_itype,
     %                 iradbnd_in_spectype,iradbnd_read_rstar,
     %                 iradbnd_read_starspec
c
      common/equatordisk/eqd_mmstar,eqd_mmdot,eqd_rinn,eqd_rout,
     %                   eqd_xi,eqd_fracnonevap
      doubleprecision eqd_mmstar,eqd_mmdot,eqd_rinn,eqd_rout
      doubleprecision eqd_xi,eqd_fracnonevap
      common/iequatordisk/ieqd_active,ieqd_simple,ieqd_type
      integer ieqd_active,ieqd_simple,ieqd_type
c
      common/radbndbackup/radbnd_cont_starspec,
     %         radbnd_cont_interstellfield
      doubleprecision radbnd_cont_starspec(FREQ_START:FRSIZE_FREQ)
      doubleprecision radbnd_cont_interstellfield
     %                   (FREQ_START:FRSIZE_FREQ)
c
