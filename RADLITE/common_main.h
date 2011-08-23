      common/mainroutine/main_anginf,freq_nu_start,freq_nu_end,
     %     iter_convcrit,main_passbwidth,iter_fluxcons,main_distance,
     %     main_globapert
      doubleprecision main_anginf,freq_nu_start,freq_nu_end
      doubleprecision iter_convcrit,main_passbwidth,iter_fluxcons
      doubleprecision main_distance,main_globapert
c
      common/imainroutine/telesc_command,main_nrphiinf,main_nrrayextra,
     %     image_ifreq,main_passbnfr,telesc_dbdr,radical_quit,iter_trans,
     %     imseq_inu0,imseq_inu1,radical_calcflux,
     %     telesc_nincl,telesc_writespec,
     %     telesc_inostar,telesc_aperture,main_imethod,main_nrref

      integer radical_quit
      integer iter_start,iter_trans,telesc_command,main_nrphiinf
      integer main_nrrayextra,image_ifreq,main_passbnfr,telesc_dbdr
      integer imseq_inu0,imseq_inu1,radical_calcflux
      integer telesc_nincl,telesc_writespec
      integer telesc_inostar,telesc_aperture
      integer main_imethod,main_nrref
c
      common/lmainroutine/iter_converged
      logical iter_converged
c
      common/convergence/cnv_src_error,cnv_temp_error
      doubleprecision cnv_src_error,cnv_temp_error
c
