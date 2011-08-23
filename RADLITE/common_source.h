c     --------------------------------------------------------------
c                THE GLOBAL COMMONS FOR THE SOURCE TERMS 
c                 INCLUDING THE EMISSION AND SCATTERING
c     --------------------------------------------------------------
#ifndef JSRC_NO_ARRAYS
      common/emission/slte_src,slte_alpha
      doubleprecision slte_src,slte_alpha
      dimension slte_src(1:FRSIZE_FREQ,0:FRSIZE_Y,0:FRSIZE_X)
      dimension slte_alpha(1:FRSIZE_FREQ,0:FRSIZE_Y,0:FRSIZE_X)
c
      common/scatteriso/scati_src,scati_alpha
      doubleprecision   scati_src,scati_alpha
      dimension scati_src(1:FRSIZE_FREQ,0:FRSIZE_Y,0:FRSIZE_X)
      dimension scati_alpha(1:FRSIZE_FREQ,0:FRSIZE_Y,0:FRSIZE_X)
c
      common/scatibk/scatibk_src
      doubleprecision scatibk_src
      dimension scatibk_src(1:FRSIZE_FREQ,0:FRSIZE_Y,
     %                      0:FRSIZE_X,0:FRSIZE_SRC_BK)
c
      common/scatisocum/scaticum_src,scaticum_alpha
      doubleprecision   scaticum_src,scaticum_alpha
      dimension scaticum_src(1:FRSIZE_FREQ,0:FRSIZE_Y,0:FRSIZE_X)
      dimension scaticum_alpha(1:FRSIZE_FREQ,0:FRSIZE_Y,0:FRSIZE_X)
c
      common/aliaccel/ali_accell,ali_approx_lambda
      doubleprecision ali_accell,ali_approx_lambda
      dimension ali_accell(1:FRSIZE_FREQ,0:FRSIZE_Y,0:FRSIZE_X)
      dimension ali_approx_lambda(1:FRSIZE_FREQ,0:FRSIZE_Y,0:FRSIZE_X)
#endif
c
      common/synchwarn/isynch_warn_lowtemp,isynch_warn_lownunub
      integer isynch_warn_lowtemp,isynch_warn_lownunub
c
      common/coolwarn/icool_warn_guessvertheight
      integer icool_warn_guessvertheight
c
      common/scatterisokap/scati_kappa,scaticum_kappa
      doubleprecision scati_kappa,scaticum_kappa
c
      common/ipntsrc/ipsrc_ir,ipsrc_itheta,ipsrc_imu,ipsrc_iphi
      integer ipsrc_ir,ipsrc_itheta,ipsrc_imu,ipsrc_iphi
c
      common/iscatmode/iscat_mode,iscat_donewarn,iscat_warn_int2,
     %                 iscat_warn_src2,iscat_done_nonlte,
     %                 iscat_done_temp,iscat_update_dust_src,
     %                 iscat_no_radtrans
      integer iscat_mode,iscat_donewarn,iscat_done_nonlte
      integer iscat_warn_int2,iscat_warn_src2,iscat_done_temp
      integer iscat_update_dust_src,iscat_no_radtrans
c
      common/alpswitch/asw_lte,asw_scat
      integer asw_lte,asw_scat
c
      common/srcswitch/ssw_lte,ssw_scat
      integer ssw_lte,ssw_scat
c
      common/radprocesses/iradproc_dum_emis,iradproc_brems,
     %     iradproc_synch,iradproc_compt,iradproc_dum_scat,
     %     iradproc_no_emis,iradproc_line,
     %     iradproc_cmptst,iradproc_pseudo_disk,iradproc_dust,
     %     iradproc_line_dust,iradproc_solve_temp,
     %     iradproc_line_starpump
      integer iradproc_dum_emis,iradproc_brems,iradproc_synch
      integer iradproc_compt,iradproc_dum_scat,iradproc_no_emis
      integer iradproc_line,iradproc_cmptst
      integer iradproc_pseudo_disk,iradproc_dust
      integer iradproc_line_dust,iradproc_line_starpump
      integer iradproc_solve_temp
c
      common/thetrans/trns_length_scale
      doubleprecision trns_length_scale(FRSIZE_Y,FRSIZE_X)
      common/ithetrans/trns_escape
      integer trns_escape
c
#ifndef NO_TAU_TEST
      common/tauerror/tauerr_r_taustep,tauerr_t_taustep
      integer tauerr_r_taustep(FRSIZE_FREQ,FRSIZE_Y_SMALL)
      integer tauerr_t_taustep(FRSIZE_FREQ,FRSIZE_X)
      common/itauerror/tauerr_r_ir,tauerr_t_it,
     %                 tauerr_r_err,tauerr_t_err
      integer tauerr_r_ir(FRSIZE_FREQ,FRSIZE_Y_SMALL),tauerr_r_err
      integer tauerr_t_it(FRSIZE_FREQ,FRSIZE_X),tauerr_t_err
#endif
c
#ifdef INCLUDE_COMPTON
      common/isocsk/angsiga,gam1a,wei1a
      real*4 gam1a,wei1a,angsiga
      dimension gam1a(NGAM1,1:FRSIZE_FREQ,1:NRTAUMAX)
      dimension wei1a(NGAM1,1:FRSIZE_FREQ,1:NRTAUMAX)
      dimension angsiga(NGAM1,1:FRSIZE_FREQ,1:NRTAUMAX)
      common/isocskd/isocsk_tau
      doubleprecision isocsk_tau
      dimension isocsk_tau(1:NRTAUMAX)
      common/iisocsk/isocsk_ntau,isocsk_gaindone
      integer isocsk_ntau,isocsk_gaindone
c
      common/warncompt/compt_warn_compt1
      integer compt_warn_compt1
c
      common/comptsolvet/cst_coolingformula
      integer cst_coolingformula
c
      common/comptgain/isocsk_gain
      doubleprecision isocsk_gain(1:FRSIZE_FREQ,1:NRTAUMAX)
#endif
c
      common/qcool/qc_qcool
      doubleprecision qc_qcool(FRSIZE_Y,FRSIZE_X)
c
      common/solvetemp/st_temp_bk
      doubleprecision st_temp_bk(FRSIZE_Y_SMALL,FRSIZE_X)
c
