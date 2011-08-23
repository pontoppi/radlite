c     --------------------------------------------------------------
c               THE VARIABLES SPECIFYING THE MEDIUM 
c     --------------------------------------------------------------
c     
      common/setup/setup_dens_read,setup_temp_read,
     %             setup_dustdens_read,setup_dusttemp_read,
     %             setup_bmag_read,setup_ion_temp_read
      integer setup_dens_read,setup_temp_read,
     %             setup_dustdens_read,setup_dusttemp_read,
     %             setup_bmag_read,setup_ion_temp_read
c
      common/windpar/wp_a,wp_b,wp_pli,wp_rho0,wp_r0,wp_rho_co
      doubleprecision wp_a,wp_b,wp_pli,wp_rho0,wp_r0,wp_rho_co
      common/iwindpar/wp_type,wp_ir0,wp_ir1,wp_it0,wp_it1
      integer wp_type,wp_ir0,wp_ir1,wp_it0,wp_it1
c
      common/medium/medium_rho,medium_el_temp,medium_bmag,
     %              medium_velocity
      doubleprecision medium_rho,medium_el_temp,medium_bmag
      doubleprecision medium_velocity(1:3)
c
      common/mediumarr/medium_arr_rho,medium_arr_el_temp,
     %                 medium_arr_bmag,medium_arr_velocity,
     %                 medium_arr_ion_temp
      doubleprecision medium_arr_rho(0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision medium_arr_el_temp(0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision medium_arr_bmag(0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision medium_arr_velocity(1:3,
     %                            0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision medium_arr_ion_temp(0:FRSIZE_Y_SMALL,0:FRSIZE_X)
c
      common/emisdum/emisdum_kappa
      doubleprecision emisdum_kappa
c
      common/eggnebula/eggneb_belt_kappa_a,eggneb_belt_kappa_s,
     %    eggneb_belt_bb,eggneb_belt_rho,eggneb_blob_kappa_a,
     %    eggneb_blob_kappa_s,eggneb_blob_bb,eggneb_blob_rho,
     %    eggneb_medium_kappa_a,eggneb_medium_kappa_s,
     %    eggneb_medium_rho0
      doubleprecision eggneb_belt_kappa_a,eggneb_belt_kappa_s,
     %    eggneb_belt_bb,eggneb_belt_rho,eggneb_blob_kappa_a,
     %    eggneb_blob_kappa_s,eggneb_blob_bb,eggneb_blob_rho,
     %    eggneb_medium_kappa_a,eggneb_medium_kappa_s,
     %    eggneb_medium_rho0
      common/ieggnebula/eggneb_belt_ir0,eggneb_belt_ir1,
     %    eggneb_belt_it0,eggneb_belt_it1,eggneb_blob_ir0,
     %    eggneb_blob_ir1,eggneb_blob_it0,eggneb_blob_it1
      integer eggneb_belt_ir0,eggneb_belt_ir1,eggneb_belt_it0,
     %    eggneb_belt_it1,eggneb_blob_ir0,eggneb_blob_ir1,
     %    eggneb_blob_it0,eggneb_blob_it1
c
      common/smeggnebula/smeggneb_belt_kappa_a,smeggneb_belt_kappa_s,
     %    smeggneb_belt_a,smeggneb_belt_b,smeggneb_belt_plaw,
     %    smeggneb_belt_co,smeggneb_belt_rco,smeggneb_belt_bb,
     %    smeggneb_belt_rho0,smeggneb_blob_kappa_a,
     %    smeggneb_blob_kappa_s,smeggneb_blob_a,smeggneb_blob_b,
     %    smeggneb_blob_plaw,smeggneb_blob_co,smeggneb_blob_rco,
     %    smeggneb_blob_bb,smeggneb_blob_rho0,smeggneb_medium_kappa_a,
     %    smeggneb_medium_kappa_s,smeggneb_medium_rho0,
     %    smeggneb_medium_plaw,smeggneb_r0,smeggneb_wall_kappa_a,
     %    smeggneb_wall_kappa_s,smeggneb_wall_a,smeggneb_wall_b,
     %    smeggneb_wall_plaw,smeggneb_wall_co,smeggneb_wall_rco,
     %    smeggneb_wall_bb,smeggneb_wall_rho0,smeggneb_wall_theta0,
     %    smeggneb_medium_co,smeggneb_medium_rco,
     %    smeggneb_belt_dust_abun,smeggneb_blob_dust_abun,
     %    smeggneb_wall_dust_abun,smeggneb_medium_dust_abun,
     %    smeggneb_belt_theta0,smeggneb_blob_theta0
      doubleprecision smeggneb_belt_kappa_a,smeggneb_belt_kappa_s,
     %    smeggneb_belt_a,smeggneb_belt_b,smeggneb_belt_plaw,
     %    smeggneb_belt_co,smeggneb_belt_rco,smeggneb_belt_bb,
     %    smeggneb_belt_rho0,smeggneb_blob_kappa_a,
     %    smeggneb_blob_kappa_s,smeggneb_blob_a,smeggneb_blob_b,
     %    smeggneb_blob_plaw,smeggneb_blob_co,smeggneb_blob_rco,
     %    smeggneb_blob_bb,smeggneb_blob_rho0,smeggneb_medium_kappa_a,
     %    smeggneb_medium_kappa_s,smeggneb_medium_rho0,
     %    smeggneb_medium_plaw,smeggneb_r0,smeggneb_wall_kappa_a,
     %    smeggneb_wall_kappa_s,smeggneb_wall_a,smeggneb_wall_b,
     %    smeggneb_wall_plaw,smeggneb_wall_co,smeggneb_wall_rco,
     %    smeggneb_wall_bb,smeggneb_wall_rho0,smeggneb_wall_theta0,
     %    smeggneb_medium_co,smeggneb_medium_rco,
     %    smeggneb_belt_theta0,smeggneb_blob_theta0
#ifdef INCLUDE_DUST
      doubleprecision smeggneb_belt_dust_abun(DUST_SPECIES_MAX)
      doubleprecision smeggneb_blob_dust_abun(DUST_SPECIES_MAX)
      doubleprecision smeggneb_wall_dust_abun(DUST_SPECIES_MAX)
      doubleprecision smeggneb_medium_dust_abun(DUST_SPECIES_MAX)
#else
      doubleprecision 
     %    smeggneb_belt_dust_abun,smeggneb_blob_dust_abun,
     %    smeggneb_wall_dust_abun,smeggneb_medium_dust_abun
#endif
      common/ismeggnebula/smeggneb_solve_temp
      integer smeggneb_solve_temp
c
      common/adaf/adaf_rho_a,adaf_rho_b,adaf_rad_rho_plaw,
     %    adaf_rho0,adaf_rho_co,
     %    adaf_temp_a,adaf_temp_b,adaf_rad_temp_plaw,
     %    adaf_temp0,adaf_temp_co,
     %    adaf_r0,adaf_rho_rco,adaf_temp_rco,adaf_magcoef
      doubleprecision adaf_rho_a,adaf_rho_b,adaf_rad_rho_plaw
      doubleprecision adaf_rho0,adaf_rho_co
      doubleprecision adaf_temp_a,adaf_temp_b,adaf_rad_temp_plaw
      doubleprecision adaf_temp0,adaf_temp_co
      doubleprecision adaf_r0,adaf_rho_rco
      doubleprecision adaf_temp_rco,adaf_magcoef
c
      common/ppatmosphere/ppatm_rho0,ppatm_rho_plaw,ppatm_t0,
     %          ppatm_t1,ppatm_r0,ppatm_r1,ppatm_tcorona
      doubleprecision ppatm_rho0,ppatm_rho_plaw,ppatm_t0,
     %          ppatm_t1,ppatm_r0,ppatm_r1,ppatm_tcorona
c
      common/gassblob/gasblob_rho0,gasblob_rho_dr,gasblob_t0,
     %                gasblob_t_dr,gasblob_t1,gasblob_t_r1,
     %                gasblob_t_smooth,gasblob_t_backgr
      doubleprecision gasblob_rho0,gasblob_rho_dr,gasblob_t0,
     %                gasblob_t_dr,gasblob_t1,gasblob_t_r1,
     %                gasblob_t_smooth,gasblob_t_backgr
c
      common/leunglist/leunglist_rpc,leunglist_tkel,
     %                 leunglist_nh2,leunglist_outflowang
      doubleprecision  leunglist_rpc,leunglist_tkel,
     %                 leunglist_nh2,leunglist_outflowang
c
      common/stellenv/stellenv_r0,stellenv_n0,stellenv_npl,
     %                stellenv_t0,stellenv_tpl,
     %                stellenv_v0,stellenv_vpl,
     %                stellenv_rho0,stellenv_rhopl
      doubleprecision stellenv_r0,stellenv_n0,stellenv_npl,
     %                stellenv_t0,stellenv_tpl,
     %                stellenv_v0,stellenv_vpl,
     %                stellenv_rho0,stellenv_rhopl
c
      common/comptest/ctc_rmax,ctc_taues,ctc_temp,ctc_embg,
     %                ctc_empeak,ctc_ev0,ctc_ev1
      doubleprecision ctc_rmax,ctc_taues,ctc_temp,ctc_embg,
     %                ctc_empeak,ctc_ev0,ctc_ev1
c
      common/ssdadaf/ssdadaf_mdot,ssdadaf_rinn,ssdadaf_rout,
     %               ssdadaf_alpha,ssdadaf_eltemp,ssdadaf_h,
     %               ssdadaf_r0,ssdadaf_Afudge,ssdadaf_Bfudge,
     %               ssdadaf_xi,ssdadaf_magcoef,ssdadaf_gamma,
     %               ssdadaf_f,ssdadaf_iont_vir
      doubleprecision ssdadaf_mdot,ssdadaf_rinn,ssdadaf_rout,
     %               ssdadaf_alpha,ssdadaf_eltemp,ssdadaf_h,
     %               ssdadaf_r0,ssdadaf_Afudge,ssdadaf_Bfudge,
     %               ssdadaf_xi,ssdadaf_magcoef,ssdadaf_gamma,
     %               ssdadaf_f,ssdadaf_iont_vir
      common/issdadaf/ssdadaf_simple
      integer ssdadaf_simple
c
      common/importrt/import_file
      character*80 import_file
c
      common/modust1/modust1_rhofact,modust1_dcontr,modust1_dexp,
     %               modust1_rho_r,modust1_rho_theta
      doubleprecision modust1_rhofact,modust1_dcontr,modust1_dexp,
     %               modust1_rho_r,modust1_rho_theta
      dimension modust1_rho_r(1:FRSIZE_X)
      dimension modust1_rho_theta(1:FRSIZE_Y)
c
      common/accrcloud/accr_cloud_nh2_a,accr_cloud_nh2_b,
     %                 accr_cloud_vr,accr_cloud_vphi,
     %                 accr_cloud_temp0,accr_cloud_temppl,
     %                 accr_cloud_nh20,accr_cloud_nh2pl,
     %                 accr_cloud_mstar,accr_cloud_r0,
     %                 accr_cloud_dustabun,accr_cloud_r1,
     %                 accr_cloud_nh2pl1,accr_cloud_temppl1
      doubleprecision  accr_cloud_nh2_a,accr_cloud_nh2_b,
     %                 accr_cloud_vr,accr_cloud_vphi,
     %                 accr_cloud_temp0,accr_cloud_temppl,
     %                 accr_cloud_nh20,accr_cloud_nh2pl,
     %                 accr_cloud_mstar,accr_cloud_r0,
     %                 accr_cloud_dustabun,accr_cloud_r1,
     %                 accr_cloud_nh2pl1,accr_cloud_temppl1
      common/iaccrcloud/accr_cloud_linedust
      integer accr_cloud_linedust
c
      common/shuconst/shu_asound,shu_rinf,shu_temp_r0,
     %                shu_temp0,shu_temppl
      doubleprecision shu_asound,shu_rinf,shu_temp_r0,
     %                shu_temp0,shu_temppl
c
cc      common/ossenkopf/ossenkopf_r,ossenkopf_nh2,ossenkopf_nh2_pl,
cc     %                 ossenkopf_temp,ossenkopf_temp_pl,
cc     %                 ossenkopf_vrad,ossenkopf_vrad_pl
cc      doubleprecision ossenkopf_r(0:4),ossenkopf_nh2(4),
cc     %                 ossenkopf_nh2_pl(4),
cc     %                 ossenkopf_temp(4),ossenkopf_temp_pl(4),
cc     %                 ossenkopf_vrad(4),ossenkopf_vrad_pl(4)
c
      common/harcalb/ hcb_mdot,hcb_mstar,hcb_tgas,hcb_vturb,
     %                hcb_eta,hcb_rc,hcb_dustabun,hcb_disktemp,
     %                hcb_oned_mu,hcb_racc,hcb_envrhopl
      doubleprecision hcb_mdot,hcb_mstar,hcb_tgas,hcb_vturb,
     %                hcb_eta,hcb_rc,hcb_dustabun,hcb_disktemp,
     %                hcb_oned_mu,hcb_racc,hcb_envrhopl
      common/iharcalb/hcb_linedust,hcb_envextr
      integer hcb_linedust,hcb_envextr
c
      common/flaredisk/fldisk_rin,fldisk_rout,fldisk_hrin,fldisk_hrout,
     %                fldisk_sigin,fldisk_sigout,fldisk_tin,fldisk_tout
      doubleprecision fldisk_rin,fldisk_rout,fldisk_hrin,fldisk_hrout,
     %                fldisk_sigin,fldisk_sigout,fldisk_tin,fldisk_tout
c
      common/zeusadaf/zeusadaf_rho_mult,zeusadaf_b_mult,zeusadaf_t_fix
      doubleprecision zeusadaf_rho_mult,zeusadaf_b_mult,zeusadaf_t_fix
c
