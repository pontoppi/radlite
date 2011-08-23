c     --------------------------------------------------------------
c                      THE SHORT CHARACTERISTICS
c     --------------------------------------------------------------
c
#     ifndef NO_SHORT_CHARS
c
#     ifdef SHORCHAR_COMPACT_ARRAYS
      integer FRSIZE_LONGESC
      parameter(FRSIZE_LONGESC=LONGESC_FRACTION*2*FRSIZE_MU_HALF*
     %                          FRSIZE_PHI*FRSIZE_Y*(FRSIZE_DRR+1))
#     endif
c
      common/shortchar/sc_ds,sc_dr,sc_dtheta,sc_phi,sc_mu,
     %                 sc_dphi,sc_dmu,sc_drrp,sc_drrn
      doubleprecision sc_ds,sc_dr,sc_dtheta,sc_phi,sc_mu,
     %                 sc_dphi,sc_dmu,sc_drrp,sc_drrn
      dimension sc_drrp(0:FRSIZE_DRR),sc_drrn(0:FRSIZE_DRR)
      dimension sc_ds(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension sc_dr(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension sc_dtheta(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension sc_phi(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension sc_mu(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension sc_dphi(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension sc_dmu(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      common/ishortchar/isc_icross,isc_nel,isc_longesc,
     %           isc_iphi,isc_imu,isc_ismax,isc_ready,
     %           isc_iangset,isc_longesc_current,isc_maxlen,
     %           isc_shortsc_current,isc_error_longesc,isc_error
      integer isc_icross,isc_nel,isc_iphi,isc_imu,isc_ismax,isc_ready
      integer isc_iangset,isc_longesc,isc_longesc_current
      integer isc_shortsc_current,isc_error_longesc
      integer isc_error,isc_maxlen
      dimension isc_icross(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension isc_nel(0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension isc_iphi(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension isc_imu(0:FRSIZE_CHR1,0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension isc_longesc(0:FRSIZE_PHI,
     %           -FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_Y_SMALL,0:FRSIZE_DRR)
      dimension isc_ismax(0:FRSIZE_DRR)
      dimension isc_iangset(0:FRSIZE_DRR)
c
#     ifdef SHORCHAR_COMPACT_ARRAYS
      common/longesc/esc_ds,esc_dr,esc_dtheta,esc_phi,esc_mu,
     %                 esc_dphi,esc_dmu
      doubleprecision esc_ds,esc_dr,esc_dtheta,esc_phi,esc_mu,
     %                 esc_dphi,esc_dmu
      dimension esc_ds(0:FRSIZE_CHAR,FRSIZE_LONGESC)
      dimension esc_dr(0:FRSIZE_CHAR,FRSIZE_LONGESC)
      dimension esc_dtheta(0:FRSIZE_CHAR,FRSIZE_LONGESC)
      dimension esc_phi(0:FRSIZE_CHAR,FRSIZE_LONGESC)
      dimension esc_mu(0:FRSIZE_CHAR,FRSIZE_LONGESC)
      dimension esc_dphi(0:FRSIZE_CHAR,FRSIZE_LONGESC)
      dimension esc_dmu(0:FRSIZE_CHAR,FRSIZE_LONGESC)
      common/ilongesc/iesc_icross,iesc_nel,iesc_iphi,iesc_imu
      integer iesc_icross,iesc_nel,iesc_iphi,iesc_imu
      dimension iesc_icross(0:FRSIZE_CHAR,FRSIZE_LONGESC)
      dimension iesc_nel(FRSIZE_LONGESC)
      dimension iesc_iphi(0:FRSIZE_CHAR,FRSIZE_LONGESC)
      dimension iesc_imu(0:FRSIZE_CHAR,FRSIZE_LONGESC)
#     endif
c
      common/onesc/osc_ds,osc_dr,osc_dtheta,osc_phi,osc_mu,
     %      osc_dphi,osc_dmu
      doubleprecision osc_ds,osc_dr,osc_dtheta,osc_phi,osc_mu,
     %      osc_dphi,osc_dmu
      dimension osc_ds(0:FRSIZE_CHAR_BIG)
      dimension osc_dr(0:FRSIZE_CHAR_BIG)
      dimension osc_dtheta(0:FRSIZE_CHAR_BIG)
      dimension osc_phi(0:FRSIZE_CHAR_BIG)
      dimension osc_mu(0:FRSIZE_CHAR_BIG)
      dimension osc_dphi(0:FRSIZE_CHAR_BIG)
      dimension osc_dmu(0:FRSIZE_CHAR_BIG)
      common/ionesc/iosc_icross,iosc_nel,iosc_iphi,iosc_imu 
      integer iosc_icross,iosc_nel,iosc_iphi,iosc_imu 
      dimension iosc_icross(0:FRSIZE_CHAR_BIG)
      dimension iosc_iphi(0:FRSIZE_CHAR_BIG)
      dimension iosc_imu(0:FRSIZE_CHAR_BIG)
c
      common/scset/sc_rratiobig
      doubleprecision sc_rratiobig
      common/iscset/sc_ndrr,sc_nstep
cccccccc ERROR (14-05-01)     integer sc_ndrr,sc_nstep(0:FRSIZE_CHAR)
      integer sc_ndrr,sc_nstep(0:FRSIZE_DRR)
c
      common/dbschist/dbsc_hist_nel
      integer dbsc_hist_nel(FRSIZE_CHAR_BIG+2)
c
#     endif /* ifndef NO_SHORT_CHARS */ 
c
