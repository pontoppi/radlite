c     --------------------------------------------------------------
c         THE GLOBAL SPATIAL AND ANGULAR GRIDDING FOR TRANSFER
c     --------------------------------------------------------------
c
c     >>>>  Spatial grid  <<<<
c
      common/irsi/irsi_place,irsi_frsizex,irsi_frsizey,
     %           irsi_inp,irsi_ibnp
      integer irsi_place,irsi_frsizex,irsi_frsizey
      integer irsi_inp,irsi_ibnp
      dimension irsi_place(0:3)
c
      common/rsi/rsi_x_c,rsi_dx_c,rsi_x_i,rsi_dx_i,
     %          rsi_dx_cil,rsi_dx_cir
      doubleprecision rsi_x_c,rsi_dx_c,rsi_x_i,rsi_dx_i
      doubleprecision rsi_dx_cil,rsi_dx_cir
      dimension rsi_x_c(-1:FRSIZE_MAX+2,1:2)
      dimension rsi_dx_c(-1:FRSIZE_MAX+2,1:2)
      dimension rsi_x_i(-1:FRSIZE_MAX+2,1:2)
      dimension rsi_dx_i(-1:FRSIZE_MAX+2,1:2)
      dimension rsi_dx_cil(-1:FRSIZE_MAX+2,1:2)
      dimension rsi_dx_cir(-1:FRSIZE_MAX+2,1:2)
c
c     >>>>  Frequency grid <<<<
c
      common/frequencies/freq_nu
      doubleprecision freq_nu
      dimension freq_nu(1:FRSIZE_FREQ)
      common/ifrequencies/freq_nr,freq_grid_type,freq_done
      integer freq_nr,freq_grid_type,freq_done
c
c     >>>>  Backup frequency grid in case line overwrites dust <<<<
c
      common/icontfr/icont_freq_nr
      integer icont_freq_nr
      common/contfr/cont_freq_nu
      doubleprecision cont_freq_nu(FRSIZE_FREQ)
c
c     >>>>  Global spacegrid info  <<<<
c
      common/spacegrid/spacegrid_metric,
     %    spacegrid_type_x,spacegrid_type_y
      character spacegrid_metric
      character spacegrid_type_x,spacegrid_type_y
      common/cspacegrid/spacegrid_dxi,spacegrid_xi,
     %       spacegrid_xo,spacegrid_equator_eps,
     %       spacegrid_refine_a,spacegrid_refine_thetar,
     %       spacegrid_dxixi
      doubleprecision spacegrid_dxi,spacegrid_xi,
     %       spacegrid_xo,spacegrid_equator_eps,
     %       spacegrid_refine_a,spacegrid_refine_thetar,
     %       spacegrid_dxixi
      common/ispacegrid/spacegrid_refine_n,spacegrid_radius_read,
     %                  spacegrid_theta_read
      integer spacegrid_refine_n,spacegrid_radius_read,
     %                  spacegrid_theta_read
c
c     >>>>  Angular grid  <<<<
c
      common/anggrid/rmu,rphi,dmuvdr,rmu_i,rphi_i
      doubleprecision rmu,rphi,dmuvdr,rmu_i,rphi_i
      dimension rmu(-FRSIZE_MU_HALF:FRSIZE_MU_HALF,0:FRSIZE_DRR)
      dimension rphi(-1:FRSIZE_PHI+2,0:FRSIZE_DRR)
      dimension rmu_i(-FRSIZE_MU_HALF-1:FRSIZE_MU_HALF+1,0:FRSIZE_DRR)
      dimension rphi_i(-1:FRSIZE_PHI+3,0:FRSIZE_DRR)
      common/ianggrid/nrmu,nrphi,iextrmu,iang_warn_mureslow
      integer nrmu,nrphi,iextrmu,iang_warn_mureslow
      dimension nrmu(0:FRSIZE_DRR),nrphi(0:FRSIZE_DRR)
      dimension iextrmu(0:FRSIZE_DRR)
c
      common/anggridinfo/anggrid_drr,anggrid_muerr_max
      doubleprecision anggrid_drr,anggrid_muerr_max
      common/ianggridinfo/anggrid_frsizemu,anggrid_frsizephi,
     %       anggrid_mu_type,anggrid_mu_zero
      integer anggrid_frsizemu,anggrid_frsizephi
      integer anggrid_mu_type,anggrid_mu_zero
c
      common/reindex/ridx_ip,ridx_it,ridx_ready
      integer ridx_ip,ridx_it,ridx_ready
      dimension ridx_ip(-FRSIZE_PHI-4:2*FRSIZE_PHI+4,
     %                    -4:FRSIZE_Y+4)
      dimension ridx_it(-4:FRSIZE_Y+4)
c
#     ifndef NO_SHORT_CHARS
      common/radscindex/ir_iscset_p,ir_iscset_n
      integer ir_iscset_p(1:FRSIZE_X),ir_iscset_n(1:FRSIZE_X)
#     endif /* ifndef NO_SHORT_CHARS */
c
