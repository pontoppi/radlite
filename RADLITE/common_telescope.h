c     --------------------------------------------------------------
c                        THE TELESCOPE VARIABLES
c     --------------------------------------------------------------
#define RAYEXPT 100
#define RAYADPT 4
#define RAYRNPT 4
#define RAYSIZE  2*(2*FRSIZE_X+FRSIZE_Y+RAYEXPT) 
#define RAYSMN 2*(2*FRSIZE_X+FRSIZE_Y)
c
c (be sure that RAYEXPT > 2 + 4*RAYADPT*RAYRNPT)
c
      common/raysparam/rp_x0,rp_z0,rp_theta0,rp_s0,rp_b,
     %                 rays_sizepix_x,rays_sizepix_y,
     %                 rays_r,minvel,maxvel
      doubleprecision rp_x0,rp_z0,rp_theta0,rp_s0,rp_b
      doubleprecision rays_sizepix_x,rays_sizepix_y,rays_r
      real minvel, maxvel
      dimension rp_x0(0:MAXRAYNR),rp_z0(0:MAXRAYNR)
      dimension rp_theta0(0:MAXRAYNR),rp_s0(0:MAXRAYNR)
      dimension rp_b(0:MAXRAYNR)
      dimension rays_r(0:FRSIZE_INF_R)
      dimension minvel(0:MAXRAYNR), maxvel(0:MAXRAYNR)
      common/iraysparam/rays_amount,rays_nrphi,rays_nrr,
     %                  rays_ready,rays_nrx,rays_nry
      integer rays_amount,rays_nrphi,rays_ready,rays_nrr
      integer rays_nrx,rays_nry
      common/iraysinfo/rp_iserth0,rp_iserrad,rp_nrrayextra,
     %        rp_irmin,rp_irmax,rp_irstep,rp_dbdr,rp_nrref
      integer rp_iserth0,rp_iserrad,rp_nrrayextra,rp_dbdr
      integer rp_irmin,rp_irmax,rp_irstep,rp_nrref
      dimension rp_iserth0(0:MAXRAYNR),rp_iserrad(0:MAXRAYNR)
      common/raysinfo/rp_thetainf
      doubleprecision rp_thetainf
c
      common/imagecir/imcir_r,imcir_ri,imcir_int,imcir_cont,
     %                imcir_convfunc,imcir_convreduction,
     %                imcir_cmask
      doubleprecision imcir_r,imcir_ri,imcir_int,imcir_cont,
     %                imcir_convfunc,imcir_convreduction
      integer         imcir_cmask
      dimension imcir_r(0:FRSIZE_INF_R)
      dimension imcir_ri(0:FRSIZE_INF_R+1)
      dimension imcir_int(1:FRSIZE_FREQ,
     %         1:FRSIZE_INF_PHI_SMALL,0:FRSIZE_INF_R)
      dimension imcir_cont(1:FRSIZE_INF_PHI_SMALL,0:FRSIZE_INF_R)
      dimension imcir_cmask(1:FRSIZE_FREQ,
     %         1:FRSIZE_INF_PHI_SMALL,0:FRSIZE_INF_R)
      dimension imcir_convfunc(0:FRSIZE_INF_R)
      common/iimageinf/imcir_nrphi,imcir_nrr
      integer imcir_nrphi,imcir_nrr
c
      common/charinttau/char_tau,char_emis,
     %          char_tau_center,char_tau_linecenter
      doubleprecision char_tau,char_emis
      doubleprecision char_tau_center,char_tau_linecenter
#ifdef INCLUDE_LINES
      dimension char_tau_linecenter(SZ_NLINES)
#endif

#ifndef RADGRID_ONEDIM
      common/imagerec/imrec_sizepix_x,imrec_sizepix_y,imrec_int
      doubleprecision imrec_sizepix_x,imrec_sizepix_y,imrec_int
      dimension imrec_int(1:FRSIZE_FREQ,
     %         1:FRSIZE_INF_X,1:FRSIZE_INF_Y)
      common/iimageinf/imrec_nrx,imrec_nry,imrec_addstar,
     %                 imrec_starunres
      integer imrec_nrx,imrec_nry,imrec_addstar,imrec_starunres
#ifndef TELESCOPE_NOIMTAU
      common/imrectau/imrec_tau,imrec_emis
      doubleprecision imrec_tau(1:FRSIZE_FREQ,
     %         1:FRSIZE_INF_X,1:FRSIZE_INF_Y)
      doubleprecision imrec_emis(1:FRSIZE_FREQ,
     %         1:FRSIZE_INF_X,1:FRSIZE_INF_Y)
#endif
#endif
c
      common/spectrum/spec_freq,spec_flux_observer,
     %                spec_beamsize,spec_distance,
     %                spec_lum,spec_aperture
      doubleprecision spec_freq(1:FRSIZE_FREQ),spec_lum
      doubleprecision spec_flux_observer(1:FRSIZE_FREQ)
      
      doubleprecision spec_beamsize,spec_distance
      doubleprecision spec_aperture(1:FRSIZE_FREQ)
      common/ispectrum/ispec_nrfreq,ispec_currentline
      integer ispec_nrfreq,ispec_currentline
c
c      common/intresult/intinf
c      doubleprecision intinf
c      dimension intinf(1:FRSIZE_FREQ,0:FRSIZE_Y)
c
      common/onetraject/tr_radius,tr_theta,tr_s,el_ds,
     %                  tr_mu,tr_phi,tr_b,tr_th0
      doubleprecision tr_radius,tr_theta,tr_s,el_ds,
     %                  tr_mu,tr_phi,tr_b,tr_th0
      dimension tr_radius(0:RAYSIZE)
      dimension tr_theta(0:RAYSIZE)
      dimension tr_s(0:RAYSIZE)
      dimension el_ds(0:RAYSIZE)
      dimension tr_mu(0:RAYSIZE)
      dimension tr_phi(0:RAYSIZE)
      common/ionetraject/tr_amount,tr_icross,tr_iradius,
     %       tr_itheta,el_amount
      integer tr_amount,tr_icross,tr_iradius
      integer tr_itheta,el_amount
      dimension tr_icross(0:RAYSIZE)
      dimension tr_iradius(0:RAYSIZE)
      dimension tr_itheta(0:RAYSIZE)
c
#ifdef TELESCOPE_TRAJARRAY
#if MAXRAYNR>10000 
#error : Too much memory required for rays. Better undefine TELESCOPE_TRAJARRAY
#endif
      common/alltraject/atr_radius,atr_theta,atr_s,ael_ds,
     %                  atr_mu,atr_phi,atr_b,atr_th0
      doubleprecision atr_radius,atr_theta,atr_s,ael_ds,
     %                  atr_mu,atr_phi,atr_b,atr_th0
      dimension atr_radius(0:RAYSIZE,0:MAXRAYNR)
      dimension atr_theta(0:RAYSIZE,0:MAXRAYNR)
      dimension atr_s(0:RAYSIZE,0:MAXRAYNR)
      dimension atr_mu(0:RAYSIZE,0:MAXRAYNR)
      dimension atr_phi(0:RAYSIZE,0:MAXRAYNR)
      dimension ael_ds(0:RAYSIZE,0:MAXRAYNR)
      dimension atr_b(0:MAXRAYNR)
      dimension atr_th0(0:MAXRAYNR)
      common/ialltraject/atr_amount,atr_icross,atr_iradius,
     %       atr_itheta,ael_amount,atr_ready
      integer atr_amount,atr_icross,atr_iradius
      integer atr_itheta,ael_amount,atr_ready
      dimension atr_amount(0:MAXRAYNR)
      dimension atr_icross(0:RAYSIZE,0:MAXRAYNR)
      dimension atr_iradius(0:RAYSIZE,0:MAXRAYNR)
      dimension atr_itheta(0:RAYSIZE,0:MAXRAYNR)
      dimension ael_amount(0:MAXRAYNR)
#endif
c
#ifdef TELESCOPE_CAMERA
      common/camanggrid/cam_rmu,cam_rphi,cam_dmuvdr,
     %                  cam_rmu_i,cam_rphi_i,cam_r,
     %                  cam_theta,cam_brefine,cam_refdmuvdr
      doubleprecision cam_rmu,cam_rphi,cam_dmuvdr
      doubleprecision cam_rmu_i,cam_rphi_i,cam_r,cam_theta
      doubleprecision cam_brefine,cam_refdmuvdr
      dimension cam_rmu(-FRSIZE_CAM_MU:FRSIZE_CAM_MU)
      dimension cam_rphi(-1:FRSIZE_CAM_PHI+2)
      dimension cam_rmu_i(-FRSIZE_CAM_MU-1:FRSIZE_CAM_MU+1)
      dimension cam_rphi_i(-1:FRSIZE_CAM_PHI+3)
      common/camianggrid/cam_nrmu,cam_nrphi,cam_gridtype,
     %            cam_iextrmu,cam_ir,cam_it,cam_ready
      integer cam_nrmu,cam_nrphi,cam_iextrmu,cam_gridtype
      integer cam_ir,cam_it,cam_ready
c
      common/camimage/camim_int
      doubleprecision camim_int
      dimension camim_int(1:FRSIZE_FREQ,1:FRSIZE_CAM_PHI,
     %                    -FRSIZE_CAM_MU:FRSIZE_CAM_MU)
c
#endif
c
      common/tellinespec/linespec_beamsize,linespec_vmax,
     %                   linespec_dv,linespec_distance,
     %                   linespec_incang,linespec_radvelo,
     %                   linespec_dishdiameter,
     %                   linespec_beamradian
      doubleprecision linespec_beamsize,linespec_vmax,
     %                linespec_dv,linespec_distance,
     %                linespec_incang,linespec_radvelo,
     %                linespec_dishdiameter,
     %                linespec_beamradian(FRSIZE_FREQ)
      common/itellinespec/linespec_beamshape,linespec_nlines,
     %                    linespec_ilinestart,linespec_command,
     %                    linespec_style,linespec_flag1
      integer linespec_beamshape,linespec_nlines,
     %                    linespec_ilinestart,linespec_command,
     %                    linespec_style,linespec_flag1
      common/stellinespec/linespec_filename
      character*80 linespec_filename
c
      common/imr/imr_spx,imr_spy,imr_phioff,imr_xoff,imr_yoff,
     %           imr_szimx,imr_szimy
      doubleprecision imr_spx,imr_spy,imr_phioff,imr_xoff,imr_yoff
      doubleprecision imr_szimx,imr_szimy
      common/iimr/imr_nx,imr_ny
      integer imr_nx,imr_ny
c
