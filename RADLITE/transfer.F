c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c        ( MODULE: PERFORM THE INTEGRATION BY SWEEPING THE GRID ) 
c
c                            Leiden, Aug 1998
c                             C.P. Dullemond   
c
c     =============================================================


#define INTEGR_F_VERSION "Trans_05oct98"
#include "main.h"

c#ifndef NOTRANS
c#ifdef ESC_MINIMAL_EXTENSION
c# error: not allowed in this version
c#endif

#ifndef COORD_SPHERICAL
#ifndef COORD_CARTESIAN
#ifndef COORD_CYLINDRICAL
#   ERROR: Must specify either COORD_SPHERICAL, COORD_CARTESIAN or COORD_CYLINDRICAL
#endif
#endif
#endif

c     --------------------------------------------------------------
c                    MAIN SWITCHBOARD FOR DO TRANS()
c     --------------------------------------------------------------
      subroutine do_trans(inu)
      implicit none
      integer inu
c
#include "common_source.h"
c
cc########################
c      trns_escape=1
cc########################
      if(trns_escape.eq.0) then
          call do_trans_char(inu)
      else
          call do_trans_escapeprob(inu)
      endif
c
      end



c     --------------------------------------------------------------
c                     ESCAPE-PROBABILITY TRANSFER
c     --------------------------------------------------------------
      subroutine do_trans_escapeprob(inu)
      implicit none
      integer inu
c
#     ifndef NO_INTENSITY_STORAGE
#     ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
c
      integer ir,it,itmax,iangset,iphimax,imu,iphi
      doubleprecision tau,alpha,src,meanint
c
c     Send stern warning, to make sure that this mode is not
c     accidently used.
c
      if(inu.eq.1) then
          write(*,*) '--- ESCAPE PROBABILITY TRANSFER ACTING ---'
      endif
      iangset  = 1
c
c     Some computations
c
      iangset = 1
      iphimax = nrphi(iangset)
      itmax   = irsi_frsizey
c
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#     endif
#     ifdef MIRROR_PHI
      iphimax = nrphi(iangset) / 2
#     endif
#     ifdef RADGRID_ONEDIM
      iphimax = 1
      itmax = 1
#     endif
c
c     First a loop over all the grid positions
c
      do ir=1,irsi_frsizex
          do it=1,itmax
c
c             Find the intensity by escape probability
c
              alpha   = scati_alpha(inu,it,ir) + slte_alpha(inu,it,ir)
              src     = scati_src(inu,it,ir) + slte_alpha(inu,it,ir)
              tau     = trns_length_scale(it,ir) * alpha
              if(tau.le.0.d0) then
                  write(*,*) 'ERROR: Escape probability method needs'
                  write(*,*) '       a length scale array'
                  stop 13
              endif
              meanint = (src/alpha) * ( exp(-tau) - 1.d0 + tau) / tau
c
c             Put the mean intensity into the full intensity, since
c             this is what this routine is supposed to deliver.
c
              if(anggrid_mu_zero.eq.1) then
                  do iphi=1,iphimax
                      intens(iphi,0,it,ir)    = meanint
                  enddo
              endif
              do imu=1,nrmu(iangset)
                  do iphi=1,iphimax
                      intens(iphi,imu,it,ir)  = meanint
                      intens(iphi,-imu,it,ir) = meanint
                  enddo
              enddo
          enddo
      enddo
c
#     endif /* ifndef NO_ARRAYS */
c
#     else /* ifndef NO_INTENSITY_STORAGE */
      write(*,*) 'ERROR: do_trans_escapeprob() has been called.'
      write(*,*) '       This function needs the intensity array.'
      write(*,*) '       But this array is not present, because'
      write(*,*) '       NO_INTENSITY_STORAGE is defined.'
      write(*,*) '       Tip: undefine NO_INTENSITY_STORAGE in'
      write(*,*) '       configure.h'
      stop 13
#     endif /* ifndef NO_INTENSITY_STORAGE */
c
      end


c     --------------------------------------------------------------
c          PERFORM FORMAL TRANSFER USING SHORT CHARACTERISTICS
c
c     PURPOSE:
c     This is the main routine for the formal transfer in 2D using
c     short characteristics. Calling this routine for a specific
c     frequency (given by the frequency index inu) will fill the 
c     array intensity(*,*,*,*) with the entire radiation field at
c     all gridpoints and all angles at this specific frequency.
c        Once this is done, this routine returns. Then the user
c     should call the routine store_radfield(), which will copy the
c     data to another array with lower angular resolution, and with
c     frequency dependence. The lower angular resolution is required
c     because otherwise storage becomes a problem for many frequencies.
c     The high resolution for the transfer is necessary because 
c     otherwise the spatial dependence of the radiation field will
c     be smeared out by numerical diffusion.
c        After store_radfield(), the user can call this subroutine
c     again to continue radiative transfer with a new frequency. After
c     repeating this cycle until all frequencies are done, the 
c     scattering integrals and atomic/molecular populations can be 
c     computed. This is of course where the short characteristic method
c     is designed for in the first place.
c
c     HOW TO USE:
c      - At the start of the program (and definitely before calling this
c        function) the short characteristics should be set up. This is
c        done by calling the subroutine create_shortchar(). This routine
c        will compute and prepare the short characteristics for all the
c        points on the grid, so that do_trans() only needs to look up
c        this information in a table, rather than computing the photon
c        trajectories on-the-fly (which would be CPU costly). 
c      - Before calling do_trans() be sure that in source.F the functions 
c        calc_src_pnt() and calc_alp_pnt() give the correct source j resp. 
c        opacity alpha that should be included. It is in these routines 
c        that the user can add source terms, and for example add the 
c        scattering emissivity. The do_trans() routine here doesn't bother
c        at all about the source terms: it simply calls the calc_src_pnt()
c        and calc_alp_pnt() routines, and leaves it to these routine to 
c        decide which emission/absorption processes are included. The
c        do_trans() routine uses these values and performs the transfer.
c        This way the do_trans() routine is independent on the physics of
c        the medium. This medium may be anything from molecular clouds to
c        X-ray emitting hot plasmas. For more information about the source
c        terms: see the header of these two functions calc_src_pnt() and 
c        calc_alp_pnt().
c      - Next, call do_trans(inu), where the index inu refers to the 
c        frequency found in the array nu_ar(*). The source terms and 
c        absorptivity in the transfer will be for this frequency.
c     
c     RESULT:
c     The routine produces the array intensity(iphi,imu,itheta,iradius) 
c     which is the radiation field for all angles at all gridpoints.
c     
c     --------------------------------------------------------------
      subroutine do_trans_char(inu)
      implicit none
      integer inu
c
#     ifndef NO_INTENSITY_STORAGE
#     ifndef ONLY_RAY_TRACING
c
      doubleprecision find_int_t
c
#     ifdef LC_INTEGRATION
      doubleprecision integrate_LC
#     else
      doubleprecision integrate_ESC
#     endif
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
c#include "common_units.h"
#include "common_boundary.h"
#include "common_shortchar.h"
c
      doubleprecision int,dr,dt,dmu,dphi,rcrp2,dummy
      doubleprecision mui,mus,mured,mu,phi,r0,rscale
      integer ir,it,imu,iphi,is,iphimax,iphimin,imumax,nel,icr,idrr
      integer imus,iphis,ithetamax,iphioff,idrr2,imu0
      integer iphihalf,iangset,i,idrrdw
c
c     Now check some short characteristics stuff
c
#     ifndef NO_SHORT_CHARS
c
c     Clear ESC histogram (debugging purposes only)
c
      do i=1,FRSIZE_CHAR 
          dbsc_hist_nel(i) = 0
      enddo
c
c     Check whether the Extended Short Characteristics have in fact
c     been rendered before this subroutine was called
c
      if(isc_ready.ne.1) then
          write(*,*) 'HOUSTON: While trying to start transfer:' 
          write(*,*) '         The Extended Short Characteristics '
          write(*,*) '         have not been redered yet. This is '
          write(*,*) '         probably the result of the fact that '
          write(*,*) '         the Extended Short Characteristics '
          write(*,*) '         were longer than allowed by the '
          write(*,*) '         pre-compiled array boundaries '
          write(*,*) '         (Check error message above!) '
          write(*,*) '         Recompile with larger FRSIZE_CHAR!'
          stop 13
      endif
c
c     First a hack... For the moment the idrr must be fixed.
c
#     ifndef SHORTCHAR_NEWSTYLE
c
c     Old style (pre 09-00)
c
      idrr    = 1
      idrr2   = 2
      idrrdw  = 1
      iangset = 1
c
#     else
c
c     New style (post 09-00)
c
      idrr    = 1
      idrr2   = 0
      idrrdw  = 1
      iangset = 1
c
#     endif
c
#     else /* ifndef NO_SHORT_CHARS */
c
c     If long characteristics, then set the following:
c
      idrr    = 1
      idrr2   = 0
      idrrdw  = 1
      iangset = 1
c
#     endif /* ifndef NO_SHORT_CHARS */
c
c     Make the index arrays for intens(). This index is used for
c     the mirror symmetries. 
c
      call make_index(iangset)
c
c     Get some important integer values...
c
      ithetamax = irsi_frsizey
      imumax    = nrmu(iangset)
      iphimin   = 1
      iphimax   = nrphi(iangset)
      iphihalf  = nrphi(iangset) / 2
      iphioff   = 0
c
c     The number of phi angles must be divisible by 4, for the
c     way the transfer is done here
c
#     ifndef RADGRID_ONEDIM
      if(4*(iphimax/4).ne.iphimax) then
          write(*,*) 'ERROR: Sorry, the number of phi ',
     %               'angles must be divisible by 4'
          stop 13
      endif
#     endif
c
c     It is important to know whether we have a mu-gripoint
c     at the value mu=0. Find that out.
c
      if(anggrid_mu_zero.eq.1) then
          imu0 = 0
      else
          imu0 = 1
      endif
c
c     Now the mirror things...
c
#ifdef MIRROR_THETA
      ithetamax = (irsi_frsizey+1)/2
#endif
#ifdef MIRROR_PHI
      iphimin = nrphi(iangset) / 4 + 1
      iphimax = 3 * nrphi(iangset) / 4 
      iphioff = - nrphi(iangset) / 4
#endif
#ifdef RADGRID_ONEDIM
      iphimin = 1
      iphimax = 1
      iphioff = 0
      ithetamax = 1
#endif
c
c     First we fill (if necessary) the inwards pointing intensity
c     at the R=R_max. This could be useful for microwave background
c     or in a 1-D plane--parallel calculation where the photosphere
c     is on the R=R_0 side, instead of R=R_max.
c
      if(iradbnd_out_itype.eq.1) then
c
c         Outer equilibrium boundary condition
c
          call outer_equilibrium(inu)
c
      elseif(iradbnd_out_itype.eq.2) then
c
c         Microwave background radiation outer BC
c
          call outer_microwave_backgr(inu)
c
      elseif(iradbnd_out_itype.eq.3) then
c
c         Interstellar radiation field outer BC
c
          call outer_interstellar_radfield(inu)
c
      elseif(iradbnd_out_itype.eq.0) then
c
c         Vacuum outer BC
c
          call outer_vacuum(inu)
c
      else
          write(*,*) 'Dont know this type of outer bc'
          stop 13
      endif
c
c     First sweep from R=infty to R=Rin. This is the first of the two 
c     main loops. We start at R=infty, and sweep the grid along 
c     all mu's phi's and Theta's. We follow the short characteristics
c     from each radial shell R_i up to the radial shell R_{i+1},
c     find the values of the intensity there, integrate the source
c     term along the short characteristics and thus update the 
c     intensity at the shell R_i everywhere and at every angle.
c     We repeat this until we reach the innermost radius.
c
      do ir=irsi_frsizex-1,1,-1
c
c         Find the R, and the scale factor by which the dimensionless
c         sc_ds should be multiplied in order to find the real ds 
c         along the elements of the extended short characteristics.
c     
c            Spherical coordinates:   Use R itself as a scaling factor
c            Cartesian coordinates:   Use R_{i+1}-R_i as a scaling factor
c
#ifdef COORD_SPHERICAL
          r0     = rsi_x_c(ir,1)
          rscale = r0
#else
          r0     = rsi_x_c(ir,1)
          rscale = rsi_x_c(ir+1,1) - rsi_x_c(ir,1)
          stop 842
#endif
c
c         In the new style of SC we may have varying dr/r, so we
c         retrieve the index of the SC-set appropriate for this
c         radius from an array
c
#         ifndef NO_SHORT_CHARS
#         ifdef SHORTCHAR_NEWSTYLE
          idrr   = ir_iscset_p(ir)
          if(ir.gt.1) then
              idrrdw = ir_iscset_n(ir)
          else
              idrrdw = idrr
          endif
#         endif
#         endif
c
c         Check if the iangset of the idrr is properly set
c
#         ifndef NO_SHORT_CHARS
          if(isc_iangset(idrr).eq.0) then
              write(*,*) 'Major internal error. Warn author.'
              stop 111
          endif
#         endif
c
c         Since we integrate towards smaller R, we need not worry
c         about mu-binning. Therefore put mured=1.d0
c
          mured  = 1.d0
c
c         Do a loop along theta from North Pole to either
c         the South Pole or the Equator (dependent on mirror 
c         the existence of symmetry).
c
c         NOTE: For ESC this order is not important, but for
c               the MESC algorithm it is essential!!
c
          do it=1,ithetamax
c
c             Then a loop over mu
c
              do imu=-imumax,-imu0
c
c                 Then a loop over phi, only for those angles that
c                 point towards positive Theta
c
                  do iphi=iphihalf+1,iphimax
c
c                     Perform the integration along the short/long char
c
#                     ifndef LC_INTEGRATION
                      intens(iphi+iphioff,imu,it,ir) = 
     %                      integrate_ESC(ir,it,imu,iphi,imumax,idrr,
     %                              idrrdw,inu,rscale,mured) 
#                     else
                      intens(iphi+iphioff,imu,it,ir) = 
     %                      integrate_LC(ir,it,imu,iphi,imumax,idrr,
     %                              inu,mured) 
#                     endif
c
c                     Done... Continue loop
c
                  enddo
c
c                 Now address the problem of the 3-rd and 4-th order
c                 interpolation in phi, where one of the points may
c                 accidently lie outside of the safe domain for the
c                 MESC algorithm. We solve this by applying a filthy
c                 trick: copy nearest neighbors. This is not a neat
c                 solution, but not wrong either, since the problem
c                 arises solely due to the fact that interpolation
c                 is a non--local effect in the phi-space, possibly
c                 involving points not causally connected to the
c                 phi at hand.
c             
#ifdef ESC_MINIMAL_EXTENSION
                  intens(iphihalf+iphioff,imu,it,ir) =
     %               intens(iphihalf+1+iphioff,imu,it,ir)
                  intens(iphihalf-1+iphioff,imu,it,ir) =
     %               intens(iphihalf+1+iphioff,imu,it,ir)
#ifndef MIRROR_PHI
                  intens(1,imu,it,ir) =
     %               intens(iphimax,imu,it,ir)
                  intens(2,imu,it,ir) =
     %               intens(iphimax,imu,it,ir)
#endif
#endif
c
              enddo
          enddo
c
c         Now do a loop along the theta from the Equator or
c         the South pole to the North Pole
c
          do it=ithetamax,1,-1
c
c             Then a loop over mu
c
              do imu=-imumax,-imu0
c
c                 Then a loop over phi, only for those angles that
c                 point towards negative Theta
c
                  do iphi=iphimin,iphihalf
c
c                     Perform the integration along the short/long char
c
#                     ifndef LC_INTEGRATION
                      intens(iphi+iphioff,imu,it,ir) = 
     %                      integrate_ESC(ir,it,imu,iphi,imumax,idrr,
     %                              idrrdw,inu,rscale,mured) 
#                     else
                      intens(iphi+iphioff,imu,it,ir) = 
     %                      integrate_LC(ir,it,imu,iphi,imumax,idrr,
     %                              inu,mured)
#                     endif
c
c                     Done... Continue loop
c
                  enddo
              enddo
c
c             If we have a disk at the equator, then we must 
c             do some special things overhere. The disk can
c             be activated by putting ieqd_active to 1. Be sure
c             to have set the hard-switch EQUATOR_DISK. If
c             this switch is not defined, the disk module is
c             not even compiled. This is for safety reasons:
c             If at some time mysterious bugs appear, parts of
c             the code can be completely switched off.
c
#             ifdef EQUATOR_DISK
              if((ieqd_active.eq.1).and.(it.eq.ithetamax)) then
#                 ifdef RADGRID_ONEDIM
                  write(*,*) 'Sorry: The equatorial disk option ',
     %            'in the transfer.F module'
                  write(*,*) 'is only available for 2-D'
                  stop 13
#                 endif
#                 ifdef ESC_COMPLETE
                  write(*,*) 'Sorry: Equatorial disk is not yet ',
     %                   'compatible with the ESC integration'
                  write(*,*) 'Instead use the MESC integration by ',
     %                   'setting the ESC_MINIMAL_EXTENSION '
                  write(*,*) 'hard-switch in the Makefile.'
                  stop 13
#                 endif
#                 ifndef MIRROR_THETA
                  write(*,*) 'Sorry: The equatorial disk option ',
     %            'in the transfer.F module'
                  write(*,*) 'is currently only available ',
     %            'when the MIRROR_THETA is active' 
                  stop 13
#                 endif
                  call equator_disk(ir,inu)
              endif
#             else
              if(ieqd_active.ne.0) then
                  write(*,*) 'Equatorial disk model can only be used'
                  write(*,*) 'when RADICAL is compiled with ',
     %                        'EQUATOR_DISK set'
                  write(*,*) 'in the configure.h file.'
                  stop 13
              endif
#             endif
          enddo
      enddo
c
c     Now add the possibility to emit radiation from the inner
c     boundary. This is handy when a star is present in the center.
c     If the inner boundary type is 2, then place the central source
c     not on the inner boundary of the radiative transfer domain, but
c     on a separate stellar surface (the intens array with ir=0).
c
      if(iradbnd_in_itype.eq.0) then
          write(*,*) 'STOPPED: Inner BC type 0 is temporarily disabled'
          write(*,*) '         Use type 2 instead, and choose R* tiny.'
          stop 13
      elseif(iradbnd_in_itype.eq.1) then
          call central_source(1,inu)
      elseif(iradbnd_in_itype.eq.2) then
          call central_source(0,inu)
      elseif(iradbnd_in_itype.eq.3) then
          call central_mirror(inu)
      else   !if(.not.(iradbnd_in_itype.eq.0)) then
          write(*,*) 'Dont know this type of inner bc'
          stop 13
      endif
c
c     If the inner boundary is transparent, and envelops a vacuum with 
c     a central star in it (iradbnd_in_itype.eq.2) then we must transfer 
c     the radiation from the inner boundary at R_in to itself, through the
c     vacuum. Some of that radiation can be intercepted by the central
c     star, so that must be eliminated, and replaced by the emission of
c     the star. The rest of the radiation simply flows freely through this
c     inner vaccum hole and reappears into the grid with opposite mu.
c
c     The reduction factor of the stellar light is given by the ratio
c     of the solid angle of the star divided by the solid angle of the
c     smallest mu-angle bin. It is given by the formula:
c
c              1 - mu_star   1 - sqrt( 1 - (R_s/R_1)^2 )
c      mured = ----------- = ---------------------------
c              1 - mu_imax           1 - mu_imax
c
c     where mu_imax is the largest mu in the mu-grid (more precise:
c     the interface of the mu=+1 bin), and mu_star is the mu which the
c     stellar surface has when projected on the sky of the observer
c     standing at R=R_1. 
c
      if((iradbnd_in_itype.eq.2).or.(iradbnd_in_itype.eq.0)) then
#ifdef COORD_SPHERICAL
          r0     = rsi_x_c(1,1)
          rscale = r0
          do it=1,ithetamax
              do imu=1,imumax
                  if(imu.eq.imumax) then
                      mui   = rmu_i(imu,iangset)   ! iangset is okay here
c
cccc BUG                     mured = min(1.d0,radbnd_rstar/
cccc BUG    %                          (rscale*sqrt(1.d0-mui*mui)))
c
cccc BUG2                    mured = min(1.d0,(1.d0-
cccc BUG2   %                   sqrt(1.d0-(radbnd_rstar/rscale)**2))
cccc BUG2   %                          /(1.d0-mui))
c
                      if(radbnd_rstar.gt.1.d-3*rscale) then
                          mured = min(1.d0,(1.d0-
     %                         sqrt(1.d0-(radbnd_rstar/rscale)**2))
     %                         /(1.d0-mui))
                      else
                          mured = min(1.d0,0.5d0*
     %                         ((radbnd_rstar/rscale)**2)/(1.d0-mui))
                      endif
                      if(iradbnd_in_itype.eq.0) then
                          mured = 0.d0
                      endif
c
c                     If the angular set resolves the star, then we must
c                     be sure that the angular spacing is fine enough!
c                     (added 05-01-01)
c
                      if(mured.eq.1.d0) then
                          dummy = (1.d0-rmu(imu-1,iangset)**2) /
     %                      (1.d0-rmu(imu,iangset)**2)
                          if(dummy.gt.1.21d0) then
                              write(*,*) 'Warning: star badly sampled',
     %                            ' (warning printed only once).'
                              transerr_starresolution = 1
                          endif
c##########################
                          if(dummy.lt.1.d0) stop 999
c##########################
                      endif
                  else
                      mured = 1.d0
                  endif
                  do iphi=iphimin,iphimax
c
#                     ifndef LC_INTEGRATION
c
c                     Find the start of the local characteristic.
c                     In the absence of a central object, all of the
c                     rays (except the mu=1 ray) will start at the
c                     same radial inner shell as it ends, namely the
c                     ir=1 shell. We just allow radiation to move
c                     freely (vacuum) along this characteristic.
c
#                     ifndef SHORCHAR_COMPACT_ARRAYS
                      nel   = isc_nel(iphi,imu,it,idrr2)
                      icr   = isc_icross(nel,iphi,imu,it,idrr2)
                      dr    = sc_dr(nel,iphi,imu,it,idrr2)
                      dt    = sc_dtheta(nel,iphi,imu,it,idrr2)
                      dmu   = sc_dmu(nel,iphi,imu,it,idrr2)
                      dphi  = sc_dphi(nel,iphi,imu,it,idrr2)
                      imus  = isc_imu(nel,iphi,imu,it,idrr2)
                      iphis = isc_iphi(nel,iphi,imu,it,idrr2)
#                     else
                      nel   = isc_nel(iphi,imu,it,idrr2)
                      call get_shortchar(nel,iphi,imu,it,idrr2,
     %                       icr,dr,dt,mu,phi,dmu,dphi,imus,iphis)
#                     endif
c
c                     The intensity at the start
c
                      int   = find_int_t(icr,1,it,dr,dt,
     %                                  imus,iphis,dmu,dphi)
c
#                     else /* ifndef LC_INTEGRATION */
c
c                     If LC integration, then call integrate_LC
c
                      int   = integrate_LC(1,it,imu,iphi,10000,idrr,
     %                                     inu,1.d0)
c
#                     endif /* ifndef LC_INTEGRATION */
c
c                     The radial ray (mu=1) is a special one, since
c                     it always hits the central star, even though
c                     most of the mu-bin it represents may pass the
c                     star if the star is tiny (and therefore unresolved
c                     in mu). So we must correct for that, if indeed
c                     the star is unresolved. If so, then part of 
c                     the radiation is allowed to move `through' the
c                     star (actually past) and part will be the 
c                     contribution of the star itself.
c
c                     WARNING: I now see that the (1.d0-mured) * 
c                              intens(1,-imu,it,1) is wrong in case
c                              of no theta-mirror-symmetry. Should
c                              (I think) be irsi_frsizey-it+1 index.
c
#                     ifndef RADGRID_ONEDIM 
#                     ifndef MIRROR_THETA
c                      write(*,*) 'WARNING: possible error in do_trans()'
                      warn_bug_centrbeam = 1
#                     endif
#                     endif
c
                      if(imu.eq.imumax) then
                          int = mured * int + (1.d0-mured) * 
     %                          intens(1,-imu,it,1) 
                      endif
c
c                     So now store it
c
                      intens(iphi+iphioff,imu,it,1) = int 
                  enddo
              enddo
          enddo
#else
#ifdef COORD_CARTESIAN
c
c         In case of Cartesian coordinates
c         Remark: in principle this type of inner bc should not
c                 have been selected for 1-D, but let's not be
c                 too picky. We'll correct for this minor mistake here:
c
          do it=1,ithetamax
              do imu=1,imumax
                  do iphi=iphimin,iphimax
                      intens(iphi+iphioff,imu,it,1) = 
     %                    intens(iphi+iphioff,imu,it,0)
                  enddo
              enddo
          enddo
#else
          write(*,*) 'inner bc not ready for cylindrical coordinates'
          stop 13
#endif
#endif
      endif
c
c     Now back to R=infty. This is the second of the two 
c     main loops. We start at R=R_in, and sweep the grid along 
c     all mu's phi's and Theta's. We follow the short characteristics
c     from each radial shell R_i down to the radial shell R_{i-1},
c     find the values of the intensity there, integrate the source
c     term along the short characteristics and thus update the 
c     intensity at the shell R_i everywhere and at every angle.
c     We repeat this until we are back at R=infinity, where we started
c     at the beginning of this subroutine.
c
      do ir=2,irsi_frsizex
c
c         Find the R, and the scale factor by which the dimensionless
c         sc_ds should be multiplied in order to find the real ds 
c         along the elements of the extended short characteristics. 
c     
c            Spherical coordinates:   Use R itself as a scaling factor
c            Cartesian coordinates:   Use R_i-R_{i-1} as a scaling factor
c
c         Also compute the mu-binning factor, used in spherical 
c         coordinates.
c     
#ifdef COORD_SPHERICAL
          r0     = rsi_x_c(ir,1)
          rscale = r0
          rcrp2  = rsi_x_c(ir-1,1) / rsi_x_c(ir,1)
          rcrp2  = rcrp2 * rcrp2
#else
          r0     = rsi_x_c(ir,1)
          rscale = rsi_x_c(ir,1) - rsi_x_c(ir-1,1)
          rcrp2  = 1.d0
#endif
c
c         In the new style of SC we may have varying dr/r, so we
c         retrieve the index of the SC-set appropriate for this
c         radius from an array
c
#         ifndef NO_SHORT_CHARS
#         ifdef SHORTCHAR_NEWSTYLE
          idrr   = ir_iscset_n(ir)
          if(ir.lt.irsi_frsizex) then
              idrrdw = ir_iscset_p(ir)
          else
              idrrdw = idrr
          endif
#         endif
#         endif
c
c         Check if the iangset of the idrr is properly set
c
#         ifndef NO_SHORT_CHARS
          if(isc_iangset(idrr).eq.0) then
              write(*,*) 'Major internal error. Warn author.'
              stop 111
          endif
#         endif
c
c         Do a loop along theta from North Pole to either
c         the South Pole or the Equator (dependent on mirror 
c         the existence of symmetry)
c
          do it=1,ithetamax
c
c             Then a loop over mu
c
              do imu=1,imumax
c
c                 If imu=imumax, then there will be cumulative binning.
c                 One could explain cumulative binning by comparing
c                 it with the human eye: details that are too far away
c                 will average out over one eye-pixel. So these are
c                 binned cumulatively onto one pixel. The one pixel in
c                 our case here is the pixel with mu=+1. So I wish to 
c                 know the averaging factor between the mu=+1 bin and 
c                 previous mu bin.
c                 NOTE: I assume that a maximum of one previous bin can 
c                 be absorbed into the cumulative bin.
c
                  if(imu.eq.imumax) then
                      mured = rcrp2
                  else
                      mured = 1.d0
                  endif
c
c                 Then a loop over phi, only for those angles that
c                 point towards positive Theta
c
                  do iphi=iphihalf+1,iphimax
c
c                     Perform transfer along ESC:
c
#                     ifndef LC_INTEGRATION
                      intens(iphi+iphioff,imu,it,ir) = 
     %                      integrate_ESC(ir,it,imu,iphi,imumax,idrr,
     %                            idrrdw,inu,rscale,mured) 
#                     else
                      intens(iphi+iphioff,imu,it,ir) = 
     %                      integrate_LC(ir,it,imu,iphi,imumax,idrr,
     %                            inu,mured)
#                     endif
c
c                     Done... Continue with loop
c
                  enddo
c
c                 Now address the problem of the 3-rd and 4-th order
c                 interpolation in phi. See above (where exactly the
c                 same has been done).
c             
#ifdef ESC_MINIMAL_EXTENSION
                  intens(iphihalf+iphioff,imu,it,ir) =
     %               intens(iphihalf+1+iphioff,imu,it,ir)
                  intens(iphihalf-1+iphioff,imu,it,ir) =
     %               intens(iphihalf+1+iphioff,imu,it,ir)
#ifndef MIRROR_PHI
                  intens(1,imu,it,ir) =
     %               intens(iphimax,imu,it,ir)
                  intens(2,imu,it,ir) =
     %               intens(iphimax,imu,it,ir)
#endif
#endif
              enddo
          enddo
c
c         Now do a loop along the theta from the Equator or
c         the South pole to the North Pole
c
          do it=ithetamax,1,-1
c
c             Then a loop over mu
c
              do imu=1,imumax
c
c                 Again, the mu-binning stuff:
c
                  if(imu.eq.imumax) then
                      mured = rcrp2
                  else
                      mured = 1.d0
                  endif
c
c                 Then a loop over phi, only for those angles that
c                 point towards positive Theta
c
                  do iphi=iphimin,iphihalf
c
c                     Perform transfer along ESC:
c
#                     ifndef LC_INTEGRATION
                      intens(iphi+iphioff,imu,it,ir) = 
     %                      integrate_ESC(ir,it,imu,iphi,imumax,idrr,
     %                            idrrdw,inu,rscale,mured)
#                     else
                      intens(iphi+iphioff,imu,it,ir) = 
     %                      integrate_LC(ir,it,imu,iphi,imumax,idrr,
     %                            inu,mured)
#                     endif
c
c                     Done... Continue with loop
c
                  enddo
              enddo
c
c             If we have a disk at the equator, then we must 
c             do some special things overhere. The disk can
c             be activated by putting ieqd_active to 1. Be sure
c             to have set the hard-switch EQUATOR_DISK. If
c             this switch is not defined, the disk module is
c             not even compiled. This is for safety reasons:
c             If at some time mysterious bugs appear, parts of
c             the code can be completely switched off.
c
#             ifdef EQUATOR_DISK
              if((ieqd_active.eq.1).and.(it.eq.ithetamax)) then
#                 ifndef MIRROR_THETA
                  write(*,*) 'Sorry: The equatorial disk option ',
     %            'in the transfer.F module'
                  write(*,*) 'is currently only available ',
     %            'when the MIRROR_THETA is active' 
                  stop 13
#                 endif
#                 ifdef RADGRID_ONEDIM
                  write(*,*) 'Sorry: The equatorial disk option ',
     %            'in the transfer.F module'
                  write(*,*) 'is only available for 2-D'
                  stop 13
#                 endif
                  call equator_disk(ir,inu)
              endif
#             endif
          enddo
c
c         Now the transfer for this ir is done.
c     
c         If mu=+1 and if not forbidden, then average out all
c         phi bins. 
c         
c           WARNING: Assume the phi bins to be equally spaced.
c           WARNING: The quadrature is LINEAR here 
c
#ifndef NO_MUONE_AVER
          do it=1,ithetamax
              imu = imumax
              int = 0.d0
              do iphi=iphimin,iphimax
                  int = int + intens(iphi+iphioff,imu,it,ir)
              enddo
              int = int / ( iphimax - iphimin + 1 )
              do iphi=iphimin,iphimax
                  intens(iphi+iphioff,imu,it,ir) = int 
              enddo                      
          enddo
#endif
      enddo
c
#     else /* ifndef ONLY_RAY_TRACING */ 
      write(*,*) 'ERROR: do_trans_char() has been called'
      write(*,*) '       while the ONLY_RAY_TRACING is active'
      stop 13
#     endif /* ifndef ONLY_RAY_TRACING */ 
c
#     else /* ifndef NO_INTENSITY_STORAGE */
      write(*,*) 'ERROR: do_trans_char() has been called.'
      write(*,*) '       This function needs the intensity array.'
      write(*,*) '       But this array is not present, because'
      write(*,*) '       NO_INTENSITY_STORAGE is defined.'
      write(*,*) '       Tip: undefine NO_INTENSITY_STORAGE in'
      write(*,*) '       configure.h'
      stop 13
#     endif /* ifndef NO_INTENSITY_STORAGE */
c
      end



c     --------------------------------------------------------------
c             DO TRANSFER ALONG ONE SHORT CHARACTERISTICS
c
c     Integrate along one single Extended Short Characteristic.
c     This routine is called by do_trans(), for each (R,Theta,mu,phi)
c     point. It follows the (extendedn) short characteristics from
c     upwind to downwind and performs the transfer. It autmatically
c     finds the I_0 at the start of the ESC, and then performs the
c     integral until the point P is reached.
c     --------------------------------------------------------------
#     ifndef NO_SHORT_CHARS
      function integrate_ESC(ir,it,imu,iphi,imumax,idrr,idrrdw,
     %                       inu,rscale,mured)
      implicit none
c
      doubleprecision integrate_ESC,rscale,mured
      integer ir,it,imu,iphi,idrr,imumax,inu,idrrdw
c
#ifdef SRCQDR_PNT_2
      doubleprecision qdr_src_2
#else
      doubleprecision qdr_src_3
#endif
      doubleprecision find_int_t,find_src_t,find_alp_t
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
c#include "common_units.h"
#include "common_boundary.h"
#include "common_shortchar.h"
c
      doubleprecision int,src,src1,src2,src3,alp1,alp2,alp3
      doubleprecision dr,dt,ds,ds1,ds2
      doubleprecision mu,phi,dmu,dphi
      integer is,nel,icr
      integer imus,iphis
c
c     Now let's find the start of the local characteristic
c
#     ifndef SHORCHAR_COMPACT_ARRAYS
      nel   = isc_nel(iphi,imu,it,idrr)
      icr   = isc_icross(nel,iphi,imu,it,idrr)
      dr    = sc_dr(nel,iphi,imu,it,idrr)
      dt    = sc_dtheta(nel,iphi,imu,it,idrr)
      mu    = sc_mu(nel,iphi,imu,it,idrr)
      phi   = sc_phi(nel,iphi,imu,it,idrr)
      dmu   = sc_dmu(nel,iphi,imu,it,idrr)
      dphi  = sc_dphi(nel,iphi,imu,it,idrr)
      imus  = isc_imu(nel,iphi,imu,it,idrr)
      iphis = isc_iphi(nel,iphi,imu,it,idrr)
#     else
      nel   = isc_nel(iphi,imu,it,idrr)
      call get_shortchar(nel,iphi,imu,it,idrr,
     %                   icr,dr,dt,mu,phi,dmu,dphi,imus,iphis)
#     endif
c
c     For debugging purposes
c
#ifndef NO_ESC_HISTOGRAM
      if(nel.gt.FRSIZE_CHAR) then
          write(*,*) 'Huh? transfer.F/integrate_ESC(): nel=',nel
          stop 778
      endif
      dbsc_hist_nel(nel) = dbsc_hist_nel(nel)+1
#endif
c
c     Next find the intensity between gridpoints,
c     by performing an interpolation procedure
c
      int = find_int_t(icr,ir,it,dr,dt,
     %    imus,iphis,dmu,dphi)
c     
c     If mu=+1, then also include the previous mu
c     bin, because the mu=+1 bin is a cumulative bin
c     (only for spherical coordinates)
c     
#ifdef COORD_SPHERICAL
      if(imu.eq.imumax) then
          int = mured * int + (1.d0-mured) * 
     %        find_int_t(icr,ir,it,dr,dt,
     %        imumax-1,iphi,0.d0,0.d0)
          if(icr.ne.1) then
              warn_mubin_theta = 1
          endif
      endif
#endif
c
c     Now we have the choice between 2-point integration
c     and 3-point integration (lin vs quadr). These are
c     switched by a define.
c
#ifndef NO_SOURCE
#ifdef SRCQDR_PNT_2
c
c     This is the 2-point source integrator
c
ccccccccc      mu   = sc_mu(nel,iphi,imu,it,idrr)
ccccccccc      phi  = sc_phi(nel,iphi,imu,it,idrr)
c
c     Now find the total source and the total opacity 
c     at this start of the characteristic
c     
      src1 = find_src_t(inu,icr,ir,it,dr,dt,mu,phi)
      alp1 = find_alp_t(inu,icr,ir,it,dr,dt,mu,phi)
c     
c     Loop in downwind direction towards point 0
c     
      do is=nel-1,0,-1
c     
c         Find the information of this point on the 
c         local characteristic
c     
          if(is.eq.0) then
              icr  = 1
              dr   = 0.d0
              dt   = 0.d0
              mu   = rmu(imu,isc_iangset(idrr))
              phi  = rphi(iphi,isc_iangset(idrr))
          else
#             ifndef SHORCHAR_COMPACT_ARRAYS
              icr  = isc_icross(is,iphi,imu,it,idrr)
              dr   = sc_dr(is,iphi,imu,it,idrr)
              dt   = sc_dtheta(is,iphi,imu,it,idrr)
              mu   = sc_mu(is,iphi,imu,it,idrr)
              phi  = sc_phi(is,iphi,imu,it,idrr)
#             else
              call get_shortchar_concise(is,iphi,imu,it,idrr,
     %                   icr,dr,dt,mu,phi)
#             endif
          endif
c     
c         Find the ds
c     
#         ifndef SHORCHAR_COMPACT_ARRAYS
          ds = rscale * sc_ds(is+1,iphi,imu,it,idrr)
#         else
          call get_shortchar_ds(is+1,iphi,imu,it,idrr,ds)
          ds = rscale * ds
#         endif
c     
c         Find the new src and alpha
c     
          src2 = find_src_t(inu,icr,ir,it,dr,dt,mu,phi)
          alp2 = find_alp_t(inu,icr,ir,it,dr,dt,mu,phi)
c     
c         Now we can actually perform the formal
c         integral over this small piece of local
c         characteristic
c         
          int = qdr_src_2(int,src1,alp1,src2,alp2,ds)
c     
c         New is old, for the next step
c     
          src1 = src2
          alp1 = alp2
c     
      enddo
c     
#else /* ifdef SRCQDR_PNT_2 */
c     
c     This is the 3-point source integrator
c     
ccccccc      mu   = sc_mu(nel,iphi,imu,it,idrr)
ccccccc      phi  = sc_phi(nel,iphi,imu,it,idrr)
c
c     Now find the total source and the total opacity 
c     at this start of the characteristic...
c     
      src1 = find_src_t(inu,icr,ir,it,dr,dt,mu,phi)
      alp1 = find_alp_t(inu,icr,ir,it,dr,dt,mu,phi)
c     
c     ... and the next point (this may be point 0 
c     already or another point)
c     
      is=nel-1
      if(is.eq.0) then
          icr  = 1
          dr   = 0.d0
          dt   = 0.d0
          mu   = rmu(imu,isc_iangset(idrr))
          phi  = rphi(iphi,isc_iangset(idrr))
      else
#         ifndef SHORCHAR_COMPACT_ARRAYS
          icr  = isc_icross(is,iphi,imu,it,idrr)
          dr   = sc_dr(is,iphi,imu,it,idrr)
          dt   = sc_dtheta(is,iphi,imu,it,idrr)
          mu   = sc_mu(is,iphi,imu,it,idrr)
          phi  = sc_phi(is,iphi,imu,it,idrr)
#         else
          call get_shortchar_concise(is,iphi,imu,it,idrr,
     %         icr,dr,dt,mu,phi)
#         endif
      endif
c     
c     Find the ds
c     
#     ifndef SHORCHAR_COMPACT_ARRAYS
      ds1 = rscale * sc_ds(is+1,iphi,imu,it,idrr)
#     else
      call get_shortchar_ds(is+1,iphi,imu,it,idrr,ds1)
      ds1 = rscale * ds1
#     endif
c     
c     Find the new src and alpha
c     
      src2 = find_src_t(inu,icr,ir,it,dr,dt,mu,phi)
      alp2 = find_alp_t(inu,icr,ir,it,dr,dt,mu,phi)
c     
c     Loop in downwind direction towards point 0
c     
      do is=nel-1,0,-1
c     
c         Find the information for the next point on the 
c         local characteristic (2nd order integration!)
c     
          if(is-1.eq.0) then
              icr  = 1
              dr   = 0.d0
              dt   = 0.d0
              mu   = rmu(imu,isc_iangset(idrr))
              phi  = rphi(iphi,isc_iangset(idrr))
          elseif(is-1.eq.-1) then
ccc BUGFIX: idrr --> idrrdw (21-06-01)
              icr  = isc_icross(0,iphi,imu,it,idrrdw)
              dr   = sc_dr(0,iphi,imu,it,idrrdw)
              dt   = sc_dtheta(0,iphi,imu,it,idrrdw)
              mu   = sc_mu(0,iphi,imu,it,idrrdw)
              phi  = sc_phi(0,iphi,imu,it,idrrdw)
          else
#             ifndef SHORCHAR_COMPACT_ARRAYS
              icr  = isc_icross(is-1,iphi,imu,it,idrr)
              dr   = sc_dr(is-1,iphi,imu,it,idrr)
              dt   = sc_dtheta(is-1,iphi,imu,it,idrr)
              mu   = sc_mu(is-1,iphi,imu,it,idrr)
              phi  = sc_phi(is-1,iphi,imu,it,idrr)
#             else
              call get_shortchar_concise(is-1,iphi,imu,it,idrr,
     %                   icr,dr,dt,mu,phi)
#             endif
          endif
c
c         Get the new ds2
c
          if(is.gt.0) then
#             ifndef SHORCHAR_COMPACT_ARRAYS
              ds2 = rscale * sc_ds(is,iphi,imu,it,idrr)
#             else
              call get_shortchar_ds(is,iphi,imu,it,idrr,ds2)
              ds2 = rscale * ds2
#             endif
          else
ccc BUGFIX: idrr --> idrrdw (21-06-01)
#             ifndef SHORCHAR_COMPACT_ARRAYS
              ds2 = rscale * sc_ds(0,iphi,imu,it,idrrdw)
#             else
              call get_shortchar_ds(0,iphi,imu,it,idrrdw,ds2)
              ds2 = rscale * ds2
#             endif
          endif
c     
c         Find the src and alpha
c     
ccc BUGFIX: For ir=1 mu<0 and ir=frsizex mu>0 no D point (21-06-01)
c          
          if(((ir.eq.1).and.(imu.lt.0)).or.
     %            ((ir.eq.irsi_frsizex).and.(imu.gt.0))) then
c
c             No downwind point present
c
              ds2  = ds1*1d-4
              src3 = src2 + (src2-src1)*1.d-4
              alp3 = alp2
          else
c
c             Yes, there is a D point
c
              src3 = find_src_t(inu,icr,ir,it,dr,dt,mu,phi)
              alp3 = find_alp_t(inu,icr,ir,it,dr,dt,mu,phi)
          endif
c     
c         Now we can actually perform the formal
c         integral over this small piece of local
c         characteristic
c     
          int = qdr_src_3(int,src1,alp1,
     %        src2,alp2,src3,alp3,ds1,ds2)
c     
c         New is old, for the next step
c     
ccc BUGFIX: I forgot to copy ds2--->ds1 (21-06-01)
          ds1  = ds2
          src1 = src2
          alp1 = alp2
          src2 = src3
          alp2 = alp3
c     
      enddo
c     
#endif /* ifdef SRCQDR_PNT_2 */
#endif /* ifdef NO_SOURCE */      
c
c     Return the answer
c     
      integrate_ESC = int
      return
c
      end
#     endif /* ifndef NO_SHORT_CHARS */ 




#ifdef LC_INTEGRATION
c     --------------------------------------------------------------
c                    LONG CHARACTERISTIC INTEGRATION
c
c     The ESC algorithm has the intrinsic weakness of being 
c     approximative, even in vacuum. The method of Long 
c     Characteristics is much slower, but also more reliable.
c     This method is implemented here in a rather blunt way, so
c     if this method is going to be used for production purposes
c     then one might consider rearranging things so that the
c     trajectories are computed once per gridpoint, and used
c     by all frequencies. But for the moment the trajectories
c     are calculated for each gridpoint separately. 
c
c     This is the integration routine for long characteristics. 
c     It has the same entries as the integrate_ESC() routine, which
c     is the default for RADICAL. But instead of using the ESC,
c     this routine uses the trajectory routines of the telescope.F
c     module to integrate the transfer equation. 
c
c     NOTE: This routine assumes that R* is small enough
c     --------------------------------------------------------------
      function integrate_LC(ir,it,imu,iphi,imumax,idrr,inu,mured)
      implicit none
      integer ir,it,imu,iphi,idrr,inu,imumax
      doubleprecision integrate_LC,mured
c
#include "common_grid.h"
#include "common_boundary.h"
#include "common_telescope.h"
c
      doubleprecision epsrad,epss
      parameter(epsrad=1.d2*TELESC_EPS)
      parameter(epss=5.d1*TELESC_EPS)
c
      doubleprecision rorig,torig,r,theta,phi,mu,b,costinf,sintinf2
      doubleprecision z00,s00,z0,s0,x0,dum,charint,thetainf
      doubleprecision src1,src2,alp1,alp2,ds,intc,inte
      doubleprecision find_int_t,find_src_t,find_alp_t,qdr_src_2
      integer tram,i
c
#ifndef COORD_SPHERICAL
      write(*,*) 'Long Characteristics only for spherical coords'
      stop 13
#endif
c
c     Check if we are at mu=1 or not, since the mu=1 case is 
c     treated in a special way (mu-binning).
c
      if((imu.lt.imumax).or.(ir.le.1)) then
c
c         Get the coordinates
c
          rorig = rsi_x_c(ir,1)
          torig = rsi_x_c(it,2)
          mu    = rmu(imu,idrr)
          phi   = rphi(iphi,idrr)
c
c         In order not to wake-up possible bugs in the telescope
c         routine, I shift the r or the Theta slightly
c
          if(mu.eq.0.d0) then
              r     = rorig * (1.d0+epsrad)
              theta = torig
          else
              r     = rorig
              if(phi.gt.3.14159265359d0) then
                  theta = torig - epsrad
              else
                  theta = torig + epsrad
              endif
          endif
c
c         Find the Thetainf, b and z0 of the trajectory going through
c         this point at this angle, and the s=s0 where it reaches this
c         point.
c
          b        = r * sqrt(1.d0-mu*mu)
          costinf  = mu * cos(theta) + sqrt(1.d0-mu*mu) * 
     %                 sin(theta) * sin(phi) 
          z00      = r * ( (1.d0-mu*mu) * cos(theta) - 
     %                 mu * sqrt(1.d0-mu*mu) * sin(theta) * sin(phi) ) 
          s00      = r * mu
          thetainf = acos(costinf)
c
c         For safety, make s00 a tiny bit larger, so that we are sure that
c         the trajectory in fact ends at P, and not accidently at one 
c         node earlier.
c
          s00      = s00 + epss*abs(r)
c         
c         Stupidity check:
c         
          if((thetainf.lt.0.d0).or.(thetainf.gt.3.15)) then
              write(*,*) 'ERROR in integrate_LC(): Thetainf wrong...'
              stop 445
          endif
c         
c         For historic reasons the parameterization of the trajectories
c         for the telescope.F routine was done differently from the
c         ESCs. 
c         
c           - telescope.F:     z0,x0,theta0,s
c           - ESC/MESC:        z00,b,theta0,s00
c         
c             b^2 = x0^2 + z0^2*(sin(thetainf))^2
c             s00 = s + z0*cos(thetainf)
c             z00 = z0*(sin(thetainf))^2
c         
          sintinf2 = 1.d0 - costinf**2
          if(sintinf2.eq.0.d0) then
              write(*,*) 'ERROR in integrate_LC(): sin(Thetainf)=0'
              stop 13
          endif
          z0       = z00 / sintinf2
          s0       = s00 - z0*costinf
          dum      = b*b - z0*z0*sintinf2
          if(dum.lt.0.d0) then
              write(*,*) 'ERROR in integrate_LC(): x0 imaginary'
              stop 13
          endif
          x0       = sqrt(dum)
c         
c         Now make the long characteristic. The istar needs not be 1, 
c         since the charint will take care of the star itself, by simple
c         impact parameter arguments. Also, by specifying s0 as not being
c         1d-30 (as in the telescope), but s0, the long characteristic is 
c         stopped at the point P.
c         
          call make_trajectory_c(x0,z0,thetainf,s0,0)
c
c##########################
          if(mu.ne.0.d0) then
              if(tr_icross(tr_amount).ne.1) then
                  write(*,*) 'INTERNAL ERROR: long char ends wrongly'
                  stop 2961
              endif
          endif
          if(abs(tr_s(tr_amount)-s0)/rorig.gt.1d-5) then
              write(*,*) 'INTERNAL ERROR: long char send-s0...'
              stop 2962
          endif
c##########################
c         
c         Now perform the transfer and return
c         
          integrate_LC = charint(inu,1,0.d0)
c
      else
c         
c         If mu=+1, then also include the previous mu
c         bin, because the mu=+1 bin is a cumulative bin
c     
c         WARNING: If one day I will make the LC method more
c                  efficient by doing all angles and frequencies 
c                  at once for each gridpoint, then I must make
c                  a separate array for the intensity at 
c                  imu=imumax and imu=imumax-1, and do something
c                  to the code below.
c          
#         ifdef SRCQDR_PNT_3
          write(*,*) 'Second order integration not yet compatible'
          write(*,*) 'with long characteristics method.'
          stop 13
#         endif
c
c         We must still do one short char integration here, since
c         we must transport this radiation from the previous radius
c         to the current.
c
          src1 = find_src_t(inu,1,ir-1,it,0.d0,0.d0,1.0,0.0)
          alp1 = find_alp_t(inu,1,ir-1,it,0.d0,0.d0,1.0,0.0)
          src2 = find_src_t(inu,1,ir,it,0.d0,0.d0,1.0,0.0)
          alp2 = find_alp_t(inu,1,ir,it,0.d0,0.d0,1.0,0.0)
          ds   = rsi_x_c(ir,1)-rsi_x_c(ir-1,1)
          intc = find_int_t(1,ir-1,it,0.d0,0.d0,
     %                   imumax,iphi,0.d0,0.d0)
          inte = find_int_t(1,ir-1,it,0.d0,0.d0,
     %                   imumax-1,iphi,0.d0,0.d0)
          integrate_LC = mured * intc + (1.d0-mured) * inte
          integrate_LC = qdr_src_2(integrate_LC,
     %                    src1,alp1,src2,alp2,ds)
      endif
c
c     Done...
c      
      return
      end
#endif /* ifdef LC_INTEGRATION */



c     --------------------------------------------------------------
c      INTEGRATE THE FORMAL SOLUTION ALONG ONE ELEMENT OF SHORT CHAR
c
c     This function computes the formal integral over a segment of
c     a local characteristic. This will be done first order here.
c
c     This linear integration is not suitable for optically thick
c     systems. In that case the 3-point quadratic integrator,
c     qdr_src_3() should be used. 
c     --------------------------------------------------------------
      function qdr_src_2(int,js1,alp1,js2,alp2,ds)
      implicit none
      doubleprecision qdr_src_2
      doubleprecision int,js1,alp1,js2,alp2,ds
c
      doubleprecision e0,e1,a,b,dtau1,xp,src1,src2,theomax
c
#     ifndef SRCQDR_NO_ERRCHECK
#include "common_radfield.h"      
#     endif
c
cc#ifdef SRCQDR_IMPLICIT
cc      qdr_src_2 = ( int + 0.5d0 * ds * (js1+js2) ) /
cc     %            ( 1.d0 + 0.5d0 * ds * (alp1+alp2) )
cc      return
cc#endif
cc#ifdef SRCQDR_EPOW
      dtau1     = 0.5d0 * ( alp1 + alp2 ) * ds
      theomax   = 0.5d0 * ( js1 + js2 ) * ds
      if(dtau1.gt.1.d-6) then
          xp        = exp(-dtau1)
          e0        = 1.d0 - xp
          e1        = dtau1 - e0
          b         = e1 / dtau1
          a         = e0 - b
      else
          a         = 0.5d0 * dtau1
          b         = 0.5d0 * dtau1
          xp        = 1.d0 - dtau1 
      endif
      if(alp1.gt.0.d0) then
          src1 = js1 / alp1
      elseif(alp2.gt.0.d0) then
          src1 = js2 / alp2
      else
          src1 = 0.d0
      endif
      if(alp2.gt.0.d0) then
          src2 = js2 / alp2
      elseif(alp1.gt.0.d0) then
          src2 = js1 / alp1
      else
          src2 = 0.d0
      endif
      if(dtau1.gt.1e-9) then
         qdr_src_2 = a * src1 + b * src2
      else
         qdr_src_2 = 0.5d0 * ( js1 + js2 ) * ds
      endif
c
c     Verify if serious overshoot occur (say: factor 2 or more)
c
#     ifndef SRCQDR_NO_ERRCHECK
      if(qdr_src_2.gt.2.d0*theomax) then
          transerr_overshoot =  transerr_overshoot + 1
      endif
#     endif
c
c     NEW: Bracket the total integrated emission between the
c          theoretical maximum and 0.
c
#     ifndef SRCQDR_NO_BRACKET
      qdr_src_2 = min(qdr_src_2,theomax)
#     endif
c
      qdr_src_2 = int * xp + qdr_src_2
      return
c
cc#endif
c
      write(*,*) 'Ehh, something wrong in dqr_src_2'
      stop 13
      return
      end




c     --------------------------------------------------------------
c      INTEGRATE THE FORMAL SOLUTION ALONG ONE ELEMENT OF SHORT CHAR
c
c     This function computes the formal integral over a segment of
c     a local characteristic. This will be done second order here.
c     This is necessary for optically thick systems, but hardy so
c     for optically thin systems.
c
c     WARNING: Be careful for optical depths close to 1e-6, because
c              this regime might be uncertain.....
c
c     --------------------------------------------------------------
      function qdr_src_3(int,js1,alp1,js2,alp2,js3,alp3,ds1,ds2)
      implicit none
      doubleprecision qdr_src_3
      doubleprecision int,js1,alp1,js2,alp2,js3,alp3,ds1,ds2
c
      doubleprecision e0,e1,e2,a,b,c,dtau1,dtau2,dtsum,xp
      doubleprecision src1,src2,src3,theomax,theomin
c
#     ifndef SRCQDR_NO_ERRCHECK
#include "common_radfield.h"      
#     endif
#include "common_check.h"
c
cc#     ifdef SRCQDR_IMPLICIT
cc      write(*,*) 'Implicit quadratic integration not possible.'
cc      write(*,*) 'Undefine SRCQDR_PNT_3.'
cc      stop 13
cc#     endif
c
      dtau1     = 0.5d0 * ( alp1 + alp2 ) * ds1
      dtau2     = 0.5d0 * ( alp2 + alp3 ) * ds2
      dtsum     = dtau1 + dtau2
c
c     If the safety checks are on, then check something out
c
#     ifdef SAFETY_CHECKS_ACTIVE
      if(dtau1.gt.selfch_taumax) then
          selfch_taumax = dtau1
      endif
#     endif
c
c     Now compute the theoretical maximum and minimum of the
c     quadrature. These are used to bracket the Olson-Kunasz formula
c     which may sometimes go wrong.
c     I did this previously in an erronous way.  But since I've used
c     that (wrong) formula for a long time (until 6 July 1999) I rather
c     not remove it completely.  I still keep it an option, just in case
c     old runs may be incompatible with this new thing.
c
#     ifndef SRCQDR_OLDTHEOMAX
c
c     The new one
c
      theomax   = max(js1,js2) * ds1
      theomin   = 0.d0        ! min(js1,js2) * ds1 * exp(-dtau1)
#     else
c
c     The old one ( < 6 july 1999 )
c
      theomax   = 0.5d0 * ( js1 + js2 ) * ds1
      theomin   = 0.d0
#     endif
c
      if(dtau1.gt.1.d-6) then
          xp        = exp(-dtau1)
          e0        = 1.d0 - xp
          e1        = dtau1 - e0
          e2        = dtau1 * dtau1 - 2.d0 * e1
          if(dtau2.gt.1.d-6) then
              b  = ( dtsum * e1 - e2 ) / ( dtau1 * dtau2 )
              a  = e0 + ( e2 - ( dtsum + dtau1 ) * e1 ) 
     %                    / ( dtau1 * dtsum )
              c  = ( e2 - dtau1 * e1 ) / ( dtau2 * dtsum )
          else  
              b  = e1 / dtau1
              a  = e0 - b
              c  = 0.d0
          endif
          if(alp1.gt.0.d0) then
             src1 = js1 / alp1
          elseif(alp2.gt.0.d0) then
             src1 = js2 / alp2
          else
             src1 = 0.d0
          endif
          if(alp2.gt.0.d0) then
             src2 = js2 / alp2
          elseif(alp1.gt.0.d0) then
             src2 = js1 / alp1
          elseif(alp3.gt.0.d0) then
             src2 = js3 / alp3
          else
             src2 = 0.d0
          endif
          if(alp3.gt.0.d0) then
             src3 = js3 / alp3
          elseif(alp2.gt.0.d0) then
             src3 = js2 / alp2
          else
             src3 = 0.d0
          endif
          qdr_src_3 = a * src1 + b * src2 + c * src3
c
c         Verify if serious overshoot occur (say: factor 2 or more)
c
#         ifndef SRCQDR_NO_ERRCHECK
          if(qdr_src_3.gt.2.d0*theomax) then
              transerr_overshoot =  transerr_overshoot + 1
          endif
#         endif
c
c         NEW: Bracket the total integrated emission between the
c              theoretical maximum and minimum.
c
#         ifndef SRCQDR_NO_BRACKET
          qdr_src_3 = min(qdr_src_3,theomax)
          qdr_src_3 = max(qdr_src_3,theomin)
#         endif
c
          qdr_src_3 = int * xp + qdr_src_3
c
      else
c          a         = 0.d0
c          b         = 0.d0
c          c         = 0.d0
c          xp        = 1.d0
         qdr_src_3 = int + 0.5d0 * ( js1 + js2 ) * ds1
      endif
      return
c
      write(*,*) 'Ehh, something wrong in dqr_src_3'
      stop 13
      return
      end






c     --------------------------------------------------------------
c                    WRITE 4-D INTENSITY ARRAY TO FILE
c     --------------------------------------------------------------
      subroutine write_intens()
      implicit none
c
#     ifndef NO_INTENSITY_STORAGE
c
#include "common_grid.h"
#include "common_radfield.h"
c     
      integer ir,it,imu,iphi,iphimax,ithmax,imirth,imirphi,iangset
c     
      iangset = 1
c
      iphimax = nrphi(iangset)
      ithmax = irsi_frsizey
      imirth = 0
      imirphi = 0
c
#ifndef MIRROR_SAVE_ALL
#ifdef MIRROR_THETA
      ithmax = (irsi_frsizey+1)/2
      imirth = 1
#endif
#ifdef MIRROR_PHI
      iphimax = nrphi(iangset)/2
      imirphi = 1
#endif
#ifdef RADGRID_ONEDIM
      iphimax = 1
      ithmax = 1
#endif
      open(unit=1,file='intensity.dat',status='unknown')
      write(1,*) irsi_frsizex,ithmax,
     %            anggrid_frsizemu,iphimax,1,imirth,imirphi
      write(1,*) ' '
      do ir=1,irsi_frsizex
          do it=1,ithmax
              do imu=-nrmu(iangset),nrmu(iangset)
                  if((imu.eq.0).and.(anggrid_mu_zero.eq.0)) then
                      goto 23
                  endif
                  write(1,*) (intens(iphi,imu,it,ir),
     %                        iphi=1,iphimax)
   23             continue
              enddo
              write(1,*) ' '
          enddo
      enddo      
#else
      write(*,*) 'SAVE ALL OPTION DISABLED!'
      stop 1
#endif
      close(1)
c     
#     else /* ifndef NO_INTENSITY_STORAGE */
      write(*,*) 'ERROR: transfer.F/write_intens()'
      write(*,*) '       This function needs the intensity array.'
      write(*,*) '       But this array is not present, because'
      write(*,*) '       NO_INTENSITY_STORAGE is defined.'
      write(*,*) '       Tip: undefine NO_INTENSITY_STORAGE in'
      write(*,*) '       configure.h'
      stop 13
#     endif /* ifndef NO_INTENSITY_STORAGE */
c
      end



c     --------------------------------------------------------------
c                           A DEBUGGING TOOL
c     --------------------------------------------------------------
#ifndef NO_INTENSITY_STORAGE
#ifdef DEBUG_FILLINTENS
      subroutine debug_fill_intens()
      implicit none
#include "common_grid.h"
c#include "common_source.h"
#include "common_radfield.h"
c#include "common_units.h"
c#include "common_boundary.h"
c#include "common_shortchar.h"
c
      doubleprecision int,dr,dt,dmu,dphi,rcrp2
      doubleprecision mui,mus,mured,mu,phi,r0,rscale
      integer ir,it,imu,iphi,is,iphimax,iphimin,imumax,nel,icr,iangset
      integer imus,iphis,ithetamax,iphioff,imu0
      integer iphihalf
c
c     First a hack... For the moment the iangset must be fixed.
c
      iangset  = 1
c
c     Get some important integer values...
c
      ithetamax = irsi_frsizey
      imumax    = nrmu(iangset)
      iphimin   = 1
      iphimax   = nrphi(iangset)
      iphihalf  = nrphi(iangset) / 2
      iphioff   = 0
c
c     Now the mirror things...
c
#ifdef MIRROR_THETA
      ithetamax = (irsi_frsizey+1)/2
#endif
#ifdef MIRROR_PHI
      iphimin = nrphi(iangset) / 4 + 1
      iphimax = 3 * nrphi(iangset) / 4 
      iphioff = - nrphi(iangset) / 4
#endif
#ifdef RADGRID_ONEDIM
      iphimin = 1
      iphimax = 1
      iphioff = 0
      ithetamax = 1
#endif
      do ir=1,irsi_frsizex
          do it=1,ithetamax
              do imu=-imumax,imumax
                  do iphi=iphimin,iphimax
                      intens(iphi+iphioff,imu,it,ir) = 1.d30
                  enddo
              enddo
          enddo
      enddo
c
      end
#endif /* ifdef DEBUG_FILLINTENS */
#endif /* ifndef NO_INTENSITY_STORAGE */



