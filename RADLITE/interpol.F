c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c        ( MODULE: FIND I,S,ALPHA BETWEEN GRIDPOINTS BY INTERPOLATION ) 
c
c                            Leiden, Aug 1998
c                             C.P. Dullemond   
c
c     =============================================================


#define INTERPOL_F_VERSION "Trans_05oct98"
#include "main.h"
#define SPL_IOFF 3
#define FRSIZE_PHI_SPLINE (FRSIZE_PHI+2*SPL_IOFF+1)

c#ifndef NOTRANS



c     --------------------------------------------------------------
c                MAKE THE INDEX ARRAYS FOR INTENSITY DATA
c
c     When mirror symmetries are imposed, the indexation of the 
c     intensity array goes rather different from normal. In order 
c     not to make the program unreadable and slow, we set up two
c     re-indexation arrays. One for iphi and one for itheta.
c
c     NOTE: The iangset argument is actually irrelevant. It could just
c           as well be 1. Note that the angle-grid (mu,phi) is (and
c           must be) the same at all ir! The reason why I have
c           a separate storage for the angle-grid for each iangset is
c           mainly for historic reason, and for a possible future
c           option of alternating mu-gridding. 
c     --------------------------------------------------------------
      subroutine make_index(iangset)
      implicit none
      integer iangset
c
#include "common_grid.h"
c
      integer it,ip,iphihalf,iphiquadr
      integer imu,ii
      doubleprecision drdivr,muu,mu1,mu2
c
#ifndef ONLY_RAY_TRACING
      iphiquadr = nrphi(iangset) / 4
      iphihalf = nrphi(iangset) / 2
      do ip=-nrphi(iangset)-4,2*nrphi(iangset)+4
          ridx_ip(ip,-4) = ip
   11     continue
          if(ridx_ip(ip,-4).lt.1) then
              ridx_ip(ip,-4) = ridx_ip(ip,-4) + nrphi(iangset)
              goto 11
          endif
   12     continue
          if(ridx_ip(ip,-4).gt.nrphi(iangset)) then
              ridx_ip(ip,-4) = ridx_ip(ip,-4) - nrphi(iangset)
              goto 12
          endif
#ifdef MIRROR_PHI
          if(ridx_ip(ip,-4).le.iphiquadr) then
              ridx_ip(ip,-4) = iphihalf + 1 - ridx_ip(ip,-4) 
          endif
          if(ridx_ip(ip,-4).gt.3*iphiquadr) then
              ridx_ip(ip,-4) = 3*iphihalf + 1 - ridx_ip(ip,-4) 
          endif
#endif    
      enddo
#endif /* ONLY_RAY_TRACING */
#ifndef RADGRID_ONEDIM
#ifndef ONLY_RAY_TRACING
      do it=-4,irsi_frsizey+4
          do ip=-nrphi(iangset)-4,2*nrphi(iangset)+4
              ridx_ip(ip,it) = ridx_ip(ip,-4)
          enddo
#ifdef MIRROR_THETA
          if(it.gt.irsi_frsizey/2) then
              do ip=-nrphi(iangset)-4,2*nrphi(iangset)+4
                  ridx_ip(ip,it) = nrphi(iangset) + 1 - ridx_ip(ip,it)
              enddo
          endif
#endif
      enddo
#endif /* ONLY_RAY_TRACING */
      do it=-4,irsi_frsizey+4
          ridx_it(it) = it
          if(ridx_it(it).lt.1) then
              ridx_it(it) = 1 - ridx_it(it)
          endif
          if(ridx_it(it).gt.irsi_frsizey) then
              ridx_it(it) = 2*irsi_frsizey + 1 - ridx_it(it)
          endif
#ifdef MIRROR_THETA
          if(ridx_it(it).gt.irsi_frsizey/2) then
              ridx_it(it) = irsi_frsizey + 1 - ridx_it(it)
          endif
#endif
      enddo
#else /* RADGRID_ONEDIM */
c
c     The case of 1-D Transfer
c
      do it=-4,irsi_frsizey+4
          ridx_it(it) = 1
#ifndef ONLY_RAY_TRACING
          do ip=-nrphi(iangset)-4,2*nrphi(iangset)+4
              ridx_ip(ip,it) = 1
          enddo
#endif
      enddo
#endif
#ifndef ONLY_RAY_TRACING
#ifdef MIRROR_PHI
      do it=-4,irsi_frsizey+4
          do ip=-nrphi(iangset)-4,2*nrphi(iangset)+4
              ridx_ip(ip,it) = ridx_ip(ip,it) - iphiquadr
          enddo
      enddo
#endif
#endif
c
c     Some deleted stuff from old times. Should not be compiled.
c     Or in other words: the INTERPOL_RADIUS_MU_CONSERV should
c     not be defined.
c
#ifdef INTERPOL_RADIUS_MU_CONSERV
#ifdef ESC_MINIMAL_EXTENSION 
      write(*,*) 'INTERNAL ERROR! WARN AUTHOR!'
      write(*,*) 'dullemon@mpa-garching.mpg.de'
      stop 55
#endif
#endif      
c
c     Set the flag that the reindexation has been performed
c
      ridx_ready = 321
c
      return
      end



c     --------------------------------------------------------------
c           FIND THE MIDPOINT INTENSITY BY INTERPOLATION
c
c     This function is able to find the intensity at a given point
c     on the spatial grid, as long as it is on one of the gridlines.
c     It performs a linear interpolation between the two gridpoints.
c     Be sure that the intensity is already known at that point. 
c     If you sweep the grid in the correct way, then this should 
c     indeed be guaranteed. 
c
c     Since the algorithm used in this code requires the intensity
c     only to be interpolated along lines of R=const, the icr is 
c     there only for verification. If this routine finds that icr=2,
c     (meaning that the point lies on a Theta=const line), then it
c     will give an error. 
c
c     This function will do tree kinds of linear interpolations:
c       1) Between gridpoints in Theta coordinate
c       2) Between gridpoint in mu-coordinate (local angle)
c       3) Between gridpoint in phi-coordinate (local angle)
c     Since the algorithm is such that in principle the mu's should
c     connect smoothly to eachother, the mu averaging should not
c     be necessary except for special cases. But for the moment I
c     don't take the risk, and I just leave it in. It will be slower,
c     but for the moment: who cares.
c
c     Also this routine accepts the mu and phi as doubles. In the
c     future the index values and the dmu and dphi should both be
c     specified instead. But also here I just skip that stuff for
c     the moment.
c
c     REMARK: If find_int_t() is used with the MESC algorithm,
c             care should be taken that the value of the intensity
c             at upwind gridpoints is known for ALL phi. For the
c             MESC method this is not guaranteed. While formally
c             they will be known, due to the 3-point or 4-point
c             interpolations in phi, some phi-gridpoints will be
c             included that have not yet been swept by the MESC
c             algorithm before. This is a subtle point, which
c             must be dealt with carefully. This find_int_t()
c             routine always presumes values at upwind gridpoints
c             to be known.
c     
c     --------------------------------------------------------------
      function find_int_t(icr,ir,it,dr,dt,imu,iphi,dmu,dphi)
      implicit none
      doubleprecision find_int_t
      integer icr,ir,it,imu,iphi
      doubleprecision dr,dt,dmu,dphi
c
#     ifndef NO_INTENSITY_STORAGE
c
#include "common_grid.h"
#include "common_radfield.h"
c#include "common_units.h"
c
c#####################
      common/dbirit/db_ir,db_it,db_im,db_ip
      integer db_ir,db_it,db_im,db_ip
      integer indextold
c#####################
      doubleprecision dtr,phdum,int0,int1,int2,int3
      doubleprecision int00,int01,int10,int11
      doubleprecision int20,int21,int30,int31
      doubleprecision mu,phi,dum,dmu1,radius
      doubleprecision thetaclip,phiclip,intclip,dy,theta,rclip
      doubleprecision ddmu,ddmu1,drabs
      dimension thetaclip(0:3),phiclip(0:3),intclip(0:3),rclip(0:2)
      dimension ddmu(0:2),ddmu1(0:2)
      integer idt,it0,it1,idr,ir0,imu0,imu1,iphi0,iphi1,i
      integer iphi00,iphi01,iphi10,iphi11,iangset
      integer indexp,indext,indexr,iymir,iphiquadr,iphihalf,ithhalf
      integer indexmu
      dimension indexp(0:3,0:3),indext(0:3),indexr(0:2)
      dimension indexmu(0:1,0:2)
c
c     Check the initialization of the reindexation.
c     
cc####################
c      write(*,*) '->',db_ir,db_it,db_im,db_ip
cc####################
      if(ridx_ready.ne.321) then
          write(*,*) 'ERROR, find_int_t: Not initialized ',
     %               'the reindexation yet!'
          stop 13
      endif
c
c     For the moment we only allow one global angular set
c
      iangset = 1
c
c     Along which gridline should we interpolate?
c       icr = 1 : Along R=const gridline (interpolate in Theta)
c       icr = 2 : Along Theta=const gridline (interpolate in R)
c
      if(icr.eq.1) then
c
c     ==========================================================
c                  -->  INTERPOLATION IN THETA  <--
c     ==========================================================c
c     Now find the index of the radial shell. This one is 
c     simple, because this function always assumes that we are
c     on a radial shell anyway. No interpolation needed!
c
      idr   = dr
      ir0   = ir + idr
c
c     The mu indices. 
c
      imu0 = imu
      imu1 = imu+1
c
c     Compute the mu itself from the information imu and dmu
c     
      mu   = rmu(imu0,iangset) + dmu * ( rmu(imu1,iangset) 
     %            - rmu(imu0,iangset) )
      dmu1 = 1.d0-dmu
c
c     From here on things depend strongly on the type of interpolation
c     we apply. These can be done to first order or to third order. They
c     can be switched by the defines.
c
c
c     If 1-D then we're done quickly (first order interpolation).
c
#ifdef RADGRID_ONEDIM
      find_int_t = dmu1*intens(1,imu0,1,ir0) + dmu*intens(1,imu1,1,ir0)
      return
#endif      
c
c     For 2-D, there are a lot of different possibilities......
c
#ifdef INTERPOL_ORDER_1
c
c     --------------------------------------------
c     --- FIRST ORDER INTERPOLATION EVERYWHERE ---
c     --------------------------------------------
c     This is a special case: things can be done particularly efficiently.
c     Normally the procedure sets up various auxiliary variables (see
c     the higher order interpolation stuff below), but this time the
c     computational cost is minimized for linear interpolation everywhere.
c
c     Find the the index of the Theta_0 from the dt
c
      idt   = dt                ! Trick: from dbl to int automatic trunc
      it0   = it + idt          ! Note: do not say it0 = it + dt
c
c     Find the real dt (the dt from Theta_1 to Theta, instead
c     of from Theta(it) to Theta)
c
      dtr   = abs(dt - idt)
c
c     Now find the index of the Theta_1. This is always one further
c     in the direction of dtr, except when we are near the z-axis.
c     If we are near the z-axis, I choose the index the same as it0.
c
      if(dt.gt.0.d0) then
          if(it0.lt.irsi_frsizey) then 
              it1 = it0 + 1
          else
              it1 = it0
          endif
      else
          if(it0.gt.1) then           
              it1 = it0 - 1
          else
              it1 = it0
          endif
      endif
c
c     Now the mirror symmetry clausule.
c     
#ifndef MIRROR_PHI
c
c     Find iphi0 and iphi1 and phi for NO mirror symmetry in phi
c
      iphi00 = iphi
      iphi10 = iphi + 1
      phi = rphi(iphi00,iangset) + dphi * 
     %             (rphi(iphi10,iangset)-rphi(iphi00,iangset))
      if(iphi10.gt.nrphi(iangset)) then
          iphi10 = 1
      endif
      if(phi.gt.6.28318530718d0) then
          phi = phi - 6.28318530718d0
      endif
#else
c
c     Find iphi0 and iphi1 and phi, still unshifted. Later, the iphi's
c     are going to be shifted so that the range starts at 1 again. For
c     not let's do as if we don't have mirror symmetry, but simply
c     confine ourselves to the range phi=[0.5\pi,1.5\pi].
c
c     NOTE: nrphi(iangset) should be divisible by 4!!!!
c
      iphiquadr = nrphi(iangset) / 4
      iphi00 = iphi
      iphi10 = iphi + 1
      phi = rphi(iphi00,iangset) + 
     %            dphi * (rphi(iphi10,iangset)-rphi(iphi00,iangset))
      if(iphi10.gt.3*iphiquadr) then
          iphi10 = 3*iphiquadr 
      endif
      if(iphi00.le.iphiquadr) then
          iphi00 = iphiquadr + 1
      endif
#endif
c
      iphi01 = iphi00
      iphi11 = iphi10
c
      if(irsi_frsizey.gt.1) then
#ifdef MIRROR_THETA
          iymir = (irsi_frsizey+1)/2
          if(it0.gt.iymir) then
              it0 = irsi_frsizey - it0 + 1
              iphi00 = nrphi(iangset) + 1 - iphi00
              iphi10 = nrphi(iangset) + 1 - iphi10
          endif
          if(it1.gt.iymir) then
              it1 = irsi_frsizey - it1 + 1
              iphi01 = nrphi(iangset) + 1 - iphi01
              iphi11 = nrphi(iangset) + 1 - iphi11
          endif
#else 
          it0 = it0
#endif
      else
          it0 = 1
          it1 = 1
          dtr = 0.d0
      endif
#ifdef MIRROR_PHI
      iphi00 = iphi00 - iphiquadr
      iphi10 = iphi10 - iphiquadr
      iphi01 = iphi01 - iphiquadr
      iphi11 = iphi11 - iphiquadr
#endif
c
c     Now the averaging over phi 
c
      int00 = (1.d0-dphi)*intens(iphi00,imu0,it0,ir0) + 
     %               dphi*intens(iphi10,imu0,it0,ir0)
      int10 = (1.d0-dphi)*intens(iphi00,imu1,it0,ir0) + 
     %               dphi*intens(iphi10,imu1,it0,ir0)
      int11 = (1.d0-dphi)*intens(iphi01,imu1,it1,ir0) + 
     %               dphi*intens(iphi11,imu1,it1,ir0)
      int01 = (1.d0-dphi)*intens(iphi01,imu0,it1,ir0) + 
     %               dphi*intens(iphi11,imu0,it1,ir0)
c
c     Averaging over mu
c
      int0  = dmu1*int00 + dmu*int10
      int1  = dmu1*int01 + dmu*int11
c
c     Averaging over Theta
c
      find_int_t = (1.d0-dtr)*int0 + dtr*int1
c
#else 
c
c     --------------------------------------------
c     ---  HIGHER ORDER INTERPOLATION METHOD   ---
c     --------------------------------------------
c
c     First make the phi-index stencil. These are simply
c     4 points around iphi.
c
      indexp(0,0) = iphi-1
      indexp(1,0) = iphi
      indexp(2,0) = iphi+1
      indexp(3,0) = iphi+2
c
c     Now find the phis. Since we have ghost-gridpoints on each side of the
c     rphi grid array, we need not worry about periodicity and boundaries
c     here. See the create_anggrid() subroutine in gridchar.F .
c
      phiclip(0) = rphi(indexp(0,0),iangset)
      phiclip(1) = rphi(indexp(1,0),iangset)
      phiclip(2) = rphi(indexp(2,0),iangset)
      phiclip(3) = rphi(indexp(3,0),iangset)
c
c     Reconstruct the phi itself from iphi and dphi.
c
      phi = rphi(indexp(1,0),iangset) + dphi * ( 
     %        rphi(indexp(2,0),iangset) - rphi(indexp(1,0),iangset) )
c
c     Find the the index of the Theta_0 from the dt. The dt is a real
c     number that is composed as:  
c                                    Theta(P_1)-Theta(P_A)
c           dt = (it(P_A)-it(P_0)) + ---------------------
c                                    Theta(P_B)-Theta(P_A)
c
c     We wish to reconstruct the index it(P_A) (which we call indext(1) here)
c     and the precise value of Theta(P_1) (which we call theta).
c
      idt       = dt                ! Trick: from dbl to int automatic trunc
      indext(1) = it + idt          ! Note: do not say indext(1) = it + dt
c
c     Find the real dt (the dt from Theta_1 to Theta, instead
c     of from Theta(it) to Theta)
c
      dtr   = abs(dt - idt)
c
c     Now complete the 4-point Theta index stencil.
c
      if(dt.gt.0.d0) then
          indext(0) = indext(1) - 1
          indext(2) = indext(1) + 1
          indext(3) = indext(1) + 2
      else
          indext(0) = indext(1) + 1
          indext(2) = indext(1) - 1
          indext(3) = indext(1) - 2
      endif
c
c     Find the corresponding values of theta on this stencil, except for
c     the thetaclip(3) because that might go wrong, since the grid rsi_x_c
c     only has one ghost-gridpoints on each border, while indext(3) could
c     be 2 gridpoints off the grid under certain circumstances.
c
      thetaclip(0) = rsi_x_c(indext(0),2)
      thetaclip(1) = rsi_x_c(indext(1),2)
      thetaclip(2) = rsi_x_c(indext(2),2)
c
c     Find these circumstances, and find the thetaclip(3). 
c
      if(indext(2).lt.1) then
          thetaclip(3) = -rsi_x_c(2,2)
      elseif (indext(2).gt.irsi_frsizey) then
          thetaclip(3) = 6.28318530718d0 - rsi_x_c(irsi_frsizey-1,2)
      else
          thetaclip(3) = rsi_x_c(indext(3),2)
      endif      
c
c     Reconstruct the Theta itself from the above
c
      theta = thetaclip(1) + dtr * ( thetaclip(2) - thetaclip(1) )
c
c     Now comes the task of finding the real indices indexp() and indext().
c     So far we just added and subtracted a few counts to find neighboring
c     points, but we did not take into account the boundaries of the
c     domain, which have mirror reflection (for Theta) and periodicity
c     (for phi). Also, when the options MIRROR_PHI and/or MIRROR_THETA
c     are enabled, there are extra mirror reflection criteria to be 
c     handled. The routine make_index() has precalculated these real
c     indices, for each `virtual' index. We merely have to look them up
c     in that table.
c
      indexp(0,1) = ridx_ip(indexp(0,0),indext(1))
      indexp(1,1) = ridx_ip(indexp(1,0),indext(1))
      indexp(2,1) = ridx_ip(indexp(2,0),indext(1))
      indexp(3,1) = ridx_ip(indexp(3,0),indext(1))
      indexp(0,2) = ridx_ip(indexp(0,0),indext(2))
      indexp(1,2) = ridx_ip(indexp(1,0),indext(2))
      indexp(2,2) = ridx_ip(indexp(2,0),indext(2))
      indexp(3,2) = ridx_ip(indexp(3,0),indext(2))
      indexp(0,3) = ridx_ip(indexp(0,0),indext(3))
      indexp(1,3) = ridx_ip(indexp(1,0),indext(3))
      indexp(2,3) = ridx_ip(indexp(2,0),indext(3))
      indexp(3,3) = ridx_ip(indexp(3,0),indext(3))
      indexp(0,0) = ridx_ip(indexp(0,0),indext(0))
      indexp(1,0) = ridx_ip(indexp(1,0),indext(0))
      indexp(2,0) = ridx_ip(indexp(2,0),indext(0))
      indexp(3,0) = ridx_ip(indexp(3,0),indext(0))
      indext(0)   = ridx_it(indext(0))
      indext(1)   = ridx_it(indext(1))
      indext(2)   = ridx_it(indext(2))
      indext(3)   = ridx_it(indext(3))
c
c     ---------------------------------------------------------------
c
c     Here now come the various interpolation methods. They can be
c     switched by use of a simple define in the Makefile
c
c     ---------------------------------------------------------------
c
#ifdef INTERPOL_PHI_3
#ifdef INTERPOL_THETA_3
c
c           ---------------------------------------------
c           |       Theta: 3-rd  ,  phi: 3-rd           |
c           ---------------------------------------------
c
c
c     First averaging over mu
c
      intclip(0) = dmu1*intens(indexp(0,IDX0),imu0,indext(0),ir0)
     %            + dmu*intens(indexp(0,IDX0),imu1,indext(0),ir0) 
      intclip(1) = dmu1*intens(indexp(0,IDX1),imu0,indext(1),ir0)
     %            + dmu*intens(indexp(0,IDX1),imu1,indext(1),ir0) 
      intclip(2) = dmu1*intens(indexp(0,IDX2),imu0,indext(2),ir0)
     %            + dmu*intens(indexp(0,IDX2),imu1,indext(2),ir0) 
      intclip(3) = dmu1*intens(indexp(0,IDX3),imu0,indext(3),ir0)
     %            + dmu*intens(indexp(0,IDX3),imu1,indext(3),ir0) 
      call polint4(thetaclip,intclip,theta,int0)
      intclip(0) = dmu1*intens(indexp(1,IDX0),imu0,indext(0),ir0)
     %            + dmu*intens(indexp(1,IDX0),imu1,indext(0),ir0) 
      intclip(1) = dmu1*intens(indexp(1,IDX1),imu0,indext(1),ir0)
     %            + dmu*intens(indexp(1,IDX1),imu1,indext(1),ir0) 
      intclip(2) = dmu1*intens(indexp(1,IDX2),imu0,indext(2),ir0)
     %            + dmu*intens(indexp(1,IDX2),imu1,indext(2),ir0) 
      intclip(3) = dmu1*intens(indexp(1,IDX3),imu0,indext(3),ir0)
     %            + dmu*intens(indexp(1,IDX3),imu1,indext(3),ir0) 
      call polint4(thetaclip,intclip,theta,int1)
      intclip(0) = dmu1*intens(indexp(2,IDX0),imu0,indext(0),ir0)
     %            + dmu*intens(indexp(2,IDX0),imu1,indext(0),ir0) 
      intclip(1) = dmu1*intens(indexp(2,IDX1),imu0,indext(1),ir0)
     %            + dmu*intens(indexp(2,IDX1),imu1,indext(1),ir0) 
      intclip(2) = dmu1*intens(indexp(2,IDX2),imu0,indext(2),ir0)
     %            + dmu*intens(indexp(2,IDX2),imu1,indext(2),ir0) 
      intclip(3) = dmu1*intens(indexp(2,IDX3),imu0,indext(3),ir0)
     %            + dmu*intens(indexp(2,IDX3),imu1,indext(3),ir0) 
      call polint4(thetaclip,intclip,theta,int2)
      intclip(0) = dmu1*intens(indexp(3,IDX0),imu0,indext(0),ir0)
     %            + dmu*intens(indexp(3,IDX0),imu1,indext(0),ir0) 
      intclip(1) = dmu1*intens(indexp(3,IDX1),imu0,indext(1),ir0)
     %            + dmu*intens(indexp(3,IDX1),imu1,indext(1),ir0) 
      intclip(2) = dmu1*intens(indexp(3,IDX2),imu0,indext(2),ir0)
     %            + dmu*intens(indexp(3,IDX2),imu1,indext(2),ir0) 
      intclip(3) = dmu1*intens(indexp(3,IDX3),imu0,indext(3),ir0)
     %            + dmu*intens(indexp(3,IDX3),imu1,indext(3),ir0) 
      call polint4(thetaclip,intclip,theta,int3)
c
c     Averaging over Theta
c
      intclip(0) = int0
      intclip(1) = int1
      intclip(2) = int2
      intclip(3) = int3
c
c     Now the 3-rd order interpolation over phi, with abs()
c
      call polint4(phiclip,intclip,phi,find_int_t)
      find_int_t = 0.5 * ( find_int_t + abs(find_int_t) )
c
c     Done...
c
#else /* ifdef INTERPOL_THETA_3 */
#ifdef INTERPOL_THETA_2
c
c
c           ---------------------------------------------
c           |       Theta: 2-nd  ,  phi: 3-rd           |
c           ---------------------------------------------
c
c
c     NOTE: We use points (1,2,3) of the clipping (0,1,2,3) in Theta
c           Therefore we shift them from (*,1,2,3) -> (1,2,3,*)
c
      indext(0)    = indext(1)
      indext(1)    = indext(2)
      indext(2)    = indext(3)
      thetaclip(0) = thetaclip(1)
      thetaclip(1) = thetaclip(2)
      thetaclip(2) = thetaclip(3)
c
c     First averaging over mu
c
      intclip(0) = dmu1*intens(indexp(0,IDX0),imu0,indext(0),ir0)
     %            + dmu*intens(indexp(0,IDX0),imu1,indext(0),ir0) 
      intclip(1) = dmu1*intens(indexp(0,IDX1),imu0,indext(1),ir0)
     %            + dmu*intens(indexp(0,IDX1),imu1,indext(1),ir0) 
      intclip(2) = dmu1*intens(indexp(0,IDX2),imu0,indext(2),ir0)
     %            + dmu*intens(indexp(0,IDX2),imu1,indext(2),ir0) 
      call polint3(thetaclip,intclip,theta,int0)
      intclip(0) = dmu1*intens(indexp(1,IDX0),imu0,indext(0),ir0)
     %            + dmu*intens(indexp(1,IDX0),imu1,indext(0),ir0) 
      intclip(1) = dmu1*intens(indexp(1,IDX1),imu0,indext(1),ir0)
     %            + dmu*intens(indexp(1,IDX1),imu1,indext(1),ir0) 
      intclip(2) = dmu1*intens(indexp(1,IDX2),imu0,indext(2),ir0)
     %            + dmu*intens(indexp(1,IDX2),imu1,indext(2),ir0) 
      call polint3(thetaclip,intclip,theta,int1)
      intclip(0) = dmu1*intens(indexp(2,IDX0),imu0,indext(0),ir0)
     %            + dmu*intens(indexp(2,IDX0),imu1,indext(0),ir0) 
      intclip(1) = dmu1*intens(indexp(2,IDX1),imu0,indext(1),ir0)
     %            + dmu*intens(indexp(2,IDX1),imu1,indext(1),ir0) 
      intclip(2) = dmu1*intens(indexp(2,IDX2),imu0,indext(2),ir0)
     %            + dmu*intens(indexp(2,IDX2),imu1,indext(2),ir0) 
      call polint3(thetaclip,intclip,theta,int2)
      intclip(0) = dmu1*intens(indexp(3,IDX0),imu0,indext(0),ir0)
     %            + dmu*intens(indexp(3,IDX0),imu1,indext(0),ir0) 
      intclip(1) = dmu1*intens(indexp(3,IDX1),imu0,indext(1),ir0)
     %            + dmu*intens(indexp(3,IDX1),imu1,indext(1),ir0) 
      intclip(2) = dmu1*intens(indexp(3,IDX2),imu0,indext(2),ir0)
     %            + dmu*intens(indexp(3,IDX2),imu1,indext(2),ir0) 
      call polint3(thetaclip,intclip,theta,int3)
c
c     Averaging over Theta
c
      intclip(0) = int0
      intclip(1) = int1
      intclip(2) = int2
      intclip(3) = int3
c
c     Now the 3-rd order interpolation over phi, with abs()
c
      call polint4(phiclip,intclip,phi,find_int_t)
      find_int_t = 0.5 * ( find_int_t + abs(find_int_t) )
c
c     Done...
c
#else  /* ifdef INTERPOL_THETA_2 */
c
c
c           ---------------------------------------------
c           |       Theta: 1-st  ,  phi: 3-rd           |
c           ---------------------------------------------
c
c
c     First averaging over mu
c
      int00  = dmu1*intens(indexp(0,IDX1),imu0,indext(1),ir0) 
     %        + dmu*intens(indexp(0,IDX1),imu1,indext(1),ir0) 
      int10  = dmu1*intens(indexp(1,IDX1),imu0,indext(1),ir0) 
     %        + dmu*intens(indexp(1,IDX1),imu1,indext(1),ir0) 
      int20  = dmu1*intens(indexp(2,IDX1),imu0,indext(1),ir0) 
     %        + dmu*intens(indexp(2,IDX1),imu1,indext(1),ir0) 
      int30  = dmu1*intens(indexp(3,IDX1),imu0,indext(1),ir0) 
     %        + dmu*intens(indexp(3,IDX1),imu1,indext(1),ir0) 
      int01  = dmu1*intens(indexp(0,IDX2),imu0,indext(2),ir0) 
     %        + dmu*intens(indexp(0,IDX2),imu1,indext(2),ir0) 
      int11  = dmu1*intens(indexp(1,IDX2),imu0,indext(2),ir0) 
     %        + dmu*intens(indexp(1,IDX2),imu1,indext(2),ir0) 
      int21  = dmu1*intens(indexp(2,IDX2),imu0,indext(2),ir0) 
     %        + dmu*intens(indexp(2,IDX2),imu1,indext(2),ir0) 
      int31  = dmu1*intens(indexp(3,IDX2),imu0,indext(2),ir0) 
     %        + dmu*intens(indexp(3,IDX2),imu1,indext(2),ir0) 
c
c     Averaging over Theta
c
      intclip(0) = (1.d0-dtr)*int00 + dtr*int01
      intclip(1) = (1.d0-dtr)*int10 + dtr*int11
      intclip(2) = (1.d0-dtr)*int20 + dtr*int21
      intclip(3) = (1.d0-dtr)*int30 + dtr*int31
c
c     Now the 3-rd order interpolation over phi, with abs()
c
      call polint4(phiclip,intclip,phi,find_int_t)
      find_int_t = 0.5 * ( find_int_t + abs(find_int_t) )
c
c     Done...
c
#endif  /* ifdef else INTERPOL_THETA_2 */
#endif  /* ifdef else INTERPOL_THETA_3 */
#else   /* ifdef INTERPOL_PHI_3 */
#ifndef INTERPOL_THETA_3
#ifndef INTERPOL_THETA_2
c
c
c           ---------------------------------------------
c           |       Theta: 1-st  ,  phi: 1-rd           |
c           ---------------------------------------------
c
      int00 = (1.d0-dphi)*intens(indexp(1,IDX1),imu0,indext(1),ir0) 
     %             + dphi*intens(indexp(2,IDX1),imu0,indext(1),ir0)
      int10 = (1.d0-dphi)*intens(indexp(1,IDX1),imu1,indext(1),ir0)  
     %             + dphi*intens(indexp(2,IDX1),imu1,indext(1),ir0)
      int11 = (1.d0-dphi)*intens(indexp(1,IDX2),imu1,indext(2),ir0) 
     %             + dphi*intens(indexp(2,IDX2),imu1,indext(2),ir0)
      int01 = (1.d0-dphi)*intens(indexp(1,IDX2),imu0,indext(2),ir0)  
     %             + dphi*intens(indexp(2,IDX2),imu0,indext(2),ir0)      
c
c     Averaging over mu
c
      int0  = dmu1*int00 + dmu*int10
      int1  = dmu1*int01 + dmu*int11
c
c     Averaging over Theta
c
      find_int_t = (1.d0-dtr)*int0 + dtr*int1
c
c     Done...
c
#else  /* ifndef INTERPOL_THETA_2 */
      write(*,*) 'This combination (Theta : 2-nd,phi : 1-st)'
      write(*,*) '  is not implemented...'
      write(*,*) '  Recompile with other combination of'
      write(*,*) '  INTERPOL_THETA_2 and INTERPOL_PHI_3'
      stop 13
#endif /* ifndef INTERPOL_THETA_2 */
#else /* ifndef INTERPOL_THETA_3 */
c
c
c           ---------------------------------------------
c           |       Theta: 3-rd  ,  phi: 1-st           |
c           ---------------------------------------------
c
c
c     First averaging over mu
c
      int00  = dmu1*intens(indexp(1,IDX0),imu0,indext(0),ir0) 
     %        + dmu*intens(indexp(1,IDX0),imu1,indext(0),ir0) 
      int10  = dmu1*intens(indexp(1,IDX1),imu0,indext(1),ir0) 
     %        + dmu*intens(indexp(1,IDX1),imu1,indext(1),ir0) 
      int20  = dmu1*intens(indexp(1,IDX2),imu0,indext(2),ir0) 
     %        + dmu*intens(indexp(1,IDX2),imu1,indext(2),ir0) 
      int30  = dmu1*intens(indexp(1,IDX3),imu0,indext(3),ir0) 
     %        + dmu*intens(indexp(1,IDX3),imu1,indext(3),ir0) 
      int01  = dmu1*intens(indexp(2,IDX0),imu0,indext(0),ir0) 
     %        + dmu*intens(indexp(2,IDX0),imu1,indext(0),ir0) 
      int11  = dmu1*intens(indexp(2,IDX1),imu0,indext(1),ir0) 
     %        + dmu*intens(indexp(2,IDX1),imu1,indext(1),ir0) 
      int21  = dmu1*intens(indexp(2,IDX2),imu0,indext(2),ir0) 
     %        + dmu*intens(indexp(2,IDX2),imu1,indext(2),ir0) 
      int31  = dmu1*intens(indexp(2,IDX3),imu0,indext(3),ir0) 
     %        + dmu*intens(indexp(2,IDX3),imu1,indext(3),ir0) 
c
c     Averaging over Phi
c
      intclip(0) = (1.d0-dphi)*int00 + dphi*int01
      intclip(1) = (1.d0-dphi)*int10 + dphi*int11
      intclip(2) = (1.d0-dphi)*int20 + dphi*int21
      intclip(3) = (1.d0-dphi)*int30 + dphi*int31
c
c     Now the 3-rd order interpolation over Theta, with abs()
c
      call polint4(thetaclip,intclip,theta,find_int_t)
      find_int_t = 0.5 * ( find_int_t + abs(find_int_t) )
c
c     Done...
c
#endif  /* ifndef else INTERPOL_THETA_3 */
#endif  /* ifdef else INTERPOL_PHI_3 */
c
#endif
      else
c     ==========================================================
c                  -->  INTERPOLATION IN RADIUS  <--
c     ==========================================================
#ifdef RADGRID_ONEDIM
      write(*,*) 'In 1-D the icr must only be 1'
      stop 13
#endif
#ifndef ESC_MINIMAL_EXTENSION
      write(*,*) 'When the full ESC method is used, the ',
     %           'find_int_t() may never be called with icr.ne.1'
      stop 13
#endif
c
c     Find the index of the theta coordinate
c
      idt       = dt                ! Trick: from dbl to int automatic trunc
      indext(0) = it + idt          ! Note: do not say indext(0) = it + dt
c
c     Find the index of the radial shell. 
c
      ir0   = ir
      drabs = abs(dr)
c
c     The radius stencil. If dr.lt.0 then the 3 radial points
c     are assumed to be (ir0-2,ir0-1,ir0) otherwise they are
c     assumed to be (ir0,ir0+1,ir0+2). 
c     
      if(dr.lt.0.d0) then
          indexr(0) = ir0
          indexr(1) = ir0-1
          indexr(2) = ir0-2
      else
          indexr(0) = ir0
          indexr(1) = ir0+1
          indexr(2) = ir0+2          
      endif
c
c     Now find the radii
c
      rclip(0) = rsi_x_c(indexr(0),1)
      rclip(1) = rsi_x_c(indexr(1),1)
      rclip(2) = rsi_x_c(indexr(2),1)
c
c     The real radius
c
      radius = rclip(0) + drabs * ( rclip(1) - rclip(0) )
c
c     Check if we are not going over array boundaries
c
      if(indexr(2).lt.1) then
          indexr(2) = 1
      endif
      if(indexr(2).gt.irsi_frsizex) then
          indexr(2) = irsi_frsizex
      endif
c
c     The mu indices
c
      imu0 = imu
      imu1 = imu+1
c
c     Some deleted stuff
c
#ifdef INTERPOL_RADIUS_MU_CONSERV
      write(*,*) 'INTERNAL ERROR! WARN AUTHOR!'
      write(*,*) 'dullemon@mpa-garching.mpg.de'
      stop 55
#endif
c
c     Compute the mu itself from the information imu and dmu
c     
      mu   = rmu(imu0,iangset) + dmu * ( rmu(imu1,iangset) 
     %                 - rmu(imu0,iangset) )
      dmu1 = 1.d0-dmu
c
c     Make the phi-index stencil. These are simply
c     4 points around iphi.
c
      indexp(0,0) = iphi-1
      indexp(1,0) = iphi
      indexp(2,0) = iphi+1
      indexp(3,0) = iphi+2
c
c     Now find the phis. Since we have ghost-gridpoints on each side of the
c     rphi grid array, we need not worry about periodicity and boundaries
c     here. See the create_anggrid() subroutine in gridchar.F .
c
      phiclip(0) = rphi(indexp(0,0),iangset)
      phiclip(1) = rphi(indexp(1,0),iangset)
      phiclip(2) = rphi(indexp(2,0),iangset)
      phiclip(3) = rphi(indexp(3,0),iangset)
c
c     Reconstruct the phi itself from iphi and dphi.
c
      phi = rphi(indexp(1,0),iangset) + dphi * ( 
     %        rphi(indexp(2,0),iangset) - rphi(indexp(1,0),iangset) )
c
c     Now comes the task of finding the real indices indexp() and indext().
c     So far we just added and subtracted a few counts to find neighboring
c     points, but we did not take into account the boundaries of the
c     domain, which have mirror reflection (for Theta) and periodicity
c     (for phi). Also, when the options MIRROR_PHI and/or MIRROR_THETA
c     are enabled, there are extra mirror reflection criteria to be 
c     handled. The routine make_index() has precalculated these real
c     indices, for each `virtual' index. We merely have to look them up
c     in that table.
c
      indexp(0,0) = ridx_ip(indexp(0,0),indext(0))
      indexp(1,0) = ridx_ip(indexp(1,0),indext(0))
      indexp(2,0) = ridx_ip(indexp(2,0),indext(0))
      indexp(3,0) = ridx_ip(indexp(3,0),indext(0))
      indextold   = indext(0)
      indext(0)   = ridx_it(indext(0))
c
c     ---------------------------------------------------------------
c
c     Here now come the various interpolation methods. They can be
c     switched by use of a simple define in the Makefile
c
c     ---------------------------------------------------------------
c
#ifdef INTERPOL_PHI_3
#ifdef INTERPOL_RADIUS_2
c
c           ---------------------------------------------
c           |       Radius: 2-nd  ,  phi: 3-rd          |
c           ---------------------------------------------
c
#ifdef INTERPOL_RADIUS_MU_CONSERV
      write(*,*) 'FLUX-CONS RADIAL INTERPOL NOT WORKING'
      write(*,*) '2ND ORDER IN R AND 3-RD IN PHI'
#endif
c
c
c     First averaging over mu
c
      intclip(0) = dmu1*intens(indexp(0,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(0,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(0,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(0,0),imu1,indext(0),indexr(1)) 
      intclip(2) = dmu1*intens(indexp(0,0),imu0,indext(0),indexr(2))
     %            + dmu*intens(indexp(0,0),imu1,indext(0),indexr(2)) 
      call polint3(rclip,intclip,radius,int0)
      intclip(0) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(1)) 
      intclip(2) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(2))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(2)) 
      call polint3(rclip,intclip,radius,int1)
      intclip(0) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(1)) 
      intclip(2) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(2))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(2)) 
      call polint3(rclip,intclip,radius,int2)
      intclip(0) = dmu1*intens(indexp(3,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(3,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(3,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(3,0),imu1,indext(0),indexr(1)) 
      intclip(2) = dmu1*intens(indexp(3,0),imu0,indext(0),indexr(2))
     %            + dmu*intens(indexp(3,0),imu1,indext(0),indexr(2)) 
      call polint3(rclip,intclip,radius,int3)
c
      intclip(0) = int0
      intclip(1) = int1
      intclip(2) = int2
      intclip(3) = int3
c
c     Now the 3-rd order interpolation over phi, with abs()
c
      call polint4(phiclip,intclip,phi,find_int_t)
      find_int_t = 0.5 * ( find_int_t + abs(find_int_t) )
#else /* ifdef INTERPOL_RADIUS_2 */
c
c
c           ---------------------------------------------
c           |       Radius: 1-st  ,  phi: 3-rd          |
c           ---------------------------------------------
c
#ifdef INTERPOL_RADIUS_MU_CONSERV
      write(*,*) 'FLUX-CONS RADIAL INTERPOL NOT WORKING'
      write(*,*) '3ND ORDER IN PHI'
#endif
c
c     First averaging over mu
c
      intclip(0) = dmu1*intens(indexp(0,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(0,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(0,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(0,0),imu1,indext(0),indexr(1)) 
      int0 = (1.d0-drabs) * intclip(0) + drabs * intclip(1)
      intclip(0) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(1)) 
      int1 = (1.d0-drabs) * intclip(0) + drabs * intclip(1)
      intclip(0) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(1)) 
      int2 = (1.d0-drabs) * intclip(0) + drabs * intclip(1)
      intclip(0) = dmu1*intens(indexp(3,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(3,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(3,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(3,0),imu1,indext(0),indexr(1)) 
      int3 = (1.d0-drabs) * intclip(0) + drabs * intclip(1)
c
      intclip(0) = int0
      intclip(1) = int1
      intclip(2) = int2
      intclip(3) = int3
c
c     Now the 3-rd order interpolation over phi, with abs()
c
      call polint4(phiclip,intclip,phi,find_int_t)
      find_int_t = 0.5 * ( find_int_t + abs(find_int_t) )
c
#endif  /* ifdef INTERPOL_RADIUS_2 else */
#else   /* ifdef INTERPOL_PHI_3 */
#ifdef INTERPOL_RADIUS_2       
c
c           ---------------------------------------------
c           |       Radius: 2-nd  ,  phi: 1-st          |
c           ---------------------------------------------
c
#ifdef INTERPOL_RADIUS_MU_CONSERV
      write(*,*) 'FLUX-CONS RADIAL INTERPOL NOT WORKING'
      write(*,*) '2ND ORDER IN R'
#endif
c
c     First averaging over mu
c
      intclip(0) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(1)) 
      intclip(2) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(2))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(2)) 
      call polint3(rclip,intclip,radius,int1)
      intclip(0) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(1)) 
      intclip(2) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(2))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(2)) 
      call polint3(rclip,intclip,radius,int2)
c
c     Now the 1-st order interpolation over phi, with abs()
c
      find_int_t = (1.d0-dphi)*int1 + dphi*int2
      find_int_t = 0.5 * ( find_int_t + abs(find_int_t) )
c
#else  /* ifdef INTERPOL_PHI_3 */
#ifdef INTERPOL_RADIUS_2
      write(*,*) 'NO 2ND RADIUS, 1ST PHI INTERPOLATION'
      stop 13
#else
c
c           ---------------------------------------------
c           |       Radius: 1-st  ,  phi: 1-st          |
c           ---------------------------------------------
c
c     First averaging over mu
c
#ifndef INTERPOL_RADIUS_MU_CONSERV
      intclip(0) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(0))
      intclip(1) = dmu1*intens(indexp(1,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(1,0),imu1,indext(0),indexr(1)) 
      int1 = (1.d0-drabs) * intclip(0) + drabs * intclip(1)
      intclip(2) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(0))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(0))
      intclip(3) = dmu1*intens(indexp(2,0),imu0,indext(0),indexr(1))
     %            + dmu*intens(indexp(2,0),imu1,indext(0),indexr(1)) 
      int2 = (1.d0-drabs) * intclip(2) + drabs * intclip(3)
#else
      write(*,*) 'INTERNAL ERROR! WARN AUTHOR!'
      write(*,*) 'dullemon@mpa-garching.mpg.de'
      stop 55
#endif
c
c     Now the 1-st order interpolation over phi, with abs()
c
      find_int_t = (1.d0-dphi)*int1 + dphi*int2
      find_int_t = 0.5 * ( find_int_t + abs(find_int_t) )
c
#endif /* ifdef INTERPOL_RADIUS_2 */
#endif /* ifdef INTERPOL_PHI_3 else */
#endif
      endif ! if(icr.eq.1) .. else 
c
      return
c
#     else /* ifndef NO_INTENSITY_STORAGE */
      write(*,*) 'ERROR: interpol.F/find_int_t()'
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
c               FIND THE MIDPOINT SOURCE BY INTERPOLATION
c
c     This routine is able to find the value of the source term at
c     any given position 'P' at any given angle. It does so by 
c     interpolation. The ir and it are the position of the point '0'
c     (which can be miles away from the point 'P' in principle).
c     The dr and dt contain the offset of 'P' from '0'. Its integer
c     root is the amount of gridpoints away from '0', while its residu
c     is the position within the grid interval. Example: ir=4, dr=3.5
c     means that '0' is at ir=4, so that 'P' is right in the middle
c     between ir=7 and ir=8. 
c
c     ARGUMENTS:
c       inu          = Index of frequency
c       icr          = 1 : Crossing with R=const grid line
c                    = 2 : Crossing with Theta=const grid line
c                    = 0 : Point lies within gridcell, not on grid line
c       ir           = Index of nearest R=const grid line in negative dir
c       it           = Index of nearest Theta=const grid line in negative dir
c       dr           = [0..1] : 0 = at R=R_ir
c                               1 = at R=R_{ir+1}
c       dt           = [0..1] : 0 = at Theta=Theta_ir , 
c                               1 = at Theta=Theta_{ir+1}
c       mu           = Angle mu
c       phi          = Angle phi
c
c     REQUIRES:
c       This routine requires the re-indexation table to have been
c       initialized. This is done by calling make_index(). This table
c       incorporates all mirror symmetries and boundaries in the 
c       spatial coordinate Theta and angular coordinate phi.
c
c     --------------------------------------------------------------
      function find_src_t(inu,icr,ir,it,dr,dt,mu,phi)
      implicit none
      doubleprecision find_src_t,calc_src_pnt
      integer icr,ir,it,inu
      doubleprecision dr,dt,mu,phi
c
#include "common_grid.h"
c
      integer idr,idt,ir0,it0,iymir
      doubleprecision drr,dtr,src,dy
      dimension src(0:3)
      integer indexr,indext
      doubleprecision thetaclip,theta,rclip,r
      dimension thetaclip(0:3),rclip(0:3)
      dimension indexr(0:3),indext(0:3)
c
c     Check the initialization of the reindexation.
c     
      if(ridx_ready.ne.321) then
          write(*,*) 'ERROR, find_src_t: Not initialized ',
     %               'the reindexation yet!'
          stop 13
      endif
c
c     If 1-D then we're done quickly, since it must be on a gridpoint
c
#ifdef RADGRID_ONEDIM
      if(icr.eq.2) then
          write(*,*) 'Inconsistency in 1-D option find_src_t()'
          stop 1
      endif
      if(icr.eq.1) then
          idr   = dr            ! Trick: from dbl to int automatic trunc
          ir0   = ir + idr      ! Note: do not say ir0 = ir + dr        
          find_src_t = calc_src_pnt(inu,ir0,1,mu,phi)
          return
      else
#ifdef SCONEDIM_ADDPOINT
c
c         If icr.eq.0 in 1-D, then follow normal procedure
c     
          if((dt.ne.0.d0).or.(it.ne.1)) then
              write(*,*) 'Internal error with SC 1-D addpoint'
              stop 63453
          endif
#else
c
c         If icr.eq.0 in 1-D, but SCONEDIM_ADDPOINT not set,
c         then error
c
          write(*,*) 'PROBLEM: in 1-D the icr=3 option is not '
          write(*,*) '         yet available for find_src_t()'
          write(*,*) '         If this problem occured while '
          write(*,*) '         rendering spectrum or image in 1-D'
          write(*,*) '         then it is advised to undefine '
          write(*,*) '         LONGCHAR_EXTRA_POINTS'
          stop 1
#endif
      endif
#endif
c
c     For 2-D we proceed doing the more complicated interpolation stuff
c
#ifndef INTERPOL_SRC_2
#define INTERPOL_SRC_WIDE
#endif
#ifdef INTERPOL_SRC_0
#define INTERPOL_SRC_PNT_2
#endif
#ifdef INTERPOL_SRC_1
#define INTERPOL_SRC_PNT_2
#endif
c
c     Find the index of the radial shell from the dr 
c
      idr   = dr                ! Trick: from dbl to int automatic trunc
      ir0   = ir + idr          ! Note: do not say ir0 = ir + dr        
c
c     Find the real dr (the dr from R_1 to R, instead of from R_0 to R)
c
      drr   = abs(dr - idr)
c
c     Now find the index of the R_1. This is always one further
c     in the direction of drr, except when we are at R_out or R_in.
c     If so, I choose the index the same as ir0.
c
      if(dr.gt.0.d0) then
#ifdef INTERPOL_SRC_WIDE
          indexr(0) = ir0
          indexr(1) = ir0 + 1
          indexr(2) = ir0 + 2
          rclip(0)  = rsi_x_c(indexr(0),1)
          rclip(1)  = rsi_x_c(indexr(1),1)
          rclip(2)  = rsi_x_c(indexr(2),1)
          if(indexr(1).gt.irsi_frsizex) then
              indexr(1) = irsi_frsizex
          endif
          if(indexr(2).gt.irsi_frsizex) then
              indexr(2) = irsi_frsizex
          endif
#else
          indexr(0) = ir0 - 1
          indexr(1) = ir0 
          indexr(2) = ir0 + 1
          rclip(0)  = rsi_x_c(indexr(0),1)
          rclip(1)  = rsi_x_c(indexr(1),1)
          rclip(2)  = rsi_x_c(indexr(2),1)
          if(indexr(0).lt.1) then
              indexr(0) = 1
          endif
          if(indexr(2).gt.irsi_frsizex) then
              indexr(2) = irsi_frsizex
          endif
#endif          
      else
#ifdef INTERPOL_SRC_WIDE
          indexr(0) = ir0
          indexr(1) = ir0 - 1
          indexr(2) = ir0 - 2
          rclip(0)  = rsi_x_c(indexr(0),1)
          rclip(1)  = rsi_x_c(indexr(1),1)
          rclip(2)  = rsi_x_c(indexr(2),1)
          if(indexr(1).lt.1) then
              indexr(1) = 1
          endif
          if(indexr(2).lt.1) then
              indexr(2) = 1
          endif
#else
          indexr(0) = ir0 + 1
          indexr(1) = ir0 
          indexr(2) = ir0 - 1
          rclip(0)  = rsi_x_c(indexr(0),1)
          rclip(1)  = rsi_x_c(indexr(1),1)
          rclip(2)  = rsi_x_c(indexr(2),1)
          if(indexr(0).gt.irsi_frsizex) then
              indexr(0) = irsi_frsizex
          endif
          if(indexr(2).lt.1) then
              indexr(2) = 1
          endif
#endif          
      endif
c
c     Reconstruct the R itself from the above
c
#ifdef INTERPOL_SRC_WIDE
      r = rclip(0) + drr * ( rclip(1) - rclip(0) )
#else
      r = rclip(1) + drr * ( rclip(2) - rclip(1) )
#endif
c
c     Only in 2-D:
c
#ifndef RADGRID_ONEDIM
c
c     Find the the index of the Theta_0 from the dt
c
      idt   = dt                ! Trick: from dbl to int automatic trunc
      it0   = it + idt          ! Note: do not say it0 = it + dt        
c
c     Find the real dt (the dt from Theta_1 to Theta, instead
c     of from Theta(it) to Theta)
c
      dtr   = abs(dt - idt)
c
c     Now complete the 3-point Theta index stencil.
c
      if(dt.gt.0.d0) then
#ifdef INTERPOL_SRC_WIDE
          indext(0) = it0 
          indext(1) = it0 + 1
          indext(2) = it0 + 2
#else
          indext(0) = it0 - 1
          indext(1) = it0 
          indext(2) = it0 + 1
#endif
      else            
#ifdef INTERPOL_SRC_WIDE
          indext(0) = it0 
          indext(1) = it0 - 1
          indext(2) = it0 - 2
#else
          indext(0) = it0 + 1
          indext(1) = it0 
          indext(2) = it0 - 1
#endif
      endif
c
c     Find the corresponding values of theta on this stencil, except for
c     the thetaclip(2) because that might go wrong, since the grid rsi_x_c
c     only has one ghost-gridpoints on each border, while indext(2) could
c     be 2 gridpoints off the grid under certain circumstances.
c
      thetaclip(0) = rsi_x_c(indext(0),2)
      thetaclip(1) = rsi_x_c(indext(1),2)
c
c     Find these circumstances, and find the thetaclip(2). 
c
      if(indext(2).lt.0) then
          thetaclip(2) = -rsi_x_c(2,2)
      elseif (indext(2).gt.irsi_frsizey+1) then
          thetaclip(2) = 6.28318530718d0 - rsi_x_c(irsi_frsizey-1,2)
      else
          thetaclip(2) = rsi_x_c(indext(2),2)
      endif      
c
c     Reconstruct the Theta itself from the above
c
#ifdef INTERPOL_SRC_WIDE
      theta = thetaclip(0) + dtr * ( thetaclip(1) - thetaclip(0) )
#else
      theta = thetaclip(1) + dtr * ( thetaclip(2) - thetaclip(1) )
#endif
c
c     Remap the theta indices in case of boundaries crossed, or a 
c     mirror line is crossed.
c
      indext(0)   = ridx_it(indext(0))
      indext(1)   = ridx_it(indext(1))
      indext(2)   = ridx_it(indext(2))
c
#else /* ifndef RADGRID_ONEDIM */
c
c     In case of 1-D, but midpoint interpolation
c
      indext(0)    = 1
      indext(1)    = 1
      thetaclip(0) = rsi_x_c(indext(0),2)
      thetaclip(1) = rsi_x_c(indext(1),2)
      theta        = thetaclip(0)
c
#endif /* ifndef RADGRID_ONEDIM */
c     
c     In case of 0-th order interpolation, do something about 
c     the drr and dtr
c
#ifdef INTERPOL_SRC_0
      if(drr.gt.0.5d0) then 
          drr=1.d0 
      else 
          drr=0.d0 
      endif
      if(dtr.gt.0.5d0) then 
          dtr=1.d0 
      else 
          dtr=0.d0 
      endif
#endif
c
c     Now that we have bracketed the point at which we wish to find
c     the source term, let us make distinction between three types 
c     of interpolation. icr=1 means that the point lies on a R=const
c     grid line, icr=2 means that the point lies on a Theta=const grid
c     line, and icr=0 means that the point really lies somewhere in
c     between both the Theta and R gridlines.
c
#ifdef INTERPOL_SRC_PNT_2
      if(icr.eq.1) then
c
c         The point 'P' lies on an R=const gridline
c
          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
          src(1) = calc_src_pnt(inu,indexr(0),indext(1),mu,phi)
          find_src_t = (1.d0-dtr)*src(0) + dtr*src(1)
c
      elseif(icr.eq.2) then
c
c         The point 'P' lies on an Theta=const gridline
c
          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
          src(1) = calc_src_pnt(inu,indexr(1),indext(0),mu,phi)
          find_src_t = (1.d0-drr)*src(0) + drr*src(1)
c
      else
c
c         The point 'P' lies somewhere in between all gridlines
c
          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
          src(1) = calc_src_pnt(inu,indexr(0),indext(1),mu,phi)
          src(2) = calc_src_pnt(inu,indexr(1),indext(0),mu,phi)
          src(3) = calc_src_pnt(inu,indexr(1),indext(1),mu,phi)
          find_src_t = (1.d0-drr) * ((1.d0-dtr)*src(0) + dtr*src(1)) 
     %                +       drr * ((1.d0-dtr)*src(2) + dtr*src(3)) 
      endif
#else
      if(icr.eq.1) then
c
c         The point 'P' lies on an R=const gridline
c
          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
          src(1) = calc_src_pnt(inu,indexr(0),indext(1),mu,phi)
          src(2) = calc_src_pnt(inu,indexr(0),indext(2),mu,phi)
          call polint3(thetaclip,src,theta,find_src_t)
c
      elseif(icr.eq.2) then
c
c         The point 'P' lies on an Theta=const gridline
c
          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
          src(1) = calc_src_pnt(inu,indexr(1),indext(0),mu,phi)
          src(2) = calc_src_pnt(inu,indexr(2),indext(0),mu,phi)
          call polint3(rclip,src,r,find_src_t)
c
      else
c
c         The point 'P' lies somewhere in between all gridlines
c
#ifdef RADGRID_ONEDIM
c
c         In 1-D this can be done like if it were a Theta=const line
c
          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
          src(1) = calc_src_pnt(inu,indexr(1),indext(0),mu,phi)
          src(2) = calc_src_pnt(inu,indexr(2),indext(0),mu,phi)
          call polint3(rclip,src,r,find_src_t)
#else
          write(*,*) 'ERROR: interpol.F/find_src_t()'
          write(*,*) '       icr=0: This option is not yet ',
     %      'available for 2nd order interpolation in 2-D'
          stop 14
#endif
      endif
c
c     Make sure that the answer is positive definite
c
      find_src_t = 0.5d0 * ( find_src_t + abs(find_src_t) )
#endif
c
      return
      end



c     --------------------------------------------------------------
c               FIND THE MIDPOINT ALPHA BY INTERPOLATION
c
c     This routine is able to find the value of alpha at
c     any given position 'P' at any given angle. See the similar
c     routine find_src_t() for details. 
c     --------------------------------------------------------------
      function find_alp_t(inu,icr,ir,it,dr,dt,mu,phi)
      implicit none
      doubleprecision find_alp_t,calc_alp_pnt
      integer icr,ir,it,inu
      doubleprecision dr,dt,mu,phi
c
#include "common_grid.h"
c
      integer idr,idt,ir0,ir1,it0,it1,iymir
      doubleprecision drr,dtr,alp
      dimension alp(0:3)
c
c     If 1-D then we're done quickly, since it must be on a gridpoint
c
#ifdef RADGRID_ONEDIM
      if(icr.eq.2) then
          write(*,*) 'Inconsistency in 1-D option find_alp_t()'
          stop 1
      endif
      if(icr.eq.1) then
          idr   = dr            ! Trick: from dbl to int automatic trunc
          ir0   = ir + idr      ! Note: do not say ir0 = ir + dr        
          find_alp_t = calc_alp_pnt(inu,ir0,1,mu,phi)
          return
      else
          if((dt.ne.0.d0).or.(it.ne.1)) then
              write(*,*) 'Internal error with SC 1-D addpoint'
              stop 63455
          endif
      endif
#endif
c
c     For 2-D we proceed doing the more complicated interpolation stuff
c
c
c     Find the index of the radial shell from the dr 
c
      idr   = dr                ! Trick: from dbl to int automatic trunc
      ir0   = ir + idr          ! Note: do not say ir0 = ir + dr        
c
c     Find the real dr (the dr from R_1 to R, instead of from R_0 to R)
c
      drr   = abs(dr - idr)
c
c     Now find the index of the R_1. This is always one further
c     in the direction of drr, except when we are at R_out or R_in.
c     If so, I choose the index the same as ir0.
c
      if(dr.gt.0.d0) then
          if(ir0.lt.irsi_frsizex) then 
              ir1 = ir0 + 1
          else
              ir1 = ir0
          endif
      else
          if(ir0.gt.1) then           
              ir1 = ir0 - 1
          else
              ir1 = ir0
          endif
      endif
c
c     Only in 2-D:
c
#ifndef RADGRID_ONEDIM
c
c     Find the the index of the Theta_0 from the dt
c
      idt   = dt                ! Trick: from dbl to int automatic trunc
      it0   = it + idt          ! Note: do not say it0 = it + dt        
c
c     Find the real dt (the dt from Theta_1 to Theta, instead
c     of from Theta(it) to Theta)
c
      dtr   = abs(dt - idt)
c
c     Now find the index of the Theta_1. This is always one further
c     in the direction of dtr, except when we are near the z-axis.
c     If we are near the z-axis, I choose the index the same as it0.
c
      if(dt.gt.0.d0) then
          if(it0.lt.irsi_frsizey) then 
              it1 = it0 + 1
          else
              it1 = it0
          endif
      else
          if(it0.gt.1) then           
              it1 = it0 - 1
          else
              it1 = it0
          endif
      endif
c
c     Now the mirror stuff. Only the Theta mirror needs to be
c     installed for this isotropic source term.
c
#ifdef MIRROR_THETA
      iymir = (irsi_frsizey+1)/2
      if(it0.gt.iymir) then
          it0 = irsi_frsizey - it0 + 1
      endif
      if(it1.gt.iymir) then
          it1 = irsi_frsizey - it1 + 1
      endif
#endif
c
#else /* ifndef RADGRID_ONEDIM */
      it0 = 1
      it1 = 1
      dtr = 0.d0
#endif /* ifndef RADGRID_ONEDIM */
c
c     For 0-th order interpolation we do an if-statement
c
c#######################
c#######################
#ifdef INTERPOL_SRC_0
      if(drr.gt.0.5d0) then 
          drr=1.d0 
      else 
          drr=0.d0 
      endif
      if(dtr.gt.0.5d0) then 
          dtr=1.d0 
      else 
          dtr=0.d0 
      endif
#endif
c
c     Now that we have bracketed the point at which we wish to find
c     the source term, let us make distinction between three types 
c     of interpolation. icr=1 means that the point lies on a R=const
c     grid line, icr=2 means that the point lies on a Theta=const grid
c     line, and icr=0 means that the point really lies somewhere in
c     between both the Theta and R gridlines.
c
      if(icr.eq.1) then
c
c         The point 'P' lies on an R=const gridline
c
          alp(0) = calc_alp_pnt(inu,ir0,it0,mu,phi)
c         NaN here
          alp(1) = calc_alp_pnt(inu,ir0,it1,mu,phi)
          find_alp_t = (1.d0-dtr)*alp(0) + dtr*alp(1)
c
      elseif(icr.eq.2) then
c
c         The point 'P' lies on an Theta=const gridline
c
          alp(0) = calc_alp_pnt(inu,ir0,it0,mu,phi)
          alp(1) = calc_alp_pnt(inu,ir1,it0,mu,phi)
          find_alp_t = (1.d0-drr)*alp(0) + drr*alp(1)
c
      else
c
c         The point 'P' lies somewhere in between all gridlines
c
          alp(0) = calc_alp_pnt(inu,ir0,it0,mu,phi)
          alp(1) = calc_alp_pnt(inu,ir0,it1,mu,phi)
          alp(2) = calc_alp_pnt(inu,ir1,it0,mu,phi)
          alp(3) = calc_alp_pnt(inu,ir1,it1,mu,phi)
          find_alp_t = (1.d0-drr) * ((1.d0-dtr)*alp(0) + dtr*alp(1)) 
     %                +       drr * ((1.d0-dtr)*alp(2) + dtr*alp(3)) 
      endif
c
      return
      end



#ifdef INCLUDE_LINES
#ifdef LINE_SRC_INTERPOL_LEVPOP
c     --------------------------------------------------------------
c               FIND THE MIDPOINT LEVELPOPS, DENSITY AND TEMP
c
c     This routine performs the interpolation for the fluid-frame
c     quanties N,n_i,T_kin. It is important to do the interpolation
c     on the fluid frame quantities for line transfer, because 
c     the direction of the velocity with respect to the characteristic
c     may change very quickly.
c
c     ARGUMENTS:
c       inu          = Index of frequency
c       icr          = 1 : Crossing with R=const grid line
c                    = 2 : Crossing with Theta=const grid line
c                    = 0 : Point lies within gridcell, not on grid line
c       ir           = Index of nearest R=const grid line in negative dir
c       it           = Index of nearest Theta=const grid line in negative dir
c       dr           = [0..1] : 0 = at R=R_ir
c                               1 = at R=R_{ir+1}
c       dt           = [0..1] : 0 = at Theta=Theta_ir , 
c                               1 = at Theta=Theta_{ir+1}
c
c     REQUIRES:
c       This routine requires the re-indexation table to have been
c       initialized. This is done by calling make_index(). This table
c       incorporates all mirror symmetries and boundaries in the 
c       spatial coordinate Theta and angular coordinate phi.
c
c     NOTE: 
c       This routine was make for a future implementation. Right now
c       I already solved the problem, which prompted this, by another
c       means. But still, in principle, in the future I should indeed
c       interpolate the sources using this routine, instead of the
c       direct interpolation (13-May-1999)
c
c     --------------------------------------------------------------
      subroutine find_levpop_and_temp(iline,icr,ir,it,dr,dt,
     %             nup,ndown,temp)
      implicit none
      doubleprecision nup,ndown,temp
      integer icr,ir,it,iline
      doubleprecision dr,dt
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_lines.h"
c
      integer idr,idt,ir0,it0,iymir
      doubleprecision drr,dtr,dy
      doubleprecision t_nup(0:3),t_ndown(0:3),t_temp(0:3)
      integer indexr,indext
      doubleprecision thetaclip,theta,rclip,r
      dimension thetaclip(0:3),rclip(0:3)
      dimension indexr(0:3),indext(0:3)
c
c     Check the initialization of the reindexation.
c     
      if(ridx_ready.ne.321) then
          write(*,*) 'ERROR, find_levpop_and_temp: Not initialized ',
     %               'the reindexation yet!'
          stop 13
      endif
c
c     If 1-D then we're done quickly, since it must be on a gridpoint
c
#ifdef RADGRID_ONEDIM
      if(icr.eq.2) then
          write(*,*) 'Inconsistency in 1-D find_levpop_and_temp()'
          stop 1
      endif
      if(icr.eq.1) then
          idr   = dr            ! Trick: from dbl to int automatic trunc
          ir0   = ir + idr      ! Note: do not say ir0 = ir + dr        
          nup   = 
          call get_levpop_temp(iline,ir0,1,nup,ndown,temp)
      else
          write(*,*) 'PROBLEM: in 1-D the icr=3 option is not '
          write(*,*) '         yet available for find_src_t()'
          write(*,*) '         If this problem occured while '
          write(*,*) '         rendering spectrum or image in 1-D'
          write(*,*) '         then it is advised to undefine '
          write(*,*) '         LONGCHAR_EXTRA_POINTS'
          stop 1
      endif
      return
#endif
c
c     For 2-D we proceed doing the more complicated interpolation stuff
c
#ifndef INTERPOL_SRC_2
#define INTERPOL_SRC_WIDE
#endif
#ifdef INTERPOL_SRC_0
#define INTERPOL_SRC_PNT_2
#endif
#ifdef INTERPOL_SRC_1
#define INTERPOL_SRC_PNT_2
#endif
c
c     Find the index of the radial shell from the dr 
c
      idr   = dr                ! Trick: from dbl to int automatic trunc
      ir0   = ir + idr          ! Note: do not say ir0 = ir + dr        
c
c     Find the real dr (the dr from R_1 to R, instead of from R_0 to R)
c
      drr   = abs(dr - idr)
c
c     Now find the index of the R_1. This is always one further
c     in the direction of drr, except when we are at R_out or R_in.
c     If so, I choose the index the same as ir0.
c
      if(dr.gt.0.d0) then
#ifdef INTERPOL_SRC_WIDE
          indexr(0) = ir0
          indexr(1) = ir0 + 1
          indexr(2) = ir0 + 2
          rclip(0)  = rsi_x_c(indexr(0),1)
          rclip(1)  = rsi_x_c(indexr(1),1)
          rclip(2)  = rsi_x_c(indexr(2),1)
          if(indexr(1).gt.irsi_frsizex) then
              indexr(1) = irsi_frsizex
          endif
          if(indexr(2).gt.irsi_frsizex) then
              indexr(2) = irsi_frsizex
          endif
#else
          indexr(0) = ir0 - 1
          indexr(1) = ir0 
          indexr(2) = ir0 + 1
          rclip(0)  = rsi_x_c(indexr(0),1)
          rclip(1)  = rsi_x_c(indexr(1),1)
          rclip(2)  = rsi_x_c(indexr(2),1)
          if(indexr(0).lt.1) then
              indexr(0) = 1
          endif
          if(indexr(2).gt.irsi_frsizex) then
              indexr(2) = irsi_frsizex
          endif
#endif          
      else
#ifdef INTERPOL_SRC_WIDE
          indexr(0) = ir0
          indexr(1) = ir0 - 1
          indexr(2) = ir0 - 2
          rclip(0)  = rsi_x_c(indexr(0),1)
          rclip(1)  = rsi_x_c(indexr(1),1)
          rclip(2)  = rsi_x_c(indexr(2),1)
          if(indexr(1).lt.1) then
              indexr(1) = 1
          endif
          if(indexr(2).lt.1) then
              indexr(2) = 1
          endif
#else
          indexr(0) = ir0 + 1
          indexr(1) = ir0 
          indexr(2) = ir0 - 1
          rclip(0)  = rsi_x_c(indexr(0),1)
          rclip(1)  = rsi_x_c(indexr(1),1)
          rclip(2)  = rsi_x_c(indexr(2),1)
          if(indexr(0).gt.irsi_frsizex) then
              indexr(0) = irsi_frsizex
          endif
          if(indexr(2).lt.1) then
              indexr(2) = 1
          endif
#endif          
      endif
c
c     Reconstruct the R itself from the above
c
#ifdef INTERPOL_SRC_WIDE
      r = rclip(0) + drr * ( rclip(1) - rclip(0) )
#else
      r = rclip(1) + drr * ( rclip(2) - rclip(1) )
#endif
c
c     Find the the index of the Theta_0 from the dt
c
      idt   = dt                ! Trick: from dbl to int automatic trunc
      it0   = it + idt          ! Note: do not say it0 = it + dt        
c
c     Find the real dt (the dt from Theta_1 to Theta, instead
c     of from Theta(it) to Theta)
c
      dtr   = abs(dt - idt)
c
c     Now complete the 3-point Theta index stencil.
c
      if(dt.gt.0.d0) then
#ifdef INTERPOL_SRC_WIDE
          indext(0) = it0 
          indext(1) = it0 + 1
          indext(2) = it0 + 2
#else
          indext(0) = it0 - 1
          indext(1) = it0 
          indext(2) = it0 + 1
#endif
      else            
#ifdef INTERPOL_SRC_WIDE
          indext(0) = it0 
          indext(1) = it0 - 1
          indext(2) = it0 - 2
#else
          indext(0) = it0 + 1
          indext(1) = it0 
          indext(2) = it0 - 1
#endif
      endif
c
c     Find the corresponding values of theta on this stencil, except for
c     the thetaclip(2) because that might go wrong, since the grid rsi_x_c
c     only has one ghost-gridpoints on each border, while indext(2) could
c     be 2 gridpoints off the grid under certain circumstances.
c
      thetaclip(0) = rsi_x_c(indext(0),2)
      thetaclip(1) = rsi_x_c(indext(1),2)
c
c     Find these circumstances, and find the thetaclip(2). 
c
      if(indext(2).lt.0) then
          thetaclip(2) = -rsi_x_c(2,2)
      elseif (indext(2).gt.irsi_frsizey+1) then
          thetaclip(2) = 6.28318530718d0 - rsi_x_c(irsi_frsizey-1,2)
      else
          thetaclip(2) = rsi_x_c(indext(2),2)
      endif      
c
c     Reconstruct the Theta itself from the above
c
#ifdef INTERPOL_SRC_WIDE
      theta = thetaclip(0) + dtr * ( thetaclip(1) - thetaclip(0) )
#else
      theta = thetaclip(1) + dtr * ( thetaclip(2) - thetaclip(1) )
#endif
c
c     Remap the theta indices in case of boundaries crossed, or a 
c     mirror line is crossed.
c
      indext(0)   = ridx_it(indext(0))
      indext(1)   = ridx_it(indext(1))
      indext(2)   = ridx_it(indext(2))
c     
c     In case of 0-th order interpolation, do something about 
c     the drr and dtr
c
#ifdef INTERPOL_SRC_0
      if(drr.gt.0.5d0) then 
          drr=1.d0 
      else 
          drr=0.d0 
      endif
      if(dtr.gt.0.5d0) then 
          dtr=1.d0 
      else 
          dtr=0.d0 
      endif
#endif
c
c     Now that we have bracketed the point at which we wish to find
c     the source term, let us make distinction between three types 
c     of interpolation. icr=1 means that the point lies on a R=const
c     grid line, icr=2 means that the point lies on a Theta=const grid
c     line, and icr=0 means that the point really lies somewhere in
c     between both the Theta and R gridlines.
c
#ifdef INTERPOL_SRC_PNT_2
      if(icr.eq.1) then
c
c         The point 'P' lies on an R=const gridline
c
          call get_levpop_temp(iline,indexr(0),indext(0),
     %                         t_nup(0),t_ndown(0),t_temp(0))
          call get_levpop_temp(iline,indexr(0),indext(1),
     %                         t_nup(1),t_ndown(1),t_temp(1))
          nup   = (1.d0-dtr)*t_nup(0)   + dtr*t_nup(1)
          ndown = (1.d0-dtr)*t_ndown(0) + dtr*t_ndown(1)
          temp  = (1.d0-dtr)*t_temp(0)  + dtr*t_temp(1)
c
ccc          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
ccc          src(1) = calc_src_pnt(inu,indexr(0),indext(1),mu,phi)
ccc          find_src_t = (1.d0-dtr)*src(0) + dtr*src(1)
c
      elseif(icr.eq.2) then
c
c         The point 'P' lies on an Theta=const gridline
c
          call get_levpop_temp(iline,indexr(0),indext(0),
     %                         t_nup(0),t_ndown(0),t_temp(0))
          call get_levpop_temp(iline,indexr(1),indext(0),
     %                         t_nup(1),t_ndown(1),t_temp(1))
          nup   = (1.d0-dtr)*t_nup(0)   + dtr*t_nup(1)
          ndown = (1.d0-dtr)*t_ndown(0) + dtr*t_ndown(1)
          temp  = (1.d0-dtr)*t_temp(0)  + dtr*t_temp(1)
c
ccc          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
ccc          src(1) = calc_src_pnt(inu,indexr(1),indext(0),mu,phi)
ccc          find_src_t = (1.d0-drr)*src(0) + drr*src(1)
c
      else
c
c         The point 'P' lies somewhere in between all gridlines
c
          call get_levpop_temp(iline,indexr(0),indext(0),
     %                         t_nup(0),t_ndown(0),t_temp(0))
          call get_levpop_temp(iline,indexr(0),indext(1),
     %                         t_nup(1),t_ndown(1),t_temp(1))
          call get_levpop_temp(iline,indexr(1),indext(0),
     %                         t_nup(2),t_ndown(2),t_temp(2))
          call get_levpop_temp(iline,indexr(1),indext(1),
     %                         t_nup(3),t_ndown(3),t_temp(3))
          nup   = (1.d0-drr) * ((1.d0-dtr)*t_nup(0) + dtr*t_nup(1)) 
     %               +   drr * ((1.d0-dtr)*t_nup(2) + dtr*t_nup(3)) 
          ndown = (1.d0-drr) * ((1.d0-dtr)*t_ndown(0) + dtr*t_ndown(1)) 
     %               +   drr * ((1.d0-dtr)*t_ndown(2) + dtr*t_ndown(3)) 
          temp  = (1.d0-drr) * ((1.d0-dtr)*t_temp(0) + dtr*t_temp(1)) 
     %               +   drr * ((1.d0-dtr)*t_temp(2) + dtr*t_temp(3)) 
c
ccc          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
ccc          src(1) = calc_src_pnt(inu,indexr(0),indext(1),mu,phi)
ccc          src(2) = calc_src_pnt(inu,indexr(1),indext(0),mu,phi)
ccc          src(3) = calc_src_pnt(inu,indexr(1),indext(1),mu,phi)
ccc          find_src_t = (1.d0-drr) * ((1.d0-dtr)*src(0) + dtr*src(1)) 
ccc     %                +       drr * ((1.d0-dtr)*src(2) + dtr*src(3)) 
      endif
#else
      if(icr.eq.1) then
c
c         The point 'P' lies on an R=const gridline
c
ccc          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
ccc          src(1) = calc_src_pnt(inu,indexr(0),indext(1),mu,phi)
ccc          src(2) = calc_src_pnt(inu,indexr(0),indext(2),mu,phi)
ccc          call polint3(thetaclip,src,theta,find_src_t)
          stop 438
c
      elseif(icr.eq.2) then
c
c         The point 'P' lies on an Theta=const gridline
c
ccc          src(0) = calc_src_pnt(inu,indexr(0),indext(0),mu,phi)
ccc          src(1) = calc_src_pnt(inu,indexr(1),indext(0),mu,phi)
ccc          src(2) = calc_src_pnt(inu,indexr(2),indext(0),mu,phi)
ccc          call polint3(rclip,src,r,find_src_t)
          stop 438
c
      else
c
c         The point 'P' lies somewhere in between all gridlines
c
ccc          write(*,*) 'ERROR: interpol.F/find_src_t()'
ccc          write(*,*) '       icr=0: This option is not yet ',
ccc     %      'available for 2nd order interpolation'
ccc          stop 14
          stop 438
      endif
c
c     Make sure that the answer is positive definite
c
ccc      find_src_t = 0.5d0 * ( find_src_t + abs(find_src_t) )
#endif
c
      return
      end
#endif /* LINE_SRC_INTERPOL_LEVPOP */
#endif /* INCLUDE_LINES */




c#endif (was this the ifndef NOTRANS? Removed on 31-07-00)
