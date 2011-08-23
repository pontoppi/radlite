c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c                    ( MODULE: SHORT CHARACTERISTICS )
c
c                            Leiden, June 1998
c                             C.P. Dullemond   
c
c     =============================================================


#define SHORTCHAR_F_VERSION "Trans_31Jul00"
#include "main.h"

#ifndef NO_SHORT_CHARS

#ifndef SHORTCHAR_EPS
#define SHORTCHAR_EPS 1.d-6
#endif


c     --------------------------------------------------------------
c                  SET UP ALL SHORT CHARACTERISTICS
c
c     This routine will first find out which types of SC's are 
c     necessary, for the radial grid used. This is done by calling
c     the routine inspect_radial_grid(). Then it will create each
c     of these sets of SC's by calling (for each of them) the 
c     routine create_shortchar(drrp,drrn,idrr,itcross). This routine
c     assumes that the radial, Theta and angular grids have been
c     set already. Also make sure that the iradbnd_*** arrays have
c     already been set properly.
c     --------------------------------------------------------------
      subroutine setup_all_shortchar_sets()
      implicit none
c
#include "common_grid.h"
#include "common_boundary.h"
#include "common_shortchar.h"
c
      doubleprecision drr
      integer i
c
c     First analyze the grid to find the grid spacings.
c
      call inspect_radial_grid()
c
c     Then reset the "isc_longesc_current" pointer to 1
c     This is for the smart (compact) storage of the ESCs
c
#     ifdef SHORCHAR_COMPACT_ARRAYS
      isc_longesc_current = 0
      isc_shortsc_current = 0
      isc_error_longesc   = 0
ccc      write(*,*) 'MAX nr of ESC ',FRSIZE_LONGESC
#     endif
c
c     Reset the error
c
      isc_error  = 0
      isc_maxlen = 0
c
c     If necessary, then first create the SC of inner vacuum
c
      if((iradbnd_in_itype.eq.2).or.(iradbnd_in_itype.eq.0)) then
          write(*,*) 'Creating the short characteristics for ',
     %               'the inner vacuum'
          if(radbnd_rstar.eq.0.d0) then
cc              write(*,*) 'ERROR: With central star as inner BC'
cc              write(*,*) '       the stellar radius must be'
cc              write(*,*) '       non-zero'
cc              stop 13
              drr = 1.d3  ! i.e. very small rstar
          else
              drr = rsi_x_c(1,1)/radbnd_rstar-1.d0
          endif
          call create_shortchar(drr,drr,0,0)
      endif
c
c     Now create the short characteristics for the logarithmic
c     grid with the possible 1/2^n refinement steps
c
      write(*,*) 'Creating the short characteristics for ',
     %     'the logaritmic grid'
      drr = sc_rratiobig-1.d0
      call create_shortchar(drr,drr,1,1)
c
c     Refinements, if present
c
      if(sc_ndrr.gt.1) then
          do i=2,sc_ndrr
              drr=sc_rratiobig**(0.5d0**sc_nstep(i))-1.d0
              write(*,*) '    Refinement ',i-1,':  n=',sc_nstep(i),
     %                   ' dr/r=',drr
              call create_shortchar(drr,drr,i,1)
          enddo
      endif
c
c     If something went wrong with SC setup, then stop here with
c     a message
c
      if(isc_error.ne.0) then
          write(*,*) '========================================'
          write(*,*) '  Extended Short Characteristic has '
          write(*,*) '      too many segments for this    '
          write(*,*) '            compilation!            '
          write(*,*) '  Aborting the construction of the  '
          write(*,*) '   Extended Short Characteristics.  '
          write(*,*) ' Recompile with FRSIZE_CHAR ',isc_maxlen
          write(*,*) '========================================'
          stop 13
      endif
c
c     Message
c
#     ifdef SHORCHAR_COMPACT_ARRAYS
      write(*,*) '    Nr of short SC    : ',isc_shortsc_current
      write(*,*) '    Nr of extended SC : ',isc_longesc_current
      write(*,*) '    Frac extended SC  : ',(isc_longesc_current+1.0)/
     %                                 isc_shortsc_current
#     endif
c
c     If something went wrong during the SC storage, then stop
c     here with a message
c
#     ifdef SHORCHAR_COMPACT_ARRAYS
      if(isc_error_longesc.ne.0) then
          write(*,*) 'ERROR during ESC buildup: LONGESC_FRACTION'
          write(*,*) '   is too small: ',LONGESC_FRACTION
          write(*,*) '   Put it in configure.h to at'
          write(*,*) '   least ',(isc_longesc_current+1.0)/
     %                                 isc_shortsc_current
          stop 13
      endif
#     endif
c
      end



c     --------------------------------------------------------------
c             INSPECT RADIAL GRID FOR SHORT CHAR SPACINGS
c
c     The short characteristic trajectories will be precomputed in
c     RADICAL. In order to save memory we require the radial grid
c     to be logarithmic in general, and if necessary to have grid 
c     refinements which are discrete factors of 2^n finer. This
c     way we need to compute only a few (minimal 2, but perhaps 
c     at most 6 or 10) sets of SC's, even when the radial grid is
c     40 or 60 long. This saves lots of memory, and is crucial for
c     making RADICAL run on small desktop/laptop computers.
c
c     The present routine inspects the radial grid to find out:
c      1) Which is the largest dr/r spacing?
c      2) Are all smaller dr/r spacings exactly factors 1/2^n times
c         smaller?
c      3) If so, then which spacings are there (i.e. which values
c         of 'n')?
c     --------------------------------------------------------------
      subroutine inspect_radial_grid()
      implicit none
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      doubleprecision rratio(FRSIZE_X),safety,rratn,rrat
      integer nstep(0:50)
      integer ir,n,ndrr,i,nold
c
c     Message
c
      write(*,*) 'Inspecting the radial grid'
c
      safety = 0.01d0
c
c     Get the spacings and sort
c
      do ir=1,irsi_frsizex-1
         rratio(ir) = rsi_x_c(ir+1,1) / rsi_x_c(ir,1)
         if(rratio(ir).lt.1.d0) then
             write(*,*) 'While sorting radial grid: negative spacing!'
             stop 13
         endif
      enddo
      call sort(irsi_frsizex-1,rratio)
c
c     Check whether the spacings are multiples of 1/2 times the 
c     largest spacing
c
      n            = 0
      ndrr         = 1
      nstep(1)     = 0
      sc_rratiobig = rratio(irsi_frsizex-1)
      do ir=irsi_frsizex-2,1,-1
          if(rratio(ir).lt.1.d0+1.d-10) then
              write(*,*) 'ERROR: Radial grid spacing is < 1e-10'
              stop 13
          endif
          nold=n
          do while(sc_rratiobig**(0.5d0**n)-rratio(ir).gt.
     %               (rratio(ir)-1.d0)*safety) 
              n=n+1
          enddo 
          if(sc_rratiobig**(0.5d0**n)-rratio(ir).lt.
     %               (1.d0-rratio(ir))*safety) then
              write(*,*) 'FATAL PROBLEM:'
              write(*,*) 'For the Short Characteristics algorithm'
              write(*,*) 'the radial grid should be log-spaced, i.e.'
              write(*,*) 'r(i+1)/r(i) must be a constant. OR it'
              write(*,*) 'must be a power 1/2^n (where n is an'
              write(*,*) 'arbitrary integer >=0) of the largest '
              write(*,*) 'r-ratio. This latter generalization of the '
              write(*,*) 'log-spacing allows you to refine the grid'
              write(*,*) 'where necessary, still keeping the general'
              write(*,*) 'logarithmic grid philosophy intact.'
              write(*,*) 'NOTE: If you use RADICAL only for ray-'
              write(*,*) 'tracing, then switch on the ONLY_RAY_TRACING'
              write(*,*) 'macro, and then any R-grid is allowed!'
              stop 13
          endif
          if(n.ne.nold) then
              ndrr = ndrr + 1
              if(ndrr.gt.50) then
                  write(*,*) 'WEIRD PROBLEM: More than 50 ',
     %                       'refinement steps!'
                  write(*,*) 'Something must be wrong.......'
                  write(*,*) 'Warn author.'
                  stop 13
              endif
              nstep(ndrr) = n
          endif
      enddo
c
c     Create the set of refinement n-values
c
      if(ndrr.gt.FRSIZE_DRR) then
          write(*,*) 'PROBLEM: The number of grid-refinement'
          write(*,*) 'steps in R exceeds FRSIZE_DRR, which is'
          write(*,*) 'the compiled limit of the number of sets'
          write(*,*) 'of short characteristics. Please increase '
          write(*,*) 'FRSIZE_DRR to ',ndrr
          stop 13
      endif
      sc_nstep(0) = 0
      do i=1,ndrr
          sc_nstep(i)=nstep(i)
      enddo
      sc_ndrr = ndrr
c
c     Message if FRSIZE_DRR larger than necessary
c
      if(sc_ndrr.lt.FRSIZE_DRR) then
          write(*,*) '    Note: FRSIZE_DRR = ',FRSIZE_DRR,
     %      '  whereas ',ndrr,'  would suffice'
      endif
c
c     Creat the ir_iscset_p and ir_iscset_n
c
      do ir=2,irsi_frsizex
          rrat = rsi_x_c(ir,1) / rsi_x_c(ir-1,1)
          do i=1,ndrr
              rratn = sc_rratiobig**(0.5d0**nstep(i))
              if(abs(rratn-rrat).le.safety*(rrat-1.d0)) then
                  ir_iscset_p(ir-1) = i
                  ir_iscset_n(ir)   = i
                  goto 50
              endif
          enddo
          write(*,*) 'INTERNAL ERROR: shortchar.F'
          write(*,*) 'Warn the author.'
          stop 13
 50       continue
      enddo
      ir_iscset_n(1) = 0
      ir_iscset_p(irsi_frsizex) = 10000  ! Safety nr: must crash if used
c
c     Message
c
      write(*,*) '    Largest DR/R = ',sc_rratiobig-1.d0
c
cccc#################################
ccc      do ir=1,irsi_frsizex
ccc          write(*,*) ir,ir_iscset_n(ir),ir_iscset_p(ir),
ccc     %            sc_rratiobig**(0.5d0**nstep(ir_iscset_n(ir)))
ccc      enddo
cccc#################################
c
      end




c     --------------------------------------------------------------
c               CREATE THE ARRAY OF SHORT CHARACTERISTICS
c
c     This subroutine finds the half-long characteristic starting
c     at some point and ending at the next radial shell.
c
c     REMARK: The angular set is always taken iangset=1. 
c     --------------------------------------------------------------
      subroutine create_shortchar(drrp,drrn,idrr,itcross)
      implicit none
      integer idrr,itcross
      doubleprecision drrp,drrn
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      integer it,imu,iphi,ismax,iangset,itmax
c
c     Take the only one angular set that is allowed right now.
c
      iangset = 1
c
c     First check
c
#ifndef COORD_CARTESIAN
      if(drrp.le.0.d0) then
          write(*,*) 'In Create_shortchar: drrp<=0'
          stop 13
      endif
      if(drrn.le.0.d0) then
          write(*,*) 'In Create_shortchar: drrn<=0'
          stop 13
      endif
#endif
      if(idrr.gt.FRSIZE_DRR) then
          write(*,*) 'idrr out of range... Redefine FRSIZE_DRR!'
          stop 13
      endif
      if((nrmu(iangset).eq.0).or.(nrphi(iangset).eq.0)) then
          write(*,*) 'For short characteristics ',idrr,
     %          ' the angular grid is not set...'
          stop 13
      endif
c
c     Now the mirror things...
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#else
      itmax = irsi_frsizey
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
c     Set and reset some constants
c
      isc_iangset(idrr) = iangset
      sc_drrp(idrr)     = drrp
      sc_drrn(idrr)     = drrn
      ismax             = 1
c
c     Loop over Theta
c     
cccccc (REDUNDANTLY LARGE: 14-05-01)     do it=1,irsi_frsizey 
c
      do it=1,itmax
c     
c         Loop over mu. In case anggrid_mu_zero.eq.0, then in 
c         principle the imu=0 should be skipped, but it does no
c         harm to include it here anyway.
c     
          do imu=-nrmu(iangset),nrmu(iangset)
c     
c             Loop over phi
c     
              do iphi=1,nrphi(iangset)
                  call find_shortchar(drrp,drrn,iangset,it,
     %                                imu,iphi,itcross)
                  call copy_shortchar(idrr,it,imu,iphi)
                  if(iosc_nel.gt.ismax) then
                      ismax = iosc_nel
                  endif
              enddo
          enddo
      enddo
      isc_ismax(idrr) = ismax
      if(isc_maxlen.lt.ismax) then
          isc_maxlen=ismax
      endif
      end




c     --------------------------------------------------------------
c                  FIND THE EXTENDED SHORT CHARACTERISTIC 
c
c     This subroutine finds the extended short characteristic.
c
c     The shape of the characteristics (i.e. the crossings with
c     the coordinate grid lines) depends on the gridding in Theta
c     and on the gridding in R. But for a logarithmic grid in R,
c     this shape is self-similar for R_i -> R_i+1. So instead of
c     storing the shape-information for the local characteristic
c     for each Theta AND R, we only store it for each Theta. 
c
c     Note: The ds array is of course valid for R=1. So multiply 
c           with the real R before use: real_ds = sd * r
c
c     Note: If COORD_CARTESIAN is set, then the short characteristics
c           are calculated not for spherical but for Cartesian 
c           coordinates. 
c
c     Note: The COORD_CARTESIAN option is now only working for the
c           1-D case. So it requires the option RADGRID_ONEDIM.
c
c     ARGUMENTS:
c       drrp   = dr/r from the present r_i to the next r_i+1. The 'p'
c                stands for positive step in r.
c       drrn   = dr/r from the present r_i to the previous r_i-1. The
c                'n' stands for negative step in r.
c       iangset= The index of the angular set.
c       it     = i_Theta, the index of the theta coordinate.
c       imu    = i_mu, the index for the local momentum-space mu angle
c       iphi   = i_phi, the index for the local momentum-space phi angle
c       itcross= If set to 0, then only calculate the crossings with
c                the r_i = constant gridlines. If set to 1, then do the
c                normal stuff. So for the normal 2-D spherical ESC, set
c                this to 1.
c     --------------------------------------------------------------
      subroutine find_shortchar(drrp,drrn,iangset,it,imu,iphi,itcross)
      implicit none
      doubleprecision drrp,drrn
      integer iangset,it,imu,iphi,itcross
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      doubleprecision s0,s1,st,scur,sother,sother0,ds,sxtr,sol1,sol2 
      doubleprecision r1,rt,r,theta,mu,phi,b,costinf,z0
      doubleprecision dum1,dum2,det,costh,costh1,th1
      doubleprecision eps0,eps1,eps2,epsmu,epsdet,pi,pitwice,pihalf
      integer ir0,idr,ixtr,itdir0,itdir,itcur,is,itnext,i,phiplus,flg
c
      parameter(eps0=SHORTCHAR_EPS)
      parameter(eps1=0.4d0*SHORTCHAR_EPS)
      parameter(eps2=1.0d1*SHORTCHAR_EPS)
      parameter(epsmu=1.d-4)
      parameter(epsdet=1.d-10)
      parameter(pi=PICONST)
      parameter(pitwice=2.d0*PICONST)
      parameter(pihalf=0.5d0*PICONST)
c
#ifdef COORD_SPHERICAL
c
c     Do check for reliability of this algorithm:
c
      if(min(drrp,drrn).lt.1.d4*eps0) then
          write(*,*) 'PROBLEM:'
          write(*,*) 'Trying to set up short characteristics for'
          write(*,*) 'refined grid: dr/r=',min(drrp,drrn),'.'
          write(*,*) 'The current safety tuning numbers for the'
          write(*,*) 'short characteristics are not safe (=small)'
          write(*,*) 'enough for this fine a mesh. Try reducing'
          write(*,*) 'the SHORTCHAR_EPS. But not too small since'
          write(*,*) 'then one might run into machine-precision'
          write(*,*) 'errors.'
          stop 13
      endif
c
c     Find the extended short characteristics for the spherical
c     coordinate case, both in 1-D and in 2-D. 
c
      isc_ready = 0
c
c     First find the R,Theta,mu and phi at the present point. We take
c     for convenience R=1, because for the parameterization of the
c     characteristics the R is irrelevant. Only the Delta R / R is
c     important.
c
      r       = 1.d0
      theta   = rsi_x_c(it,2)
      mu      = rmu(imu,iangset)
      phi     = rphi(iphi,iangset)
      if(mu*mu.gt.1.d0) then
          write(*,*) 'ERROR: mu^2 > 1'
          stop 131
      endif
      if((phi.gt.pitwice).or.(phi.lt.0.d0)) then
          write(*,*) 'ERROR: phi not right'
          stop 132
      endif
      if((theta.gt.pi).or.(theta.lt.0.d0)) then
          write(*,*) 'ERROR: theta out of range'
          write(*,*) ' Theta = ',theta
          stop 133
      endif
      if((phi.lt.pihalf).or.(phi.gt.pi+pihalf)) then
          phiplus = 1
      else
          phiplus = -1
      endif
      if(abs(1.d0-mu*mu).lt.2*MUTHRES) then
          if(mu.lt.0) then
              mu=-1.d0
          else
              mu=1.d0
          endif
      endif
c      if(iangset.eq.0) then
c          write(*,*) 'ERROR: iangset=0 in find_shortchar()'
c          stop 13
c      endif
c
c     Find the parameters (b,thetainf,z0) of the characteristic curve
c
      b       = r * sqrt(abs(1.d0-mu*mu))
      costinf = mu * cos(theta) + sqrt(1.d0-mu*mu) * 
     %                  sin(theta) * sin(phi) 
      z0      = r * ( (1.d0-mu*mu) * cos(theta) - 
     %             mu * sqrt(1.d0-mu*mu) * sin(theta) * sin(phi) ) 
      s0      = r * mu
c
c     Find the ds corresponding to the extended characteristic
c     in upwind direction from R till the upwind R=const
c     grid line. If the upwind direction is off the grid, then
c     simply return without result.
c
      if(mu.gt.MUTHRES) then
          r1  = r / ( 1.d0 + drrn )
          ir0 = -1
          idr = -1
      else
          r1  = r * ( 1.d0 + drrp )
          ir0 = 0
          idr = 1
      endif
c
      dum2 = r*r - b*b + eps0*r*r
      if(dum2.lt.0.d0) then
          write(*,*) 'ERROR: r<b'
          stop 140
      endif
c
c     If the ESC begins/ends at the same radius, the idr=0, and
c     the ds has special form. Otherwise use the usual form.
c
      if(r1.lt.b-eps1*b) then
          dum1 = 0.d0
          ds   = - 2.d0 * abs(sqrt(dum2))
          idr  = 0
      else
          dum1 = r1*r1 - b*b + eps0*r1*r1
          ds   = - abs(sqrt(dum1)-sqrt(dum2))
      endif
c
c     Find the s1 where the ray hits the upwind R=const grid line. 
c     For the pure ESC method, this point acts as the endpoint of
c     our local characteristic. Or to say more accurately: 
c     the starting point (since you start the integration 
c     at point s1 and end at s0). For the MESC variant this point
c     merely acts as a potential end-point, if either the SC hits
c     no Theta gridline, or causality requires the SC to be extended.
c     
      s1     = s0 + ds
c
c     BUGFIX 11-04-99: In case b=0 and drrn>>>>1 (e.g. for the 
c     inner vacuum with central star), then s0 and ds almost cancel
c     eachother. In this case, set s1 directly!
c
      if((b.eq.0.d0).and.(mu.eq.1.d0).and.(drrn.gt.1.d3)) then
          s1 = r1
      endif
c
c     If b=0, or if the itcross=0, then direct hit the R=const grid line
c
      if((b.eq.0.d0).or.(itcross.eq.0)) then
          is   = 1
          scur = s0
          goto 100
      endif
c
c     If 1-D short characteristics transfer, then we also skip
c     the search for Theta crossings.
c
#ifdef RADGRID_ONEDIM
      is   = 1
      scur = s0
      goto 100
#endif
c
c     Find the s where the maximum of theta takes place
c
      sxtr   = b * b * costinf / z0
c
c     Does this point lie upwind or downwind from point 0?
c     Meaning: ixtr=1 means extremum lies upwind from s0.
c     So if ixtr=1 we should anticipate meeting this extremum
c     during our tracing of the characteristic.
c
      if(s0.gt.sxtr) then
          ixtr = 1
      else
          ixtr = -1
      endif
c
c     Find the theta-direction of the characteristic at the
c     point 0 (i.e. towards positive or negative theta?)
c     
      if((phi.gt.0.d0).and.(phi.lt.pi)) then
          itdir0 = -1
      else
          itdir0 =  1
      endif
c
c     Get the current s and it's cousin
c
      scur   = s0
      itcur  = it
      costh  = cos(rsi_x_c(itcur,2))
      det    = z0*z0 + b*b*(costinf*costinf-
     %                costh*costh)
      if(det.lt.0.d0) then
          write(*,*) 'Weird error...'
          write(*,*) '  z0      = ',z0
          write(*,*) '  b       = ',b
          write(*,*) '  costinf = ',costinf
          write(*,*) '  costinh = ',costh
          write(*,*) '  det     = ',det  
          write(*,*) 'It could be related to having only one'
          write(*,*) 'phi-point in 2-D. Or having phi-points'
          write(*,*) 'being at phi=0.00d0 or so.'
          stop 14
      endif
      dum1  = costh * sqrt( det )
      dum2  = costh*costh - costinf*costinf
      if(abs(1.d0-mu*mu).lt.2*MUTHRES) then
          write(*,*) 'HE?'
          stop 13
      endif
      if(abs(dum2).gt.eps2) then
          sol1  = ( z0*costinf + dum1 ) / dum2
          sol2  = ( z0*costinf - dum1 ) / dum2
          if((abs(s0-sol1).gt.eps2*(abs(sol1)+abs(sol2))).and.
     %        (abs(s0-sol2).gt.eps2*(abs(sol1)+abs(sol2)))) then
              write(*,*) 'Hey! One of them must fit?'
              stop 16
          endif
          if(abs(s0-sol1).gt.eps2*(abs(sol1)+abs(sol2))) then
              sother = sol1
          else
              sother = sol2
          endif
      else
          sother = 1e30
      endif
      sother0 = sother
c
c     Do a loop along the ray in upwind direction. Each step is
c     a search for the next crossing with a theta=const gridline.
c     If the characteristic (followed upstream) crosses a R=const
c     gridline before crossing a theta=const gridline, then this
c     loop should be passed only once (a failed search for the
c     crossing with theta=const). 
c
      do is=1,irsi_frsizey
c
c         Safety check. The is parameter may not be too long. There 
c         is already a check below, but somehow there is at the date
c         12 october 1998 still some rare bug (bus error), so let's
c         put in safety checks...
c
          if(is.gt.FRSIZE_CHAR_BIG) then
              write(*,*) 'ERROR: ESC too long!!!'
              stop 13
          endif
c
c         Find the local theta direction. Are the photons moving
c         towards larger or smaller Theta? Meaning: itdir=1 means
c         towards larger Theta (i.e. downwards). 
c
          if(ixtr*(scur-sxtr).gt.0.d0) then
              itdir =  itdir0
          else
              itdir = -itdir0
          endif              
c
c         The index of the previous theta (i.e. next
c         theta in upwind direction)
c
          itnext = itcur - itdir
c
c         Try to find the crossing with the next theta coordinate
c
          if((itnext.lt.1).or.(itnext.gt.irsi_frsizey)) then
              det = -1.d0 
          else
              costh = cos(rsi_x_c(itnext,2))
              det   = z0*z0 + b*b*(costinf*costinf-costh*costh)
          endif
          if(det.le.epsdet) then
c
c             There exist no solutions at all for this theta.
c             This means that the turnover point in theta must
c             lie within this step, OR we simply reached an
c             R=const gridline. In order to find out, we check
c             the other s-solution from the previous crossing
c             with a theta=const grid line. If it lies in the
c             domain, then it means that indeed the characteristic
c             has turned back and crossed the same theta=const
c             line before the R=const line is hit.
c
c             For the MESC variant: The det.lt.0 means that the
c             either the SC should be extended, or we reached
c             an R=const gridline anyway. In any case: we need
c             not abort the procedure of finding the ESC: we
c             will continue the loop until we reach an R=const
c             gridline.
c
              if((sother.gt.scur).or.(sother.lt.s1)) then
c
c                 The other s from the previous theta step was 
c                 outside of the domain, so this means that
c                 we must have arrived at an R=const gridline. 
c
                  goto 100
              else
c
c                 The other s is inside the domain, so logically
c                 it should indicate that we are dealing with a
c                 theta extremum here. So store the information
c                 about this line element into the array, and
c                 continue the loop.
c
                  iosc_icross(is) = 2
                  st = sother
                  rt = sqrt( b*b + st*st )
                  itnext = itcur
                  osc_ds(is) = scur - st
                  if(r1.eq.r) then
                      write(*,*) '?'
                      stop 150
                  endif
                  if(r1.gt.r) then
                      osc_dr(is) = ( rt - r ) / ( r1 - r )
                  else
                      osc_dr(is) = ( rt - r ) / ( r - r1 )
                  endif
                  osc_dtheta(is)  = itnext - it
                  call calc_phi(is,iangset,iphi,b,costinf,z0,st,phiplus)
                  call calc_dmu(is,iangset,st/rt)
              endif
          else
c     
c             There are solutions for this theta, but the question
c             is: are they within the domain in s?
c
c             For MESC: The fact that the determinant.gt.0.d0 means that
c             the short characteristic should not be extended. After 
c             either a Theta=const gridline or R=const gridline we must
c             stop. Since we stop anyway at an R=const gridline, we need
c             only specially care for the Theta=const gridline.
c
c             Find the two solutions for s
c
              dum1  = costh * sqrt( det )
              dum2  = costh*costh - costinf*costinf
              if(dum2.eq.0.d0) then
                  write(*,*) 'ERROR: sol1,sol2 undefined'
              endif
              sol1  = ( z0*costinf + dum1 ) / dum2
              sol2  = ( z0*costinf - dum1 ) / dum2
c
c             Now check if they lie within the domain [s1,scur>
c             (remember we are tracing upwind towards lower s).
c
              scur = scur - eps2
              if(((sol1.gt.scur).or.(sol1.lt.s1)).and.
     %             ((sol2.gt.scur).or.(sol2.lt.s1))) then
c     
c                 Both sol1 and sol2 are outside of the domain. 
c                 This means that we must have arrived at an
c                 R=const grid line.
c
                  goto 100
              else
c
c                 At least one of them lies within the domain.
c                 Find which one is the closest from the current
c                 position
c              
                  if((sol1.gt.scur).or.(sol1.lt.s1)) then
c
c                     If not sol1 then it must be sol2 only
c
                      st     = sol2
                      sother = sol1
                      if((z0+st*costinf)*costh.lt.0.d0) goto 100
                  else if((sol2.gt.scur).or.(sol2.lt.s1)) then
c
c                     It must be sol1 only
c
                      st     = sol1
                      sother = sol2
                      if((z0+st*costinf)*costh.lt.0.d0) goto 100
                  else
c
c                     Both of them... So which one is the 
c                     closest to scur?
c
                      if((scur-sol1).lt.(scur-sol2)) then
                          st     = sol1
                          sother = sol2
                      else
                          st     = sol2
                          sother = sol1
                      endif
                      if((z0+st*costinf)*costh.lt.0.d0) then
                          dum1   = st
                          st     = sother
                          sother = dum1
                      endif
                  endif
c
c                 So now that we arrived here, we know that we have
c                 hit a Theta=const gridline and we also know the
c                 corresponding s. So we can fill the array again.
c
                  iosc_icross(is) = 2
                  rt = sqrt( b*b + st*st )
                  osc_ds(is) = scur - st
                  if(r1.eq.r) then
                      write(*,*) '?'
                      stop 150
                  endif
                  if(r1.gt.r) then
                      osc_dr(is) = ( rt - r ) / ( r1 - r )
                  else
                      osc_dr(is) = ( rt - r ) / ( r - r1 )
                  endif
                  osc_dtheta(is)  = itnext - it
                  call calc_phi(is,iangset,iphi,b,costinf,z0,st,phiplus)
                  call calc_dmu(is,iangset,st/rt)
c
c                 For MESC: once we arrived here, and is.eq.1 (i.e. we are
c                 in the first loop) we know for certain that the SC should 
c                 not be extended, and we must finish the tracing of the 
c                 short characteristic here.
c
#ifdef ESC_MINIMAL_EXTENSION
                  if(is.eq.1) then
                      iosc_nel = 1
                      goto 105
                  endif
#endif                  
              endif
          endif
c
c         New --> Old
c
          scur  = st
          itcur = itnext
c
c         Test if we are out of range of is
c
          if(is.gt.FRSIZE_CHAR_BIG) then 
              write(*,*) '===================================='
              write(*,*) '  Extended Short Characteristic has '
              write(*,*) '      too many segments, even many  '
              write(*,*) '   more than FRSIZE_CHAR_BIG!!!     '
              write(*,*) '  You must make FRSIZE_CHAR much    '
              write(*,*) '    larger!!!                       '
              write(*,*) '===================================='
              stop 13
          endif
          if(is.gt.FRSIZE_CHAR) then
c              write(*,*) '===================================='
c              write(*,*) '  Extended Short Characteristic has '
c              write(*,*) '      too many segments for this    '
c              write(*,*) '            compilation!            '
c              write(*,*) '  Aborting the construction of the  '
c              write(*,*) '   Extended Short Characteristics.  '
c              write(*,*) ' Recompile with larger FRSIZE_CHAR  '
c              write(*,*) '===================================='
cc              write(*,*) '----> For the moment, resuming and see '
cc              write(*,*) '      what we can do without the ESC. '
c              stop 13
cc              goto 159
              isc_error = 1
          endif
      enddo
c
c     Ehh, howcome we arrive here? 
c
      write(*,*) 'Somehow missed out somewhere'
      stop 20
  100 continue
c
c     Store how many line elements there are in the local 
c     characteristic. The amount of theta=const grid lines
c     crossed is this number minus 1. 
c
      iosc_nel  = is
c
c     Safety check again...
c
      if(iosc_nel.gt.FRSIZE_CHAR_BIG) then
          write(*,*) 'ERROR: iosc_nel>FRSIZE_CHAR_BIG'
          stop 13
      endif
      if(is.gt.FRSIZE_CHAR) then
          isc_error = 1
      endif
c
c     In 1-D add a point at closest approach to the star, since
c     otherwise flux will be lost in optically thick media.
c     [Added 22-06-01]
c
#ifdef RADGRID_ONEDIM
#ifdef SCONEDIM_ADDPOINT
      if((idr.eq.0).and.(itcross.eq.1)) then
c
c         This closest approach point is in the middle of r-grid points
c
          iosc_icross(is) = 0
          st              = 0.d0
          rt              = b
          osc_ds(is)      = scur - st
          if(osc_ds(is).lt.0.d0) then 
              write(*,*) 'ERROR adding point: ds<0'
              write(*,*) scur,s0,s1,mu,flg
              stop 183
          endif
          osc_dtheta(is) = 0.d0
          osc_dr(is)     = ( rt - r ) / ( r - r1 )
          call calc_phi(is,iangset,iphi,b,costinf,z0,st,phiplus)
          call calc_dmu(is,iangset,st/rt)
          scur  = st
c
c#######################################
          write(*,*) imu,rmu(imu,iangset),b,osc_dr(is)
c#######################################
c
c         Increase counter
c      
          is       = is + 1
          iosc_nel = is
c
c         Safety check again...
c
          if(iosc_nel.gt.FRSIZE_CHAR_BIG) then
              write(*,*) 'ERROR: iosc_nel>FRSIZE_CHAR_BIG'
              stop 13
          endif
          if(is.gt.FRSIZE_CHAR) then
              isc_error = 1
          endif
      endif
#endif
#endif
c     
c     We arrived at the upwind R=const grid line. This is the end
c     of the local characteristic, or as seen from the affine
c     parameter: the beginning. 
c
      iosc_icross(is)  = 1
      st         = s1
      rt         = sqrt( b*b + st*st )
      osc_ds(is) = scur - s1
      if(osc_ds(is).lt.0.d0) then 
          write(*,*) 'ERROR: ds<0'
          write(*,*) scur,s0,s1,mu,flg
          stop 14
      endif
c
#ifndef RADGRID_ONEDIM
      if(b.gt.0.d0) then
          costh1     = ( z0 + costinf*s1 ) / sqrt( b*b + s1*s1 )
          if(abs(costh1).gt.1.d0) then
              write(*,*) 'ERROR: costh1 out of range'
              stop 170
          endif
          th1        = acos(costh1)
      else
          costh1     = cos(theta)
          th1        = theta
      endif
c
c     Find the index of theta
c
      do i=1,irsi_frsizey
          if(rsi_x_c(i,2).gt.th1) goto 120
      enddo
  120 continue
      i = i - 1
      if(i.lt.it) then
          osc_dtheta(is)  = i - it + 1 +
     %        (th1-rsi_x_c(i+1,2)) / 
     %        (rsi_x_c(i+1,2)-rsi_x_c(i,2)) 
      else
          osc_dtheta(is)  = i - it +
     %        (th1-rsi_x_c(i,2)) / 
     %        (rsi_x_c(i+1,2)-rsi_x_c(i,2)) 
      endif
c
#else
      osc_dtheta(is) = 0.d0
#endif
c
c     The dr is simple. Since we are at the end of the
c     extended short characteristic, this must be idr.
c     which can have values -1,0,+1. 
c
      osc_dr(is) = idr
      call calc_phi(is,iangset,iphi,b,costinf,z0,st,phiplus)
      call calc_dmu(is,iangset,st/rt)
c
  105 continue
c     ----------------------------------------------------
c     Now find the 1 single point downwind (for the second
c     order integration).
c     ----------------------------------------------------
c
c     First find the next R=const grid line downwind
c
      if(mu.lt.-MUTHRES) then
          r1 = r / ( 1 + drrn ) 
          idr = -1
      else
          r1 = r * ( 1 + drrp ) 
          idr = 1
      endif
c
c     Now find the ds to the next R=const grid line (downwind)
c     
      if(r1.lt.b-eps1*b) then
          ds  = 2.d0 * abs(sqrt(r*r-b*b+eps0*r*r))
          idr = 0
      else
          ds = abs(sqrt(r1*r1-b*b+eps0*r1*r1)-
     %                sqrt(r*r-b*b+eps0*r*r))
      endif
c
      s1 = s0 + ds
c
c     If 1-D transfer, then don't worry about theta
c
#ifdef RADGRID_ONEDIM
      is              = 1
      iosc_icross(0)  = 1 
      osc_ds(0)       = s1 - s0
      st              = s1
      rt              = sqrt( b*b + st*st )
      osc_dtheta(0)   = 0.d0
      osc_dr(0)       = idr
      call calc_phi(0,iangset,iphi,b,costinf,z0,st,phiplus)
      call calc_dmu(0,iangset,st/rt)
c
c     If necessary: replace the above with closest approach point
c     [Added 22-06-01]
c
#     ifdef SCONEDIM_ADDPOINT
      if((idr.eq.0).and.(itcross.eq.1)) then
          is              = 1
          iosc_icross(0)  = 0 
          osc_ds(0)       = 0.d0 - s0
c#################################
          if(osc_ds(0).lt.0.d0) stop 356345
c#################################
          st              = 0.d0
          rt              = b
          osc_dtheta(0)   = 0.d0
          osc_dr(0)       = ( rt - r ) / ( r - r1 )
          call calc_phi(0,iangset,iphi,b,costinf,z0,st,phiplus)
          call calc_dmu(0,iangset,st/rt)
      endif
#     endif
c
c     Skip rest
c
      goto 390
#endif
c
c     If mu=+-1 or if irsi_frsizey.eq.1 then don't worry about theta
c
      if((b.eq.0.d0).or.(itcross.eq.0)) then
          is   = 1
          iosc_icross(0)  = 1 
          osc_ds(0)       = s1 - s0
          st              = s1
          rt              = sqrt( b*b + st*st )
          osc_dtheta(0)   = 0.d0
          osc_dr(0)       = idr
          call calc_phi(0,iangset,iphi,b,costinf,z0,st,phiplus)
          call calc_dmu(0,iangset,st/rt)
          goto 390
      endif
c
c     Now we are not sure if the characteristic will cross first
c     a Theta=const grid line or first an R=const grid line. To
c     find out, we try to solve for the crossings with Theta=const
c
      itdir = itdir0
c
c     The index of the next theta (i.e. next
c     theta in downwind direction)
c
      itnext = it + itdir
c
c     Try to find the crossing with the next theta coordinate
c     Same story as above...
c     
      if((itnext.lt.1).or.(itnext.gt.irsi_frsizey)) then
          det = -1.d0 
      else
          costh = cos(rsi_x_c(itnext,2))
          det   = z0*z0 + b*b*(costinf*costinf-costh*costh)
      endif
      if(det.le.epsdet) then
c
c         No solutions for crossing with theta=const grid line...
c
          if((sother0.lt.s0).or.(sother0.gt.s1)) then
c
c             Reached R=const grid line
c
              iosc_icross(0)  = 1 
              osc_ds(0)       = s1 - s0
              st              = s1
              rt              = sqrt( b*b + st*st )
              if(b.gt.0.d0) then
                  costh1     = ( z0 + costinf*s1 ) / sqrt(b*b+s1*s1) 
                  if(abs(costh1).gt.1.d0) then
                      write(*,*) 'ERROR: costh1 out of range'
                      stop 170
                  endif
                  th1        = acos(costh1)
              else
                  costh1     = cos(theta)
                  th1        = theta
              endif
              do i=1,irsi_frsizey
                  if(rsi_x_c(i,2).gt.th1) goto 320
              enddo
  320         continue
              i = i - 1
              if(i.lt.it) then
                  osc_dtheta(0)  = i - it + 1 +
     %                (th1-rsi_x_c(i+1,2)) / 
     %                (rsi_x_c(i+1,2)-rsi_x_c(i,2)) 
              else
                  osc_dtheta(0)  = i - it +
     %                (th1-rsi_x_c(i,2)) / 
     %                (rsi_x_c(i+1,2)-rsi_x_c(i,2)) 
              endif
              osc_dr(0) = idr
              call calc_phi(0,iangset,iphi,b,costinf,z0,st,phiplus)
              call calc_dmu(0,iangset,st/rt)
          else
c
c             There is an extremum in theta, so the characteristic
c             is crossing the same theta=const grid line again
c              
              iosc_icross(0) = 2
              st = sother0
              rt = sqrt( b*b + st*st )
              osc_ds(0) = st - s0
              if(r1.eq.r) then
                  write(*,*) '?'
                  stop 150
              endif
              if(r1.gt.r) then
                  osc_dr(0) = ( rt - r ) / ( r1 - r )
              else
                  osc_dr(0) = ( rt - r ) / ( r - r1 )
              endif
              osc_dtheta(0)  = 0
              call calc_phi(0,iangset,iphi,b,costinf,z0,st,phiplus)
              call calc_dmu(0,iangset,st/rt)
          endif
      else
c     
c         There are solutions for this theta, but the question
c         is: are they within the domain in s?
c
c         Find the two solutions for s
c
          dum1  = costh * sqrt( det )
          dum2  = costh*costh - costinf*costinf
          if(dum2.eq.0.d0) then
              write(*,*) 'ERROR: sol1,sol2 undefined'
          endif
          sol1  = ( z0*costinf + dum1 ) / dum2
          sol2  = ( z0*costinf - dum1 ) / dum2
c
c         Now check if they lie within the domain <s0,s1]
c         (remember we are tracing downwind towards higher s).
c         
          s0 = s0 + eps2
          if(((sol1.gt.s1).or.(sol1.lt.s0)).and.
     %         ((sol2.gt.s1).or.(sol2.lt.s0))) then
c         
c             Both sol1 and sol2 are outside of the domain. 
c             This means that we must have arrived at an
c             R=const grid line.
c         
              iosc_icross(0)  = 1 
              osc_ds(0)       = s1 - s0
              st              = s1
              rt              = sqrt( b*b + st*st )
              if(b.gt.0.d0) then
                  costh1     = ( z0 + costinf*s1 ) / sqrt(b*b+s1*s1) 
                  if(abs(costh1).gt.1.d0) then
                      write(*,*) 'ERROR: costh1 out of range'
                      stop 170
                  endif
                  th1        = acos(costh1)
              else
                  costh1     = cos(theta)
                  th1        = theta
              endif
              do i=1,irsi_frsizey
                  if(rsi_x_c(i,2).gt.th1) goto 350
              enddo
  350         continue
              i = i - 1
              if(i.lt.it) then
                  osc_dtheta(0)  = i - it + 1 +
     %                (th1-rsi_x_c(i+1,2)) / 
     %                (rsi_x_c(i+1,2)-rsi_x_c(i,2)) 
              else
                  osc_dtheta(0)  = i - it +
     %                (th1-rsi_x_c(i,2)) / 
     %                (rsi_x_c(i+1,2)-rsi_x_c(i,2)) 
              endif
              osc_dr(0) = idr
              call calc_phi(0,iangset,iphi,b,costinf,z0,st,phiplus)
              call calc_dmu(0,iangset,st/rt)
          else
c         
c             At least one of them lies within the domain.
c             Find which one is the closest from the current
c             position
c          
              if((sol1.gt.s1).or.(sol1.lt.s0)) then
c         
c                 If not sol1 then it must be sol2 only
c         
                  st     = sol2
                  sother = sol1
              else if((sol2.gt.s1).or.(sol2.lt.s0)) then
c         
c                 It must be sol1 only
c         
                  st     = sol1
                  sother = sol2
              else
c         
c                 Both of them... So which one is the 
c                 closest to scur?
c         
                  if((sol1-s0).lt.(sol2-s0)) then
                      st     = sol1
                      sother = sol2
                  else
                      st     = sol2
                      sother = sol1
                  endif
              endif
c         
c             Now fill the array again
c         
              iosc_icross(0) = 2
              rt = sqrt( b*b + st*st )
              osc_ds(0) = st - s0 
              if(r1.eq.r) then
                  write(*,*) '?'
                  stop 150
              endif
              if(r1.gt.r) then
                  osc_dr(0) = ( rt - r ) / ( r1 - r )
              else
                  osc_dr(0) = ( rt - r ) / ( r - r1 )
              endif
              osc_dtheta(0)  = itnext - it
              call calc_phi(0,iangset,iphi,b,costinf,z0,st,phiplus)
              call calc_dmu(0,iangset,st/rt)
          endif          
      endif
  390 continue
      isc_ready = 1
c      return
      goto 159
#endif /* ifdef COORD_SPHERICAL */
#ifdef COORD_CARTESIAN
c
c     Make the short characteristics for 1-D plan-parallel.
c     This is easy. Later we might consider implementing 
c     Cartesian 2-D. But for now we skip that possibility.
c
c     The drrp and drrn have no meaning here. In the transfer
c     routines the ds will be scaled with the Delta X at the
c     position where the transfer is done. Here in the present
c     routine the ds will be always normalized (in 1-D) to 
c     ds=1.0 for mu=1.0 or mu=-1.0.
c     
c     It is recommended to take iangset=1 always.
c
      mu      = rmu(imu,iangset)
      phi     = rphi(iphi,iangset)      
      iosc_icross(1) = 1
      osc_ds(1)      = abs ( 1.d0 / ( mu + epsmu ) )
      if(mu.gt.0.d0) then
          osc_dr(1) = -1.d0
      else
          osc_dr(1) = 1.d0
      endif
      osc_dtheta(1) = 0.d0
      osc_phi(1)    = phi
      osc_mu(1)     = mu
      osc_dphi(1)   = 0.d0
      osc_dmu(1)    = 0.d0 
      iosc_iphi(1)  = iphi 
      iosc_imu(1)   = imu 
      iosc_icross(0) = 1
      osc_ds(0)      = abs ( 1.d0 / ( mu + epsmu ) )
      if(mu.lt.0.d0) then
          osc_dr(0) = -1.d0
      else
          osc_dr(0) = 1.d0
      endif
      osc_dtheta(0) = 0.d0
      osc_phi(0)    = phi
      osc_mu(0)     = mu
      osc_dphi(0)   = 0.d0
      osc_dmu(0)    = 0.d0 
      iosc_iphi(0)  = iphi 
      iosc_imu(0)   = imu 
      iosc_nel      = 1
c
      isc_ready = 1
      goto 159
c      return
c
#endif /* ifdef COORD_CARTESIAN */
c
      write(*,*) 'No shortchars made...'
      stop 13
  159 continue
      return
      end


c     --------------------------------------------------------------
c                         COMPUTE PHI AND SUCH
c     
c     This routine computes all the relevant information about
c     the local angle phi in the photon momentum space at the
c     current position along the characteristic. This information
c     is immediately stored in the /onesc/ arrays.
c
c      Input:
c        is         = The index along the characteristic
c        iangset    = The index of the set of angles
c        b,costinf,
c          z0       = Parameters of the characteristic
c        st         = The current position along the characteristic
c        phiplus    = Since the algorithm cannot distinguish
c                     whether the photon is moving towards positive
c                     of negative Phi, this should be told by the
c                     caller.
c
c     --------------------------------------------------------------
      subroutine calc_phi(is,iangset,iphi,b,costinf,z0,st,phiplus)
      implicit none
      doubleprecision b,costinf,z0,st
      integer is,iangset,phiplus,iphi
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      doubleprecision sinphi,phitmp
      doubleprecision pi,pihalf,pitwice,epssin
      integer i
c
      parameter(epssin=1.d-5)
      parameter(pi=PICONST)
      parameter(pitwice=2.d0*PICONST)
      parameter(pihalf=0.5d0*PICONST)
c
#ifdef RADGRID_ONEDIM
      iosc_iphi(is) = iphi 
      osc_phi(is)   = rphi(iphi,iangset)
      osc_dphi(is)  = 0.d0
      goto 300
#endif
c
      if(b.gt.0.d0) then
c     
c         If b>0 then the phi is well-defined
c     
c         The sin(phi) is given as an explicit formula:
c
          sinphi = ( b*b*costinf - z0*st ) / 
     %        ( b*sqrt(b*b+st*st-(z0+costinf*st)**2) )
          if(abs(sinphi-epssin*sinphi).ge.1.d0) then
              write(*,*) 'ERROR: sinphi out of range'
              stop 160
          endif
c
c         Find the phi itself
c
          osc_phi(is) = asin(sinphi-epssin*sinphi)
c
c         Depending on the direction of the characteristic (towards
c         positive or towards negative coordinate Phi?) we 
c         either mirror the phi or not.
c
          if(phiplus.lt.0) then
              osc_phi(is) = pi - osc_phi(is)
          endif
c
c         Negative phi ---> one period advanced so that it 
c         becomes positive.
c
          if(osc_phi(is).lt.0.d0) then
              osc_phi(is) = osc_phi(is) + pitwice
          endif
c
c         Find the index of the phi on the grid rphi
c         This index is such that:
c             rphi(iosc_iphi(is)) < phi < rphi(iosc_iphi(is)+1)
c
          do i=1,nrphi(iangset)
              if(rphi(i,iangset).gt.osc_phi(is)) goto 210
          enddo
  210     continue
          iosc_iphi(is) = i - 1
c
c         Treat the region aroung phi=0 or 2pi separately.
c         When 0<phi<rphi(1) or rphi(nrphi)<phi<2pi, then
c         always choose as the index iphi=nrphi.
c
          phitmp = osc_phi(is)
          if(i.eq.1) then
              iosc_iphi(is) = nrphi(iangset)
              i             = nrphi(iangset) + 1
              phitmp        = phitmp + pitwice
          endif
c
c         Now determine the dphi
c
          osc_dphi(is) = ( phitmp - rphi(i-1,iangset) ) / 
     %        ( rphi(i,iangset) - rphi(i-1,iangset) ) 
      else
c     
c         But if b=0, then there exists no phi (degenerated)
c         So take the iphi as prescribed in the arguments.
c     
          iosc_iphi(is) = iphi 
          osc_phi(is)   = rphi(iphi,iangset)
          osc_dphi(is)  = 0.d0
      endif
  300 continue
c
      return
      end



c     --------------------------------------------------------------
c                    COMPUTE MU-INDEX AND DMU
c     
c     The mu is simple to obtain, but now we must find to which index
c     it belongs in the mu-grid, and which is the value of dmu.
c     This routine is used by the short characteristics generator.
c     --------------------------------------------------------------
      subroutine calc_dmu(is,iangset,mu)
      implicit none
      integer is,iangset
      doubleprecision mu
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      doubleprecision epsmu
      integer i
c
      parameter(epsmu=1.d-2)
c     
      osc_mu(is) = mu
c
c     Find out which index belongs to this mu. If the anggrid_mu_zero
c     is 0, then the mu(0)=0 index is non-existent (i.e. should be 
c     ignored), so we skip that one then. Else we take it into account.
c
      do i=-nrmu(iangset)+1,nrmu(iangset)-1
          if((rmu(i,iangset).ge.mu).and.((anggrid_mu_zero.eq.1)
     %        .or.(i.ne.0))) then
              goto 10
          endif
      enddo
   10 continue
c
c     Now determine the position in between the two adjacent gridpoints
c     of this mu.
c
      iosc_imu(is) = i-1
      osc_dmu(is)  = ( mu - rmu(i-1,iangset) )  / 
     %               ( rmu(i,iangset) - rmu(i-1,iangset) )
      if((abs(1.d0-osc_dmu(is)).lt.epsmu).and.
     %            (i.lt.nrmu(iangset)-1)) then
          iosc_imu(is) = iosc_imu(is) + 1
          osc_dmu(is) = 0.d0
      endif
      if(iosc_imu(is).eq.nrmu(iangset)) then
          write(*,*) '---------> ',rmu(nrmu(iangset),iangset),
     %                  mu,i,nrmu(iangset)
      endif
      return
c     
      end




c     --------------------------------------------------------------
c                    COPY THE SHORT CHARACTERISTIC INTO 
c                             THE BIG ARRAY
c     --------------------------------------------------------------
      subroutine copy_shortchar(idrr,it,imu,iphi)
      implicit none
      integer idrr,it,imu,iphi
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      integer is,iangset,itmax
c
c     First check if there has not been an error so far
c
      if(isc_error.ne.0) then
          return
      endif
c
c     Not, okay then store in array
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#else
      itmax = irsi_frsizey
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
      iangset = isc_iangset(idrr)
      if((idrr.lt.0).or.(idrr.gt.FRSIZE_DRR)) then
          write(*,*) 'ERROR: idrr out of range in copy_shortchar()'
          stop 13
      endif
ccccc (FROM NOW: ONLY COMPACT 14-05-01) if((it.lt.1).or.(it.gt.irsi_frsizey)) then
      if((it.lt.1).or.(it.gt.itmax)) then
          write(*,*) 'ERROR: it out of range in copy_shortchar()'
          stop 13
      endif
      if((imu.lt.-nrmu(iangset)).or.(imu.gt.nrmu(iangset))) then
          write(*,*) 'ERROR: imu out of range in copy_shortchar()'
          stop 13
      endif
      if((iphi.lt.1).or.(iphi.gt.nrphi(iangset))) then
          write(*,*) 'ERROR: iphi out of range in copy_shortchar()'
          stop 13
      endif
c
#     ifndef SHORCHAR_COMPACT_ARRAYS
c
c     Long arrays: dimension 'is' is always FRSIZE_CHAR long.
c     This eats a lot of unnenecssary memory... OLD METHOD,
c     IF NEW METHOD WORKS: OLD ONE WILL BE REMOVED.
c
      isc_nel(iphi,imu,it,idrr) = iosc_nel
      do is=0,iosc_nel
          sc_ds(is,iphi,imu,it,idrr)      = osc_ds(is)
          sc_dr(is,iphi,imu,it,idrr)      = osc_dr(is)
          sc_dtheta(is,iphi,imu,it,idrr)  = osc_dtheta(is)
          sc_phi(is,iphi,imu,it,idrr)     = osc_phi(is)
          sc_dphi(is,iphi,imu,it,idrr)    = osc_dphi(is)
          sc_dmu(is,iphi,imu,it,idrr)     = osc_dmu(is)
          sc_mu(is,iphi,imu,it,idrr)      = osc_mu(is)
          isc_icross(is,iphi,imu,it,idrr) = iosc_icross(is)
          isc_iphi(is,iphi,imu,it,idrr)   = iosc_iphi(is)
          isc_imu(is,iphi,imu,it,idrr)    = iosc_imu(is)
      enddo
#     else
c
c     Compact arrays: Only when an ESC is longer than 1 element,
c     will memory space be requested. IF THIS METHOD WORKS, THEN
c     THIS WILL BECOME STANDARD.
c 
      isc_nel(iphi,imu,it,idrr) = iosc_nel
      do is=0,1
          sc_ds(is,iphi,imu,it,idrr)      = osc_ds(is)
          sc_dr(is,iphi,imu,it,idrr)      = osc_dr(is)
          sc_dtheta(is,iphi,imu,it,idrr)  = osc_dtheta(is)
          sc_phi(is,iphi,imu,it,idrr)     = osc_phi(is)
          sc_dphi(is,iphi,imu,it,idrr)    = osc_dphi(is)
          sc_dmu(is,iphi,imu,it,idrr)     = osc_dmu(is)
          sc_mu(is,iphi,imu,it,idrr)      = osc_mu(is)
          isc_icross(is,iphi,imu,it,idrr) = iosc_icross(is)
          isc_iphi(is,iphi,imu,it,idrr)   = iosc_iphi(is)
          isc_imu(is,iphi,imu,it,idrr)    = iosc_imu(is)
      enddo
      if(iosc_nel.gt.1) then
          isc_longesc_current = isc_longesc_current + 1
          if(isc_longesc_current.gt.FRSIZE_LONGESC) then
              if(isc_error_longesc.eq.0) then
                  write(*,*) '---------------------------------------'
                  write(*,*) 'ERROR in creating short characteristics'
                  write(*,*) 'according to the compact storage method'
                  write(*,*) 'LONGESC_FRACTION is too small!! Change '
                  write(*,*) 'this in the configure.h file!!!'
                  write(*,*) '---------------------------------------'
c                  stop 13
              endif
              isc_error_longesc = 1
          else
              isc_longesc(iphi,imu,it,idrr)   = isc_longesc_current
              do is=0,iosc_nel
                  esc_ds(is,isc_longesc_current)      = osc_ds(is)
                  esc_dr(is,isc_longesc_current)      = osc_dr(is)
                  esc_dtheta(is,isc_longesc_current)  = osc_dtheta(is)
                  esc_phi(is,isc_longesc_current)     = osc_phi(is)
                  esc_dphi(is,isc_longesc_current)    = osc_dphi(is)
                  esc_dmu(is,isc_longesc_current)     = osc_dmu(is)
                  esc_mu(is,isc_longesc_current)      = osc_mu(is)
                  iesc_icross(is,isc_longesc_current) = iosc_icross(is)
                  iesc_iphi(is,isc_longesc_current)   = iosc_iphi(is)
                  iesc_imu(is,isc_longesc_current)    = iosc_imu(is)
                  iesc_nel(isc_longesc_current)       = iosc_nel
              enddo
          endif
      else
          isc_longesc(iphi,imu,it,idrr)   = 0
      endif
c
#     endif
c
c     Just for fun: count nr of SC
c
      isc_shortsc_current = isc_shortsc_current + 1
c
      end


c     --------------------------------------------------------------
c                WRITE SHORT CHARACTERISTICS TO FILE
c     --------------------------------------------------------------
      subroutine write_shortchar(idrr)
      implicit none
      integer idrr
c
#     ifndef SHORCHAR_COMPACT_ARRAYS
#include "common_grid.h"
#include "common_shortchar.h"
c
      integer iphi,is,imu,iy,nel,iangset
      character*80 filename,base,ext
c
      iangset = isc_iangset(idrr)
c
c     Make filename
c
      base='shortchar_'
      ext ='.dat'
      call make_indexed_filename(base,idrr,ext,filename)
c
c     Open file
c
      open(unit=1,file=filename,status='unknown')
      write(1,*) irsi_frsizey,anggrid_frsizemu,anggrid_frsizephi
      write(1,*) isc_ismax(idrr)
      write(1,*) sc_drrp(idrr),sc_drrn(idrr)
      write(1,*) ' '
c   
c     Loop over the Theta direction
c
      do iy=1,irsi_frsizey
c
c         Loop over mu, but skip the mu(0) if the mu=0.0 gridpoint
c         is non-existent (i.e. if anggrid_mu_zero.eq.0, or in other
c         words if anggrid_frsizemu is an even number).
c          
          do imu=-nrmu(iangset),nrmu(iangset)
c
c             Check if mu(0) is non-existent, if so then skip imu=0
c
              if((imu.eq.0).and.(anggrid_mu_zero.eq.0)) then 
                  goto 23
              endif
c
c             Loop over phi
c
              do iphi=1,nrphi(iangset)
                  nel = isc_nel(iphi,imu,iy,idrr)
                  write(1,*) nel
c
c                 Loop over s
c                 
                  do is=0,nel
                      write(1,11) isc_icross(is,iphi,imu,iy,idrr),
     %                      sc_ds(is,iphi,imu,iy,idrr),
     %                      sc_dr(is,iphi,imu,iy,idrr),
     %                      sc_dtheta(is,iphi,imu,iy,idrr),
     %                      sc_mu(is,iphi,imu,iy,idrr),
     %                      sc_phi(is,iphi,imu,iy,idrr)
   11                 format(I4,1X,E13.6,1X,E13.6,1X,E13.6,1X,
     %                       E13.6,1X,E13.6)
#ifdef GRIDCHAR_WRITE_IANG
                      write(1,12) sc_dmu(is,iphi,imu,iy,idrr),
     %                       sc_dphi(is,iphi,imu,iy,idrr),
     %                       isc_imu(is,iphi,imu,iy,idrr),
     %                       isc_iphi(is,iphi,imu,iy,idrr)
   12                 format(4X,E13.6,1X,E13.6,1X,I4,1X,I4)
#endif
                  enddo
                  write(1,*) ' '
              enddo
   23         continue
          enddo
      enddo
c
      close(1)
#else  /* ifndef SHORCHAR_COMPACT_ARRAYS */
      write(*,*) 'CANNOT WRITE SHORTCHAR IN COMPACT METHOD'
      stop 13
#endif /* ifndef SHORCHAR_COMPACT_ARRAYS */
c
      end


c     -------------------------------------------------------------
c                 SUBROUTINE FOR GETTING SHORTCHAR INFO
c     
c     This routine retrieves the information about the current
c     shortchar element.
c     --------------------------------------------------------------
      subroutine get_shortchar(iel,iphi,imu,it,idrr,
     %                       icr,dr,dt,mu,phi,dmu,dphi,imus,iphis)
      implicit none
      integer iel,iphi,imu,it,idrr,icr,imus,iphis
      doubleprecision dr,dt,dmu,dphi,mu,phi
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      integer ilesc
c
#     ifndef SHORCHAR_COMPACT_ARRAYS
      icr   = isc_icross(iel,iphi,imu,it,idrr)
      dr    = sc_dr(iel,iphi,imu,it,idrr)
      dt    = sc_dtheta(iel,iphi,imu,it,idrr)
      mu    = sc_mu(iel,iphi,imu,it,idrr)
      phi   = sc_phi(iel,iphi,imu,it,idrr)
      dmu   = sc_dmu(iel,iphi,imu,it,idrr)
      dphi  = sc_dphi(iel,iphi,imu,it,idrr)
      imus  = isc_imu(iel,iphi,imu,it,idrr)
      iphis = isc_iphi(iel,iphi,imu,it,idrr)
#     else
      if(iel.le.1) then
          icr   = isc_icross(iel,iphi,imu,it,idrr)
          dr    = sc_dr(iel,iphi,imu,it,idrr)
          dt    = sc_dtheta(iel,iphi,imu,it,idrr)
          mu    = sc_mu(iel,iphi,imu,it,idrr)
          phi   = sc_phi(iel,iphi,imu,it,idrr)
          dmu   = sc_dmu(iel,iphi,imu,it,idrr)
          dphi  = sc_dphi(iel,iphi,imu,it,idrr)
          imus  = isc_imu(iel,iphi,imu,it,idrr)
          iphis = isc_iphi(iel,iphi,imu,it,idrr)
      else
          ilesc = isc_longesc(iphi,imu,it,idrr)
          if(ilesc.eq.0) then
              write(*,*) 'INTERNAL ERROR IN GET_SHORTCHAR()'
              stop 4234
          endif
          if(ilesc.gt.FRSIZE_LONGESC) then
              write(*,*) 'INTERNAL ERROR 2 IN GET_SHORTCHAR()'
              stop 6234
          endif
          icr   = iesc_icross(iel,ilesc)
          dr    = esc_dr(iel,ilesc)
          dt    = esc_dtheta(iel,ilesc)
          mu    = esc_mu(iel,ilesc)
          phi   = esc_phi(iel,ilesc)
          dmu   = esc_dmu(iel,ilesc)
          dphi  = esc_dphi(iel,ilesc)
          imus  = iesc_imu(iel,ilesc)
          iphis = iesc_iphi(iel,ilesc)
      endif
#     endif
c
      end



c     -------------------------------------------------------------
c                 SUBROUTINE FOR GETTING SHORTCHAR INFO
c     
c     This routine retrieves the information about the current
c     shortchar element.
c     --------------------------------------------------------------
      subroutine get_shortchar_concise(iel,iphi,imu,it,idrr,
     %                      icr,dr,dt,mu,phi)
      implicit none
      integer iel,iphi,imu,it,idrr,icr
      doubleprecision dr,dt,mu,phi
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      integer ilesc
c
#     ifndef SHORCHAR_COMPACT_ARRAYS
      icr   = isc_icross(iel,iphi,imu,it,idrr)
      dr    = sc_dr(iel,iphi,imu,it,idrr)
      dt    = sc_dtheta(iel,iphi,imu,it,idrr)
      mu    = sc_mu(iel,iphi,imu,it,idrr)
      phi   = sc_phi(iel,iphi,imu,it,idrr)
#     else
      if(iel.le.1) then
          icr   = isc_icross(iel,iphi,imu,it,idrr)
          dr    = sc_dr(iel,iphi,imu,it,idrr)
          dt    = sc_dtheta(iel,iphi,imu,it,idrr)
          mu    = sc_mu(iel,iphi,imu,it,idrr)
          phi   = sc_phi(iel,iphi,imu,it,idrr)
      else
          ilesc = isc_longesc(iphi,imu,it,idrr)
          if(ilesc.eq.0) then
              write(*,*) 'INTERNAL ERROR IN GET_SHORTCHAR()'
              stop 4234
          endif
          if(ilesc.gt.FRSIZE_LONGESC) then
              write(*,*) 'INTERNAL ERROR 2 IN GET_SHORTCHAR()'
              stop 6234
          endif
          icr   = iesc_icross(iel,ilesc)
          dr    = esc_dr(iel,ilesc)
          dt    = esc_dtheta(iel,ilesc)
          mu    = esc_mu(iel,ilesc)
          phi   = esc_phi(iel,ilesc)
      endif
#     endif
c
      end


c     -------------------------------------------------------------
c                 SUBROUTINE FOR GETTING SHORTCHAR INFO
c     
c     This routine retrieves the information about the current
c     shortchar element.
c     --------------------------------------------------------------
      subroutine get_shortchar_ds(iel,iphi,imu,it,idrr,ds)
      implicit none
      integer iel,iphi,imu,it,idrr
      doubleprecision ds
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      integer ilesc
c
#     ifndef SHORCHAR_COMPACT_ARRAYS
      ds    = sc_ds(iel,iphi,imu,it,idrr)
#     else
      if(iel.le.1) then
          ds    = sc_ds(iel,iphi,imu,it,idrr)
      else
          ilesc = isc_longesc(iphi,imu,it,idrr)
          if(ilesc.eq.0) then
              write(*,*) 'INTERNAL ERROR IN GET_SHORTCHAR()'
              stop 4237
          endif
          if(ilesc.gt.FRSIZE_LONGESC) then
              write(*,*) 'INTERNAL ERROR 2 IN GET_SHORTCHAR()'
              stop 6234
          endif
          ds    = esc_ds(iel,ilesc)
      endif
#     endif
c
      end


c     --------------------------------------------------------------
c                          WRITE THE ESC HISTOGRAM
c     --------------------------------------------------------------
      subroutine write_esc_histogram()
      implicit none
c
#include "common_grid.h"
#include "common_shortchar.h"
c
      integer i
c
      open(unit=1,file='esc_histogram.dat',status='unknown')
      write(1,*) FRSIZE_CHAR
      do i=1,FRSIZE_CHAR
          write(1,*) dbsc_hist_nel(i)
      enddo
      close(1)
c
      end



#endif /* ifndef NO_SHORT_CHARS */ 







