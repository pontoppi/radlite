c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c        ( MODULE: ITERATION PROCEDURE FOR THE RADIATIVE TRANSFER)
c
c                            Leiden, Aug 1998
c                             C.P. Dullemond   
c
c     =============================================================


#define ITERATION_F_VERSION "Trans_15aug98"
#include "main.h"

#ifndef ONLY_RAY_TRACING

#define TAU_EPS 1.d-4

c     --------------------------------------------------------------
c                    ITERATION OF THE RADIATIVE TRANSFER
c
c     This is the main routine for the iteration of the transfer, 
c     to obtain the scattering- and non-LTE data. The iteration 
c     procedure updates the radiation field and the sources. 
c     Sometimes one also wishes to update the temperature alongc
c     with the radiation field. For the dust module this is all
c     done within the dust solve procedure, and is of no concern here.
c     But for other procecces, such as line transfer or compton
c     transfer, the temperature iteration is a separate issue. 
c     Here it is possible to specify that the temperature iteration
c     be performed.
c
c     ARGUMENTS:
c       itermax       The maximum total number of iterations
c       start         Normally this should be 1. Only if the
c                          iteration is re-started from a previous
c                          run, then this value should be equal to 
c                          the previous itermax PLUS 1, and be
c                          sure to set the itermax to a larger
c                          value this time. 
c       method        0  = LI
c                          1  = ALI
c                          2  = ALI + Ng
c                          3  = LI + Ng
c                          -1 = IS
c       guess         Type of initial guess used
c       convcrit      The confidence level for convergence
c       itertemp      0  = No temperature iteration
c                     >0 = The maximum nr of iterations before the
c                     temperature solving routine is called
c       convcrittemp  The convergence criterion for the temperature
c                     iteration.
c     
c     
c     NOTE: When this routine is finished, it is assumed that the
c           Non-LTE stuff has been done, and therefore the flag
c           iscat_done_nonlte will be set to 1. This flag can be
c           used by routines that depend on the Non-LTE data to verify
c           that indeed this data has either been read from a file OR 
c           has been computed.
c     --------------------------------------------------------------
      subroutine iteration_transfer(start,itermax,method,guess,
     %                convcrit,itertemp,convcrittemp,iterdone,
     %                converged) 
      implicit none
      integer start,itermax,method,guess,itertemp,iterdone
      doubleprecision convcrit,convcrittemp
      logical converged
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
#include "common_isave.h"
#include "common_lines.h"
#include "common_dust.h"
#include "common_check.h"
c
      integer iter,ifreq
c
c     Reset error count
c
      transerr_overshoot = 0
c
c     First check some things
c
cc#     ifdef DEBUG_TEMP_ITER
cc      write(*,*) '------------------------------------------'
cc      write(*,*) '         WARNING: DEBUGGING MODE          '
cc      write(*,*) '------------------------------------------'
cc#     endif
      if((itertemp.ne.0).and.(iradproc_solve_temp.eq.0)) then
          write(*,*) 'Problem: iteration.F/iteration_transfer()'
          write(*,*) '         Wanting to iterate on temperature,'
          write(*,*) '         but the iradproc_solve_temp is not set'
          stop 13
      endif
      if((itertemp.ne.0).and.(method.eq.-1)) then
          write(*,*) 'ERROR: Cannot use temperature iteration with'
          write(*,*) '       iterative scattering method'
          stop 13
      endif
      if(start.gt.itermax) then
          write(*,*) 'WARNING: No iterations done'
          write(*,*) '         Starting iteration number ',
     %           'is greater than amount of iterations'
          return
      endif
      if((start.gt.1).and.(guess.ne.-1)) then
          write(*,*) 'WARNING: When restarting, the ',
     %               'initial guess should be -1'
          write(*,*) '         I will do this...'
          guess = -1
      endif
      if(iradproc_dust.eq.1) then
c          if((method.ne.0).and.(method.ne.3)) then
c              write(*,*) 'Dust transfer only Lambda Iteration!'
c              method = 0
c          endif
#         ifndef INCLUDE_DUST
          write(*,*) 'Dust RT requested, while not compiled with ',
     %            'INCLUDE_DUST.'
          write(*,*) 'Recompile with -DINCLUDE_DUST'
          stop 13
#         endif
      endif
      if(iradproc_line.eq.1) then
#         ifndef INCLUDE_LINES
          write(*,*) 'Line RT requested, while not compiled with ',
     %            'INCLUDE_LINES.'
          write(*,*) 'Recompile with -DINCLUDE_LINES'
          stop 13
#         endif
      endif
c
c     Now reset the convergence error file
c
#     ifdef WRITE_CONVERGENCE
      open(unit=4,file='conv.dat',status='unknown')
      write(4,*) ' '
      close(4)
#     endif
c
c     Some warnings
c
      if((guess.eq.1).and.(iradproc_dust.ne.0)) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Do you really want to have initial guess   |'
          write(*,*) '|  type 1 for dust radiative transfer?        |'
          write(*,*) '|  It is very well OK, but experience tells   |'
          write(*,*) '|  that it is not such a good choice. Type    |'
          write(*,*) '|  0 is better.                               |'
          write(*,*) '-----------------------------------------------'
      endif
      if((guess.ne.1).and.(iradproc_line.ne.0)) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Do you really want to have initial guess   |'
          write(*,*) '|  type 0 for line ransfer?                   |'
          write(*,*) '|  It is very well OK, but experience tells   |'
          write(*,*) '|  that it is not such a good choice. Type    |'
          write(*,*) '|  1 is better.                               |'
          write(*,*) '-----------------------------------------------'
      endif
c
c     A debugging tool
c
#ifdef DEBUG_FILLINTENS
      call debug_fill_intens()
#endif
c
c     Figure out which iteration method must be used
c
      if(method.eq.0) then
          call lambda_iteration(start,itermax,0,guess,
     %         convcrit,itertemp,convcrittemp,iterdone,converged)
      elseif(method.eq.1) then
          call accel_lambda_iteration(start,itermax,
     %         0,guess,convcrit,itertemp,convcrittemp,
     %         iterdone,converged)
      elseif(method.eq.2) then
          call accel_lambda_iteration(start,itermax,
     %         1,guess,convcrit,itertemp,convcrittemp,
     %         iterdone,converged)
      elseif(method.eq.3) then
          call lambda_iteration(start,itermax,1,guess,
     %         convcrit,itertemp,convcrittemp,iterdone,converged)
      elseif(method.eq.-1) then
          call iterative_scattering(itermax,convcrit,iterdone,converged)
      else
          write(*,*) 'I dont know iteration type ',method
          stop 13
      endif
c
c     Check if there were warning flags set
c
      if(warn_bug_centrbeam.ne.0) then
          write(*,*) ' ' 
          write(*,*) '###############################################' 
          write(*,*) '#  WARNING:                                   #'
          write(*,*) '#  Uncorrected bug: the central beam bug in   #'
          write(*,*) '#  the do_trans() routine...                  #'
          write(*,*) '###############################################' 
      endif
      if(warn_mubin_theta.ne.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Mu-binnig happened on Theta-crossing.      |'
          write(*,*) '|  This is not very serious, but might be     |'
          write(*,*) '|  slightly decreasing reliability. Tip: use  |'
          write(*,*) '|  more mu-points, and this problem is solved.|'
          write(*,*) '-----------------------------------------------' 
      endif
      if(icool_warn_guessvertheight.ne.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  In the cooling formulae: I made the crude  |'
          write(*,*) '|  guess of H=R for the ADAF...               |'
          write(*,*) '-----------------------------------------------' 
      endif
#ifdef VERBOSE
      if(iscat_warn_int2.eq.1) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING: (repeat)                          |'
          write(*,*) '|  Second order source integration not set!   |'
          write(*,*) '|  Recompile with SRCQDR_PNT_3                |'
          write(*,*) '-----------------------------------------------' 
      endif
      if(iscat_warn_src2.eq.1) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING: (repeat)                          |'
          write(*,*) '|  Second order source interpolation not set! |'
          write(*,*) '|  Recompile with INTERPOL_SRC_2              |'
          write(*,*) '-----------------------------------------------' 
      endif
      if(iang_warn_mureslow.eq.1) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING: (repeat)                          |'
          write(*,*) '|  The mu-gridding is insufficient!           |'
          write(*,*) '|  Use larger NrMu!                           |'
          write(*,*) '-----------------------------------------------' 
      endif
#endif
#ifdef SAFETY_CHECKS_ACTIVE
      if(selfch_taumax.gt.SAFETY_CHECK_TAUMAX) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Found rather high one-cell optical depths  |'
          write(*,*) '|  Largest optical depth found is:            |'
          write(*,*) '|  ',selfch_taumax
          write(*,*) '-----------------------------------------------' 
      endif
#endif
#ifdef INCLUDE_LINES
      if(line_warn_cr_extrapol.ne.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Collision rates were computed outside of   |'
          write(*,*) '|  the temperature range of the tabulated     |'
          write(*,*) '|  lab-values. A blunt extrapolation was done.|'
          write(*,*) '-----------------------------------------------' 
      endif
      if(line_warn_negpop.ne.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Some level populations have become         |'
          write(*,*) '|  negative at least once during the transfer.|'
          write(*,*) '-----------------------------------------------' 
      endif
      if(line_warn_umass.ne.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  The umass_av (average molec mass) seems to |'
          write(*,*) '|  have been set previously to another value  |'
          write(*,*) '|  than is given by the line.inp file.        |'
          write(*,*) '-----------------------------------------------' 
      endif
      if(line_warn_lineout.ne.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Passband was too small to harbor the line  |'
          write(*,*) '|  profile at all width and velocity shifts.  |'
          write(*,*) '-----------------------------------------------' 
      endif
      if(line_warn_badsampling.ne.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Line profile frequency sampling is too     |'
          write(*,*) '|  coarse. Might result in errors.            |'
          write(*,*) '-----------------------------------------------' 
      endif
      if(line_warn_velostep.eq.1) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Velocity step exceeded line width.         |'
          write(*,*) '|  Might result in errors.                    |'
          write(*,*) '-----------------------------------------------' 
      endif
      if(line_warn_velostep.eq.2) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|         >>>>> RED ALTERT <<<<<<<<           |'
          write(*,*) '|  Velocity step MUCH larger than line width! |'
          write(*,*) '|  Big errors are expected!                   |'
          write(*,*) '-----------------------------------------------' 
      endif
#ifdef LINE_LOCAL_ABUNDANCE
      if(locabun_done.eq.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Abundances are taken globally (value from  |'
          write(*,*) '|  the line.dat file).                        |'
          write(*,*) '-----------------------------------------------' 
      endif
#endif
#ifndef LINE_NO_LIMIT_MASERS
      if(line_warn_opacity_neg.ne.0) then
          write(*,*) ' ' 
          write(*,*) '===============================================' 
          write(*,*) '|  MAJOR WARNING:                             |'
          write(*,*) '|  In last iteration step, the opacity became |'
          write(*,*) '|  negative at some point. I have bluntly put |'
          write(*,*) '|  this to zero!! I switched masers off!      |'
          write(*,*) '===============================================' 
      endif
#endif
      if((line_warn_opacity_neg.eq.0).and.
     %   (line_warn_neg_opacity.ne.0)) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Negative opacity was detected during the   |'
          write(*,*) '|  iteration procedure. Fortunatly not at the |'
          write(*,*) '|  final iteration step.                      |'
          write(*,*) '-----------------------------------------------' 
      endif
#endif
#ifdef INCLUDE_DUST
      if(dust_warn_few_freqs.ne.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Very few frequencies for dust temperature  |'
          write(*,*) '|  solver. Serious errors expected.           |'
          write(*,*) '-----------------------------------------------' 
      endif
      if(dust_warn_zero_temp.gt.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  In last iteration there were still points  |'
          write(*,*) '|  with zero (or very tiny) dust temperature. |'
          write(*,*) '-----------------------------------------------' 
      endif
#endif
#ifdef INCLUDE_COMPTON
      if(compt_warn_compt1.ne.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  Compton switch irradproc_compt is set      |'
          write(*,*) '|  to 1. This is only Thompson scattering.    |'
          write(*,*) '|  It is okay, but just remember that for     |'
          write(*,*) '|  isotropic CSK you should put it to 3,      |'
          write(*,*) '|  and for diffusion equation to 2.           |'
          write(*,*) '|----------------------------------------------' 
      endif
#endif
      if(isynch_warn_lowtemp.gt.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,163) isynch_warn_lowtemp,SYNCH_LOWER_T_CUTOFF
 163      format(' |  Found ',I4,' times an Te<',E9.2,' error in  |')
          write(*,*) '|  synchj() subroutine for synchrotron in the |'
          write(*,*) '|  latest iteration.                          |'
          write(*,*) '|----------------------------------------------' 
      endif
      if(isynch_warn_lownunub.gt.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,164) isynch_warn_lownunub
 164      format(' |  Found ',I4,' times an chi<10 error in the    |')
          write(*,*) '|  synchj() subroutine for synchrotron in the |'
          write(*,*) '|  latest iteration.                          |'
          write(*,*) '|----------------------------------------------' 
      endif
      if(transerr_overshoot.gt.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  During Short Characteristic integration    |'
          write(*,*) '|  there were overshoots in the quadrature.   |'
          write(*,*) '|  These have been corrected automatically,   |'
          write(*,*) '|  but they indicate poor grid resolution!    |'
          write(*,*) '-----------------------------------------------' 
      endif
      if(transerr_starresolution.gt.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  The central star is big enough to be       |'
          write(*,*) '|  resolved by the angular grid at the inner  |'
          write(*,*) '|  edge of the domain. But the mu-spacing is  |'
          write(*,*) '|  too course to sample the star properly.    |'
          write(*,*) '|  This could result in errors of more than   |'
          write(*,*) '|  10% in the stellar flux.                   |'
          write(*,*) '-----------------------------------------------' 
      endif
c
c     Set the Non-LTE / scattering flag
c
      iscat_done_nonlte = 1
c     
      end



c     --------------------------------------------------------------
c                      ORDINARY LAMBDA ITERATION
c
c     This is the slowly converging Lambda iteration method. No
c     accelleration applied.
c
c     ARGUMENTS:
c        start       The iteration count at which to (re)start
c        itermax     The maximum number of iterations
c        ng          Switch for Ng acceleration
c        init        The type of initial guess
c        convcrit    The convergence criterion for the sourcefunc
c        itertemp    =0 No iteration on temperature (default)
c                    >0 Maximum number of iterations of RT before
c                       solving temperature. =1 -> every iteration.
c        convcrtemp  The convergence criterion for the temperature
c
c     NOTE: The ssw_scat = 0 for the first iteration is silly. 
c           I'd prefer to do this differently, but for the moment
c           I don't have the guts to change thing (07-11-00).
c
c     --------------------------------------------------------------
      subroutine lambda_iteration(start,itermax,ng,init,convcrit,
     %                itertemp,convcrtemp,iterdone,converged)
      implicit none
c
      integer start,itermax,init,ng,itertemp,iterdone
      doubleprecision convcrit,convcrtemp
      doubleprecision conv_error,tempsolve_error
      logical converged
c
#include "common_grid.h"
#include "common_source.h"
#include "common_isave.h"
#include "common_lines.h"
#include "common_dust.h"
#include "common_main.h"

#ifdef PARALLEL
!for parallel execution we need the array intmom_0 
!declared in this headerfile
#include "common_radfield.h"

!parallel headerfile
#include "mpif.h"

      integer ierr         !error code
      integer MyRank       !rank of this process
      integer NrProc       !number of processes
      integer MyNrFreq     !number of frequencies for this process
      integer FreqBegin    !start of frequency domain
      integer FreqEnd      !end of frequency domain
      integer SendNrFr     !nr of freqs that are send to this process
      integer SendFrBegin  !start of frequency domain of sending process
      integer SendFrEnd    !end of frequency domain of sending process
      integer NrTheta, NrR !number of Theta and R points
      integer i,j,k,p      !counters

      !buffer for sending and receiving data
      doubleprecision SendRecvArray(1:FRSIZE_Y*FRSIZE_X)
#endif

c
      integer iter,ifreq,icnt_ng,icnt_temp,ididsolvetemp
c
#ifdef VERBOSE
      write(*,*) '===================================================='
      write(*,*) '====   Starting the Ordinary Lambda Iteration   ===='
      write(*,*) '===================================================='
#endif
c
c     If line radiation: check level populations, and compute the
c     sources for the first time.
c     
#     ifdef INCLUDE_LINES
c
      if(iradproc_line.eq.1) then
          call check_levelpopul()
#         ifndef LINE_VELOCITIES
          call calc_src_from_populations()
#         else
          call global_prepare_transitions()
c
c         NEW: 27.03.06: Removed the next, because it is now done
c                        in the read_linedata() between reading the
c                        line characteristics and setting up the
c                        passbands.
c
ccc          if(iradproc_line_dust.eq.1) then
ccc              call global_prepare_line_dust()
ccc          endif
#         endif
      endif
c
#     endif 
c     
c     If isotropic csk is selected, read the data files
c
c      if(iradproc_compt.eq.3) call read_isocsk()
c     
c     Enable the emission processes and all opacities, by setting the
c     appropriate flags
c     
      asw_lte  = 1
      asw_scat = 1
      ssw_lte  = 1
      ssw_scat = 0
c
c     If requested, make an initial guess. There is a variety of types
c     of initial guesses. They are selected by the value of init.
c
      if(init.eq.1) then
c
c         Initial guess = no scattering/non-lte (for continuum tranfer)
c                       = populations in RE without radiation field (for lines)
c
          write(*,*) 'Make initial guess'
c
c         If lines: then compute the LTE level populations
c
#         ifdef INCLUDE_LINES
c
          if(iradproc_line.eq.1) then
              call calc_re_levelpopul()
          endif
c
#         endif
c
c         Now do the transfer for each freq
c
          if(iradproc_line.eq.0) then
              asw_scat = 0
              ssw_scat = 0
              asw_lte  = 1
              ssw_lte  = 1
              do ifreq=1,freq_nr
c
c                 Perform the transfer in the array intens(iphi,imu,it,ir)
c
                  call do_trans(ifreq)
c
c                 Perform the angular integrals for the mean intensity (J)
c                 or perform the angular regridding to a lower angular 
c                 resolution (if necessary), and then store the result into 
c                 an array with a frequency index. 
c    
                  call intens_prep_store(ifreq)
c    
              enddo
              asw_scat = 1
              ssw_scat = 1
c    
c             Compute the scattering integrals, from the given J(nu) or 
c             I(inu,iphi,imu) (whichever is appropriate).
c    
              if(iradproc_line.eq.0) then
                  call solve_statequil_and_scat(0)
              endif
          endif
c
      elseif(init.eq.2) then
          if(iradproc_dust.ne.0) then
c
c             [ Dust radiative transfer ]
c
c             Init guess = read dust temperatures from file
c
#             ifdef INCLUDE_DUST
              call read_dust_temperatures(.true.)
              call scat_src_clear()
              call calc_emis_dust()
              ssw_scat = 1
#             else
              stop 676
#             endif
c
          elseif((iradproc_compt.ne.0).and.(itertemp.gt.0)) then
             write(*,*) 'Compton scattering disabled, Klaus 22/03/07'
             stop 13
          else
              write(*,*) 'ERROR: Init guess type 2 only for compton ',
     %               'with electron temperature solving.'
              stop 13
          endif
      elseif(init.eq.-1) then
c
c         Initial guess is already made, including the source
c
          ssw_scat = 1
c
      endif
c
c     If line transfer is active, then write the initial populations here
c
#     ifdef INCLUDE_LINES
      if(iradproc_line.eq.1) then
ccc          call write_levelpopul(0)
          call write_iteration_sub_result(1,0)
      endif
#     endif
c
c     If Line processes, then switch on ssw_scat!
c
#     ifdef INCLUDE_LINES
      if(iradproc_line.eq.1) then
          ssw_scat = 1
      endif
#     endif
c     
c     Reset the Ng counter. The Ng counter is meant to count to 4
c     and then reset, triggering the Ng accelleration at each 4-th
c     iteration step.
c
      icnt_ng = 0
c
c     Reset the temperature iteration counter
c
      icnt_temp = 0
c
c     Reset the temperature convergence error
c
      cnv_temp_error = 1.d0
c
c     For safety, backup the current sources (if e.g. needed for the
c     convergence criterion on sources).
c
      if(iradproc_line.eq.0) then
          call scat_copy(-1,0)
      endif

#ifdef PARALLEL
      !get rank of this process and the number of processes
      call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, NrProc, ierr)

      !set Theta and R points
      NrTheta = irsi_frsizey
      NrR = irsi_frsizex

#ifdef MIRROR_THETA
      NrTheta = (NrTheta+1)/2
#endif

      !decompose frequency domain
      MyNrFreq = freq_nr/NrProc
      if((MyRank+1).LE.(Mod(freq_nr, NrProc))) then
         MyNrFreq = MyNrFreq + 1
         FreqBegin = MyNrFreq*MyRank + 1
         FreqEnd = FreqBegin + MyNrFreq - 1
      else
         FreqEnd = freq_nr - ((NrProc-MyRank-1)*MyNrFreq)
         FreqBegin = FreqEnd - MyNrFreq + 1
      end if
#endif

c
c     Now we start with the actual Lambda iteration procedure.
c     
      do iter=start,itermax
c
c         Reset warning flags.
c
#         ifdef INCLUDE_LINES
          line_warn_opacity_neg = 0
#         endif
#         ifdef INCLUDE_DUST
          dust_warn_zero_temp = 0
#         endif
          isynch_warn_lowtemp  = 0
          isynch_warn_lownunub = 0
c     
          write(*,*) 'Iteration ',iter,'.....'
          write(*,*) '------------------------'
c
c         If requested, save LTE and Scattering source functions to file
c
          icount_source = icount_source + 1
          if((isave_source.gt.0).and.
     %         (icount_source.ge.isave_source)) then
              write(*,*) '   (Writing sources to file)'
              call write_source(iter,lastsave_source+1)
              icount_source = 0
          endif
c
c         Increase the Ng counter
c
          icnt_ng = icnt_ng + 1
c     
c         Do a loop over frequency for the radiative transfer
c         
          write(*,*) '   (Perform formal radiative transfer)'

#ifdef PARALLEL
          do ifreq=FreqBegin, FreqEnd
#else
          do ifreq=1, freq_nr
#endif
#             ifndef RADGRID_ONEDIM
#             ifdef INCLUDE_LINES
              if(line_ilifreq(ifreq).eq.1) then
                  write(*,*) ' Line ',line_iline(ifreq)
              endif
#             endif
#             endif
c     
c             Call the radiative transfer
c              
              call do_trans(ifreq)
c     
c             Store this data in compact form into arrays
c             This can be in the form of a mean intensity, or
c             in the form of a lower-angular-resolution 
c             radiation field. This is decided in the scattering
c             module
c     
              call intens_prep_store(ifreq)
c     
c             If requested by the user, and if in the last
c             iteration cycle, and at the right frequency
c             then save the entire radiation field to 
c             a file
c     
              if(((iter.eq.itermax).and.
     %           (ifreq.eq.isave_intens_inu)).or.
     %           (ifreq.eq.-isave_intens_inu)) then
                  write(*,*) '   (Writing Intensity to file for ',
     %                       'ifreq=',ifreq,')'
                  call write_intens()
              endif
c
c             Compute the luminosity at outer edge
c     
#             ifndef JSRC_NO_ARRAYS
              if(isave_fluxcons.eq.0) then 
                  call outer_luminosity_freq(ifreq)
              else
                 write(*,*) 
     %           'compute_lumfunc_freq disabled, Klaus 22/03/07'
              endif
#             endif
c
          enddo

#ifdef PARALLEL
          !broadcast and receive data
          do p=0, NrProc-1
             if(p.EQ.MyRank) then
                do i=FreqBegin, FreqEnd              
                   !put mean intensity data in send-receive array
                   do j=1, NrTheta
                      do k=1, NrR
                         SendRecvArray((j-1)*NrR+k) = intmom_0(i, j, k)
                      end do
                   end do

                   !broadcast data
                   call MPI_Bcast(SendRecvArray(1), NrTheta*NrR, 
     %                  MPI_DOUBLE_PRECISION, p, MPI_COMM_WORLD, ierr)
                end do
                    
             else
                !determine how much and which frequencies will be send
                SendNrFr = freq_nr/NrProc
                if((p+1).LE.(Mod(freq_nr, NrProc))) then
                   SendNrFr = SendNrFr + 1
                   SendFrBegin = SendNrFr*p + 1
                   SendFrEnd = SendFrBegin + SendNrFr - 1
                else
                   SendFrEnd = freq_nr - ((NrProc-p-1)*SendNrFr)
                   SendFrBegin = SendFrEnd - SendNrFr + 1
                end if

                !receive data
                do i=SendFrBegin, SendFrEnd
                   call MPI_Bcast(SendRecvArray(1), NrR*NrTheta,
     %                  MPI_DOUBLE_PRECISION, p, MPI_COMM_WORLD, ierr)

                   !put received data into mean intensity array
                   do j=1, NrTheta
                      do k=1, NrR
                         intmom_0(i, j, k) = SendRecvArray((j-1)*NrR+k)
                      end do
                   end do
                end do
             end if

          end do
#endif

c     
c         >>> ORIGINAL POSITION OF WRITE_INTENS_PREP() <<< (10apr99)
c     
c         Now that the radiation field is known, we can compute
c         the scattering integrals, the dust temperatures and/or
c         the level populations.
c     
          write(*,*) '   (Compute local non-LTE and/or scattering)'
          call solve_statequil_and_scat(0)
c     
c         Find the convergence of the radiation field
c
          cnv_src_error = conv_error()
c
c         Now, if appropriate, solve the temperature. Note that this
c         is not the temperature solving of the dust grains. In dust
c         continuum emission, the dust temperature is so fundamentally
c         coupled to the transfer that it is considered a source term.
c         It is solved in the dust.F module. But here, we solve for the
c         gas/electron temperature, for the case of Comptonization for
c         instance. 
c
          if(itertemp.gt.0) then
              icnt_temp = icnt_temp + 1
              if((icnt_temp.ge.itertemp).or.
     %               (cnv_src_error.lt.convcrit)) then
c                  write(*,*) '   (Solving the temperature)'
c                  call solve_the_temperature()
c                  ididsolvetemp = 1
c                  icnt_temp = 0
                 write(*,*) 
     %             'Solving temperature disabled, Klaus 22/03/07'
                 stop 13
              else
                  ididsolvetemp = 0
              endif
          endif
c
c         If requested by the user, save the moments of the 
c         radiation field, the possible molecular level populations,
c         and the possible dust/electron temperatures to a file
c         
          call write_iteration_sub_result(1,iter)
c
c         Write a warning message if appropriate
c
#         ifdef INCLUDE_DUST
          if(dust_warn_zero_temp.gt.0) then
              write(*,*) '  Warning: detected zero T_dust at ',
     %             dust_warn_zero_temp,' points'
          endif
#         endif
c     
c         Compare the current source with the previous one, in order
c         to find out if we have converged
c
          if(itertemp.eq.0) then
              write(*,*) '   (Error =',cnv_src_error,')'
              if(cnv_src_error.lt.convcrit) then
                  converged = .true.
                  goto 30
              endif
          else
              if(ididsolvetemp.eq.1) then
                  cnv_temp_error = tempsolve_error()
              endif
              write(*,*) '   (Error =',cnv_src_error,')',
     %         '   [ Error T =',cnv_temp_error,']'
              if((cnv_src_error.lt.convcrit).and.
     %           (cnv_temp_error.lt.convcrtemp)) then
                  converged = .true.
                  goto 30
              endif
          endif
c
c         Using the current scattering source and the memorized
c         scattering sources of previous iteration steps, we can
c         (each 4 steps) make an extrapolation of the convergence
c         process, by using Ng-accelleration.
c
          if(ng.ne.0) then
              if(iradproc_line.ne.0) then
#                 ifdef LINE_NG_ACCEL
c
c                 For line transfer, do the Ng on the Jbar
c
                  if(icnt_ng.eq.4) then
                      write(*,*) '   Do the Ng accelleration'
                      call line_ng_accel_jbar()
                      icnt_ng = 0
                  endif
c    
c                 Now memorize the source functions for Ng-accelleration 
c                 later. This is because Ng accelleration uses the memory
c                 of the previous steps to make its prediction
c
                  call line_jbar_copy(2,3)
                  call line_jbar_copy(1,2)
                  call line_jbar_copy(-1,1)
#                 else
                  write(*,*) 'WARNING: Ng for lines de-activated!'
#                 endif
              elseif(iradproc_dust.ne.0) then
#                 ifdef INCLUDE_DUST 
#                 ifndef NO_DUST_NG
c
c                 For dust transfer: do the Ng on T_dust^4. Here we
c                 recompute the emission, but we ignore the fact that
c                 we should recompute the opacity as well, if the opacities
c                 are temperature dependent. This will solve itself.
c
                  if(icnt_ng.eq.4) then
                      write(*,*) '   Do the Ng accelleration'
                      call dust_ng_accellerate()
                      call scat_src_clear()
                      call do_dust_iso_scat()
                      call calc_emis_dust()
                      icnt_ng = 0
                  endif
c    
c                 Now memorize the source functions for Ng-accelleration 
c                 later. This is because Ng accelleration uses the memory
c                 of the previous steps to make its prediction
c
                  call dusttemp_copy(2,3)
                  call dusttemp_copy(1,2)
                  call dusttemp_copy(-1,1)
#                 else
                  write(*,*) 'WARNING: Ng for dust de-activated!'
#                 endif
#                 endif
              else
c
c                 For non-line and non-dust transfer, do the 
c                 Ng on the emissivity
c
                  if(icnt_ng.eq.4) then
                      write(*,*) '   Do the Ng accelleration'
                      call ng_accellerate()
                      icnt_ng = 0
                  endif
c    
c                 Now memorize the source functions for Ng-accelleration 
c                 later. This is because Ng accelleration uses the memory
c                 of the previous steps to make its prediction
c    
c                 NOTE: The implementation here is quite blunt. Should be
c                       done more elegantly
c    
                  call scat_copy(2,3)
                  call scat_copy(1,2)
                  call scat_copy(-1,1)
              endif
          endif
c     
c         Now make sure that the scattering emissivity is
c         now activated for the rest of the cycles
c         
          ssw_scat = 1
c
      enddo
      converged = .false.
   30 continue
      iterdone  = min(iter,itermax)
c
c     Write the final mean intensity to the file 'intmean_final.dat'
c
#     ifndef JSRC_NO_ARRAYS
      if(isave_meanint.ne.0) then
          write(*,*) '   (Writing last intensity moments to file)'
          call write_mean_intensity(iterdone,-2)
      endif
#     endif
c
c     Write the final line Jbar (still odfashioned way)
c
#     ifdef INCLUDE_LINES
      if(iradproc_line.eq.1) then
          write(*,*) '   (Writing last Jbar to file)'
          call write_line_jbar(-3)
      endif
#     endif
c
c     Write the final dust temperature 
c
#     ifdef INCLUDE_DUST
#     ifndef EXCLUDE_DUST_ITERATION
      if((isave_dusttemp.ne.0).and.(iradproc_dust.eq.1)) then
          write(*,*) '   (Writing last dust temperatures to file)'
          call write_dust_temp(iterdone,-2)
      endif
#     endif
#     endif
c
c     Write the final source function
c
      if(isave_source.ne.0) then
          write(*,*) '   (Writing last source to file)'
          call write_source(iterdone,-2)
      endif
c
      end



c     --------------------------------------------------------------
c                      ACCELLERATED LAMBDA ITERATION
c
c     This is the ALI method with Olsen & Kunasz (1987) diagonal
c     approximate operator. In the present implementation this is
c     still not the exact diagonal operator, but to a high degree 
c     it is. For an excellent review on ALI methods, see Hubeny, 
c     1992, "Accellerated Lambda Iteration", in "The Atmospheres 
c     of Early-Type Stars", Lecture Notes in Physics 401, Springer,
c     pp. 377-392, eds. Heber and Jeffery
c
c     The ALI method, as implemented here, works only for isotropic
c     scattering, and for now also only for conservative scattering. 
c
c     REMARK: ng.eq.1 means apply Ng accelleration.
c     REMARK: The Ng accelleration works on each frequency separately.
c             This might be incorrect for line transfer. 
c
c     ARGUMENTS:
c        start       The iteration count at which to (re)start
c        itermax    The maximum number of iterations
c        ng          Switch for Ng acceleration
c        init        The type of initial guess
c        convcrit    The convergence criterion for the sourcefunc
c        itertemp    =0 No iteration on temperature (default)
c                    >0 Maximum number of iterations of RT before
c                       solving temperature. =1 -> every iteration.
c        convcrtemp  The convergence criterion for the temperature
c
c     --------------------------------------------------------------
      subroutine accel_lambda_iteration(start,itermax,ng,init,convcrit,
     %                   itertemp,convcrtemp,iterdone,converged)
      implicit none
c
      integer start,itermax,ng,init,itertemp,iterdone
      doubleprecision convcrit,convcrtemp
      doubleprecision conv_error
      logical converged
c
#include "common_grid.h"
#include "common_source.h"
#include "common_isave.h"
#include "common_lines.h"
#include "common_dust.h"
#include "common_main.h"
c
ccc      doubleprecision converr
      integer iter,ifreq,icnt_ng
c
      if(itertemp.ne.0) then
          write(*,*) 'SORRY: temperature iteration not yet implemented'
          write(*,*) '       in Accelerated Lambda Iteration.'
          stop 13
      endif
c
#     ifdef VERBOSE
      write(*,*) '===================================================='
      write(*,*) '==== Starting the Accellerated Lambda Iteration ===='
      write(*,*) '===================================================='
#     endif
c
c     If line radiation: check level populations, and compute the
c     sources for the first time.
c
#     ifdef INCLUDE_LINES
c     
      if(iradproc_line.eq.1) then
          call check_levelpopul()
c          call write_levelpopul(0)
#         ifndef LINE_VELOCITIES
          call calc_src_from_populations()
#         else
          call global_prepare_transitions()
c
c         NEW: 27.03.06: Removed the next, because it is now done
c                        in the read_linedata() between reading the
c                        line characteristics and setting up the
c                        passbands.
c
ccc          if(iradproc_line_dust.eq.1) then
ccc              call global_prepare_line_dust()
ccc          endif
#         endif
      endif
c
#     endif
c
c     First compute the diagonal of the Lambda operator. This can be
c     done once and for all, as long as the opacities do not change
c     in the course of the iteration. We must set the appropriate
c     flags for the enabling of the alpha_scat and alpha_abs
c
#     ifdef VERBOSE
      write(*,*) 'Computing the ALI accelleration factor'
#     endif
      asw_lte  = 1
      asw_scat = 1
      if(iradproc_line.eq.0) then
          call calc_ali_accell()
          if(isave_alioper.ne.0) then
              write(*,*) '   (Writing accelleration factor to file)'
              call write_ali_diaglamb()
          endif
      endif
c     
c     Enable the emission processes and all opacities, by setting the
c     appropriate flags.
c     
      asw_lte  = 1
      asw_scat = 1
      ssw_lte  = 1
      ssw_scat = 0
c
c     If requested, make an initial guess. There is a variety of types
c     of initial guesses. They are selected by the value of init.
c
      if(init.eq.1) then
c
c         Initial guess = no scattering (for continuum tranfer)
c                       = populations in RE without radiation field (for lines)
c
          write(*,*) 'Make initial guess'
c
c         If lines: then compute the LTE level populations
c
#         ifdef INCLUDE_LINES
c
          if(iradproc_line.eq.1) then
              call calc_re_levelpopul()
          endif
c
#         endif
c
c         Now do the transfer for each freq
c
          if(iradproc_line.eq.0) then
              asw_scat = 0
              ssw_scat = 0
              asw_lte  = 1
              ssw_lte  = 1
              do ifreq=1,freq_nr
c
c                 Perform the transfer in the array intens(iphi,imu,it,ir)
c
                  call do_trans(ifreq)
c
c                 Perform the angular integrals for the mean intensity (J)
c                 or perform the angular regridding to a lower angular 
c                 resolution (if necessary), and then store the result into 
c                 an array with a frequency index. 
c
                  call intens_prep_store(ifreq)
c     
              enddo
              asw_scat = 1
              ssw_scat = 1
c
c             Do a copy of the j_nu source, for ALI
c
              call scat_copy(-1,0)
c    
c             Compute the scattering integrals, from the given J(nu) or 
c             I(inu,iphi,imu) (whichever is appropriate).
c    
              if(iradproc_line.eq.0) then
                  call solve_statequil_and_scat(0)
              endif
          endif
c
      elseif(init.eq.2) then
          if(iradproc_dust.ne.0) then
c
c             [ Dust radiative transfer ]
c
c             Init guess = read dust temperatures from file
c
#             ifdef INCLUDE_DUST
              call read_dust_temperatures(.true.)
              call scat_src_clear()
              call calc_emis_dust()
              ssw_scat = 1
#             else
              stop 676
#             endif
          endif
c
      elseif(init.eq.-1) then
c
c         Initial guess is already made, including the source
c
          ssw_scat = 1
c
      endif
c
c     If line transfer is active, then write the initial populations here
c
#     ifdef INCLUDE_LINES
      if(iradproc_line.eq.1) then
cc          call write_levelpopul(0)
          call write_iteration_sub_result(1,0)
      endif
#     endif
c
c     If Line processes, then switch on ssw_scat!
c
#     ifdef INCLUDE_LINES
      if(iradproc_line.eq.1) then
          ssw_scat = 1
      endif
#     endif
c     
c     Reset the Ng counter. The Ng counter is meant to count to 4
c     and then reset, triggering the Ng accelleration at each 4-th
c     iteration step.
c
      icnt_ng = 0
c
c     For safety, backup the current sources (if e.g. needed for the
c     convergence criterion on sources).
c
      if(iradproc_line.eq.0) then
          call scat_copy(-1,0)
      endif
c
c     Now start the loop for the ALI iteration
c
      do iter=start,itermax
c
c         Reset warning flags
c
#         ifdef INCLUDE_LINES
          line_warn_opacity_neg = 0
#         endif
#         ifdef INCLUDE_DUST
          dust_warn_zero_temp = 0
#         endif
          isynch_warn_lowtemp  = 0
          isynch_warn_lownunub = 0
c     
          write(*,*) 'Iteration ',iter,'.....'
          write(*,*) '------------------------'
c
c         If requested, save LTE and Scattering source functions to file
c
          icount_source = icount_source + 1
          if((isave_source.gt.0).and.
     %         (icount_source.ge.isave_source)) then
              write(*,*) '   (Writing sources to file)'
              call write_source(iter,lastsave_source+1)
              icount_source = 0
          endif
c
c         Increase the Ng counter
c
          icnt_ng = icnt_ng + 1
c
c         For dust ALI, first save the mean intensity and dust temperature
c
#         ifdef DUST_ALI
          if(iradproc_dust.eq.1) then
              call copy_dust_temp_and_meanint()
          endif
#         endif
c     
c         Do a loop over frequency for the radiative transfer
c         
          do ifreq=1,freq_nr
#             ifndef RADGRID_ONEDIM
#             ifdef INCLUDE_LINES
              if(line_ilifreq(ifreq).eq.1) then
                  write(*,*) ' Line ',line_iline(ifreq)
              endif
#             endif
#             endif
c     
c             Call the radiative transfer for the current frequency
c              
              call do_trans(ifreq)
c     
c             Store this data in compact form into arrays
c             This can be in the form of a mean intensity, or
c             in the form of a lower-angular-resolution 
c             radiation field. This is decided in the scattering
c             module
c     
              call intens_prep_store(ifreq)
c     
c             If requested by the user, and if in the last
c             iteration cycle, and at the right frequency
c             then save the entire radiation field to 
c             a file
c     
              if(((iter.eq.itermax).and.
     %           (ifreq.eq.isave_intens_inu)).or.
     %           (ifreq.eq.-isave_intens_inu)) then
                  write(*,*) '   (Writing Intensity to file for ',
     %                       'ifreq=',ifreq,')'
                  call write_intens()
              endif
c
c             Compute the luminosity at outer edge
c     
#             ifndef JSRC_NO_ARRAYS
              if(isave_fluxcons.eq.0) then 
                  call outer_luminosity_freq(ifreq)
              else
                  write(*,*) 
     %         'compute_lumfunc_freq disabled, Klaus 22/03/07'
              endif
#             endif
c
          enddo
c     
c         >>> ORIGINAL POSITION OF WRITE_INTENS_PREP() <<< (10apr99)
c
c         Copy source function for ALI. Note that the conv_error
c         function may also make this copy, in order to find whether
c         we have converged. But this is no problem.
c
          if(iradproc_line.eq.0) then
#             ifdef VERBOSE
              write(*,*) '   Copy source function j'
#             endif
              call scat_copy(-1,0)
          endif
c     
c         Now that the radiation field is known, we can compute
c         the scattering source of the formal solution, i.e.
c            j^FS = alpha_s*Lambda[j^n/alpha]
c     
          if((iradproc_line.eq.0).and.(iradproc_dust.eq.0)) then
c
c             If not line and not dust, then it is either simple
c             scattering, or a process for which ALI is not yet
c             implemented.
c
              write(*,*) '   (Compute local non-LTE and/or scattering)'
              call solve_statequil_and_scat(0)
c
          else
c
c             For line and dust transfer, do the whole process of 
c             computing the level populations from a statistical 
c             equilibrium calculation and compute the source terms. 
c             Give as argument a 1, in order to tell this thing to 
c             use MALI to solve the Non-LTE lines or the DALI for
c             lines.
c
              write(*,*) '   (Compute ALI non-LTE and/or scattering)'
              call solve_statequil_and_scat(1)
c
          endif
c
c         If requested by the user, save the moments of the 
c         radiation field, the possible molecular level populations,
c         and the possible dust temperatures to a file
c         
          call write_iteration_sub_result(1,iter)
c
c         Next call the ALI accellerator. This routine computes
c         the diagonal of the Lambda operator, and uses the 
c         ALI formula to improve the scattering source term so
c         that it converges quicker. 
c
c         NOTE: This ALI accelleration acts on a single frequency
c               only. For line transfer this might be not the correct
c               procedure. Keep this in mind.
c
          if((iradproc_line.eq.0).and.(iradproc_dust.eq.0)) then
#             ifdef VERBOSE
              write(*,*) '   (do the approximate operator acceleration)'
#             endif
              call ali_accelerate()
          endif
c     
c         Compare the current source with the previous one, in order
c         to find out if we have converged
c
          cnv_src_error = conv_error()
          write(*,*) '   (Error =',cnv_src_error,')'
          if(cnv_src_error.lt.convcrit) then
              converged = .true.
              goto 30
          endif
c
c         Using the current scattering source and the memorized
c         scattering sources of previous iteration steps, we can
c         (each 4 steps) make an extrapolation of the convergence
c         process, by using Ng-accelleration.
c
          if(ng.ne.0) then
              if(iradproc_line.ne.0) then
#                 ifdef LINE_NG_ACCEL
c
c                 For line transfer, do the Ng on the Jbar
c
                  if(icnt_ng.eq.4) then
                      write(*,*) '   Do the Ng accelleration'
                      call line_ng_accel_jbar()
                      icnt_ng = 0
                  endif
c    
c                 Now memorize the source functions for Ng-accelleration 
c                 later. This is because Ng accelleration uses the memory
c                 of the previous steps to make its prediction
c
                  call line_jbar_copy(2,3)
                  call line_jbar_copy(1,2)
                  call line_jbar_copy(-1,1)
#                 else
                  write(*,*) 'WARNING: Ng for lines de-activated!'
#                 endif
              elseif(iradproc_dust.ne.0) then
#                 ifdef INCLUDE_DUST 
#                 ifndef NO_DUST_NG
c
c                 For dust transfer: do the Ng on T_dust^4. Here we
c                 recompute the emission, but we ignore the fact that
c                 we should recompute the opacity as well, if the opacities
c                 are temperature dependent. This will solve itself.
c
                  if(icnt_ng.eq.4) then
                      write(*,*) '   Do the Ng accelleration'
                      call dust_ng_accellerate()
                      call scat_src_clear()
                      call do_dust_iso_scat()
                      call calc_emis_dust()
                      icnt_ng = 0
                  endif
c    
c                 Now memorize the source functions for Ng-accelleration 
c                 later. This is because Ng accelleration uses the memory
c                 of the previous steps to make its prediction
c
                  call dusttemp_copy(2,3)
                  call dusttemp_copy(1,2)
                  call dusttemp_copy(-1,1)
#                 else
                  write(*,*) 'WARNING: Ng for dust de-activated!'
#                 endif
#                 endif
              else
c
c                 For non-line and non-dust transfer, do the 
c                 Ng on the emissivity
c
                  if(icnt_ng.eq.4) then
                      write(*,*) '   Do the Ng accelleration'
                      call ng_accellerate()
                      icnt_ng = 0
                  endif
c    
c                 Now memorize the source functions for Ng-accelleration 
c                 later. This is because Ng accelleration uses the memory
c                 of the previous steps to make its prediction
c    
c                 NOTE: The implementation here is quite blunt. Should be
c                       done more elegantly
c    
                  call scat_copy(2,3)
                  call scat_copy(1,2)
                  call scat_copy(-1,1)
              endif
          endif
c     
c         Now make sure that the scattering emissivity is
c         now activated for the rest of the cycles
c         
          ssw_scat = 1
c     
      enddo
      converged = .false.
   30 continue
      iterdone = min(iter,itermax)
c
c     Write the final mean intensity to the file 'intmean_final.dat'
c
#     ifndef JSRC_NO_ARRAYS
      if(isave_meanint.ne.0) then
          write(*,*) '   (Writing last intensity moments to file)'
          call write_mean_intensity(iterdone,-2)
      endif
#     endif
c
c     Write the final line Jbar (still odfashioned way)
c
#     ifdef INCLUDE_LINES
      if(iradproc_line.eq.1) then
          write(*,*) '   (Writing last Jbar to file)'
          call write_line_jbar(-3)
      endif
#     endif
c
c     Write the final dust temperature 
c
#     ifdef INCLUDE_DUST
#     ifndef EXCLUDE_DUST_ITERATION
      if((isave_dusttemp.ne.0).and.(iradproc_dust.eq.1)) then
          write(*,*) '   (Writing last dust temperatures to file)'
          call write_dust_temp(iterdone,-2)
      endif
#     endif
#     endif
c
c     Write the final source function
c
      if(isave_source.ne.0) then
          write(*,*) '   (Writing last source to file)'
          call write_source(iterdone,-2)
      endif
c
      end




c     --------------------------------------------------------------
c                      ALI ACCELLERATION ALGORITHM
c     
c     This routine takes as input the current scattering emissivity
c     and the previous one, and uses the diagonal of the Lambda
c     operator to improve the current scattering emissivity. 
c     --------------------------------------------------------------
      subroutine ali_accelerate()
      implicit none
c
#ifndef JSRC_NO_ARRAYS
c
#include "common_source.h"
#include "common_grid.h"
#include "common_radfield.h"
c
      integer inu,ir,it,itmax
c
c     Next do some stuff 
c
      itmax   = irsi_frsizey
c
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#     endif
#     ifdef RADGRID_ONEDIM
      itmax = 1
#     endif
c
c     Loop over all space and frequency
c     
      do ir=1,irsi_frsizex
      do it=1,itmax
      do inu=1,freq_nr
c
c         Check if the diagonal Lambda operator array is non-zero here
c
          if(ali_accell(inu,it,ir).eq.0.d0) then
              write(*,*) 'ALI accelleration factor = 0 at',inu,it,ir
              stop 13
          endif
c
c         Update the real scattering source, by applying the
c         formula:
c     
c                                       j_s^FS - j_s^n
c          j_s^n+1 = j_s^n + -------------------------------------
c                            1 - (alp_s/(alp_a+alp_s))*Lambda_diag
c     
          scati_src(inu,it,ir) = scatibk_src(inu,it,ir,0)
     %      + ( scati_src(inu,it,ir) - scatibk_src(inu,it,ir,0) )
     %      * ali_accell(inu,it,ir)
c
      enddo
      enddo
      enddo      
c
#else
      write(*,*) 'Routine ali_accelerate() only when source'
      write(*,*) 'and mean int arrays are present.'
      stop 13
#endif
c
      end


c     --------------------------------------------------------------
c       COMPUTE THE DIAGONAL OF LAMBDA AND THE ACCELLERATION FACTOR
c
c     This routine computes the accelleration factor for the ALI
c     accelleration at each point on the grid at each frequency.
c     This routine needs to be called only once, before the ALI 
c     iteration is stared, unless the opacities change per iteration,
c     which in general they don't.
c     --------------------------------------------------------------
      subroutine calc_ali_accell()
      implicit none
c
#     ifndef NO_SHORT_CHARS
c
      doubleprecision find_alp_t
c
#ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
#include "common_shortchar.h"
c
      integer ifreq,ir,it,itmax,imu,iphi,iphimin,iphimax,iphioff
      integer idrrp,idrrn,icr,iangset
      doubleprecision phiweight,diaglamb,dummy,tau,r0,rscale
      doubleprecision alp0,alp1,alp2,dr,dt,mu,phi
      doubleprecision tau1,tau2,tsum,e0,e1,e2,xp
c
c     Check whether the ESC exist
c
      if(isc_ready.ne.1) then
          write(*,*) 'ERROR: ALI called while ESC non-existent.'
          stop 13
      endif      
c
c     Check whether both scattering and emission are switched on
c
      if((asw_lte.ne.1).or.(asw_scat.ne.1)) then
          write(*,*) 'Trying to compute ALI accelleration factor'
          write(*,*) 'while NOT both scattering and absorption are'
          write(*,*) 'switched on... '
          write(*,*) 'asw_lte = ',asw_lte,' asw_scat = ',asw_scat
          write(*,*) 'Aborting...'
          stop 13
      endif
c
c     For the moment we only allow one global angular set
c
      idrrp   = 1
      idrrn   = 1
      iangset = 1
c
c     Next do some stuff 
c
      itmax   = irsi_frsizey
      iphimin = 1
      iphimax = nrphi(iangset)
      iphioff = 0
c
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#     endif
#     ifdef MIRROR_PHI
      iphimin = nrphi(iangset) / 4 + 1
      iphimax = 3 * nrphi(iangset) / 4 
      iphioff = - nrphi(iangset) / 4
#     endif
#     ifdef RADGRID_ONEDIM
      iphimin = 1
      iphimax = 1
      iphioff = 0
      itmax = 1
#     endif
      phiweight = 1.d0 / (1.d0*(iphimax-iphimin+1))
c
c     Loop over all space and frequencies
c     
      do ifreq=1,freq_nr
c
c         First set accel factor to 1 (Lambda*=0) at ir=1 and ir=frsizex
c
          do it=1,itmax
              ali_accell(ifreq,it,1) = 1.d0
              ali_accell(ifreq,it,irsi_frsizex) = 1.d0
          enddo
c
c         Now do the real thing
c
          do ir=1+1,irsi_frsizex-1
c
c             Find the R, and the scale factor by which the dimensionless
c             sc_ds should be multiplied in order to find the real ds 
c             along the elements of the extended short characteristics.
c         
c                Spherical coordinates:   Use R itself as a scaling factor
c                Cartesian coordinates:   Use R_{i+1}-R_i as a scaling factor
c
#             ifdef COORD_SPHERICAL
              r0     = rsi_x_c(ir,1)
              rscale = r0
#             else
              r0     = rsi_x_c(ir,1)
              rscale = rsi_x_c(ir+1,1) - rsi_x_c(ir,1)
#             endif
c
c             In the new style of SC we may have varying dr/r, so we
c             retrieve the index of the SC-set appropriate for this
c             radius from an array
c             
#             ifdef SHORTCHAR_NEWSTYLE
              idrrp = ir_iscset_p(ir)
              idrrn = ir_iscset_n(ir)
#             endif
c
c             Now the loop over Theta
c
              do it=1,itmax
c
c                 Find the relevant opacity terms
c
                  alp0 = scati_alpha(ifreq,it,ir) + 
     %                    slte_alpha(ifreq,it,ir) 
c
c                 Find the diagonal of the Lambda Operator
c
#ifndef ALI_PNT_3
c
c                 ( The case of 1st order integration of source )
c
                  diaglamb = 0.d0
                  dummy    = 0.d0
c
c                 If there is a mu=0 gridpoint, then do this one
c                 first
c     
                  if(anggrid_mu_zero.eq.1) then
                      imu      = 0
                      do iphi=iphimin,iphimax
                          icr   = isc_icross(1,iphi,imu,it,idrrp)
                          dr    = sc_dr(1,iphi,imu,it,idrrp)
                          dt    = sc_dtheta(1,iphi,imu,it,idrrp)
                          mu    = sc_mu(1,iphi,imu,it,idrrp)
                          phi   = sc_phi(1,iphi,imu,it,idrrp)
                          alp1  = find_alp_t(ifreq,icr,ir,it,
     %                        dr,dt,mu,phi)
                          tau   = sc_ds(1,iphi,imu,it,idrrp) * rscale * 
     %                        0.5d0 * ( alp0 + alp1 )
c     MAJOR BUG FIXED: CATCH SMALL TAU   (07-06-01)
                          if(tau.gt.TAU_EPS) then
                              dummy = dummy + (tau-1.d0+exp(-tau))/tau
                          elseif(tau.gt.0.d0) then
                              dummy = dummy + 0.5d0 * tau
                          endif
                      enddo
                      dummy = dummy * phiweight
                      diaglamb = diaglamb + dummy * 
     %                    (rmu_i(1,iangset)-rmu_i(-1,iangset))
                  endif
c
c                 Next do the other mu gridpoints
c
                  do imu=1,nrmu(iangset)
                      dummy = 0.d0
                      do iphi=iphimin,iphimax
                          icr   = isc_icross(1,iphi,imu,it,idrrp)
                          dr    = sc_dr(1,iphi,imu,it,idrrp)
                          dt    = sc_dtheta(1,iphi,imu,it,idrrp)
                          mu    = sc_mu(1,iphi,imu,it,idrrp)
                          phi   = sc_phi(1,iphi,imu,it,idrrp)
                          alp1  = find_alp_t(ifreq,icr,ir,it,
     %                              dr,dt,mu,phi)
                          tau   = sc_ds(1,iphi,imu,it,idrrp) * rscale * 
     %                            0.5d0 * ( alp0 + alp1 )
c     MAJOR BUG FIXED: CATCH SMALL TAU   (07-06-01)
                          if(tau.gt.TAU_EPS) then
                              dummy = dummy + (tau-1.d0+exp(-tau))/tau
                          elseif(tau.gt.0.d0) then
                              dummy = dummy + 0.5d0 * tau
                          endif
                          icr   = isc_icross(1,iphi,-imu,it,idrrn)
                          dr    = sc_dr(1,iphi,-imu,it,idrrn)
                          dt    = sc_dtheta(1,iphi,-imu,it,idrrn)
                          mu    = sc_mu(1,iphi,-imu,it,idrrn)
                          phi   = sc_phi(1,iphi,-imu,it,idrrn)
                          alp1  = find_alp_t(ifreq,icr,ir,it,
     %                              dr,dt,mu,phi)
c     MAJOR BUG FIXED: IMU --> -IMU     (05-04-99)
                          tau   = sc_ds(1,iphi,-imu,it,idrrn) * rscale * 
     %                            0.5d0 * ( alp0 + alp1 )
c     MAJOR BUG FIXED: CATCH SMALL TAU   (07-06-01)
                          if(tau.gt.TAU_EPS) then
                              dummy = dummy + (tau-1.d0+exp(-tau))/tau
                          elseif(tau.gt.0.d0) then
                              dummy = dummy + 0.5d0 * tau
                          endif
                      enddo
                      dummy = dummy * phiweight
                      diaglamb = diaglamb + dummy * 
     %                        (rmu_i(imu+1,iangset)-rmu_i(imu,iangset))
                  enddo
                  diaglamb = 0.5 * diaglamb
#else
c
c                 ( The case of 2nd order integration of source )
c                 ( This is an approximation to the diagonal, not the
c                   entire diagonal. But the difference with the real
c                   diagonal is not very large, since these extra terms
c                   only involve the contributions of the 2nd order 
c                   integration, and they are multiplied by exp(-tau1)
c                   so that they vanish in the optically thick limit. )
c
                  diaglamb = 0.d0
                  dummy    = 0.d0
c
c                 If there is a mu=0 gridpoint, then do this one
c                 first
c     
                  if(anggrid_mu_zero.eq.1) then
                      imu      = 0
                      do iphi=iphimin,iphimax
                          icr   = isc_icross(1,iphi,imu,it,idrrp)
                          dr    = sc_dr(1,iphi,imu,it,idrrp)
                          dt    = sc_dtheta(1,iphi,imu,it,idrrp)
                          mu    = sc_mu(1,iphi,imu,it,idrrp)
                          phi   = sc_phi(1,iphi,imu,it,idrrp)
                          alp1  = find_alp_t(ifreq,icr,ir,it,
     %                        dr,dt,mu,phi)
                          tau1   = sc_ds(1,iphi,imu,it,idrrp) * rscale * 
     %                        0.5d0 * ( alp0 + alp1 )
                          icr   = isc_icross(0,iphi,imu,it,idrrp)
                          dr    = sc_dr(0,iphi,imu,it,idrrp)
                          dt    = sc_dtheta(0,iphi,imu,it,idrrp)
                          mu    = sc_mu(0,iphi,imu,it,idrrp)
                          phi   = sc_phi(0,iphi,imu,it,idrrp)
                          alp2  = find_alp_t(ifreq,icr,ir,it,
     %                        dr,dt,mu,phi)
                          tau2   = sc_ds(0,iphi,imu,it,idrrp) * rscale * 
     %                        0.5d0 * ( alp0 + alp2 )
                          tsum   = tau1 + tau2
                          if(tau1.gt.1.d-6) then
                              xp        = exp(-tau1)
                              e0        = 1.d0 - xp
                              e1        = tau1 - e0
                              e2        = tau1 * tau1 - 2.d0 * e1
                              if(tau2.gt.1.d-6) then
                                  dummy  = dummy + ( tsum * e1 - e2 ) /
     %                                ( tau1 * tau2 )
                              else  
                                  dummy  = dummy + e1 / tau1
                              endif
                          endif
                      enddo
                      dummy = dummy * phiweight
                      diaglamb = diaglamb + dummy * 
     %                    (rmu_i(1,iangset)-rmu_i(-1,iangset))
                  endif
                  do imu=1,nrmu(iangset)
                      dummy = 0.d0
                      do iphi=iphimin,iphimax
                          icr   = isc_icross(1,iphi,imu,it,idrrp)
                          dr    = sc_dr(1,iphi,imu,it,idrrp)
                          dt    = sc_dtheta(1,iphi,imu,it,idrrp)
                          mu    = sc_mu(1,iphi,imu,it,idrrp)
                          phi   = sc_phi(1,iphi,imu,it,idrrp)
                          alp1  = find_alp_t(ifreq,icr,ir,it,
     %                              dr,dt,mu,phi)
                          tau1   = sc_ds(1,iphi,imu,it,idrrp) * rscale * 
     %                             0.5d0 * ( alp0 + alp1 )
                          icr   = isc_icross(0,iphi,imu,it,idrrp)
                          dr    = sc_dr(0,iphi,imu,it,idrrp)
                          dt    = sc_dtheta(0,iphi,imu,it,idrrp)
                          mu    = sc_mu(0,iphi,imu,it,idrrp)
                          phi   = sc_phi(0,iphi,imu,it,idrrp)
                          alp2  = find_alp_t(ifreq,icr,ir,it,
     %                              dr,dt,mu,phi)
                          tau2  = sc_ds(0,iphi,imu,it,idrrp) * rscale * 
     %                            0.5d0 * ( alp0 + alp2 )
                          tsum  = tau1 + tau2
                          if(tau1.gt.1.d-6) then
                              xp        = exp(-tau1)
                              e0        = 1.d0 - xp
                              e1        = tau1 - e0
                              e2        = tau1 * tau1 - 2.d0 * e1
                              if(tau2.gt.1.d-6) then
                                  dummy  = dummy + 
     %                           ( tsum * e1 - e2 ) / ( tau1 * tau2 )
                              else  
                                  dummy  = dummy + e1 / tau1
                              endif
                          endif
c     MAJOR BUGS FIXED: IMU --> -IMU     (05-04-99)
                          icr   = isc_icross(1,iphi,-imu,it,idrrn)
                          dr    = sc_dr(1,iphi,-imu,it,idrrn)
                          dt    = sc_dtheta(1,iphi,-imu,it,idrrn)
                          mu    = sc_mu(1,iphi,-imu,it,idrrn)
                          phi   = sc_phi(1,iphi,-imu,it,idrrn)
                          alp1  = find_alp_t(ifreq,icr,ir,it,
     %                              dr,dt,mu,phi)
                          tau1  = sc_ds(1,iphi,-imu,it,idrrn) * rscale * 
     %                             0.5d0 * ( alp0 + alp1 )
                          icr   = isc_icross(0,iphi,-imu,it,idrrn)
                          dr    = sc_dr(0,iphi,-imu,it,idrrn)
                          dt    = sc_dtheta(0,iphi,-imu,it,idrrn)
                          mu    = sc_mu(0,iphi,-imu,it,idrrn)
                          phi   = sc_phi(0,iphi,-imu,it,idrrn)
                          alp2  = find_alp_t(ifreq,icr,ir,it,
     %                              dr,dt,mu,phi)
                          tau2  = sc_ds(0,iphi,-imu,it,idrrn) * rscale * 
     %                            0.5d0 * ( alp0 + alp2 )
                          tsum  = tau1 + tau2
                          if(tau1.gt.1.d-6) then
                              xp        = exp(-tau1)
                              e0        = 1.d0 - xp
                              e1        = tau1 - e0
                              e2        = tau1 * tau1 - 2.d0 * e1
                              if(tau2.gt.1.d-6) then
                                  dummy  = dummy + 
     %                          ( tsum * e1 - e2 ) / ( tau1 * tau2 )
                              else  
                                  dummy  = dummy + e1 / tau1
                              endif
                          endif
                      enddo
                      dummy = dummy * phiweight
                      diaglamb = diaglamb + dummy * 
     %                        (rmu_i(imu+1,iangset)-rmu_i(imu,iangset))
                  enddo
                  diaglamb = 0.5 * diaglamb
#endif
c
c                 Store the result (for line transport and dust continuum)
c
                  ali_approx_lambda(ifreq,it,ir) = diaglamb
c    
c                 Now compute the accelleration factor (only for
c                 scattering)
c
c                 NOTE: if alp_s/alpha does not chance too drastic,
c                        then this should do.
c    
                  if(alp0.ne.0.d0) then
                      dummy = scati_alpha(ifreq,it,ir) / alp0
                  else
                      dummy = 1.d0
                  endif
                  ali_accell(ifreq,it,ir) = 1.d0 /
     %                        ( 1.d0 - dummy*diaglamb )
              enddo
          enddo
      enddo
c
#else
      write(*,*) 'Routine calc_ali_accell() only when source'
      write(*,*) 'and mean int arrays are present.'
      stop 13
#endif
c
#     else /* ifndef NO_SHORT_CHARS */ 
      write(*,*) 'ERROR: calc_ali_accell() has been called'
      write(*,*) '       while the short characteristics have'
      write(*,*) '       been switched off by NO_SHORT_CHARS'
      write(*,*) '       Tip: undefine NO_SHORT_CHARS in'
      write(*,*) '       configure.h'
      stop 13
#     endif /* ifndef NO_SHORT_CHARS */ 
c
      end


c     --------------------------------------------------------------
c              WRITE THE APPROXIMATE LAMBDA OPERATOR TO FILE
c                            (DEBUGGING ONLY)
c     --------------------------------------------------------------
      subroutine write_ali_diaglamb()
      implicit none
c
#ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_source.h"
c
      integer ir,it,imirt,itmax
      integer inu
c
c     Some computations
c
      itmax   = irsi_frsizey
      imirt   = 0
c
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
      imirt = 1
#     endif
#     ifdef RADGRID_ONEDIM
      itmax = 1
      imirt = 0
#     endif
c
      open(unit=1,file='alidiag.dat',status='unknown')
      write(1,*) freq_nr,irsi_frsizex,itmax,imirt
      write(1,*) ' '
      do inu=1,freq_nr
          do ir=1,irsi_frsizex
              do it=1,itmax
                  write(1,*) ali_approx_lambda(inu,it,ir)
              enddo
          enddo
          write(1,*) ' '
      enddo
      close(1)
c
#else
      write(*,*) 'Routine write_ali_diaglamb() only when source'
      write(*,*) 'and mean int arrays are present.'
      stop 13
#endif
c
      end




c     --------------------------------------------------------------
c                           NG ACCELLERATION 
c
c     This version of the Ng acceleration does the acceleration on
c     the entire solution only.
c     --------------------------------------------------------------
      subroutine ng_accellerate()
      implicit none
c
#ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
c
      doubleprecision a1,a2,b1,b2,c1,c2,a,b,c,x0,x1,x2,x3,w
      doubleprecision dum1,dum2,dum3
      integer inu,ir,it,itmax
c
      itmax   = irsi_frsizey
c
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#     endif
#     ifdef RADGRID_ONEDIM
      itmax = 1
#     endif
c
      a1 = 0.d0
      a2 = 0.d0
      b1 = 0.d0
      b2 = 0.d0
      c1 = 0.d0
      c2 = 0.d0
      do inu=1,freq_nr
          do ir=1,irsi_frsizex
              do it=1,itmax
                  x0   = scati_src(inu,it,ir)
                  x1   = scatibk_src(inu,it,ir,1)
                  x2   = scatibk_src(inu,it,ir,2)
                  x3   = scatibk_src(inu,it,ir,3)
                  w    = 1.d0 / intmom_0(inu,it,ir)
                  dum1 = x0 - 2.d0*x1 + x2
                  dum2 = x0 - x1 - x2 + x3
                  dum3 = x0 - x1
                  a1   = a1 + dum1 * dum1 * w
                  b1   = b1 + dum1 * dum2 * w
                  c1   = c1 + dum1 * dum3 * w
                  a2   = b1
                  b2   = b2 + dum2 * dum2 * w
                  c2   = c2 + dum2 * dum3 * w
              enddo
          enddo
      enddo
      a = ( c1*b2 - c2*b1 ) / ( a1*b2 - a2*b1 ) 
      b = ( c2*a1 - c1*a2 ) / ( a1*b2 - a2*b1 ) 
      c = 1.d0 - a - b
#     ifdef VERBOSE
      write(*,*) '       Ng : a = ',a,'  b = ',b
#     endif
c
c     Now perform the accelleration step
c
      do inu=1,freq_nr
          do ir=1,irsi_frsizex
              do it=1,itmax
                  if((scati_src(inu,it,ir).lt.0.d0).or.
     %               (scatibk_src(inu,it,ir,1).lt.0.d0).or.
     %               (scatibk_src(inu,it,ir,2).lt.0.d0)) then
                      write(*,*) '---> Negative source detected at Ng'
                  endif
                  scati_src(inu,it,ir) = c * scati_src(inu,it,ir) 
     %                      + a * scatibk_src(inu,it,ir,1)
     %                      + b * scatibk_src(inu,it,ir,2)
                  if(scati_src(inu,it,ir).lt.0.d0) then
                      write(*,*) 'Ng ERROR: negative source ',ir,it,inu
                      write(*,*) '   src   = ',scati_src(inu,it,ir),
     %                  scatibk_src(inu,it,ir,1), 
     %                  scatibk_src(inu,it,ir,2)
                      write(*,*) '   a,b,c = ',a,b,c
                      scati_src(inu,it,ir) = 0.d0
                  endif
              enddo
          enddo
      enddo
c
#else
      write(*,*) 'Routine ng_accellerate() only when source'
      write(*,*) 'and mean int arrays are present.'
      stop 13
#endif
c
      end



c     --------------------------------------------------------------
c               NG ACCELLERATION FOR EACH FREQ SEPARATELY
c
c     --------------------------------------------------------------
      subroutine ng_accellerate_freq()
      implicit none
c
#ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
c
      doubleprecision a1,a2,b1,b2,c1,c2,a,b,c,x0,x1,x2,x3,w
      doubleprecision dum1,dum2,dum3
      integer inu,ir,it,itmax
c
      itmax   = irsi_frsizey
c
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#     endif
#     ifdef RADGRID_ONEDIM
      itmax = 1
#     endif
c
      do inu=1,freq_nr
          a1 = 0.d0
          a2 = 0.d0
          b1 = 0.d0
          b2 = 0.d0
          c1 = 0.d0
          c2 = 0.d0
          do ir=1,irsi_frsizex
              do it=1,itmax
                  x0   = scati_src(inu,it,ir)
                  x1   = scatibk_src(inu,it,ir,1)
                  x2   = scatibk_src(inu,it,ir,2)
                  x3   = scatibk_src(inu,it,ir,3)
                  w    = 1.d0 / intmom_0(inu,it,ir)
                  dum1 = x0 - 2.d0*x1 + x2
                  dum2 = x0 - x1 - x2 + x3
                  dum3 = x0 - x1
                  a1   = a1 + dum1 * dum1 * w
                  b1   = b1 + dum1 * dum2 * w
                  c1   = c1 + dum1 * dum3 * w
                  a2   = b1
                  b2   = b2 + dum2 * dum2 * w
                  c2   = c2 + dum2 * dum3 * w
              enddo
          enddo
          a = ( c1*b2 - c2*b1 ) / ( a1*b2 - a2*b1 ) 
          b = ( c2*a1 - c1*a2 ) / ( a1*b2 - a2*b1 ) 
          c = 1.d0 - a - b
#         ifdef VERBOSE
          write(*,*) '       Ng : a = ',a,'  b = ',b
#         endif
c
c         Now perform the accelleration step
c
          do ir=1,irsi_frsizex
              do it=1,itmax
                  scati_src(inu,it,ir) = c * scati_src(inu,it,ir) 
     %                      + a * scatibk_src(inu,it,ir,1)
     %                      + b * scatibk_src(inu,it,ir,2)
              enddo
          enddo
      enddo
c
#else
      write(*,*) 'Routine ali_accelerate_freq() only when source'
      write(*,*) 'and mean int arrays are present.'
      stop 13
#endif
c
      end


c     --------------------------------------------------------------
c         COMPUTE THE LEVEL POPULATIONS GIVEN THE SOURCE TERMS
c
c     If the source terms are given, but the rest is not, then in
c     principle the level populations are known. What has to be
c     done is first the Lambda Operation, and then the local
c     statistical equilibrium equation.
c     --------------------------------------------------------------
#     ifdef INCLUDE_LINES
c
      subroutine calc_levelpop_from_source()
      implicit none
c
#include "common_grid.h"
#include "common_source.h"
c
      integer ifreq
c
#     ifdef VERBOSE
      write(*,*) 'Performing SC transfer to find ',
     %    'the levelpopulations from the given '
      write(*,*) 'source terms'
#     endif
c
c     First be sure that the line data is set up properly
c
      call global_prepare_transitions()
c
c     Do a loop over frequency for the radiative transfer
c     This must be done in order to find J
c
      do ifreq=1,freq_nr
          call do_trans(ifreq)
          call intens_prep_store(ifreq)
      enddo
c
c     Now find Jbar and solve the level populations
c
      call calc_level_populations(0)
c
      end
#     endif


c     --------------------------------------------------------------
c                      ITERATIVE SCATTERING METHOD
c
c     This is the iterative scattering method that is handy for 
c     Compton scattering in optically thin media. The intensity
c     is written as a sum of contributions of photons that have
c     scattered 0,1,2,3,4,... etc times, 
c
c       I = I_0 + I_1 + I_2 + I_3 + I_4 + I_5 + .....
c
c     We first compute I_0, then I_1 from I_0, then I_2 from I_1,
c     etc. In this procedure we must memorize only the previous 
c     I_{i-1} to calculate I_{i}. And, in order to find the final
c     answer, we must store also the sum I_0 + ... + I_{i-1}, 
c     so that we can compute the final spectrum. In practise this
c     is not stored as I_0 + ... + I_{i-1} but instead as the sum
c     of the scattering emissivities j_0 + ... + j_{i-1}, which is 
c     after all what we want to know.
c     --------------------------------------------------------------
      subroutine iterative_scattering(itermax,theconvcrit,iterdone,
     %               converged)
      implicit none
c
      integer itermax,iterdone
      doubleprecision theconvcrit
      logical converged
c
      doubleprecision conv_error_cumul_intmean
c
#ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_source.h"
#include "common_isave.h"
#include "common_boundary.h"
#include "common_radfield.h"
#include "common_main.h"
c
      integer iter,ifreq
ccc      doubleprecision converr
      doubleprecision convcrit
c
c     First parse the convcrit
c
      if((theconvcrit.gt.0.d0).and.(theconvcrit.lt.3.d-1)) then
          convcrit = theconvcrit
      else
          convcrit = 0.d0
      endif
c
c     If line transfer, then error
c
      if(iradproc_line.eq.1) then
          write(*,*) 'Iiterative scattering: no line transfer!'
          stop 13
      endif
c
c     Let's set the intmom_iterscat flag, to signal to the rest of
c     the program that the intmom_0 an iterative snapshot is, rather
c     than the total.
c
      intmom_iterscat = 1
c
#     ifdef VERBOSE
      write(*,*) '===================================================='
      write(*,*) '====      Starting the Iterative Scattering     ===='
      write(*,*) '===================================================='
#     endif
c
c     First clear the arrays of the cumulative scattering
c     
      call scat_cumul_clear()
c
#ifdef CUMULATIVE_INTMEAN
      call intmom_cumul_clear()
#endif
c     
c     Then enable the emission processes and all opacities
c     
      asw_lte  = 1
      asw_scat = 1
      ssw_lte  = 1
      ssw_scat = 0
c     
c     If isotropic csk is selected, read the data files
c
c      if(iradproc_compt.eq.3) call read_isocsk()
c
      do iter=1,itermax
c
c         Reset warning flags
c     
          isynch_warn_lowtemp  = 0
          isynch_warn_lownunub = 0
c
          write(*,*) 'Iteration ',iter,'.....'
          write(*,*) '------------------------'
c     
c         If requested, save LTE and Scattering source functions to file
c
          if(isave_source.gt.0) then
              write(*,*) '   (Writing sources to file)'
              call write_source(iter,lastsave_source+1)
          endif
c     
c         Do a loop over frequency for the radiative transfer
c         
          do ifreq=1,freq_nr
#             ifdef VERBOSE
              if(freq_nr.gt.1) then
                  write(*,*) ' Frequency nr ',ifreq
              endif
#             endif
c     
c             Call the radiative transfer
c              
#             ifdef VERBOSE
              write(*,*) '   Performing SC transfer'
#             endif
              call do_trans(ifreq)
c     
c             Store this data in compact form into arrays
c             This can be in the form of a mean intensity, or
c             in the form of a lower-angular-resolution 
c             radiation field. This is decided in the scattering
c             module
c     
#             ifdef VERBOSE
              write(*,*) '   Computing mean intensity J'
#             endif
              call intens_prep_store(ifreq)
c     
c             If requested by the user, and if in the last
c             iteration cycle, and at the right frequency
c             then save the entire radiation field to 
c             a file
c     
              if(((iter.eq.itermax).and.
     %           (ifreq.eq.isave_intens_inu)).or.
     %           (ifreq.eq.-isave_intens_inu)) then
                  write(*,*) '   (Writing Intensity to file)'
                  call write_intens()
              endif
          enddo
c     
c         >>> ORIGINAL POSITION OF WRITE_INTENS_PREP() <<< (10apr99)
c
c         For the sake of user-friendliness, it might be handy 
c         to keep track of the total cumulative mean intensity
c         so that this can be saved to file lateron.
c
#         ifdef CUMULATIVE_INTMEAN
          write(*,*) '   (Updating the cumulative intmean)'
          call intmom_cumul_add()
#         endif
c     
c         Now that the radiation field is known, we can compute
c         the scattering integrals or the level populations.
c     
          call solve_statequil_and_scat(0)
c     
c         We should make a copy of the scattering source terms, 
c         or level populations into a cumulative array.
c     
          call scat_cumul_add()
c
c         If requested by the user, save the moments of the 
c         radiation field to a file
c         
          call write_iteration_sub_result(1,iter)
c     
c         Compare the cumulative source with the added one, in order
c         to find out if we have converged
c
          cnv_src_error = conv_error_cumul_intmean()
          write(*,*) '   (Error =',cnv_src_error,')'
          if(cnv_src_error.lt.convcrit) then 
              converged = .true.
              goto 30
          endif
c
c         Now make sure that the scattering emissivity is
c         now activated for the rest of the cycles, and that the
c         LTE emission processes are switched off, and the 
c	  central source, if any, if switched off as well
c         
          ssw_scat = 1
          ssw_lte  = 0
          do ifreq=1,freq_nr
              radbnd_starspec(ifreq) = 0.d0
          enddo
      enddo
      converged = .false.
   30 continue
      iterdone = min(iter,itermax)
c
c     Now that the Iterative Scattering is finished, we insert the
c     cumulative scattering source/emission into the usual scattering
c     emission. When the telescope is then called, then the result
c     should be the intensity resulting from LTE + scattering! 
c     
      write(*,*) 'Storing the cumulative source into ',
     %   'the scattering source'
      call scat_cumul_to_src()
c
c     If the user wants, we can write the final mean intensity
c     (i.e. cumulative mean intensity) to the file 'intmean_final.dat'
c
      if(isave_meanint.ne.0) then
#     ifdef CUMULATIVE_INTMEAN
          write(*,*) '   (Writing Last Intensity Moments to file)'
          call write_cumul_mean_intensity(iterdone,-2)
#     else
          write(*,*) 'Sorry: if CUMULATIVE_INTMEAN not set, then '
          write(*,*) 'the iterative scattering method can not '
          write(*,*) 'write last cumulative mean intensity...'
#     endif          
      endif
#else
      write(*,*) 'Iterative scattering mode only works when the'
      write(*,*) 'source and mean intensity arrays are present.'
      write(*,*) 'Since the macro JSRC_NO_ARRAYS is set, these'
      write(*,*) 'are not included right now. Aborting...'
      stop 13
#endif
c
      end





c     --------------------------------------------------------------
c               WRITE THE (SUB-) RESULT OF THE ITERATION
c
c     This routine writes everything that the user may want to see
c     as the subresult of the iteration process. These things may
c     include the mean intensity, the dust temperatures, the level
c     populations etc. These data are meant only for diagnostic
c     purposes, not as the actual output of the iteration procedure.
c     For that, please use the routine write_iteration_end_result().
c     --------------------------------------------------------------
      subroutine write_iteration_sub_result(isave,iter)
      implicit none
      integer isave,iter
c
#include "common_source.h"
#include "common_main.h"
#include "common_dust.h"
#include "common_isave.h"
c
c     Write the mean intensity to file
c
#     ifndef JSRC_NO_ARRAYS
c
c     First increase the counter 
c
      icount_meanint = icount_meanint + 1
c
c     Then if required, write mean intensity
c
      if((isave_meanint.gt.0).and.
     %     (icount_meanint.ge.isave_meanint)) then
          write(*,*) '   (Writing Intensity Moments to file)'
          call write_mean_intensity(iter,lastsave_meanint+1)
#         ifdef INCLUDE_LINES
          if(iradproc_line.eq.1) then
              write(*,*) '   (Writing line Jbar and S to file)'
              call write_line_jbar(iter)
          endif
#         endif
          icount_meanint = 0 
      endif
#     endif
c
c     Write electron temperature to file
c
#     ifdef INCLUDE_COMPTON
#     ifdef COMPTON_SOLVE_TE
c
c     First increase the counter
c
      icount_electemp = icount_electemp + 1
c
c     Then if required, write the electron temperature
c
      if((iradproc_solve_temp.ne.0).and.(isave_electemp.gt.0)
     %     .and.(icount_electemp.ge.isave_electemp)) then
          call write_electron_temperature(iter)
          icount_electemp = 0
      endif
#     endif
#     endif
c
c     If lines are present, write the level populations
c
#     ifdef INCLUDE_LINES
c
c     First increase the counter
c
      icount_levpop = icount_levpop + 1
c
c     Then if required, write the level populations
c
      if((iradproc_line.eq.1).and.(isave_levpop.gt.0).and.
     %       (icount_levpop.ge.isave_levpop)) then
          write(*,*) '   (Writing level populations to file)'
          call write_levelpopul(iter)
          call write_lineinfo(iter)
          icount_levpop = 0
      endif
#     endif
c     
c     If dust is present, and if dust temperature is iterated
c     then write the dust temperatures
c
#     ifdef INCLUDE_DUST
#     ifndef EXCLUDE_DUST_ITERATION
c
c     First increase the counter 
c
      icount_dusttemp = icount_dusttemp + 1
c
c     Then, if required, write the dust temperature
c
      if((iradproc_dust.eq.1).and.(isave_dusttemp.gt.0).and.
     %     (icount_dusttemp.ge.isave_dusttemp)) then
          write(*,*) '   (Writing dust temperatures to file)'
          call write_dust_temp(iter,lastsave_dusttemp+1)
          icount_dusttemp = 0
      endif
c
#     endif
#     endif
c
c     Write the data from which the iteration can be restarted, or
c     from which a spectrum can be made. But do this only when the
c     isave_scatnonlte equals 2.
c
      if(isave_scatnonlte.eq.2) then
          if(iradproc_line.ne.0) then
#             ifdef INCLUDE_LINES
              call write_levelpopul(-3)
#             else
              write(*,*) 'LINE TRANSFER DISABLED'
              stop 13
#             endif
          else
#             ifndef JSRC_NO_ARRAYS
              call write_scattering()
#             else
              write(*,*) 'ERROR: ',
     %          'iteration.F/write_iteration_sub_result()'
              write(*,*) '       Attempt to write the nonlte.dat file'
              write(*,*) '       while the src arrays are not compiled!'
              stop 13
#             endif
          endif
#         ifndef JSRC_NO_ARRAYS
          call write_nonlte_info(iter,.false.,.false.)
#         endif
      endif
c     
      return
      end




c     --------------------------------------------------------------
c                 MAIN ROUTINE FOR THE CONVERGENCE ERROR
c
c     This routine will call the appropriate convergence error
c     computation subroutines, and return the answer. It can 
c     return the convergence error on the mean intensity or on
c     the source function. This can be selected by the flag
c     iter_convtype.
c     --------------------------------------------------------------
      function conv_error()
      implicit none
      doubleprecision conv_error,conv_error_jbar,conv_error_intmean
      doubleprecision conv_error_sources
      logical fconv_exists
c
#include "common_source.h"
#include "common_main.h"
c
      if(iradproc_line.ne.0) then
#         ifdef INCLUDE_LINES
          if(iter_convtype.eq.0) then
              conv_error = conv_error_jbar()
          elseif(iter_convtype.eq.1) then
              write(*,*) 'Error: for line transfer is the ',
     %             'convergence on src not yet implemented'
              stop 13
          else
              stop 693
          endif             
#         endif
      else
#         ifndef JSRC_NO_ARRAYS
          if(iter_convtype.eq.0) then
              conv_error = conv_error_intmean()
          elseif(iter_convtype.eq.1) then
              conv_error = conv_error_sources()
          else
              stop 693
          endif
#         endif
      endif
#     ifdef WRITE_CONVERGENCE
      inquire(file='conv.dat',EXIST=fconv_exists)
      if(.not.fconv_exists) then
          open(unit=4,file='conv.dat',status='unknown')
          write(4,*) ' '
      else
          open(unit=4,file='conv.dat',status='old',OPENACC='append')
      endif
      write(4,*) conv_error
      close(4)
#     endif
c
      end


c     --------------------------------------------------------------
c          MAIN ROUTINE FOR THE TEMPERATURE CONVERGENCE ERROR
c     --------------------------------------------------------------
      function tempsolve_error()
      implicit none
      doubleprecision tempsolve_error
      logical fconv_exists
c
#include "common_grid.h"
#include "common_source.h"
#include "common_setup.h"
c
      doubleprecision error
      integer ir,it,itmax
c
      itmax   = irsi_frsizey
c
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#     endif
#     ifdef RADGRID_ONEDIM
      itmax = 1
#     endif
c
c     Find error
c
      tempsolve_error = 0.d0
      do ir=1,irsi_frsizex
          do it=1,itmax
              error = abs ( medium_arr_el_temp(it,ir) - 
     %                      st_temp_bk(it,ir) ) / 
     %                (medium_arr_el_temp(it,ir) + 1.d-30 ) 
                  if(error.gt.tempsolve_error) then
                      tempsolve_error = error
                  endif
          enddo
      enddo
c
c     Copy the current to the old
c
      do ir=1,irsi_frsizex
          do it=1,itmax
              st_temp_bk(it,ir) = medium_arr_el_temp(it,ir)
          enddo
      enddo
      
#     ifdef WRITE_CONVERGENCE
      inquire(file='convtemp.dat',EXIST=fconv_exists)
      if(.not.fconv_exists) then
          open(unit=4,file='convtemp.dat',status='unknown')
          write(4,*) ' '
      else
          open(unit=4,file='convtemp.dat',status='old',
     %         OPENACC='append')
      endif
      write(4,*) tempsolve_error
      close(4)
#     endif
c
      end

      
#     endif /* ifndef ONLY_RAY_TRACING */ 
