c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c               ( MODULE: MAIN MODULE FOR THE SOURCE TERMS)
c
c                            Leiden, Aug 1998
c                             C.P. Dullemond   
c
c     BUGFIX 25.09.06: Catched it.eq.0 case. Did not cause troubles
c                      for dust module because it would return 0.
c                      But it is incorrect and for lines it can
c                      cause troubles.
c     =============================================================


#define SOURCE_F_VERSION "Trans_05oct98"
#include "main.h"


c     --------------------------------------------------------------
c                CALCULATE THE SOURCE TERM AT A GRID NODE 
c
c     This is the main routine for the evaluation of the total 
c     source term at a given spatial grid point (ir,it) during the
c     integration of the formal solution in 2D by the routine 
c     do_trans(). 
c
c     The routine may be arbitrarily modified by the user. But as 
c     it is implemented now, the routine makes distinction between
c     LTE emission/absorption processes and non-LTE processes. The
c     LTE processes can be computed in advance, and depend only on 
c     position in space. Indeed, as it is implemented now, the user
c     should call the routine make_lte_src() (which is also user
c     adaptable) to prepare the arrays containing these LTE emission
c     absorption processes for each gridpoint in space. This is done
c     before calling do_trans() for the first time. The non-LTE 
c     processes depend of course on the photon distribution function
c     of the previous do_trans() run. They are also computed before
c     each call to do_trans(), but should be repeated after each 
c     do_trans(). This is done in make_nonslte_src(). 
c
c     BUGFIX 25.09.06: Catched it.eq.0 case. Did not cause troubles
c                      for dust module because it would return 0.
c                      But it is incorrect and for lines it can
c                      cause troubles.
c     --------------------------------------------------------------
      function calc_src_pnt(inu,ir,it,mu,phi)
      implicit none
      doubleprecision calc_src_pnt
      integer inu,ir,it
      integer itt
      doubleprecision mu,phi,srcdum,dum
c
#include "common_source.h"
#include "common_grid.h"
c
      calc_src_pnt = 0.d0
c
c     Check for it out of bound (which can happen when the
c     crossings with the radial grid happen inward of the
c     smallest theta-grid-point). Exceptional, but it happens.
c     Part of the BUGFIX 25.09.06
c
      itt = it
      if(itt.le.0) then
          itt = 1
      endif
      if(itt.gt.irsi_frsizey) then
          itt = irsi_frsizey
      endif
c
c     First the usual emissions: copy from array
c     This only works if these arrays have been compiled
c     along. For 2-D problems involving huge amounts of
c     frequencies (such as line transfer computations) and
c     which compute their sources on-the-fly (for line
c     transfer with LINE_VELOCITIES set!), these arrays
c     are not used anyway and may save lots of memory.
c
c     cccNOTE: JSRC_NO_ARRAYS should never be set by the user.
c     ccc      It is automatically set, for instance, when
c     ccc      LINE_VELOCITIES and
c
c     NEW: 30.03.06 Added #ifndef LINE_VELOCTIES
c
c     NOTE: If no LINE_VELOCITIES is set, then the dust source
c           is added into the source and alp arrays. If it
c           IS set, then below the dust src and alp is added
c           by calling line_dust_src_alp().
c
#     ifndef JSRC_NO_ARRAYS
#     ifndef LINE_VELOCITIES      
      calc_src_pnt = calc_src_pnt + 
     %               ssw_lte * slte_src(inu,itt,ir) + 
     %               ssw_scat * scati_src(inu,itt,ir) 
#     else
      if(iradproc_line.eq.0) then
          calc_src_pnt = calc_src_pnt + 
     %                   ssw_lte * slte_src(inu,itt,ir) + 
     %                   ssw_scat * scati_src(inu,itt,ir) 
      endif
#     endif
#     endif
c
c     Include (if necessary) the angle-dependent line emission
c     and, possibly, line dust emission
c
#     ifdef LINE_VELOCITIES
      if(iradproc_line.eq.1) then
          call line_velo_src_alp(inu,ir,itt,mu,phi,srcdum,dum)
          calc_src_pnt = calc_src_pnt + srcdum
      endif
      if(iradproc_line_dust.eq.1) then
          call line_dust_src_alp(inu,ir,itt,srcdum,dum)
          calc_src_pnt = calc_src_pnt + srcdum              
      endif
#     endif
c
c     Include the user-defined source
c
#     ifdef USER_DEFINED_SOURCE
      calc_src_pnt = calc_src_pnt +
     %      ssw_lte * USER_DEFINED_SOURCE(inu,itt,ir,mu,phi)
#     endif
c
      return
      end



c     --------------------------------------------------------------
c             COMPUTE THE ABSORPTIVITY ALPHA AT A GRID POINT
c     
c     Similar to calc_src_pnt(), but the scattering opacity now needs
c     only be calculated once and for all, for the entire iterative
c     transfer procedure.
c
c     BUGFIX 25.09.06: Catched it.eq.0 case. Did not cause troubles
c                      for dust module because it would return 0.
c                      But it is incorrect and for lines it can
c                      cause troubles.
c     --------------------------------------------------------------
      function calc_alp_pnt(inu,ir,it,mu,phi)
      implicit none
      doubleprecision calc_alp_pnt
      integer inu,ir,it
      integer itt
      doubleprecision mu,phi,dum,alpdum
c
#include "common_source.h"
#include "common_grid.h"
c
      calc_alp_pnt = 0.d0
c
c     Check for it out of bound (which can happen when the
c     crossings with the radial grid happen inward of the
c     smallest theta-grid-point). Exceptional, but it happens.
c     Part of the BUGFIX 25.09.06
c
      itt = it
      if(itt.le.0) then
          itt = 1
      endif
      if(itt.gt.irsi_frsizey) then
          itt = irsi_frsizey
      endif
c
c     First the usual absorptions: copy from array
c     This only works if these arrays have been compiled
c     along. For 2-D problems involving huge amounts of
c     frequencies (such as line transfer computations) and
c     which compute their sources on-the-fly (for line
c     transfer with LINE_VELOCITIES set!), these arrays
c     are not used anyway and may save lots of memory.
c
c     NEW: 30.03.06 Added #ifndef LINE_VELOCTIES
c
c     NOTE: If no LINE_VELOCITIES is set, then the dust source
c           is added into the source and alp arrays. If it
c           IS set, then below the dust src and alp is added
c           by calling line_dust_src_alp().
c
#     ifndef JSRC_NO_ARRAYS
#     ifndef LINE_VELOCITIES      
      calc_alp_pnt = calc_alp_pnt +
     %               asw_lte * slte_alpha(inu,itt,ir) +
     %               asw_scat * scati_alpha(inu,itt,ir)
#     else
      if(iradproc_line.eq.0) then
          calc_alp_pnt = calc_alp_pnt +
     %               asw_lte * slte_alpha(inu,itt,ir) +
     %               asw_scat * scati_alpha(inu,itt,ir)
      endif
#     endif
#     endif
c
c     Include (if necessary) the angle-dependent line opacity
c
#     ifdef LINE_VELOCITIES
      if(iradproc_line.eq.1) then
          call line_velo_src_alp(inu,ir,itt,mu,phi,dum,alpdum)
          calc_alp_pnt = calc_alp_pnt + alpdum
      endif
      if(iradproc_line_dust.eq.1) then
          call line_dust_src_alp(inu,ir,itt,dum,alpdum)
          calc_alp_pnt = calc_alp_pnt + alpdum              
      endif
#     endif
c
      return
      end



c     --------------------------------------------------------------
c                  COMPUTE THE COOLING/HEATING RATE
c
c     This subroutine computes:
c
c       Q_cool = 4*pi*(j_th+j_scat) - 4*pi*(alpha_th+alpha_scat)*J
c
c     For this it needs the evaluated arrays of slte_* and scati_*.
c     So, at least one iteration must have been done. The unit is
c     erg/sec/cm^3.
c     --------------------------------------------------------------
      subroutine compute_cooling_rate()
      implicit none
#ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_radfield.h"
#include "common_source.h"
c
      doubleprecision dummy,srcfnc,alpha
      integer ir,it,inu,itmax,imirt
c
      itmax = irsi_frsizey
      imirt = 0
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
      imirt = 1
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
      do ir=1,irsi_frsizex
          do it=1,itmax
              dummy = 0.d0
              do inu=2,freq_nr-1
                  alpha  = slte_alpha(inu,it,ir) +
     %                         scati_alpha(inu,it,ir)
                  srcfnc = ( slte_src(inu,it,ir) + 
     %                 scati_src(inu,it,ir) ) / alpha
                  dummy = dummy + 12.5663706144d0 * alpha *
     %                 (srcfnc-intmom_0(inu,it,ir)) *
     %                   (log(freq_nu(inu+1))-log(freq_nu(inu-1))) *
     %                  0.5d0 * freq_nu(inu) 
              enddo
              inu = 1
              alpha  = slte_alpha(inu,it,ir) +
     %             scati_alpha(inu,it,ir)
              srcfnc = ( slte_src(inu,it,ir) + 
     %             scati_src(inu,it,ir) ) / alpha
              dummy = dummy + 12.5663706144d0 * alpha *
     %             (srcfnc-intmom_0(inu,it,ir)) *
     %             (log(freq_nu(inu+1))-log(freq_nu(inu))) *
     %             0.5d0 * freq_nu(inu) 
              inu = freq_nr
              alpha  = slte_alpha(inu,it,ir) +
     %             scati_alpha(inu,it,ir)
              srcfnc = ( slte_src(inu,it,ir) + 
     %             scati_src(inu,it,ir) ) / alpha
              dummy = dummy + 12.5663706144d0 * alpha *
     %             (srcfnc-intmom_0(inu,it,ir)) *
     %             (log(freq_nu(inu))-log(freq_nu(inu-1))) *
     %             0.5d0 * freq_nu(inu) 
              qc_qcool(it,ir) = dummy              
          enddo
      enddo
c
#endif
      end



c     --------------------------------------------------------------
c                WRITE THE SOURCE/SINK TERMS TO FILE
c
c     This routine is just for debugging purposes, so that you can
c     view in e.g. IDL how the source is distributed. 
c
c     ARGUMENTS:
c       iter:   The iteration to which this source file belongs
c       filenr: -1 = Save to intmean.dat
c                0,1,2... = Save to intmean_0.dat etc...
c     --------------------------------------------------------------
      subroutine write_source(iter,filenr)
      implicit none
      integer ifr,nrfr,filenr,iter
c
#ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_isave.h"
#include "common_source.h"
c
      integer ir,it,imirt,itmax
      character*80 filename,base,ext
c
c     Some computations
c
      nrfr    = freq_nr 
      itmax   = irsi_frsizey
      imirt   = 0
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
      imirt = 1
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
c     Make filename
c
      if(filenr.eq.-1) then
          write(*,*) 'No append for source.dat'
          stop 13
      elseif(filenr.eq.-2) then
          filename = 'source_final.dat'
      else
          base='source_'
          ext ='.dat'
          call make_indexed_filename(base,filenr,ext,filename)
      endif    
c
c     Open file
c
      open(unit=1,file=filename,status='unknown')
      write(1,*) nrfr,irsi_frsizex,itmax,imirt
      write(1,*) ' '
      do ifr=1,freq_nr
          do ir=1,irsi_frsizex
              do it=1,itmax
                  write(1,*) slte_src(ifr,it,ir),slte_alpha(ifr,it,ir),
     %                 scati_src(ifr,it,ir),scati_alpha(ifr,it,ir)
              enddo
              if(itmax.gt.1) then
                  write(1,*) ' '
              endif
          enddo
          write(1,*) ' '
      enddo
      close(1)
c
c     If this save was succesful, and if it was not the final save,
c     then put filenr into lastsave
c
      if(filenr.gt.0) then
          lastsave_source = filenr
      endif
c
c     After this has worked, we write a small info file for
c     the IDL software, unless filenr.eq.-3, which means that
c     the source file was requested externally, and therefore
c     should not count as a 'real' source.
c
      if(filenr.ne.-3) then
          open(unit=2,file='source.info',status='unknown')
c
c         The number of the latest file (-1 means: 'source.dat', 
c         -2 means: 'source_final.dat')
c
          write(2,*) filenr
c
c         The last written source file which is not source_final.dat
c
          write(2,*) lastsave_source
c
c         The iteration nr to which the last file belongs
c
          write(2,*) iter
c
          close(2)
c
      endif
c
#else
c
c     Tell the source.info that no sources have been written
c
      open(unit=2,file='source.info',status='unknown')
      write(2,*) -100
      write(2,*) 
      write(2,*) "No sources have been written because"
      write(2,*) "the source arrays are not active."
      write(2,*) "Sources are now computed on-the-fly "
      write(2,*) "during the transfer computation."
      close(2)      
#endif
c
      end



c     --------------------------------------------------------------
c                 COMPUTE TAU(NU,THETA) AT INFINITY
c     --------------------------------------------------------------
      function tau_infinity(inu,it)
      implicit none
      doubleprecision tau_infinity
      integer inu,it,ir
c
#include "common_grid.h"
#include "common_source.h"
c
      doubleprecision calc_alp_pnt,mu,phi,alp,alpold,ds,tauinf
c
      mu     = 1.d0
      phi    = 0.d0
      tauinf = 0.d0
      alp    = calc_alp_pnt(inu,1,it,mu,phi)
      do ir=2,irsi_frsizex
          alpold = alp
          alp    = calc_alp_pnt(inu,ir,it,mu,phi)
          ds     = rsi_x_c(ir,1) - rsi_x_c(ir-1,1)
          tauinf = tauinf + 0.5 * ds * ( alp + alpold )
      enddo
c
      tau_infinity = tauinf
      return
      end


c     --------------------------------------------------------------
c              WRITE FILE WITH TAU(NU,THETA) AT INFINITY
c     --------------------------------------------------------------
      subroutine write_tau_at_infinity()
      implicit none
c
#include "common_grid.h"
#include "common_source.h"
c
      integer inu,ir,it,imirt,itmax
      doubleprecision tau,tau_infinity
      integer alte,ascat
c
      alte   = asw_lte
      ascat  = asw_scat
      asw_lte  = 1
      asw_scat = 1
c
c     Some computations
c
      itmax   = irsi_frsizey
      imirt   = 0
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
      imirt = 1
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
      open(unit=1,file='tauinf.dat',status='unknown')
      write(1,*) freq_nr,itmax,imirt
      write(1,*) 
      do it=1,itmax
          do inu=1,freq_nr
              tau = tau_infinity(inu,it)
              write(1,10) tau
 10           format(E13.6)
          enddo
          write(1,*) 
      enddo
c
      close(1)
c
      asw_lte  = alte 
      asw_scat = ascat
c
      end


c     --------------------------------------------------------------
c                     READ THE MC SCATTERING SOURCE
c
c     Added January 2005: Also read quantsource.dat
c     --------------------------------------------------------------
      subroutine read_scatsource()
      implicit none
#ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_source.h"
#include "common_dust.h"
c
      doubleprecision dum
      integer it,ir,itmax,ifreq,imirt,i,j,k,l
      logical fex
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
      imirt = 1
#else
      itmax = irsi_frsizey
      imirt = 0
#endif
c
      inquire(file='scatsource.dat',exist=fex)
      if(fex) then
          write(*,*) 'Found file scatsource.dat! Loading...'
          open(unit=1,file='scatsource.dat',status='old')
          read(1,*) i,j,k,l
          do ifreq=1,freq_nr
              do it=1,itmax      
                  do ir=1,irsi_frsizex
                      read(1,*) dum
                      scati_src(ifreq,it,ir) = 
     %                      scati_src(ifreq,it,ir) + dum
                  enddo
              enddo
          enddo
          close(1)
      endif
c
      inquire(file='quantsource.dat',exist=fex)
      if(fex) then
          write(*,*) 'Found file quantsource.dat! Loading...'
          open(unit=1,file='quantsource.dat',status='old')
          read(1,*) i,j,k,l
          do ifreq=1,freq_nr
              do it=1,itmax      
                  do ir=1,irsi_frsizex
                      read(1,*) dum
                      scati_src(ifreq,it,ir) = 
     %                      scati_src(ifreq,it,ir) + dum
                  enddo
              enddo
          enddo
          close(1)
      endif
#endif /* JSRC_NO_ARRAYS */
      end

c     ==============================================================
c                      NOW THE USER-DEFINED SOURCE
c
c     Since it is not very useful to make the user source a separate
c     module, it is better to include it into this module.
c     ==============================================================

#ifdef USER_DEFINED_SOURCE
#include userdefsource.F
#endif

