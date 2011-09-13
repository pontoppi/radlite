c     =============================================================
c
c                RADLITE - a raytracer for polar coordinates and 
c                axial symmetry. Derived from the Radical package by
c                C.P. Dullemond. 
c
c                  K. M. Pontoppidan, Caltech, November 2008 
c
c                                  Based on:
c
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c                            Leiden, June 1998
c                             C.P. Dullemond
c
c     =============================================================

#define MAIN_F_VERSION "22-03-07"
#include "main.h"


c     --------------------------------------------------------------
c                          PROGRAM RADLITE
c     --------------------------------------------------------------
      program radicalprog
      implicit none

      logical fexists
      integer ierr    !error code

#include "common_main.h"
c
c     Default stuff
c

      main_imethod              = 1
      main_nrref                = 10      

c
c     Welcome banner
c
      write(*,*) '=============================================='
      write(*,*) '                 RADLite                      '
      write(*,*) '         (C) 2007-2011    V. 1.2              '
      write(*,*) '     K. M. Pontoppidan &  C.P. Dullemond      '
      write(*,*) '=============================================='
c
c     Messages
c
c     Check if input file "radlite.inp" exists
c
      inquire(file='radlite.inp',EXIST=fexists)
      if(.not.fexists) then
         write(*,*) 'ERROR: Could not find the main input file: ',
     %        'radlite.inp'
         write(*,*) '       This is the most crucial input file, ',
     %        'so you will have to provide one.'
         stop 13
      endif
c
c         If "radlite.inp" exists, then do the interpretation and 
c         run the code.
c
      call main_routine()
c
c     If we ended up here, no major problems have been encountered
c
      open(unit=1,file='radlite.success',status='unknown')
      write(1,*) '1'
      close(1)
c
      end


c     --------------------------------------------------------------
c                             MAIN ROUTINE
c     --------------------------------------------------------------
      subroutine main_routine()
      implicit none
c
#include "common_grid.h"
#include "common_source.h"
#include "common_setup.h"
#include "common_isave.h"
#include "common_boundary.h"
#include "common_telescope.h"
#include "common_lines.h"
#include "common_dust.h"
#include "common_main.h"

c
      doubleprecision dustmass,total_dust_mass
      integer inu
      logical fexists
c
c     Permanently set the code to do line raytracing only without
c     the now completely disabled non-lte iteration:
c

      iradproc_line = 1
      iradproc_line_dust = 1

c
c     Defaults
c
      main_distance = 3.08572d18
c
c     Now copy the radlite.inp file to temp.inp, but with the comment
c     lines removed.
c
      open(unit=1,file='radlite.inp')
      open(unit=2,file='temp.inp')
      call strip_comments(1,2)
      close(2)
      close(1)
c
c     Open the stripped input file
c
      open(unit=2,file='temp.inp',status='old',err=1000)
c
c     Call the interpreter routine for radlite.inp . This routine
c     will set all internal variables according to the wish of the
c     user, but it will not yet call the routines for setting up
c     the grid variables.
c
      call read_radical_input(2)
c
c     Now set up the frequency/spatial/angular grids
c
      call create_all_grids()
c
c     Now put the problem setup onto the grid: the geometry of the density
c     profile, and the type of radiative processes to solve etc. This
c     is done by the medium.F module. Also it is decided here whether
c     to do line transfer or continuum etc. This module reads some data
c     from the input file, so therefore it must be read here.
c
      call read_compute_medium_lean(2)
c
c     Close unit 2
c
      close(2)
c
c     Remove the temp.inp file
c
      call system('rm -f temp.inp')
c
c     NEW: 27.03.06: Moved the read_dustdata() to here (from much further
c     down in the radical_start() routine).
c
c     Read the dust data
c
#     ifdef INCLUDE_DUST
      if((iradproc_dust.ne.0).or.(iradproc_line_dust.ne.0)) then
          call read_dustdata()
          dustmass = total_dust_mass(-1,0)
          write(*,'(a,e9.4,a)') 'Total dust mass in setup: ',dustmass*
     %               5.02512562814e-34,' M_sun (cell-centered grid)'
          dustmass = total_dust_mass(-1,1)
          write(*,'(a,e9.4,a)') 'Total dust mass in setup: ',dustmass*
     %               5.02512562814e-34,' M_sun (corner-centered grid)'
      endif
#     endif
c
c     If line is set and line_dust is also set, then read the dust
c     scattering sources in (= scattering + PAH).
c

c      if(iradproc_line_dust.ne.0) then
c          call read_scatsource()
c      endif

c
c     If line radiation is activated then read the line information,
c     set up the frequency passbands, read the radiative and collisional
c     rates, set up the line profiles etc. 
c
c$$$      if(iradproc_line.eq.1) then
c$$$#         ifdef INCLUDE_LINES
c$$$          call read_linedata(0)
c$$$          call write_global_lineprofile()
c$$$#         else
c$$$          write(*,*) 'PROBLEM: Line processes switched off'
c$$$          write(*,*) '         Recompile with appropriate switch'
c$$$          stop 1
c$$$#         endif
c$$$      endif
c
c     Do some more initialization stuff
c
      imr_nx     = 50
      imr_ny     = 50
      imr_spx    = rsi_x_c(irsi_frsizex,1) / (0.5d0*imr_nx)
      imr_spy    = rsi_x_c(irsi_frsizex,1) / (0.5d0*imr_ny)
      imr_phioff = 0.d0
      imr_xoff   = 0.d0
      imr_yoff   = 0.d0
      imr_szimx  = -1.d0
      imr_szimy  = -1.d0

      call read_telesc_linespec()
      call read_linedata(0)
 
      linespec_ilinestart = 1
      linespec_distance = 1.d0  ! pc
      main_passbnfr   = line_nfreq(linespec_ilinestart)
      main_passbwidth = line_rangewidth(linespec_ilinestart)
      main_passbnfr      = 2.d0*linespec_vmax / linespec_dv + 1
      linespec_vmax = 0.5d0 * (main_passbnfr-1) * linespec_dv
      main_passbwidth    = linespec_vmax
      main_anginf        = linespec_incang
      linespec_distance = linespec_distance * 3.08572d18
c      
c     permanently disable all non-lte remnants
      iter_trans = -1
c
c     Now that we know all we need, lets do some action:
c
      call radical_start()
c
c     And now we're done...
c
      goto 2000
 1000 continue
      write(*,*) 'Could not open file temp.inp'
      stop 13
 2000 continue
      return
      end


c     --------------------------------------------------------------
c                      READ THE APERTURE FILE
c
c     This is for the spectrum, if an aperture as a function of 
c     wavelength is requested. This aperture should be given 
c     in arcsecond.
c     --------------------------------------------------------------
      subroutine read_aperture_file()
      implicit none
c
#include "common_grid.h"
#include "common_main.h"
#include "common_telescope.h"
c
      integer nf,inu
      logical fex
c
      inquire(file='aperture.inp',exist=fex)
      if(fex) then
          write(*,*) '***********************'
          write(*,*) 'Reading aperture.inp...'
          write(*,*) '***********************'
          open(unit=1,file='aperture.inp',status='old')
          read(1,*) nf
          if(nf.ne.freq_nr) then
              write(*,*) 'ERROR: Nr of frequency points in aperture.inp'
              write(*,*) '       is inconsistent with rest of the setup'
              stop
          endif
          do inu=1,freq_nr
              read(1,*) spec_aperture(inu)
          enddo
          close(1)
          telesc_aperture = 1
      endif
c
      end


c     --------------------------------------------------------------
c                    MAIN ROUTINE FOR CREATING GRIDS
c
c     --------------------------------------------------------------
      subroutine create_all_grids()
      implicit none
c
#include "common_grid.h"
#include "common_main.h"
#include "common_boundary.h"
c
      integer i
c
c     Do some checking on the dimensions
c
      if((irsi_frsizex.gt.FRSIZE_X).and.(spacegrid_type_x.ne.'-')) then
          write(*,*) 'x-size too large!'
          write(*,*) 'Requested size = ',irsi_frsizex,
     %               '; Max size = ',FRSIZE_X
          write(*,*) 'Recompile with larger array size FRSIZE_X'
          stop 1
      endif
      if((irsi_frsizey.gt.FRSIZE_Y).and.(spacegrid_type_y.ne.'-')) then
          write(*,*) 'y-size too large!'
          write(*,*) 'Requested size = ',irsi_frsizey,
     %               '; Max size = ',FRSIZE_Y
          write(*,*) 'Recompile with larger array size FRSIZE_Y'
          stop 1
      endif
      if(anggrid_frsizemu.gt.FRSIZE_MU) then
          write(*,*) 'mu-size too large!'
          write(*,*) 'Requested size = ',anggrid_frsizemu,
     %               '; Max size = ',FRSIZE_MU
          write(*,*) 'Recompile with larger array size FRSIZE_MU'
          stop 1
      endif
      if(anggrid_frsizephi.gt.FRSIZE_PHI) then
          write(*,*) 'phi-size too large!'
          write(*,*) 'Requested size = ',anggrid_frsizephi,
     %               '; Max size = ',FRSIZE_PHI
          write(*,*) 'Recompile with larger array size FRSIZE_PHI'
          stop 1
      endif
      if((freq_nr.gt.FRSIZE_FREQ).and.(freq_grid_type.ne.-1)) then
          write(*,*) 'freq-size too large!'
          write(*,*) 'Requested size = ',freq_nr,
     %               '; Max size = ',FRSIZE_FREQ
          write(*,*) 'Recompile with larger array size FRSIZE_FREQ'
          stop 1
      endif
c
c     Set up the frequency array
c
      call create_freqgrid(freq_nu_start,freq_nu_end,
     %                     freq_nr,freq_grid_type,0)
      if(freq_done.ne.0) then 
          write(*,*) 'Creating frequency grid'
      else
          write(*,*) 'No global frequency grid created...'
      endif
c 
c     Setup the spatial grid (R_i,Theta_j)
c
      write(*,*) 'Creating spatial grid'
      call create_spacegrid()
c
c     Set up the angular grid (mu_k,phi_l)
c
#     ifndef ONLY_RAY_TRACING
      write(*,*) 'Creating angular grid'
c
c     First the main grid (largest deltaR/R)
c
ccccc (CHANGED 18-09-00)  anggrid_drr = spacegrid_dxi/spacegrid_xi
      call find_largest_rgrid_drr(anggrid_drr)
      call create_anggrid(anggrid_mu_type,1)
c
c     For the ESC/MESC method this grid is automatically copied to
c     the rest. For LC the copy anggrid is not done automatically, 
c     so do it here
c
#     ifdef LC_INTEGRATION
      call copy_anggrid(1,0)
      do i=1,FRSIZE_DRR 
          call copy_anggrid(1,i)
      enddo
#     endif
c
#     endif      
c
c     Check whether star radius and grid radius are not conflicting
c
      if(radbnd_rstar.ge.spacegrid_xi) then
          write(*,*) 'Sorry, star is bigger than inner grid radius'
          write(*,*) 'Rin = ',spacegrid_xi
          write(*,*) 'R*  = ',radbnd_rstar
          stop 13
      endif
c
      end





c     --------------------------------------------------------------
c                        READ THE INPUT FILE
c     
c     This routine reads the "radlite.inp" file, which is a fixed-
c     format file which the user should provide, to define the 
c     problem to be solved. 
c
c     The file begins with a format-version-number. During the
c     future development of RADLITE, this number is raised each 
c     time a new fixed-format-line is added to the input structure, 
c     or if the input structure is modified altogether. This is 
c     meant to allow for future upgrades and simplicifations of 
c     the code. 
c
c     The format of the file is NOT ONLY determined by the format-
c     version-number at the first line. This number only determines
c     the format of the first part of the input file (which provides
c     information on the grid setup, the number of iterations etc). 
c     The second part is the part that specifies the particular
c     physical problem to solve. This part starts with a line that
c     specifies the index of the internally-generated-setup (see the
c     medium.F file, which the user may modify and add her/his own
c     internally generated setups). The remaining lines in the
c     "radlite.inp" file are in a format completely up to the
c     author of the internally-generated-setup style at hand. 
c     Usually again a fixed format style, but the format depends
c     on the setup. 
c     
c     NOTE: For some types of simulations (such as line transfer)
c           other input files may be necessary. For line transfer
c           the line.inp file is required. 
c
c     --------------------------------------------------------------
      subroutine read_radical_input(ifile)
      implicit none
c
      integer ifile,input_style
      character*80 string
c
#include "common_main.h"      
c     
c     Write warnings to screen
c
#     ifdef INTERPOL_SRC_0
      write(*,*) '===> WARNING: Source has 0-th ',
     %              'order interpolation'
#     endif
#     ifdef ONLY_RAY_TRACING
      write(*,*) '===> WARNING: Short Characteristics disabled'
      write(*,*) '===>          All non-lte/scattering iterations ',
     %           'are forbidden'
#     endif
c
c     Some inits
c
      image_ifreq = 1
      iter_start  = 1
c
c     Read the input style format version
c
      read(ifile,*,err=9001,end=9002) input_style
c
c     Switch input styles
c
      if(input_style.ge.100) goto 2
c
      write(*,*) 'ERROR while reading radlite.inp:'
      write(*,*) 'First line states that the format style of '
      write(*,*) 'radlite.inp is ',input_style,'.'
      write(*,*) 'This format version is not supported.'
      stop 13
c
 1    continue
c      call radical_input_style_one(ifile)
      goto 9999
c
 2    continue
      call radical_input_style_short(ifile,input_style)
      goto 9999
c
 9001 continue
      write(*,*) 'ERROR while reading first two lines of ',
     %     'radlite.inp'
      stop 13
 9002 continue
      write(*,*) 'ERROR while reading first two lines of ',
     %     'radlite.inp: end of file...'
      stop 13
c
 9999 continue
c
      end





c     --------------------------------------------------------------
c                        INPUT STYLE NUMBER 100
c
c     The input style number 100 represents a new philosophy in the
c     input style of RADICAL. In this mode, which will be the future
c     choice, all input will be done in separate input files. No
c     internally generated grids will be made, nor will the setup
c     be internally generated. Everything will be loaded. Only a very
c     few things will be done internally, such as the angles of the
c     radiation. 
c
c     CURRENT FORMAT NR: 102
c
c     --------------------------------------------------------------
      subroutine radical_input_style_short(ifile,iformat)
      implicit none
c
      integer ifile,iformat
      character*80 string
c
#include "common_grid.h"
#include "common_main.h"
#include "common_isave.h"
#include "common_boundary.h"
c
      doubleprecision inclination,ddummy
      integer idummy,intstore,tempstore
c
c     Make a set of defaults
c
      irsi_frsizex            = 4
      irsi_frsizey            = 4
      spacegrid_metric        = 'S'
      spacegrid_type_x        = '-'
      spacegrid_type_y        = '-' 
      spacegrid_xi            = 1.0
      spacegrid_xo            = 0.0
      spacegrid_dxixi         = 0.1
      spacegrid_refine_thetar = 1.0
      spacegrid_refine_a      = 0.0
      spacegrid_refine_n      = 0
      spacegrid_dxi           = 0.1
      freq_nr                 = 3
      freq_grid_type          = -1
      freq_nu_start           = 1.0
      freq_nu_end             = 2.0
      anggrid_drr             = 1.d-30 ! Safety value
      isave_short_char        = 0
      isave_long_char         = 0
      iradbnd_in_spectype     = -2
c
c     Conditional defines
c
#     ifdef RADGRID_ONEDIM
      nrphi(1) = 1
      main_nrphiinf = 1
#     endif
c
c     Parameters for the iteration procedure (if done)
c
c      read(ifile,*,err=1,end=2) iter_trans           ! Maximum allowed number of iterations 
      read(ifile,*,err=1,end=2) idummy           ! Maximum allowed number of iterations 
c      read(ifile,*,err=1,end=2) iter_method          ! Method of iteration [ 0=LI,1=ALI,2=ALI+Ng,3=LI=Ng,-1=IS ]
      read(ifile,*,err=1,end=2) idummy          ! Method of iteration [ 0=LI,1=ALI,2=ALI+Ng,3=LI=Ng,-1=IS ]
      if(iformat.ge.103) then
c          read(ifile,*,err=1,end=2) iter_fluxcons    ! Enable flux cons tricks
          read(ifile,*,err=1,end=2) idummy    ! Enable flux cons tricks
c      else
c          iter_fluxcons = 0.d0
      endif
c      read(ifile,*,err=1,end=2) iter_convcrit        ! Convergence criterion (tolerance)
c      read(ifile,*,err=1,end=2) iter_convtype        ! Convergence style [ 0=J_nu, 1=j_nu ]
c      read(ifile,*,err=1,end=2) iter_guess           ! Initial guess style [ 0=none,1=noscat/radequil ]
      read(ifile,*,err=1,end=2) ddummy        ! Convergence criterion (tolerance)
      read(ifile,*,err=1,end=2) idummy        ! Convergence style [ 0=J_nu, 1=j_nu ]
      read(ifile,*,err=1,end=2) idummy           ! Initial guess style [ 0=none,1=noscat/radequil ]
c
c     Then read the setup of the angular grid for short chars
c
      read(ifile,*,err=1,end=2) anggrid_frsizemu
      read(ifile,*,err=1,end=2) anggrid_frsizephi
      read(ifile,*,err=1,end=2) anggrid_mu_type
      read(ifile,*,err=1,end=2) anggrid_muerr_max
      read(ifile,*,err=1,end=2) dmuvdr
      read(ifile,*,err=1,end=2) iextrmu(1)
c
c     Now read info about the inner and outer boundary
c
      read(ifile,*,err=1,end=2) iradbnd_out_itype    ! Outer boundary condition
      read(ifile,*,err=1,end=2) iradbnd_in_itype     ! Inner boundary condition
      read(ifile,*,err=1,end=2) idummy               ! Equatorial bndry (future)
c
c     Next read the info about the distant observer
c      (even when telesc_command.eq.0)
c
      read(ifile,*,err=1,end=2) inclination
      read(ifile,*,err=1,end=2) main_nrphiinf
      read(ifile,*,err=1,end=2) telesc_dbdr  
      read(ifile,*,err=1,end=2) main_nrrayextra
c
c     Then read which things to write to a file during
c     the transfer computation. Most of these things are
c     only for the debugging and testing phase.
c     
      read(ifile,*,err=1,end=2) isave_scatnonlte
      read(ifile,*,err=1,end=2) isave_intens_inu
      read(ifile,*,err=1,end=2) isave_source
      read(ifile,*,err=1,end=2) intstore
      if(iformat.ge.101) then
          read(ifile,*,err=1,end=2) tempstore
      else
          tempstore=intstore
      endif
      if(iformat.ge.102) then
          read(ifile,*,err=1,end=2) isave_alioper
      else
          isave_alioper = 0
      endif
      read(ifile,*,err=1,end=2) isave_physvar
      if(iformat.ge.104) then
          read(ifile,*,err=1,end=2) isave_fluxcons
      else
          isave_fluxcons = 0
      endif
c
c     Do some interpretation  (23-11-00)
c
      if(iformat.lt.100) then
          if(intstore.ne.0) then
              isave_meanint  = -1
              isave_dusttemp = -1  
              isave_levpop   = -1
              isave_electemp = -1
          endif
      else
          isave_meanint  = intstore
          isave_dusttemp = tempstore
          isave_levpop   = tempstore
          isave_electemp = tempstore
      endif
c
c     Some computations
c
      main_anginf = inclination * 0.01745329252d0
c
c     Do a check
c
      if(anggrid_mu_type.ne.2) then
          write(*,*) '--------------------------------------'
          write(*,*) '  ONLY MUGRID TYPE 2 ALLOWED FOR NOW  '
          write(*,*) '  Switching automatically to type 2...'
          write(*,*) '--------------------------------------'
          anggrid_mu_type = 2
      endif
      if(main_nrrayextra.eq.0) then
          write(*,*) '**************************************'
          write(*,*) 'Must set nr of rays inward of Rin to a'
          write(*,*) 'value > 0... Otherwise flux might be'
          write(*,*) 'lost...'
          write(*,*) '**************************************'
          stop
      endif
      if(main_nrrayextra.gt.0) then
          write(*,*) '============================================'
          write(*,*) 'WARNING: Still using old way of distributing'
          write(*,*) '         rays inward of Rin (for spectrum).'
          write(*,*) '         Could be dangerous! Better use a'
          write(*,*) '         negative value of nrrayextra, which'
          write(*,*) '         means that it distributes the rays'
          write(*,*) '         in a more dynamic way.'
          write(*,*) '============================================'
      endif
c
c     Done
c
      goto 10
c
 1    continue
c
c     Error: could not read line in input file
c     
      read(ifile,*,err=1003,end=1003) string
      write(*,*) 'ERROR in reading input file. Last line:'
      write(*,*) string
      write(*,*) 'Next few lines:'
      read(ifile,*,err=1004,end=1004) string
      write(*,*) string
      read(ifile,*,err=1004,end=1004) string
      write(*,*) string
      read(ifile,*,err=1004,end=1004) string
      write(*,*) string      
      read(ifile,*,err=1004,end=1004) string
      write(*,*) string      
      read(ifile,*,err=1004,end=1004) string
      write(*,*) string      
      read(ifile,*,err=1004,end=1004) string
      write(*,*) string      
      read(ifile,*,err=1004,end=1004) string
      write(*,*) string      
      stop 133
 1003 continue
      write(*,*) 'FATAL ERROR in reading input file'
      stop 133
 1004 continue
      stop 133
c
 2    continue
      write(*,*) 'ERROR in reading input file. EOF detected'
      stop 133
c
 10   continue
c
      end






c     --------------------------------------------------------------
c                         START THE COMPUTATION
c
c     This routine is basically the main switchboard for doing the
c     various computations, and making sure that the appropriate
c     output files will be written. 
c
c     The telesc_command has the following meanings:
c       
c       1      Image       (circular, i.e. only for internal use)
c       2      Spectrum    (continuum)
c       3      Image       (normal, i.e. rectangular)
c       4      Spectrum    (line, i.e. only when INCLUDE_LINES is set)
c       5      Spectrum    (lines, i.e. more than a single line)
c       6      P-V         (line, i.e. position-velocity diagram)
c       7      Img Series  (circular images at all frequencies)
c       8      Img Series  (rectangular images at all frequencies)
c       9      Spectra     (dust, spectrum at all theta)
c
c     --------------------------------------------------------------
      subroutine radical_start()
      implicit none
c
#include "common_grid.h"
#include "common_source.h"
#include "common_setup.h"
#include "common_isave.h"
#include "common_boundary.h"
#include "common_telescope.h"
#include "common_lines.h"
#include "common_main.h"
c
c      doubleprecision dustmass,total_dust_mass
      doubleprecision luminosity,outer_luminosity_total
      doubleprecision vet_fraccov
      integer iline,iln,itertrans_orig,ifreq
      logical fexists,fex
c
      if(radical_quit.ne.0) then
          goto 999
      endif
c
c     Some checks
c
      if((iradproc_line.eq.0).and.(iradbnd_out_itype.eq.2)) then
          write(*,*) 'INCONSISTENCY:'
          write(*,*) '  Outer BC type 2 only for line processes'
          stop 1
      endif
c
c     First read the stellar radius, if necessary
c
      call read_starinfo()
c
c     Render the short characteristics, if iteration will be done
c

#     ifndef NO_SHORT_CHARS
      if((iter_trans.gt.0).or.(radical_calcflux.ne.0)) then
#         ifndef SHORTCHAR_NEWSTYLE
c
c         Old style (pre 09-00)
c
          if((iradbnd_in_itype.eq.2).or.(iradbnd_in_itype.eq.0)) then
              write(*,*) 'Creating the short characteristics for ',
     %                   'the inner vacuum'
              if(radbnd_rstar.le.0.d0) then
                  write(*,*) 'ERROR: With central star as inner BC'
                  write(*,*) '       the stellar radius must be'
                  write(*,*) '       non-zero'
                  stop 13
              endif
              call copy_anggrid(1,2)
              call create_shortchar(rsi_x_c(1,1)/radbnd_rstar-1.d0,
     %                     rsi_x_c(1,1)/radbnd_rstar-1.d0,2,0)
          endif
          write(*,*) 'Creating the short characteristics for ',
     %               'the logaritmic grid'
          call create_shortchar(spacegrid_dxi/spacegrid_xi,
     %                     spacegrid_dxi/spacegrid_xi,1,1)
c
#         else
c
c         New style (post 09-00)
c
          call setup_all_shortchar_sets()
c
#         endif
      endif
#     endif
c
c     Render (or at least set up) the long characteristics, if 
c     spectrum or image is going to be made. 
c
#     ifndef RADGRID_CARTESIAN
      if((telesc_command.eq.1).or.(telesc_command.eq.2).or.
     %   (telesc_command.eq.4).or.(telesc_command.eq.5).or.
     %   (telesc_command.eq.7)) then
          write(*,*) 'Setting up the long characteristics'
          call setup_rays_circular(1,irsi_frsizex,1,main_anginf,
     %              main_nrphiinf,main_nrrayextra,telesc_dbdr,
     %              radbnd_rstar,main_imethod,main_nrref)
      endif
      if((telesc_command.eq.3).or.(telesc_command.eq.6).or.
     %   (telesc_command.eq.8)) then
#         ifndef RADGRID_ONEDIM
          write(*,*) 'Setting up the long characteristics'
          call setup_rays_rectang(imr_nx,imr_ny,imr_spx,imr_spy,
     %                  main_anginf,imr_phioff,imr_xoff,imr_yoff)
#         else
          write(*,*) 'No rectangular image for 1-D'
          stop 1
#         endif
      endif
#     endif
c
c     Now prepare the LTE source function arrays. This works only
c     in case the source function arrays are compiled along, which
c     they are by default. In some cases (for instance, for line
c     transfer) the arrays tend to become so big that they constrain
c     the 2-D transfer too much. For line transfer this is solved
c     by simply removing them. The source terms for line transfer
c     are computed on-the-fly anyway, so the arrays are no longer
c     necessary.
c
#     ifndef LINE_DUSTCONT   /* NEW CONDITION 27-03-06 */
#     ifndef JSRC_NO_ARRAYS
      if(iradproc_line.eq.0) then
          write(*,*) 'Making the emission and absorption functions'
          call make_lte_src()
      endif
#     endif
#     endif
c
c     Now read the dust information, if the dust modules are
c     activated by their hard-switch
c
c     NEW: 27.03.06: MOVED read_dustdata() to earlier (before
c     start of this routine).
c
c#     ifdef INCLUDE_DUST
c      if((iradproc_dust.ne.0).or.(iradproc_line_dust.ne.0)) then
c          call read_dustdata()
c          dustmass = total_dust_mass(-1,0)
c          write(*,*) 'Total dust mass in setup: ',dustmass *
c     %               5.02512562814e-34,' M_sun (cell-centered grid)'
c          dustmass = total_dust_mass(-1,1)
c          write(*,*) 'Total dust mass in setup: ',dustmass *
c     %               5.02512562814e-34,' M_sun (corner-centered grid)'
c      endif
c#     endif
c
c     Now prepare the Scattering 
c       (of course the souce terms not yet)
c       (also possible non-LTE stuff not, because the opacity
c        could then depend on the level populations)
c       (if no scattering is specified, then it returns without
c        doing anything)
c
c$$$#     ifndef JSRC_NO_ARRAYS
c$$$      if(iradproc_line.eq.0) then
c$$$          write(*,*) 'Making the scattering/nonlte opacity function'
c$$$          call make_scat_alpha(1,freq_nr)
c$$$      endif
c$$$#     endif
c
c     Always write the spatial grid to file
c
      call write_grid()
c
c     If requested, save the short characteristics to file
c
#     ifndef NO_SHORT_CHARS
c#     ifndef SHORTCHAR_NEWSTYLE
      if((isave_short_char.eq.1)) then
          write(*,*) 'Writing short characteristics to file'
          call write_shortchar(1)
          if((iradbnd_in_itype.eq.2).or.(iradbnd_in_itype.eq.0)) then
              write(*,*) 'Writing inner short characteristics to file'
#             ifdef SHORTCHAR_NEWSTYLE
              call write_shortchar(0)
#             else
              call write_shortchar(2)
#             endif
          endif
      endif
c#     endif
#     endif
c
c     If requested, save the density and temperature (and possible
c     other physical variables like the magnetic field) to file
c
      if(isave_physvar.eq.1) then
         write(*,*) 'Writing physical variables of medium to file'
         call write_phys_var()
         if(ieqd_active.eq.1) then
c             call write_disk_model()
             write(*,*) 'writing disk model disabled, Klaus 22/03/07'
             stop 13
         endif
      endif
c
c     If line processes activated, AND the LINE_VELOCITIES is set,
c     then we must prepare all the local line widths. Also we must 
c     copy the global abundances into the local arrays if the local 
c     arrays are empty.
c
#     ifdef LINE_VELOCITIES
      if(radical_calcflux.ne.0) then
c
c         ...Iteration, therefore: with check of passband size/resolution
c
          call line_velo_prep_profiles_abun(1)
      else
c
c         ...Spectrum/Image, therefore: without check of passband size/resolution
c
          call line_velo_prep_profiles_abun(0)
      endif
#     endif
c
c     Now check if the gas velocity sampling is not too coarse
c
#     ifdef LINE_VELOCITIES
      call check_velo_step_linewidth()
#     endif
c
c     Now fill the stellar spectrum
c
      if(iradproc_dust.ne.0) then
          call setup_starspec(0)
      endif
      if(iradproc_line_starpump.ne.0) then
          call setup_starspec(1)
      endif
c
c     If appropriate, read the interstellar radiation field
c     Also if lines are active. It will then below be replaced
c     by the boundary conditions at the line frequencies.
c
      if(iradbnd_out_itype.eq.3) then
          if(iradproc_dust.ne.0) then
              call read_interstellfield(0)
          elseif(iradproc_line.ne.0) then
              call read_interstellfield(1)
          else
              stop 83192
          endif
      endif
c
c     If line_dust and/or line_starpump or iradbnd_out_itype=3 then do 
c     interpolation of quantities onto the spectral lines
c
c$$$      if(iradproc_line.ne.0) then 
c$$$c
c$$$c         If star is allowed to pump the lines, compute these
c$$$c
c$$$          if(iradproc_line_starpump.ne.0) then
c$$$              call line_init_starbc()
c$$$          endif
c$$$c
c$$$c         If outer boundary has interstellar radiation field
c$$$c         NOTE: In my thesis I used outer boundary 2 for the
c$$$c               cosmic microwave background radiation. Here
c$$$c               I read in an entire spectrum (see above) and
c$$$c               convert this into the boundary conditions for
c$$$c               the lines.           
c$$$c
c$$$          if(iradbnd_out_itype.eq.3) then
c$$$              call line_init_outerbc()
c$$$          endif
c$$$      endif
c
c     Now all the preparations are finished, and we must go to
c     business. First read in the level populations
c
c
      call read_nonlte_scat()

c
c     Now that the sources and opacities have been
c     read in, we dump a file containing the optical depths at radial
c     infinity. This is only for convenience to the user.
c
      call write_tau_at_infinity()
c
c     --------------------------------------------------------------
c             PRODUCE THE IMAGE AT INFINITY    (IF REQUESTED) 
c     --------------------------------------------------------------
c
c     Write the start log message to radlite.log
c
      call write_logmessage_start('spectrum/image.')
c
c     If requested, produce an image/spectrum using the telescope
c
      if(telesc_command.gt.0) then
          if(iscat_done_nonlte.ne.1) then
              write(*,*) '--- WARNING: Non-LTE stuff not yet ',
     %             'read/computed! ---'
          endif

c         Enable the emission processes and all opacities, by setting the
c         appropriate flags
c     
          asw_lte    = 1
          asw_scat   = 1
          ssw_lte    = 1
          ssw_scat   = 1
c
c         Now make the image at earth, or the corresponding spectrum
c
c$$$          if(telesc_command.eq.4) then
c$$$c
c$$$c             The spectrum for a single line  (linespec_command=0)
c$$$c
c$$$#             ifdef INCLUDE_LINES
c$$$              write(*,*) 'Rendering spectrum of line ',
c$$$     %                    linespec_ilinestart
c$$$              write(*,*) '   Inclination = ',
c$$$     %                     main_anginf*57.2957795129d0,' degrees'
c$$$              call calc_line_spectrum(linespec_ilinestart,
c$$$     %               linespec_distance,main_passbwidth,main_passbnfr)
c$$$              write(*,*) 'Writing the spectrum'
c$$$              call header_line_spectrum(1,main_passbnfr,
c$$$     %               linespec_distance,linespec_radvelo,
c$$$     %               main_anginf,linespec_style)
c$$$              call write_line_spectrum(linespec_style)
c$$$c
c$$$c             Display line-center optical depths
c$$$c
c$$$              iline=linespec_ilinestart
c$$$              write(*,'(a,i,a,e9.4)') 'Line ',iline,'  tau = ',
c$$$     %                 char_tau_linecenter(iline)
c$$$#             else
c$$$              write(*,*) 'Line transfer is not active'
c$$$              stop 13
c$$$#             endif
c$$$          endif
          if(telesc_command.eq.5) then
c
c             The spectrum for a series of lines (linespec_command=0)
c
#             ifdef INCLUDE_LINES
              call header_line_spectrum(linespec_nlines,
     %                      main_passbnfr,linespec_distance,
     %                      linespec_radvelo,main_anginf,
     %                      linespec_style)
              write(*,*) '   Inclination = ',
     %                     main_anginf*57.2957795129d0,' degrees'
              do iln=linespec_ilinestart,linespec_ilinestart+
     %                          linespec_nlines-1
                  call calc_line_spectrum(iln,linespec_distance,
     %               main_passbwidth,main_passbnfr)
                  write(*,*) 'Rendered spectrum of line ',iln
                  call write_line_spectrum(linespec_style)
              enddo
c
c             Display line-center optical depths
c
c              do iline=linespec_ilinestart,linespec_nlines
c                  write(*,*) 'Line ',iline,'  tau = ',
c     %                 char_tau_linecenter(iline)
c              enddo
#             else
              write(*,*) 'Line transfer is not active'
              stop 13
#             endif
          endif
          if(telesc_command.eq.6) then
c
c             The 3-D position-velocity data cube (linespec_command=2)
c
#             ifdef INCLUDE_LINES
              do iln=linespec_ilinestart,linespec_ilinestart+
     %                          linespec_nlines-1
                  write(*,*) 'Rendering position-velocity diagram ',
     %                       'of line ',iln
                  call calc_write_line_posvel(iln,linespec_distance,
     %                   main_passbwidth,main_passbnfr)
              enddo
#             else
              write(*,*) 'Line transfer is not active'
              stop 13
#             endif
          endif
          if(telesc_command.eq.3) then
c
c             The rectangular image
c              
#             ifndef RADGRID_ONEDIM
              write(*,*) 'Rendering the rectangular image (NEW)'
              call read_scatsource()
              call make_freq_image_rectang(image_ifreq)
              write(*,*) 'Writing the image to file'
              call write_image_rectang(image_ifreq)
#             endif
          endif
      endif
c
c     Write the finish log message to radlite.log
c
      call write_logmessage_end()
c
c     If requested, save the long characteristics to file
c
#     ifdef TELESCOPE_TRAJARRAY
      if(isave_long_char.eq.1) then
          write(*,*) 'Writing long characteristics to file'
          call write_trajects()
      endif
#     endif
c
  999 continue
c
      end




c     --------------------------------------------------------------
c                   WRITE A LOG MESSAGE TO A FILE
c
c     This includes some messages concerning the settings of the
c     defines, and the time of start and finish of a run.
c     All the messages are being appended to any previously 
c     existing file named radlite.log
c     --------------------------------------------------------------
      subroutine write_logmessage_start(action)
      implicit none
c
#include "common_grid.h"
#include "common_main.h"
#include "common_boundary.h"
#include "common_shortchar.h"
#include "common_isave.h"
c
      character*24 hpdate
      integer PGtime(3)
      character*80 dd,tt,action
      character*12 hh, mm, ss
c
      doubleprecision dummy,error
      integer iangset
c
#     ifdef COMPILER_HPUX
        call FDATE(hpdate)      
#     endif
#     ifdef COMPILER_GNU
        call date_and_time(dd,tt)
#     endif
#     ifdef COMPILER_INTEL
        call date_and_time(dd,tt)
#     endif
#     ifdef COMPILER_PG
        call date(dd)
        call itime(PGtime)
        write(hh,*) PGtime(1)
        write(mm,*) PGtime(2)
        write(ss,*) PGtime(3)        
#     endif
c
      open(unit=1,file='radlite.log',status='unknown',OPENACC='append')
      write(1,*)
      write(1,*)
      write(1,*)
      write(1,*) '=========== RUN OF RADICAL ==========='
      write(1,*) 'Action       : ',action(1:index(action,'.')-1)
#     ifdef COMPILER_HPUX
        write(1,*) 'Starting Date: ',hpdate(1:10),hpdate(20:24)
        write(1,*) '         Time: ',hpdate(12:19)
#     endif
#     ifdef COMPILER_GNU
        write(1,*) 'Starting Date: ',dd(7:8),'-',dd(5:6),'-',dd(1:4)
        write(1,*) '         Time: ',tt(1:2),':',tt(3:4),':',tt(5:6),
     %                '.',tt(8:10)
#     endif
#     ifdef COMPILER_INTEL
        write(1,*) 'Starting Date: ',dd(7:8),'-',dd(5:6),'-',dd(1:4)
        write(1,*) '         Time: ',tt(1:2),':',tt(3:4),':',tt(5:6),
     %                '.',tt(8:10)
#     endif
#     ifdef COMPILER_PG
        write(1,*) 'Starting Date: ',dd(1:20)
        write(1,*) '         Time: ',hh(11:12),':',mm(11:12),':',
     %              ss(11:12)
#     endif
      write(1,*) '......................................'
      write(1,*) '      COMPILATION OF CODE:            '
#ifdef RADGRID_ONEDIM     
      write(1,*) 'Dimension = 1              (RADGRID_ONEDIM)'
#else
      write(1,*) 'Dimension = 2              (RADGRID_TWODIM)'
#endif
#ifdef ONLY_RAY_TRACING
      write(1,*) 'Only raytracing active (no ESC/MESC!)'
#endif
#ifdef INCLUDE_LINES
      write(1,*) 'Include lines              (INCLUDE_LINES)'
#endif
#ifdef INCLUDE_DUST
      write(1,*) 'Include dust continuum     (INCLUDE_DUST)'
#endif
#ifdef INCLUDE_COMPTON
      write(1,*) 'Include comptonization     (INCLUDE_COMPTON)'
#endif
#ifdef ESC_MINIMAL_EXTENSION
      write(1,*) 'Transfer: MESC             (ESC_MINIMAL_EXTENSION)'
#endif
#ifdef ESC_COMPLETE
      write(1,*) 'Transfer: ESC              (ESC_COMPLETE)'
#endif
#ifdef LONGCHAR_EXTRA_POINTS
      write(1,*) 'Long char: refined         (LONGCHAR_EXTRA_POINTS)'
#else
      write(1,*) 'Long char: not refined'
#endif
cc#ifdef SMALL_MEMORY
cc      write(1,*) 'Memory usage : small       (SMALL_MEMORY)'
cc#else
cc      write(1,*) 'Memory usage : big'
cc#endif
#ifdef SHORCHAR_COMPACT_ARRAYS
      write(1,*) 'Compact ESC arrays         (SHORCHAR_COMPACT_ARRAYS)'
      write(1,*) '   Fraction of long ESC: ',LONGESC_FRACTION
#else
      write(1,*) 'Maximum ESC arrays'
#endif
#ifdef RADGRID_ONEDIM
#ifdef SCONEDIM_ADDPOINT
      write(1,*) 'SC: Add Pt at r=b (s=0)    (SCONEDIM_ADDPOINT)'
#else
      write(1,*) 'SC: No Pt added to SC'
#endif 
#endif 
#ifdef MIRROR_THETA
      write(1,*) 'Mirror in Theta : yes      (MIRROR_THETA)'
#else
      write(1,*) 'Mirror in Theta : no'
#endif
#ifdef MIRROR_PHI
      write(1,*) 'Mirror in phi : yes        (MIRROR_PHI)'
#else
      write(1,*) 'Mirror in phi : no'
#endif
#ifdef EQUATOR_DISK
      write(1,*) 'Equatorial disk : yes      (EQUATOR_DISK)'
#endif
#ifdef EQUATOR_DISK_INCLMU1
      write(1,*) 'Equatorial disk : incl mu1 (EQUATOR_DISK_INCLMU1)'
#endif
#ifdef INTERPOL_PHI_3
      write(1,*) 'Interpol                   (INTERPOL_PHI_3)'
#endif
#ifdef INTERPOL_THETA_3
      write(1,*) 'Interpol                   (INTERPOL_THETA_3)'
#endif
#ifdef INTERPOL_THETA_2
      write(1,*) 'Interpol                   (INTERPOL_THETA_2)'
#endif
#ifdef INTERPOL_RADIUS_2
      write(1,*) 'Interpol                   (INTERPOL_RADIUS_2)'
#endif
#ifdef INTERPOL_ORDER_1
      write(1,*) 'Interpol                   (INTERPOL_ORDER_1)'
#endif
#ifdef INTERPOL_SRC_0
      write(1,*) 'Interpol                   (INTERPOL_SRC_0)'
#endif
#ifdef INTERPOL_SRC_1
      write(1,*) 'Interpol                   (INTERPOL_SRC_1)'
#endif
#ifdef INTERPOL_SRC_2
      write(1,*) 'Interpol                   (INTERPOL_SRC_2)'
#endif
#ifdef INTERPOL_SRC_WIDE
      write(1,*) 'Interpol                   (INTERPOL_SRC_WIDE)'
#endif
#ifdef INTERPOL_FREQ_1
      write(1,*) 'Interpol                   (INTERPOL_FREQ_1)'
#endif
#ifdef INTERPOL_FREQ_3
      write(1,*) 'Interpol                   (INTERPOL_FREQ_3)'
#endif
#ifdef SRCQDR_PNT_2
      write(1,*) 'Quadrature                 (SRCQDR_PNT_2)'
#endif
#ifdef SRCQDR_PNT_3
      write(1,*) 'Quadrature                 (SRCQDR_PNT_3)'
#endif
#ifdef ALI_PNT_2
      write(1,*) 'ALI                        (ALI_PNT_2)'
#endif
#ifdef ALI_PNT_3
      write(1,*) 'ALI                        (ALI_PNT_3)'
#endif
      write(1,*) 
      write(1,*) 'FRSIZE_X       = ',FRSIZE_X    
      write(1,*) 'FRSIZE_Y       = ',FRSIZE_Y    
      write(1,*) 'FRSIZE_FREQ    = ',FRSIZE_FREQ 
      write(1,*) 'FRSIZE_MU      = ',FRSIZE_MU   
      write(1,*) 'FRSIZE_PHI     = ',FRSIZE_PHI  
      write(1,*) 'FRSIZE_CHAR    = ',FRSIZE_CHAR 
      write(1,*) 'FRSIZE_INF_PHI = ',FRSIZE_INF_PHI  
      write(1,*) 'FRSIZE_INF_R   = ',FRSIZE_INF_R    
      write(1,*) 'FRSIZE_INF_X   = ',FRSIZE_INF_X    
      write(1,*) 'FRSIZE_INF_Y   = ',FRSIZE_INF_Y    
      write(1,*) 
#ifdef INCLUDE_LINES
      write(1,*) 'SZ_LINEPROFILE       = ',SZ_LINEPROFILE       
      write(1,*) 'SZ_LINEPROFILE_SMALL = ',SZ_LINEPROFILE_SMALL 
      write(1,*) 'SZ_NLINES            = ',SZ_NLINES            
      write(1,*) 'SZ_NLEVELS           = ',SZ_NLEVELS           
      write(1,*) 'SZ_NCOLTRANS         = ',SZ_NCOLTRANS         
      write(1,*) 'SZ_NCTTEMP           = ',SZ_NCTTEMP           
      write(1,*) 
#endif
#ifdef INCLUDE_DUST
      write(1,*) 'DUST_SPECIES_MAX     = ',DUST_SPECIES_MAX
      write(1,*) 'DUST_SIZE_MAX        = ',DUST_SIZE_MAX    
      write(1,*) 'DUST_TRANGE_MAX      = ',DUST_TRANGE_MAX
      write(1,*) 'DUST_TEMP_MIN        = ',DUST_TEMP_MIN
      write(1,*) 'DUST_TEMP_MAX        = ',DUST_TEMP_MAX
      write(1,*) 'DUST_TEMP_ACCUR      = ',DUST_TEMP_ACCUR
#ifdef DUST_ALI
      write(1,*) 'DUST ALI             : installed'
#else
      write(1,*) 'DUST ALI             : disabled'
#endif
#endif
#ifdef INCLUDE_COMPTON
      write(1,*) 'ISOCSK_NGAM          = ',ISOCSK_NGAM1
#ifdef CSK_TAU_INTERPOL
      write(1,*) 'CSK_TAU_INTERPOL     = ',CSK_TAU_INTERPOL
#endif
#endif
      write(1,*) ''
      write(1,*) '......................................'
      write(1,*) '          SETUP OF PROBLEM:           '
c      write(1,*) 'max nr iter        = ',iter_trans
c      write(1,*) 'iter method        = ',iter_method
c      write(1,*) 'converg crit       = ',iter_convcrit
c      write(1,*) 'converg crit style = ',iter_convtype
c      write(1,*) 'initial guess type = ',iter_guess
      write(1,*) 'metric             = ',spacegrid_metric
      write(1,*) 'frsize radius      = ',irsi_frsizex 
      write(1,*) 'frsize theta       = ',irsi_frsizey
      write(1,*) 'nr mu              = ',anggrid_frsizemu
      write(1,*) 'nr phi             = ',anggrid_frsizephi
      write(1,*) 'mu-grid type       = ',anggrid_mu_type
      write(1,*) 'mu-grid allow err  = ',anggrid_muerr_max
      write(1,*) 'mu-grid coarseness = ',dmuvdr
      write(1,*) 'mu-grid extra mu=0 = ',iextrmu(1)
      write(1,*) 'bound cond inner   = ',iradbnd_in_itype
      write(1,*) 'bound cond outer   = ',iradbnd_out_itype
      write(1,*) 'telescope nphi     = ',main_nrphiinf
      write(1,*) 'telescope dbdr     = ',telesc_dbdr
      write(1,*) 'telescope extra b  = ',main_nrrayextra
      write(1,*) 'isave physvar      = ',isave_physvar
      write(1,*) 'isave scatnonlte   = ',isave_scatnonlte
      write(1,*) 'isave intensity    = ',isave_intens_inu
      write(1,*) 'isave source       = ',isave_source
      write(1,*) 'isave meanint      = ',isave_meanint
      write(1,*) 'isave dusttemp     = ',isave_dusttemp
      write(1,*) 'isave levpop       = ',isave_levpop
      write(1,*) 'isave electemp     = ',isave_electemp
      write(1,*) '......................................'
      write(1,*) '            SOME REMARKS:             '
      iangset=1
      dummy = rmu(nrmu(iangset)-1,iangset)
      error = abs(2.d0*abs(1.d0-dummy)/(1.d0-(dummy**2))-1.d0) 
      write(1,*) '    Max |mu|       : ',
     %     rmu(nrmu(iangset)-1,iangset),
     %     '  --> Error = ',error
      write(1,*) '    Min |mu|       : ',rmu(1,iangset)
#     ifdef SHORCHAR_COMPACT_ARRAYS
      write(1,*) 'Nr of short SC     : ',isc_shortsc_current
      write(1,*) 'Nr of extended SC  : ',isc_longesc_current
      write(1,*) 'Frac extended SC   : ',(isc_longesc_current+1.0)/
     %                                 isc_shortsc_current
#     endif
#     ifndef NO_SHORT_CHARS
      write(1,*) 'Greatest DR/R      : ',sc_rratiobig-1.d0
      write(1,*) 'Nr of R-refine     : ',sc_ndrr-1
#     endif
      write(1,*) '' 
      close(1)
c
      end


c     --------------------------------------------------------------
c                   WRITE A LOG MESSAGE TO A FILE
c
c     This includes some messages concerning the settings of the
c     defines, and the time of start and finish of a run.
c     All the messages are being appended to any previously 
c     existing file named radlite.log
c     --------------------------------------------------------------
      subroutine write_logmessage_end()
      implicit none
      character*24 hpdate
      integer PGtime(3)
      character*80 dd,tt
      character*12 hh, mm, ss
      logical yeswarnings
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
#include "common_isave.h"
#include "common_check.h"
#include "common_lines.h"
#include "common_dust.h"
#include "common_main.h"
c
      yeswarnings = .false.
c
#     ifdef COMPILER_HPUX
        call FDATE(hpdate)      
#     endif
#     ifdef COMPILER_GNU
        call date_and_time(dd,tt)
#     endif
#     ifdef COMPILER_INTEL
        call date_and_time(dd,tt)
#     endif
#     ifdef COMPILER_PG
        call date(dd)
        call itime(PGtime)
        write(hh,*) PGtime(1)
        write(mm,*) PGtime(2)
        write(ss,*) PGtime(3)        
#     endif

c
      open(unit=1,file='radlite.log',status='unknown',OPENACC='append')
      write(1,*) '......................................'
      write(1,*) 'Warnings:'
      if(warn_bug_centrbeam.ne.0) then
          write(1,*) ' * Uncorrected bug: the central beam bug in'
          write(1,*) '   the do_trans() routine...'
          yeswarnings = .true.
      endif
      if(warn_mubin_theta.ne.0) then
          write(1,*) ' * Mu-binning at Theta-crossing. Increase '
          write(1,*) '   the nr of mu-points to avoid this.'
          yeswarnings = .true.
      endif
      if(iscat_warn_int2.eq.1) then
          write(1,*) ' * Second order source integration not set!'
          write(1,*) '   Better recompile with SRCQDR_PNT_3 '
          yeswarnings = .true.
      endif
      if(iscat_warn_src2.eq.1) then
          write(1,*) ' * Second order source interpolation not set!'
          write(1,*) '   Better recompile with INTERPOL_SRC_2'
          yeswarnings = .true.
      endif
      if(iang_warn_mureslow.eq.1) then
          write(1,*) ' * The mu-gridding is insufficient!'
          write(1,*) '   Use larger NrMu!'
          yeswarnings = .true.
      endif
#     ifdef SAFETY_CHECKS_ACTIVE
      if(selfch_taumax.gt.SAFETY_CHECK_TAUMAX) then
          write(1,*) ' * Found rather high one-cell optical depths'
          write(1,*) '   Largest optical depth found is:',
     %         selfch_taumax
          yeswarnings = .true.
      endif
#     endif
#     ifdef INCLUDE_LINES
      if(iradproc_line.ne.0) then
      if(line_warn_cr_extrapol.ne.0) then
          write(1,*) ' * Collision rates were computed outside of'
          write(1,*) '   the temperature range of the tabulated     '
          write(1,*) '   lab-values. A blunt extrapolation was done.'
          yeswarnings = .true.
      endif
      if(line_warn_negpop.ne.0) then
          write(1,*) ' * Some level populations have become         '
          write(1,*) '   negative at least once during the transfer.'
          yeswarnings = .true.
      endif
      if(line_warn_umass.ne.0) then
          write(1,*) ' * The umass_av (average molec mass) seems to  '
          write(1,*) '   have been set previously to another value   '
          write(1,*) '   than is given by the line.inp file.         '
          yeswarnings = .true.
      endif
      if(line_warn_lineout.ne.0) then
          write(1,*) ' * Passband was too small to harbor the line   '
          write(1,*) '   profile at all width and velocity shifts.   '
          yeswarnings = .true.
      endif
      if(line_warn_badsampling.ne.0) then
          write(1,*) ' * Line profile frequency sampling is too      '
          write(1,*) '   coarse. Might result in errors.             '
          yeswarnings = .true.
      endif
      if(line_warn_velostep.eq.1) then
          write(1,*) ' * Velocity step exceeded line width.          '
          write(1,*) '   Might result in errors.                     '
          yeswarnings = .true.
      endif
      if(line_warn_velostep.eq.2) then
          write(1,*) ' *        >>>>> RED ALTERT <<<<<<<<            '
          write(1,*) '   Velocity step MUCH larger than line width!  '
          write(1,*) '   Big errors are expected!                    '
          yeswarnings = .true.
      endif
#ifdef LINE_LOCAL_ABUNDANCE
      if(locabun_done.eq.0) then
          write(1,*) ' * Abundances are taken globally (value from   '
          write(1,*) '   the line.dat file).                         '
          yeswarnings = .true.
      endif
#endif
#ifndef LINE_NO_LIMIT_MASERS
      if(line_warn_opacity_neg.ne.0) then
          write(1,*) ' * MAJOR WARNING:                              '
          write(1,*) '   In last iteration step, the opacity became  '
          write(1,*) '   negative at some point. I have bluntly put  '
          write(1,*) '   this to zero!! I switched masers off!       '
          yeswarnings = .true.
      endif
#endif
      if((line_warn_opacity_neg.eq.0).and.
     %   (line_warn_neg_opacity.ne.0)) then
          write(1,*) ' * Negative opacity was detected during the    '
          write(1,*) '   iteration procedure. Fortunatly not at the  '
          write(1,*) '   final iteration step.                       '
          yeswarnings = .true.
      endif
      endif 
#endif
#ifdef INCLUDE_DUST
      if(dust_warn_few_freqs.ne.0) then
          write(1,*) ' * Very few frequencies for dust temperature   '
          write(1,*) '   solver. Serious errors expected.            '
          yeswarnings = .true.
      endif
#endif
#ifdef INCLUDE_COMPTON
      if(compt_warn_compt1.ne.0) then
          write(1,*) ' * Compton switch irradproc_compt is set       '
          write(1,*) '   to 1. This is only Thompson scattering.     '
          write(1,*) '   It is okay, but just remember that for      '
          write(1,*) '   isotropic CSK you should put it to 3,       '
          write(1,*) '   and for diffusion equation to 2.            '
          yeswarnings = .true.
      endif
#endif
      if(transerr_overshoot.gt.0) then
          write(1,*) ' * During Short Characteristic integration     '
          write(1,*) '   there were overshoots in the quadrature.    '
          write(1,*) '   These have been corrected automatically,    '
          write(1,*) '   but they indicate poor grid resolution!     '
          yeswarnings = .true.
      endif
      if(.not.yeswarnings) then
          write(1,*) ' <none>'
      endif
      write(1,*) '......................................'
c$$$      if(iter_trans.gt.0) then
c$$$          write(1,*) 'Number of iterations done: ',iter_done
c$$$          if(iter_converged) then
c$$$              write(1,*) 'Converged:                 yes'
c$$$          else
c$$$              write(1,*) 'Converged:                 no'
c$$$          endif
c$$$      endif
#     ifdef COMPILER_HPUX
        write(1,*) 'Finishing Date: ',hpdate(1:10),hpdate(20:24)
        write(1,*) '          Time: ',hpdate(12:19)
#     endif
#     ifdef COMPILER_GNU
        write(1,*) 'Finishing Date: ',dd(7:8),'-',dd(5:6),'-',dd(1:4)
        write(1,*) '          Time: ',tt(1:2),':',tt(3:4),':',tt(5:6),
     %                '.',tt(8:10)
#     endif
#     ifdef COMPILER_INTEL
        write(1,*) 'Finishing Date: ',dd(7:8),'-',dd(5:6),'-',dd(1:4)
        write(1,*) '          Time: ',tt(1:2),':',tt(3:4),':',tt(5:6),
     %                '.',tt(8:10)
#     endif
#     ifdef COMPILER_PG
        write(1,*) 'Finishing Date: ',dd(1:20)
        write(1,*) '         Time: ',hh(11:12),':',mm(11:12),':',
     %              ss(11:12)
#     endif
      write(1,*) '========= END OF RUN RADLITE ========='
      close(1)
c
      end


