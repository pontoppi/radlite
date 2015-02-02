c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c                     ( MODULE: LINE TRANSFER PACKAGE ) 
c
c                            Leiden, Nov 1998
c                   C.P. Dullemond  /  Gerd-Jan v. Zadelhoff
c
c     Part of this file was copied from the program line.F in the
c     directory Spec2D/Line/line.F_14.11.98, but they were modified.
c     In fact, the later version of Spec2D/Line/line.F also includes
c     some of the  changes made here! So effectively this is built 
c     from: Spec2D/Line/line.F_14.11.98_B
c
c     IMPORTANT NOTE: 
c       This module does the line transfer in two different ways:
c
c       'Simple mode':
c       --------------
c       The default method is a simple method, in which the line 
c       profile and abundances are global, and any systematic 
c       velocities of the underlying medium are excluded.
c
c       'Expert mode':
c       --------------
c       The expert-method (activated when LINE_VELOCITIES is defined)
c       includes local line profiles (if defined) and allows for
c       systematic local velocities, and thereby angle--dependent
c       source terms and line profiles. This expert-method in many
c       ways  circumvents the usual philosophy of the transfer
c       algorithm, which uses angle-averaged source terms and 
c       opacities, and which stores the mean intensity only. But
c       when LINE_VELOCITIES is set, various other switches in the
c       other modules are automatically switched on as well, such
c       as the declaration of a big array for storage of the angle-
c       dependent intensity, and the computation of the angle-
c       dependent source terms and opacities. The switches affected
c       by defining LINE_VELOCITIES are:
c         - INTENSITY_ARRAY
c         - INTARRAY_FREQ
c         - PHYSVAR_VELOCITIES
c         - LINE_LOCAL_WIDTH
c         - LINE_LOCAL_ABUNDANCE
c       It should be noted that when LINE_VELOCITIES is set, the
c       program absorbs an enormous chunk of memory, because the
c       angle--dependent intensities have to be stored. So unless
c       it is absolutely necessary to include velocities and/or
c       variable line widths, it is advised NOT to activate this 
c       expert mode. 
c
c       Once LINE_VELOCITIES is set, you still have the option to
c       switch off the LINE_LOCAL_WIDTH. This may also save memory.
c       If the LINE_LOCAL_WIDTH is activated, the line profile is
c       local, and it will be computed and stored in a big array
c       before the transfer. This may also eat memory, so it may
c       be switched off. But this is a slightly less standard 
c       option, so if you wish to switch it off, then please do
c       so in the main.h file, where the LINE_VELOCITIES-dependent
c       macros are set.
c
c       Also, the LINE_VELOCITIES allows space-dependent abundances.
c       This is regulated by the macro LINE_LOCAL_ABUNDANCE, which
c       is automatically switched on when LINE_VELOCITIES is switched
c       on. 
c
c     NOTE: Whereever there is >>>> ADDED KEES/KLAUS MARCH 2007 <<<<
c           This is where surgical modifications are made for the
c           more stable line-integration routines for very narrow
c           line widths, and for getting a more accurate (i.e.
c           freq-dependent) dust continuum.
c
c     =============================================================

c#define SZ_LINEPROFILE 40
c#define SZ_NLINES 110
c#define SZ_NLEVELS 110
c#define SZ_NCOLTRANS 500
c#define SZ_NCTTEMP 20

#define TAU_EPS 1.d-4

#define LINE_F_VERSION "06-11-05"
#include "main.h"


#ifdef INCLUDE_LINES  /* Compile only if necessary */

c     --------------------------------------------------------------
c                     READ THE LINE DATA FILE
c
c     This routine reads the information about the lines and the
c     levels. It reads the molecular data, sets up the line profiles
c     (if fixed for the global space) and assigns transfer freq-bins
c     to the line freq-bins.
c
c     In the file line.inp one finds the necessary information
c     about the lines. It contains the following data (roughly):
c
c        How many species of molecules are there
c        Average atomic weight (2.d0 if H2)
c        =======================================
c        File name of species level data
c        Lower level to take into account
c        Upper level to take into account
c        Abundance of molecule (n(molec)/n(h2))   Now disabled, KMP 18/8/07
c        Abundance of collision partner (n(molec)/n(h2))  Now disabled ---
c        ---------------------------------------
c        Specify info for first line? 0=no,take-default,1=yes. If 1, then:
c        Symmetry in the line (0=nosymm,1=symm)
c        Number of frequency points in line
c        Width of line range from center to edge (both sides) in km/s
c        Line profile (0=block,1=Gauss,2=Voigt)
c        Additional information: (0=nothing,1=fixed-line-width) if >0 then:
c        Fixed line width (km/s)
c        Additional information: (0=nothing) if =0 then stop, next line
c        ---------------------------------------
c        Specify info for second line? 0=no,take-default,1=yes. If 1, then:
c             ................etc...etc..........
c
c     This routine automatically reads the appropriate molecular 
c     data file (for instance co.dat). 
c
c     This routine can also be asked to read only the main information
c     from the line.inp file, and skip all the setting-up of the
c     frequencies / reading rate files etc. This is useful if you 
c     have a problem setup that needs to know e.g. which molecule is
c     going to be modeled. You can then call this routine with argument
c     1, wich reads the first couple of lines of the line.inp file,
c     without messing anythnig up. You may call read_linedata(1) many
c     times in a row. But only call read_linedata(0) ONCE!
c
c     WARNING: For now only one species is allowed at any given run. 
c     
c
c     ARGUMENTS:
c       ionlyhead      = 0 --> Read all data, including rate files
c                      = 1 --> Only read the general information 
c                              from the line.inp file
c
c     --------------------------------------------------------------
      subroutine read_linedata(ionlyhead)
      implicit none
c
      integer ionlyhead
c
#include "common_grid.h"
#include "common_source.h"
#include "common_lines.h"
c
      character*80 comstring,string
      integer nrspecies,gotit,iinfo
      integer ir,it,ithetamax,iformat,istyle
      integer i_ground_level,i_top_level,imore,i,k,itr_freq
      doubleprecision phi(SZ_LINEPROFILE),dnu(SZ_LINEPROFILE)
      doubleprecision wgt(SZ_LINEPROFILE)
      doubleprecision dum1,dum2,umass_orig
c
c     NEW: 30.03.06 REMOVED FREQ_NR=0 HERE
c
c     Very perhaps the umass_av is already set accidently by the
c     medium.F module, to convert number densities into mass density.
c     If so, detect this here, and warn that this value may have 
c     changed.
c
      umass_orig = umass_av
c
      open(unit=3,file='line.inp',status='unknown',err=980)
c
c     Read how many species (must be 1 for the current version!)
c     NEW: 06.04.06
c
      read(3,*) iformat
      if(iformat.eq.1) then
         nrspecies = iformat
      elseif(iformat.eq.2) then
         read(3,*) nrspecies
      else
         stop 23223
      endif
c
c     Read the average molecular mass in units of m_u
c
      read(3,*) umass_av
      if((umass_orig.ne.0.d0).and.(umass_av.ne.umass_orig)) then
          write(*,*) 'WARNING: the umass_av (average molec mass) '
          write(*,*) '         seems to have been set previously '
          write(*,*) '         to another value than is given by '
          write(*,*) '         the line.inp file. Is this okay? '
          write(*,*) 'Original umass_av:', umass_orig
          write(*,*) 'New umass_av:', umass_av
          
          line_warn_umass = 1
      endif
c
c     Read comment line
c
      read(3,*) comstring

c
c     Now we have read all the general information. If ionlyhead.eq.1
c     then we are done now. If ionlyhead.eq.0 then let us read the
c     molecular rotational/vibrational transition rate files as well
c     and set up the frequency bins etc. etc.
c
      if(ionlyhead.eq.1) then
          close(3)
          goto 999
      endif
c
c     Now read the information for the various lines
c
c      call read_species_lambda(molec_filename,nlev_orig)
      do i=1,nlines
c
c         Read comment line
c
          read(3,*) comstring
c
c         Find out if line info should be read, or default
c         settings are taken. For now we always read line info.
c
          read(3,*) iinfo
          if((iinfo.ne.1).and.(i.eq.1)) then
              write(*,*) 'Default lines not installed yet'
              stop 13
          endif
          if(iinfo.eq.0) then
              line_symm(i)          = line_symm(i-1)
              line_nfreq(i)         = line_nfreq(i-1)
              line_rangewidth(i)    = line_rangewidth(i-1)
              line_profiletype(i)   = line_profiletype(i-1)
              line_fix_linewidth(i) = line_fix_linewidth(i-1)              
              goto 30
          endif
c
c         Symmetry? (for now no symmetry)
c   
          read(3,*) line_symm(i)
          if(line_symm(i).ne.0) then
              write(*,*) 'Symmetry in lines not yet implemented'
              stop 13
          endif
c
c         Read the number of freq points in line
c
          read(3,*) line_nfreq(i)
c
c         Read width of freq domain around line center
c
          read(3,*) line_rangewidth(i)
c
c         Read profile type for this line
c
          read(3,*) line_profiletype(i)
c
c         Is there additional information?
c
  111     continue
          read(3,*) imore
          if(imore.eq.0) then
              goto 112
          endif
          if(imore.eq.1) then
c
c             Line width is fixed. Read the line width.
c
              read(3,*) line_fix_linewidth(i)
          else
              write(*,*) 'Dont know the line info type ',imore
              stop 13
          endif
          goto 111
  112     continue
c
c         Done reading
c
   30     continue
      enddo  ! NEW: 27.03.06: Cut loop here 
c
c     NEW: 27.03.06: Find the location of the line centers in the
c                    continuum frequency array (which will be overwritten
c                    later, in case of line transfer), and put the
c                    continuum alpha and source terms also already in the
c                    appropriate places...
c

      if(iradproc_line_dust.ne.0) then
          call global_prepare_line_dust()
      endif
c
c     Now make all the passbands 
c
c     IMPORTANT NOTE: From here on the dust source terms etc are 
c                     going to be  over-written!
c
c     NEW: 27.03.06: New loop here. Also put itr_freq = 1 here instead above.
c
c     Set freq_nr = 0, because now that lines are present, we
c     dictate the amount of frequencies here from the line module.
c
      freq_nr  = 0
      itr_freq = 1
      do i=1,nlines
c
c         Set up the passband for this line. This routine will 
c         set the values of the line_nu,line_dnu,line_wgt arrays
c         and it will automatically attach the line-frequency bins
c         to the global frequency bins (freq_nu() etc). The counter
c         itr_freq will be updated to the next free global bin.
c         Also it will automatically update freq_nr, i.e. the global
c         number of frequencies.
c
ccc          call line_setup_passband(i,line_rangewidth(i),
ccc     %                             line_nfreq(i),itr_freq)
c
c         Now construct the fixed line profile, if a fixed line profile
c         is demanded. 
c
#         ifndef LINE_LOCAL_WIDTH
              call make_voigt_profile(linefreq(i),0.d0,
     %                line_fix_linewidth(i),line_rangewidth(i),
     %                line_nfreq(i),dnu,phi,wgt)
c
c             Copy this info in the relevant arrays
c
              do k=1,line_nfreq(i)
ccc                 line_dnu(k,i) = dnu(k)
ccc                 line_nu(k,i)  = linefreq(i) + dnu(k) 
ccc                 line_wgt(k,i) = wgt(k)
                  line_phi(k,i) = phi(k)
              enddo
#         endif
      enddo
c
      close(3)


c
c     Save the value of Aud to field array. This is not necessary, but 
c     just for safety. Aud is a constant and is therefore global, but 
c     just for the sake of clarity, ALL Einstein coefficients are made 
c     local. And in MALI the Aud() arrays are molested, so we must be
c     sure that the Aud_field() array is correct!!
c
#ifdef MIRROR_THETA
      ithetamax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      ithetamax = 1
#endif
c
      do ir=1,irsi_frsizex
          do it=1,ithetamax
              call save_transinfo(ir,it)
          enddo
      enddo
c
c     The end
c
      goto 999
  980 continue
      write(*,*) 'ERROR: Could not open file line.inp'
      stop 13
  999 continue
      return
      end

c     --------------------------------------------------------------
c                 READ THE NON-LTE DATA (level populations)
c
c     This is the main routine that reads the results of the non-lte
c     transfer in. 
c     --------------------------------------------------------------
      subroutine read_nonlte_scat()
      implicit none
      logical fexists
      integer iter,iconv,isolvt
c
#include "common_source.h"
#include "common_main.h"
c
      if(iscat_done_nonlte.eq.0) then
          write(*,*) 'Reading the level populations from file'
          if(iradproc_line.ne.0) then
c         
c             If line transfer, then do not read in the source functions.
c             Instead, read in the level populations. Right now, line
c             transfer is incompatible with any other processes.
c              
#             ifdef INCLUDE_LINES
              call read_levelpopul(-10)
              iscat_done_nonlte = 1
#             else
              write(*,*) 'LINE TRANSFER DISABLED'
              stop 13
#             endif
          else
          endif
      endif

      end



c     --------------------------------------------------------------
c                    PRODUCE A PASSBAND FOR A LINE
c
c     This routine sets up the frequencies for the line profile
c     and attaches these frequencies to the global frequency bins
c     of the radiative transfer. It does not set up a line profile.
c     It should be noted that this routine overrides any possible 
c     original setup of the global frequencies, and the stellar
c     spectrum. It basically confiscates these global frequency
c     bins. Since the program is not meant to be doing various types
c     of RT simultaneously, this should be no problem.
c
c       iline            = Index of the line 
c       width [km/s]     = Width of the passband
c       nfr              = Amount of frequencies for this line
c       ifreq            = Starting index for this line in global 
c                          frequency array. This will be updated 
c                          to the new value.
c
c     --------------------------------------------------------------

      subroutine line_setup_passband(iline,wwwidth,nfr,ifreq)
      implicit none

      doubleprecision wwwidth
      integer iline,nfr,ifreq

      doubleprecision passb,nu1,nu2,dnu,nu0,width
      integer k
c
#include "common_grid.h"
#include "common_lines.h"
#include "common_boundary.h"
c
      width = wwwidth
      nu0 = abs(linefreq(iline))
c
c     Consistency checks
c
      if(nu0.eq.0.d0) then
          write(*,*) 'Problem in line_setup_passband(): nu0=0 !'
          stop 13
      endif
      if((nfr.gt.SZ_LINEPROFILE).or.(nfr.lt.1)) then
          write(*,*) 'Number of frequencies for this line is ',
     %          'out of range: ',nfr
          write(*,*) 'Perhaps recompile with larger SZ_LINEPROFILE'
          stop 13
      endif
      if(nfr.eq.1) then
          write(*,*) 'ERROR: Simple square line profile deactivated'
          stop 13
      endif
c     
c     First convert the line widths from km/s to Hz
c
      passb  = 3.33567d-6 * nu0 * width    ! in Hz
      nu1    = 0.d0 - passb                ! nu0 already subtracted
      nu2    = 0.d0 + passb                ! nu0 already subtracted
      dnu    = 2.d0*passb / (nfr-1.d0)
      do k=1,nfr
          line_dnu(k,iline)  = nu1 + (k-1)*dnu
          line_nu(k,iline)   = nu0 + line_dnu(k,iline) 
      enddo
c
c     Determine the weights very crudely for now: equal weights for all
c
      do k=2,nfr-1
          line_wgt(k,iline) = dnu
      enddo
      line_wgt(1,iline)   = dnu
      line_wgt(nfr,iline) = dnu
c
c     Assign linefreq to line_nu0 (don't ask me why I did this so
c     unhandy...) 
c
      line_nu0(iline)     = linefreq(iline)
c
c     Just in case, make sure that the line_rangewidth is filled.
c
      line_rangewidth(iline)  = width
c
c     Tell the line module about how many line-fre-bins:
c
      line_nfreq(iline) = nfr
c
c     Now assign a set of transfer frequency bins to the line frequency
c     bins. This means that each frequency point in the line actually
c     is represented by an actual radiation field at this frequency.
c     
      do k=1,nfr
c     
c         First check if there are still radiation fields left
c     
          if(ifreq.gt.FRSIZE_FREQ) then
              write(*,*) 'Problem: I dont have enough radiation ',
     %            'field frequency bins to harbor all the line ',
     %            'frequencies'
              write(*,*) 'line = ',iline,' ifr = ',k
              write(*,*) 'Recompile with larger FRSIZE_FREQ'
              stop 13
          endif
c     
c         Update freq_nr
c     
          freq_nr = ifreq
c     
c         Now assign this line freq to a radiation field
c     
          line_trifreq(k,iline) = ifreq
c     
c         ...and tell the radiation field to which line and to
c         which freq in that line it has been assigned.
c     
          line_iline(ifreq)   = iline
          line_ilifreq(ifreq) = k
c     
c         Copy the frequency into the freq_nu array (the global
c         frequency array)
c     
          freq_nu(ifreq) = line_nu(k,iline)
c
c         For safety reasons I set the stellar spectrum to 0,
c         since we're overriding the original frequency setup.
c
c          write(*,*) 'Warning: central star is added, '
c          write(*,*) 'but not fully tested!'
c          radbnd_starspec(ifreq) = 0
c     
c         Increase the radiation field frequency counter
c     
          ifreq = ifreq + 1
      enddo
c
c     NEW 27.03.06: For combining dust and lines
c     
      iline_freq_nr = freq_nr
c    
      return
      end





c     --------------------------------------------------------------
c                  SAVE LEVEL POPULATIONS INTO ARRAYS
c     --------------------------------------------------------------
      subroutine save_levelpop(ir,it)
      implicit none
      integer ir,it,ilev
c
#include "common_lines.h"
c
      do ilev=1,nlevels
          line_level_popul(ilev,it,ir) = npopul(ilev) 
      enddo
c
      end

c     --------------------------------------------------------------
c                  COPY LEVEL POPULATIONS INTO ARRAYS
c     --------------------------------------------------------------
      subroutine copy_levelpop(ir,it)
      implicit none
      integer ir,it,ilev
c
#include "common_lines.h"
c
      do ilev=1,nlevels
          npopul(ilev) = line_level_popul(ilev,it,ir)
      enddo
c
      end
      

c     --------------------------------------------------------------
c                    SAVE LINE/COLTRANS INFO TO ARRAYS
c
c     Since the lines are handled locally (solving the NLTE stuff
c     is done without reference to the place on the grid), we must
c     copy this information to arrays if we want to write this info
c     to a file later (for experimentation purposes). 
c     --------------------------------------------------------------
      subroutine save_transinfo(ir,it)
      implicit none
      integer ir,it,iline,ictrans
c
#include "common_lines.h"
c
      do iline=1,nlines
          Aud_field(iline,it,ir) = Aud_const(iline)
          Bud_field(iline,it,ir) = Bud(iline)
          Bdu_field(iline,it,ir) = Bdu(iline)
      enddo
      do ictrans=1,ncoltrans
          Cud_field(ictrans,it,ir) = Cud(ictrans)
          Cdu_field(ictrans,it,ir) = Cdu(ictrans)
      enddo
c
      end


c     --------------------------------------------------------------
c             COPY LINE/COLTRANS INFO TO LOCAL DATA
c     --------------------------------------------------------------
      subroutine copy_transinfo(ir,it)
      implicit none
      integer ir,it,iline,ictrans
c
#include "common_lines.h"
c
      do iline=1,nlines
          Aud(iline) = Aud_field(iline,it,ir) 
          Bud(iline) = Bud_field(iline,it,ir) 
          Bdu(iline) = Bdu_field(iline,it,ir) 
      enddo
      do ictrans=1,ncoltrans
          Cud(ictrans) = Cud_field(ictrans,it,ir)
          Cdu(ictrans) = Cdu_field(ictrans,it,ir)
      enddo
c
      end


c     --------------------------------------------------------------
c                     CHECK LEVELS FOR CONSISTENCY
c     
c     This routine can be called at the beginning of the iteration.
c     It will check whether the level populations add up to 1. If
c     not, it will warn, and insert the local LTE values.
c     --------------------------------------------------------------
      subroutine check_levelpopul()
      implicit none
c
#include "common_grid.h"
#include "common_lines.h"
c
      integer ir,it,itmax,imirt,i,iflag
      doubleprecision dummy
      character*80 filename
      character*12 ch
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
      imirt = 0
#endif
c
c     First a loop over all the grid positions
c
      iflag = 0
      do ir=1,irsi_frsizex
          do it=1,itmax
              dummy = 0.0d0
              do i=1,nlevels
                  dummy = dummy + line_level_popul(i,it,ir)
              enddo
              if(abs(1.d0-dummy).gt.1.e-2) then
                  iflag=1
              endif
          enddo
      enddo
      if(iflag.eq.1) then
          write(*,*) 'WARNING: The level populations are not normalized'
          write(*,*) '=========> Inserting LTE values instead'
          call calc_lte_levelpopul()
      endif
c
      end      


c     --------------------------------------------------------------
c               WRITE THE OCCUPATION NUMBERS TO FILE
c
c     This is the routine that dumps the level populations to a
c     file. The filenr specification denotes the file:
c 
c       filenr: -1 = Save to levelpop.dat
c               -2 = Save to levelpop_final.dat
c               -3 = Save to levelpop_<molecname>.dat
c                0,1,2... = Save to intmean_0.dat etc...
c
c     --------------------------------------------------------------
      subroutine write_levelpopul(filenr)
      implicit none
      integer filenr
c
#include "common_grid.h"
#include "common_lines.h"
c
      integer ir,it,itmax,imirt,i
      integer dum_freq_nr,dum_frsizex,dum_itmax,dum_imirt
      character*80 filename,frmt
      character*12 ch
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
      imirt = 0
#endif
c
c     Write the info file
c
      open(unit=1,file='levelpop.info',status='unknown')
      write(1,*) filenr
      if(filenr.eq.-2) then
          write(1,*) line_pop_lastsave
      endif
      if(filenr.eq.-3) then
          write(1,*) 'levelpop_'//molec_name(1:molec_name_len)//'.dat'
          write(1,*) line_pop_lastsave
      endif
      close(1)
c
c     Open the file and write header
c
      if(filenr.eq.-1) then
         filename = 'levelpop.dat'
      elseif(filenr.eq.-2) then
         filename = 'levelpop_final.dat'
      elseif(filenr.eq.-3) then
         filename = 'levelpop_'//molec_name(1:molec_name_len)//'.dat'
      else
         write(ch,*) filenr
#ifdef COMPILER_PG
         if(filenr.lt.10) then
            filename = 'levelpop_'//ch(12:12)//'.dat'
         elseif(filenr.lt.100) then
            filename = 'levelpop_'//ch(11:12)//'.dat'
         elseif(filenr.lt.1000) then
            filename = 'levelpop_'//ch(10:12)//'.dat'
         endif
#else
         if(filenr.lt.10) then
            filename = 'levelpop_'//ch(2:2)//'.dat'
         elseif(filenr.lt.100) then
            filename = 'levelpop_'//ch(2:3)//'.dat'
         elseif(filenr.lt.1000) then
            filename = 'levelpop_'//ch(2:4)//'.dat'
         endif
#endif
         line_pop_lastsave = filenr
      endif
      write(ch,*) nlevels
#ifdef COMPILER_PG
      if(nlevels.lt.10) then
          frmt = '('//ch(12:12)//'(E13.6,1X))'
      elseif(nlevels.lt.100) then
          frmt = '('//ch(11:12)//'(E13.6,1X))'
      elseif(nlevels.lt.1000) then
          frmt = '('//ch(10:12)//'(E13.6,1X))'
      endif
#else
      if(nlevels.lt.10) then
          frmt = '('//ch(2:2)//'(E13.6,1X))'
      elseif(nlevels.lt.100) then
          frmt = '('//ch(2:3)//'(E13.6,1X))'
      elseif(nlevels.lt.1000) then
          frmt = '('//ch(2:4)//'(E13.6,1X))'
      endif
#endif
      open(unit=1,file=filename,status='unknown')
      write(1,*) irsi_frsizex,itmax,nlevels,imirt
c
c     First write all the line frequencies
c
      write(1,*) (enerlevel(i),i=1,nlevels)
      write(1,*) (gdeg(i),i=1,nlevels)
      write(1,*) ' '      
c
c     First a loop over all the grid positions
c
      do ir=1,irsi_frsizex
          do it=1,itmax
              write(1,frmt) (line_level_popul(i,it,ir),i=1,nlevels)
          enddo
      enddo
c
      close(1)
c
      return
c
      end

c     --------------------------------------------------------------
c               READ THE OCCUPATION NUMBERS FROM FILE
c
c     --------------------------------------------------------------
      subroutine read_levelpopul(filenrinp)
      implicit none
      integer filenrinp
c
#include "common_grid.h"
#include "common_lines.h"
c
      integer ir,it,itmax,imirt,i,imt,itm,nr,filenr,gotit,ilen
      integer dum_freq_nr,dum_frsizex,dum_itmax,dum_imirt
      character*80 filename      !,frmt
      character*80 molnametmp,string
      character*12 ch
      doubleprecision enl(SZ_NLEVELS),gg(SZ_NLEVELS),dum
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
      imirt = 0
#endif
c
c     If filenr lt -5, then this means that the user want
c     simply to get the most recent level populations file
c
      if(filenrinp.lt.-5) then 
          open(unit=1,file='levelpop.info',status='old')
          read(1,*) filenr
          if(filenr.eq.-3) then
              read(1,*) string
              if(string(1:9).ne.'levelpop_') then
                  write(*,*) 'ERROR: while attempting to read ',
     %                'levelpop.info'
                  stop 13
              endif
              do i=10,80-4
                  if(string(i:i+3).eq.'.dat') then
                      molnametmp=string(10:i-1)
                      gotit = 1
                      ilen  = i-10+1
                  endif
              enddo
              if(gotit.eq.0) then
                  write(*,*) 'Could not make sense of string ',string
                  write(*,*) 'while reading the levelpop.info file'
                  write(*,*) 'File name must end with .dat'
                  stop 13
              endif
          endif
          close(1)
      else
          filenr=filenrinp
      endif
c
c     Open the file and read header
c
      if(filenr.eq.-1) then
         filename = 'levelpop.dat'
      elseif(filenr.eq.-2) then
         filename = 'levelpop_final.dat'
      elseif(filenr.eq.-3) then
         filename = 'levelpop_'//molec_name(1:molec_name_len)//'.dat'
      else
         write(ch,*) filenr
#ifdef COMPILER_PG
         if(filenr.lt.10) then
            filename = 'levelpop_'//ch(12:12)//'.dat'
         elseif(filenr.lt.100) then
            filename = 'levelpop_'//ch(11:12)//'.dat'
         elseif(filenr.lt.1000) then
            filename = 'levelpop_'//ch(10:12)//'.dat'
         endif
#else
         if(filenr.lt.10) then
            filename = 'levelpop_'//ch(2:2)//'.dat'
         elseif(filenr.lt.100) then
            filename = 'levelpop_'//ch(2:3)//'.dat'
         elseif(filenr.lt.1000) then
            filename = 'levelpop_'//ch(2:4)//'.dat'
         endif
#endif
      endif
      open(unit=1,file=filename,status='old',err=980)
      read(1,*) nr,itm,nlevels,imt
      if((nr.ne.irsi_frsizex).or.(itm.ne.itmax).or.
     %   (imt.ne.imirt)) then
          write(*,*) 'File ',filename,' has inconsistent dimensions'
          write(*,*) nr,irsi_frsizex,itm,itmax,imt,imirt
          stop 13
      endif
c
c     First read dummy stuff
c
      read(1,*) (enl(i),i=1,nlevels)
      read(1,*) (gg(i),i=1,nlevels)
c
c     First a loop over all the grid positions
c
      do ir=1,irsi_frsizex
          do it=1,itmax
              read(1,*) (line_level_popul(i,it,ir),i=1,nlevels)
          enddo
      enddo
c
c     NEW: 30.03.06: Check that the level populations add up to 1
c
      do ir=1,irsi_frsizex
          do it=1,itmax
              dum=0.d0
              do i=1,nlevels
                  dum = dum + line_level_popul(i,it,ir)
              enddo
c              if(abs(dum-1.d0).gt.1d-5) then
c                  write(*,*) 'ERROR: Level populations do not',
c     %             ' add up to 1.'
c                  write(*,*) 'ir,it = ',ir,it
c                  write(*,*) 'sum   = ',dum
c     Disabled this stop to account for the use of a proper partition function
c     when doing LTE calculations, KMP 16/11/06
c                  stop 91298
c              endif
          enddo
      enddo
c
      close(1)
c
      goto 999
  980 continue
      write(*,*) 'ERROR: Could not open file ',filename
      stop 13
  999 continue
      return
c
      end



c     --------------------------------------------------------------
c                             WRITE LINE INFO
c
c     This routine writes the values of various local line-related
c     quantities. Make sure that the transfer has been done and that
c     the line processes have been done as well.
c
c            (Bud*Jbar)/Aud      = Relative importance of stim-emission
c            Cud/Aud             = Relative importance of collisions 
c            epsilon             = Cul/(Aul+BulJbar+Cul)
c            bar alpha           = Average opacity
c            bar alpha * DR      = Average optical depth of the grid cell
c            T_ex                = Excitation temperature
c            T_ex/T              = How much is it in LTE?
c
c     --------------------------------------------------------------
      subroutine write_lineinfo(filenr)
      implicit none
      integer filenr
c
#include "common_grid.h"
#include "common_lines.h"
#include "common_setup.h"
c
      integer ir,it,itmax,imirt,i,ilev_up,ilev_down,iline,ifr,idum
      doubleprecision a1,a2,a3,a4,a5,a6,a7,a8,npop_up,npop_down,temp
      doubleprecision dfrkms(1:SZ_LINEPROFILE),ddfr
      character*80 filename,frmt
      character*12 ch
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
      imirt = 0
#endif
c
c     Open the file and write header
c
      if(filenr.lt.0) then
         filename = 'lineinfo.dat'
      else
#ifdef COMPILER_PG
         if(filenr.lt.10) then
            filename = 'lineinfo_'//ch(12:12)//'.dat'
         elseif(filenr.lt.100) then
            filename = 'lineinfo_'//ch(11:12)//'.dat'
         elseif(filenr.lt.1000) then
            filename = 'lineinfo_'//ch(10:12)//'.dat'
         endif
#else
         if(filenr.lt.10) then
            filename = 'lineinfo_'//ch(2:2)//'.dat'
         elseif(filenr.lt.100) then
            filename = 'lineinfo_'//ch(2:3)//'.dat'
         elseif(filenr.lt.1000) then
            filename = 'lineinfo_'//ch(2:4)//'.dat'
         endif
#endif
      endif
      open(unit=1,file=filename,status='unknown')
      write(1,*) irsi_frsizex,itmax,nlines,imirt
c
c     First write all the line frequencies
c
      frmt = '(4(E11.4,1X))'
      write(1,frmt) (linefreq(iline),iline=1,nlines)
c
c     Then write all the frequency sample points in offset of the
c     line center in km/s
c
      idum=0
      do iline=1,nlines
          if(line_nfreq(iline).gt.idum) then
              idum = line_nfreq(iline)
          endif
      enddo
      write(1,*) idum
      do iline=1,nlines
          ddfr      = 2.d0*line_rangewidth(iline)/
     %                    (line_nfreq(iline)-1.d0)
          dfrkms(1) = -line_rangewidth(iline)
          do ifr=2,line_nfreq(iline)
              dfrkms(ifr) = dfrkms(ifr-1)+ddfr
          enddo
c          write(ch,*) line_nfreq(iline)
#ifdef COMPILER_PG
c          if(line_nfreq(iline).lt.10) then
c              frmt = '('//ch(12:12)//'(E11.4,1X))'
c          elseif(line_nfreq(iline).lt.100) then
c              frmt = '('//ch(11:12)//'(E11.4,1X))'
c          elseif(line_nfreq(iline).lt.1000) then
c              frmt = '('//ch(10:12)//'(E11.4,1X))'
c          endif
#else
c          if(line_nfreq(iline).lt.10) then
c              frmt = '('//ch(2:2)//'(E11.4,1X))'
c          elseif(line_nfreq(iline).lt.100) then
c              frmt = '('//ch(2:3)//'(E11.4,1X))'
c          elseif(line_nfreq(iline).lt.1000) then
c              frmt = '('//ch(2:4)//'(E11.4,1X))'
c          endif
#endif
          frmt = '(6(E11.4,1X))'
          write(1,*) line_nfreq(iline)
          write(1,frmt) (dfrkms(ifr),ifr=1,line_nfreq(iline))
      enddo
      write(1,*) ' '      
c
c     First a loop over all the grid positions
c
#     ifdef LINE_INFO_COMPLETE
      do ir=1,irsi_frsizex
          do it=1,itmax
              call get_phys_var(ir,it)
              temp = medium_el_temp
              do iline=1,nlines
                  a1 = line_jbar_field(iline,it,ir) *
     %                 Bud_field(iline,it,ir) / Aud_field(iline,it,ir)
                  a2 = Cud_field(iline,it,ir) / Aud_field(iline,it,ir)
                  a7 = Cud_field(iline,it,ir) / 
     %                  ( Aud_field(iline,it,ir) +
     %                    line_jbar_field(iline,it,ir) *
     %                        Bud_field(iline,it,ir) +
     %                    Cud_field(iline,it,ir) )
                  a3 = 0.d0         ! For the moment...
                  a4 = 0.d0         ! For the moment...
                  ilev_up   = lev_up(iline)
                  ilev_down = lev_down(iline)
                  npop_up   = line_level_popul(ilev_up,it,ir)
                  npop_down = line_level_popul(ilev_down,it,ir)
                  a5 = 7.2427030d15 * ( enerlevel(ilev_down) - 
     %                 enerlevel(ilev_up) ) /
     %                 log((gdeg(ilev_down)*npop_up)/
     %                     (gdeg(ilev_up)*npop_down))
                  a6 = a5 / medium_arr_el_temp(it,ir)
                  a8 = 1.47455253991d-47 * ( linefreq(iline)**3 ) /
     %                   ( exp(4.79915984645d-11*linefreq(iline)/
     %                         temp) - 1.d0 )
                  a8 = line_jbar_field(iline,it,ir) / a8
                  write(1,101) a1,a2,a7,a3,a4,a5,a6,a8
  101             format(8(E11.4,1X))
              enddo
          enddo
      enddo
#     endif
c
      close(1)
c
      return
c
      end

c     --------------------------------------------------------------
c           WRITE THE LINE-AVERAGED MEAN INTENSITY TO FILE
c
c     This is the routine that dumps the \Bar J to a file. The 
c     filenr specification denotes the file:
c 
c       filenr: -1 = Save to linemeanint.dat
c               -2 = Save to linemeanint_final.dat
c               -3 = Save to linemeanint_<molecname>.dat
c                0,1,2... = Save to intmean_0.dat etc...
c
c     --------------------------------------------------------------
      subroutine write_line_jbar(filenr)
      implicit none
      integer filenr
c
#include "common_grid.h"
#include "common_lines.h"
c
      integer ir,it,itmax,imirt,i
      integer dum_freq_nr,dum_frsizex,dum_itmax,dum_imirt
      character*80 filename,frmt
      character*12 ch
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
      imirt = 0
#endif
c
c     Write the info file
c
      open(unit=1,file='linemeanint.info',status='unknown')
      write(1,*) filenr
      if(filenr.eq.-2) then
          write(1,*) line_jbar_lastsave
      endif
      if(filenr.eq.-3) then
          write(1,*) 'linemeanint_'//molec_name
     %                        (1:molec_name_len)//'.dat'
          write(1,*) line_jbar_lastsave
      endif
      close(1)
c
c     Open the file and write header
c
      if(filenr.eq.-1) then
         filename = 'linemeanint.dat'
      elseif(filenr.eq.-2) then
         filename = 'linemeanint_final.dat'
      elseif(filenr.eq.-3) then
         filename = 'linemeanint_'//molec_name
     %                  (1:molec_name_len)//'.dat'
      else
         write(ch,*) filenr
#ifdef COMPILER_PG
         if(filenr.lt.10) then
            filename = 'linemeanint_'//ch(12:12)//'.dat'
         elseif(filenr.lt.100) then
            filename = 'linemeanint_'//ch(11:12)//'.dat'
         elseif(filenr.lt.1000) then
            filename = 'linemeanint_'//ch(10:12)//'.dat'
         endif
#else
         if(filenr.lt.10) then
            filename = 'linemeanint_'//ch(2:2)//'.dat'
         elseif(filenr.lt.100) then
            filename = 'linemeanint_'//ch(2:3)//'.dat'
         elseif(filenr.lt.1000) then
            filename = 'linemeanint_'//ch(2:4)//'.dat'
         endif
#endif
         line_jbar_lastsave = filenr
      endif
      write(ch,*) nlines
#ifdef COMPILER_PG
      if(nlines.lt.10) then
          frmt = '('//ch(12:12)//'(E13.6,1X))'
      elseif(nlines.lt.100) then
          frmt = '('//ch(11:12)//'(E13.6,1X))'
      elseif(nlines.lt.1000) then
          frmt = '('//ch(10:12)//'(E13.6,1X))'
      endif
#else
      if(nlines.lt.10) then
          frmt = '('//ch(2:2)//'(E13.6,1X))'
      elseif(nlines.lt.100) then
          frmt = '('//ch(2:3)//'(E13.6,1X))'
      elseif(nlines.lt.1000) then
          frmt = '('//ch(2:4)//'(E13.6,1X))'
      endif
#endif
      open(unit=1,file=filename,status='unknown')
      write(1,*) irsi_frsizex,itmax,nlines,imirt
c
c     First write all the line frequencies
c
      write(1,*) (linefreq(i),i=1,nlines)
      write(1,*) (gratio(i),i=1,nlines)
      write(1,*) (Aud_const(i),i=1,nlines)
      write(1,*) (lev_up(i),i=1,nlines)
      write(1,*) (lev_down(i),i=1,nlines)
      write(1,*) ' '      
c
c     First a loop over all the grid positions
c
      do ir=1,irsi_frsizex
          do it=1,itmax
              write(1,frmt) (line_jbar_field(i,it,ir),i=1,nlevels)
          enddo
      enddo
      write(1,*) ' '      
      do ir=1,irsi_frsizex
          do it=1,itmax
              write(1,frmt) (line_source_func(i,it,ir),i=1,nlevels)
          enddo
      enddo
      write(1,*) ' '      
      do ir=1,irsi_frsizex
          do it=1,itmax
              write(1,frmt) (line_approx_lambdabar(i,it,ir),i=1,nlevels)
          enddo
      enddo
c
      close(1)
c
      return
c
      end

c     --------------------------------------------------------------
c                    WRITE GLOBAL LINE PROFILE 
c     --------------------------------------------------------------
      subroutine write_global_lineprofile()
      implicit none
      integer filenr
c
#include "common_grid.h"
#include "common_lines.h"
#include "common_setup.h"
c
      integer i,iline,ilfr,nfrmax
      doubleprecision a1,a2,a3,a4,a5,a6,a7
c
      nfrmax=0
      do iline=1,nlines
          if(line_nfreq(iline).gt.nfrmax) then
              nfrmax = line_nfreq(iline)
          endif
      enddo
c
      open(unit=1,file='lineprofile.dat',status='unknown')
      write(1,*) nlines,nfrmax
      write(1,*) ' '      
c
      do iline=1,nlines
          write(1,*) line_nfreq(iline),line_symm(iline)
          write(1,*) line_nu0(iline)
          do ilfr=1,line_nfreq(iline)
              write(1,*) line_dnu(ilfr,iline),line_phi(ilfr,iline),
     %                   line_wgt(ilfr,iline)
          enddo
      enddo
ccc   20 format(E21.13,1X,E14.7,1X,E14.7)
c
      end

c     --------------------------------------------------------------
c          INTEGRATE THE J OVER THE LINE PROFILE TO GET JBAR
c
c     Once the complete radiation field at all frequencies has been
c     computed, and the J(nu) has been found, we can compute the
c     Jbar for all the lines. We do this from the perspective of
c     the lines: we handle each line separately, find the radiation
c     field data corresponding to the frequencies in that line and
c     do a loop over space and integrate over the line profile at
c     each point in space. The formula we compute here is
c
c                             /
c                     bar J = | J \phi(\nu) d\nu
c                             /
c
c     When velocity fields are present (iradproc_line_velocities.eq.1)
c     then we shift the line profile according to the velocity v/c.
c
c     This routine will be called by calc_level_populations().
c     --------------------------------------------------------------
      subroutine calc_jbar()
      implicit none
c
#ifndef JSRC_NO_ARRAYS
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
#include "common_lines.h"
#include "common_setup.h"
c
      doubleprecision dummy
      integer iline, ilfr,ir,it,ithetamax
ccc
ccc     Checks
ccc      
cc      if((iscat_mode.ne.0).and.(iscat_mode.ne.1)) then
cc          write(*,*) 'iscat_mode must be either 0 or 1'
cc          stop 13
cc      endif
c
c     Now the mirror things...
c
#ifdef MIRROR_THETA
      ithetamax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      ithetamax = 1
#endif
c
c     Loop over space
c
      do ir=1,irsi_frsizex
          do it=1,ithetamax
c
c             Loop over the lines
c
              do iline=1,nlines
c
c                 Set the jbar zero 
c
                  line_jbar_field(iline,it,ir) = 0.d0
c
c                 Compute the integral over the frequencies
c                 in the line. 
c
c                             /
c                     bar J = | J \phi(\nu) d\nu
c                             /
c
c                 NOTE: The line profile is presumed global here!
c                       This is temporary. In the future it will
c                       be done more realistically.
c
                  do ilfr=1,line_nfreq(iline)
                      line_jbar_field(iline,it,ir) = 
     %                  line_jbar_field(iline,it,ir) +
     %                  line_wgt(ilfr,iline) * line_phi(ilfr,iline) *
     %                  intmom_0(line_trifreq(ilfr,iline),it,ir) 
                  enddo
#ifdef LINE_PHI_TEST
                  dummy = 0.d0
                  do ilfr=1,line_nfreq(iline)
                      dummy = dummy +
     %                  line_wgt(ilfr,iline) * line_phi(ilfr,iline)
                  enddo                  
                  write(*,*) ir,it,iline,dummy
#endif
              enddo
          enddo
      enddo
c
c     If everything went okay, we now have bar J as a function of
c     R and Theta, and for each of the lines.
c
#else
      write(*,*) 'Routine calc_jbar() only when source'
      write(*,*) 'and mean int arrays are present.'
      stop 13
#endif
c
      end





c     --------------------------------------------------------------
c                  SET THE VALUE OF ALL JBAR TO ZERO
c     --------------------------------------------------------------
      subroutine delete_jbar()
      implicit none
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
#include "common_lines.h"
#include "common_setup.h"
c
      doubleprecision dummy
      integer iline, ilfr,ir,it,ithetamax
c
c     Now the mirror things...
c
#ifdef MIRROR_THETA
      ithetamax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      ithetamax = 1
#endif
c
c     Loop over space
c
      do ir=1,irsi_frsizex
          do it=1,ithetamax
c
c             Loop over the lines
c
              do iline=1,nlines
c
c                 Set the jbar zero 
c
                  line_jbar_field(iline,it,ir) = 0.d0
c
              enddo
          enddo
      enddo
      end



c     --------------------------------------------------------------
c              CALCULATE THE LINE SOURCE FUNCTION J/ALPHA
c
c     For the MALI acceleration we need the source function S of 
c     the old populations.
c
c     WARNING: THERE IS SERIOUS DANGER FOR A BUG: You MUST first
c              call this function BEFORE calling
c
c     --------------------------------------------------------------
      subroutine calc_line_source_func()
      implicit none
c
#include "common_grid.h"
#include "common_radfield.h"
#include "common_source.h"
#include "common_lines.h"
#include "common_setup.h"
c
      integer ir,it,inu,itmax,ilevel,ilfr,iline
      doubleprecision dummy,nup,nlow,nmol,dens,temp
c
      itmax   = irsi_frsizey
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
      do ir=1,irsi_frsizex
          do it=1,itmax
c
c             Find the density and temperature
c
              call get_phys_var(ir,it)
              temp = medium_el_temp
              dens = medium_rho
c
c             Check if the level populations add up to 1
c
              dummy = 0.d0
              do ilevel=1,nlevels
                  dummy = dummy + line_level_popul(ilevel,it,ir) 
              enddo
              if(abs(1.d0-dummy).gt.1.e-2) then
                  write(*,*) 'At ir=',ir,' it=',it,': Level ',
     %                'populations add up to: ',dummy
                  stop 13
              endif
c
c             Now loop over all lines 
c
              do iline=1,nlines
c
c                 Find the upper and lower level population
c
                  nup  = line_level_popul(lev_up(iline),it,ir)
                  nlow = line_level_popul(lev_down(iline),it,ir)
c
c                 The S_kl is now:
c
c                                  n_k A_kl
c                   S_kl  =  -------------------         ( k > l )
c                            n_l B_lk - n_k B_kl
c
                  line_source_func(iline,it,ir) = 
     %               nup * Aud_field(iline,it,ir) / 
     %               ( nlow * Bdu_field(iline,it,ir) -
     %                 nup  * Bud_field(iline,it,ir) )
c
              enddo
          enddo
      enddo
c
      return
      end


c     --------------------------------------------------------------
c                   COMPUTE THE LTE LEVEL POPULATIONS
c
c     This routine will compute the LTE level populations 
c     everywhere on the grid and will insert these values into the
c     global array line_level_popul(ilevel,it,ir). So, if you
c     simply wish to compute LTE models, then call this routine
c     and you'll have all the populations there.
c     --------------------------------------------------------------
      subroutine calc_lte_levelpopul()
      implicit none
c
#include "common_grid.h"
#include "common_radfield.h"
#include "common_source.h"
#include "common_lines.h"
#include "common_setup.h"
c
      integer ir,it,inu,itmax,ilevel,ilfr,iline
      doubleprecision dummy,nup,nlow,nmol,temp
c
      itmax   = irsi_frsizey
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
      do ir=1,irsi_frsizex
          do it=1,itmax
              call get_phys_var(ir,it)
              temp = medium_el_temp
c
c             First assume that n_1=1.d0, then compute all the
c             higher levels, using
c
c               n_i   g_i
c               --- = --- exp ( - Delta E / kT )
c               n_1   g_1 
c
              line_level_popul(1,it,ir) = 1.d0
              dummy = 1.d0
              do ilevel=2,nlevels
                  line_level_popul(ilevel,it,ir) = 
     %             line_level_popul(1,it,ir) *
     %               (gdeg(ilevel)/gdeg(1)) *
     %               exp((enerlevel(1)-enerlevel(ilevel))/
     %                       (1.3807d-16*temp))
                  dummy = dummy + line_level_popul(ilevel,it,ir) 
              enddo
c
c             Now normalize to unit
c
              dummy = 1.d0 / dummy
              do ilevel=1,nlevels
                  line_level_popul(ilevel,it,ir) =
     %               line_level_popul(ilevel,it,ir) * dummy
              enddo
          enddo
      enddo
c
      return
      end


c     ==============================================================
c              THE MOLECULAR/ATOMIC PHYSICS FOLLOWS BELOW
c     ==============================================================



c     --------------------------------------------------------------
c         GLOBALLY PREPARE ALL TRANSITIONS, BOTH LINES AND COLL
c
c     If the density and temperature are known, the Bud, Bdu, Cud,
c     and Cdu can be computed. This routine does this for all points
c     on the grid. It puts this data into the appropriate arrays.
c     --------------------------------------------------------------
      subroutine global_prepare_transitions()
      implicit none
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_lines.h"
c
      doubleprecision temp,dens,ncol,nmolec
      integer ilev,iline,ir,it,ithetamax
c
c     Now the mirror things...
c
#ifdef MIRROR_THETA
      ithetamax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      ithetamax = 1
#endif
c
c     First restore the fundamental value of Aud. The Aud
c     variable may have been changed by MALI to the 
c     value Aud - Lambda*, so we must restore this. 
c     This was the bug that was so hard to find on 5-4-99.
c     
      do iline=1,nlines
          Aud(iline) = Aud_const(iline)
      enddo
c
c     Now we can solve the level populations everywhere.
c
      do ir=1,irsi_frsizex
          do it=1,ithetamax
c
c             Find density and temperature
c
              call get_phys_var(ir,it)
              temp = medium_el_temp
              dens = medium_rho
c
c             If local abundances, then also copy the local 
c             abundances in the appropriate slots.
c
#             ifdef LINE_LOCAL_ABUNDANCE
              call copy_abundances(ir,it)
#             endif
c
c             Compute the number density of the molecule itself
c
              nmolec = abund_molec * dens / (umass_av*1.6726d-24) 
c     
c             Compute the number density of the collision partner 
c     
              ncol = abund_colpartner * dens / (umass_av*1.6726d-24) 
c
c             Now prepare the lines and coll transitions here
c
              call prepare_transitions(temp,nmolec,ncol)
c
c             Copy the information into the arrays
c
              call save_transinfo(ir,it)
c     
          enddo
      enddo
c
c     Now, just for safety, tell the other routines that I have 
c     indeed done my job here. If this flag is not set, then apparently
c     I have not been called. This enables other routines to do 
c     internal consistency checks to verify if the transitions have
c     been prepared.
c     
      line_warn_preptrans_done = 1
c
      end



c     --------------------------------------------------------------
c         PREPARE ALL TRANSITIONS, BOTH LINES AND COLL (LOCAL)
c     --------------------------------------------------------------
      subroutine prepare_transitions(temp,ndens,ncol)
      implicit none
      doubleprecision temp,ndens,ncol
c
c     First compute the Cud and Cdu from the Kud values. This
c     depends on both the temperature and the density of the
c     collision partner (ncol)
c
c      call calc_collrates(temp,ncol)
c
c     Use the Einstein relations to compute Einstein coefficients
c     and compute all other information about the lines from the
c     levels, like the line frequency etc.
c
      call prepare_lines(temp,ndens)
c
      return
      end



c     --------------------------------------------------------------
c                   PREPARE ALL THE LINES AND SUCH
c
c     This routine computes the obvious things of the lines from 
c     the things we know of the levels, and of the temperature and
c     density. 
c     --------------------------------------------------------------
      subroutine prepare_lines(temp,ndens)
      implicit none
      doubleprecision temp,ndens
c
#include "common_lines.h"
c
      integer i
c
c     Check if nlines and nlevels are within range
c
      if(nlines.gt.SZ_NLINES) then
          write(*,*) 'Too many lines. Recompile with larger SZ_NLINES'
          stop 13
      endif
      if(nlevels.gt.SZ_NLEVELS) then
          write(*,*) 'Too many levels. Recompile with larger SZ_NLEVELS'
          stop 13
      endif
      if(nlines.lt.1) then
          write(*,*) 'Minimum of 1 line!'
          stop 13
      endif
      if(nlevels.lt.2) then
          write(*,*) 'Minimum of 2 levels!'
          stop 13
      endif
c
      do i=1,nlines
c
c         Check if the lev_up and lev_down are in range
c
          if((lev_up(i).lt.1).or.(lev_down(i).lt.1).or.
     %       (lev_up(i).gt.nlevels).or.(lev_down(i).gt.nlevels)) then
              write(*,*) 'At line ',i,' levels should be within the ',
     %         'range 1,',nlevels
              stop 13
          endif
c
c         Check if the lev_up > lev_down
c          
          if(lev_up(i).le.lev_down(i)) then
              write(*,*) 'The lines should be given from upper to ',
     %                   'lower level, in that order.'
              write(*,*) 'Line number ',i,' violates that!'
              stop 13
          endif
c
c         Check if the line frequencies have already been computed
c         (this should be done in read_species()).
c
          if(linefreq(i).eq.0.d0) then
              write(*,*) 'HUH?? Somehow linefreq() is 0.d0?'
              stop 13
          endif
c
c         Compute this frequency in units [K] , so that the
c         exponent exp(-h\nu/kT) will be exp(-dek(i)/temp).
c
          dek(i)  = 4.799160d-11 * linefreq(i) 
c
c         The ratio of the degeneracy parameters g_up/g_down
c
          gratio(i) = gdeg(lev_up(i)) / gdeg(lev_down(i)) 
c
c         The Bud from the Aud:
c    
c                     c^2
c             Bud = -------- Aud 
c                   2 h nu^3
c    
          Bud(i) = 6.78171833781d46 * Aud(i) / 
     %                 (linefreq(i)*linefreq(i)*linefreq(i))
c    
c         The Bdu from the Bud:
c    
c                   gu
c             Bdu = -- Bud
c                   gd
c    
          Bdu(i) = Bud(i) * gratio(i)
      enddo
c
      do i=1,ncoltrans
c
c         Compute the energy step in units of [K], so that the
c         exponent exp(-dE/kT) will be exp(-dek_ct(i)/temp).
c
          dek_ct(i) = 7.24270297675d15 *
     %                (enerlevel(lev_ct_up(i))-
     %                 enerlevel(lev_ct_down(i)))
c
c         The ratio of the degeneracy parameters g_up/g_down
c
          gratio_ct(i) = gdeg(lev_ct_up(i)) / gdeg(lev_ct_down(i)) 
c
c         The Cdu from the Cud:
c
c                   gu           /   dE \
c             Cdu = -- Cud * exp | - -- |
c                   gd           \   kT /
c
          Cdu(i) = Cud(i) * gratio_ct(i) * exp(-dek_ct(i)/temp)
c
#ifdef DEBUG
          write(*,*) 'frq = ',linefreq(i)
          write(*,*) 'dek = ',dek(i)
          write(*,*) 'Aud = ',Aud(i)
          write(*,*) 'Bud = ',Bud(i)
          write(*,*) 'Bdu = ',Bdu(i)
          write(*,*) 'Cud = ',Cud(i)
          write(*,*) 'Cdu = ',Cdu(i)
#endif
c
      enddo
c    
      end


      subroutine read_species_lambda(filename,ilevel_top)
      implicit none
      character*80 filename
      integer ilevel_top,idum
c
#include "common_lines.h"
c
      character*160 thefilename,header
      integer i
c
      write(*,*) 'Reading molecular data - lambda style'
      thefilename = 
     % MOLECULAR_DATA_DIR

      i=index(thefilename,' ')
      if(thefilename(i-1:i-1).ne.'/') then
          thefilename = thefilename(1:i-1)//'/'//filename
      else
          thefilename = thefilename(1:i-1)//filename
      endif
c      write(*,*) i,thefilename

      molec_filename=thefilename
      write(*,500) thefilename
  500 format(' Reading line-data file: ',A)
      open(unit=5,file=thefilename,status='old',err=777)
      goto 779
  777 continue
      close(5)
      thefilename = './'//filename
      open(unit=5,file=thefilename,status='old',err=778)
      goto 779
  778 continue
      write(*,*) 'ERROR: Could not open molecule file'
      stop 13
  779 continue
c
c     Read the first dummy strings
c      
      read(5,*) header
      read(5,*) header
      read(5,*) header      
c
c     Read the atomic weight
c
      read(5,'(F4.1)') umass_molec
      read(5,*) header      
c
c     Read the amount of levels are present in the file
c
      read(5,'(I6)') nlev_orig
      if(ilevel_top.ne.0) then 
         nlevels = min(ilevel_top,nlev_orig)
      else
         nlevels = nlev_orig
      endif
c
c     Another dummy line
c
      read(5,*) header
c
c     Checks
c
      if(nlev_orig.gt.SZ_NLEVELS) then
          write(*,*) 'Too many levels. Recompile with larger SZ_NLEVELS'
          stop 13
      endif
      if(nlev_orig.lt.2) then
          write(*,*) 'Minimum of 2 levels!'
          stop 13
      endif
c
c     Now do a loop over the levels and read their energies in cm^-1
c     and convert to erg
c
      do i=1,nlev_orig
          read(5,'(I5,F12.4,F7.1)') idum, enerlevcm(i), gdeg(i)
          enerlevel(i) = 1.986468498d-16 * enerlevcm(i)
      enddo
c
c     Another dummy line
c
      read(5,*) header
c
c     Read the amount of lines are present in the file
c
      read(5,'(I6)') nlin_orig
c
c
c
      if(nlin_orig.gt.SZ_NLINES) then
          write(*,*) 'Too many lines. Recompile with larger SZ_NLINES'
          stop 13
      endif
      if(nlin_orig.lt.1) then
          write(*,*) 'Minimum of 1 line!'
          stop 13
      endif
c
c     Another dummy line
c
      read(5,*) header
c
c     Now do a loop over the levels and read their indices and 
c     Einstein A's
c
      do i=1,nlin_orig
          read(5,'(I5,I5,I5,E12.3)') 
     %        idum, lev_up(i), lev_down(i), Aud_const(i)
      enddo
c
c     Now read the Kud values. This is a bit more subtle since these
c     depend on the temperature. That's why I'll call a separate 
c     reading routine to do this for me... This will read the Kud
c     values for a series of temperature samples. The real Cud values
c     must be computed from these by interpolation in temperature and
c     by multiplication with the number density of the collision 
c     partner. 
c
c
c LATER MUST INCLUDE COLLISION RATES PROPERLY!!!!
c
      call read_collisions_lambda()
c
c     Next compute the real nlines. See remark in the 
c     header about the ordering of the line and collisional 
c     transition data in the molecule file.
c
c      do i=1,nlin_orig
c          if(lev_up(i).le.nlevels) then
              nlines = nlin_orig
c          endif
c      enddo
c
c     Now we're done reading... Remember to call calc_collrates()
c     to compute the Cud and Cdu from the Kud!
c
c      write(*,*) ' '

      write(*,*) '  nlevels  = ',nlev_orig
      write(*,*) '  nlines   = ',nlin_orig
c
c     If you like, here's a level diagram in ascii format
c
      write(*,*) 'Level Diagram:'
      call pictlevels()
c
      close(5)
c
c     Compute the frequency of the transitions
c     
c            [ enerlevel ] = erg
c            [ freq ]      = Hz
c
      do i=1,nlin_orig
          linefreq(i) = 1.509160d26 * ( enerlevel(lev_up(i)) - 
     %        enerlevel(lev_down(i)) ) 
      enddo
c
      end


c     --------------------------------------------------------------
c                      READ THE COLLISION RATES
c
c     The files contain the K_ij values for the collision rates
c     at a discrete set of temperature points. The C_ij Einstein
c     coefficients are related to the K_ij by
c     
c       C_ij = K_ij * n_col
c
c     where n_col is the density of the collision partner in units
c     of 1/cm^3. 
c     
c     In this routine we merely read in the K_ij values. The 
c     computation of the C_ij and C_ji values are done in another
c     routine: calc_collisions().
c
c     NOTE: Be sure that the input file is open, and that the file
c           pointer is at the start of the data for the K_ij 
c           (this includes the first text comment line)
c
c     NOTE: Be sure that the nlevels is set
c     --------------------------------------------------------------
      subroutine read_collisions()
c
#include "common_lines.h"
c
      character*80 string
      integer itemp,itrans,i
c
c     Check
c
      if(nlevels.eq.0) then
          write(*,*) 'ERROR read_collisions: Forgot to set nlevels!'
          stop 13
      endif
c
c     Read a string
c
      read(1,*) string
c
c     First read the number of collisional transitions, the
c     number of sampled temperatures at which this data is
c     available, and the value of the temperatures at which
c     the data is available.
c
      read(1,101) ncoltr_orig,ncttemp,(temptrans(i),i=1,ncttemp)
  101 FORMAT(I3,I3,7F6.2,3F6.1)
c
c     Check if the ncoltr_orig and ncttemp are within the required
c     boundaries
c
      if((ncoltr_orig.lt.1).or.(ncoltr_orig.gt.SZ_NCOLTRANS)) then
          write(*,*) 'ncoltr_orig is not between 1 and ',SZ_NCOLTRANS
          write(*,*) 'Perhaps recompile with larger SZ_NCOLTRANS'
          stop 1
      endif
      if(ncttemp.eq.0) then 
          ncttemp=1
      endif
      if((ncttemp.lt.1).or.(ncttemp.gt.SZ_NCTTEMP)) then
          write(*,*) 'ncttemp is not between 1 and ',SZ_NCTTEMP
          write(*,*) 'Perhaps recompile with larger SZ_NCTTEMP'
          stop 1
      endif
c
c     Now read the up and down level numbers
c     
      read(1,102) (lev_ct_up(i),i=1,ncoltr_orig)
      read(1,102) (lev_ct_down(i),i=1,ncoltr_orig)
  102 format(25I3)
c     
c     Now do a loop over the temperatures
c
      do itemp=1,ncttemp
c
c         Read the K_ij for these transitions for this temperature
c
          read(1,104) (Kud(i,itemp),i=1,ncoltr_orig)
  104     format(8E8.1)
      enddo
c
c     Next compute the real nlines. See remark in the 
c     header about the ordering of the line and collisional 
c     transition data in the molecule file.
c
      do i=1,ncoltr_orig
          if(lev_ct_up(i).le.nlevels) then
              ncoltrans = i
          endif
      enddo
c
      end
c
c=======================================================
cAlternative lambda format
c=======================================================
      subroutine read_collisions_lambda()
c
#include "common_lines.h"
c
      character*80 string
      integer itemp,itrans,i
c
      ncoltr_orig=1
      ncttemp=2
      temptrans(1) = 10
      temptrans(2) = 2000

      lev_ct_up(1) = 2
      lev_ct_down(2) = 1
      
      Kud(1,1) = 1e-12
      Kud(1,2) = 1e-12
c
      end



c     --------------------------------------------------------------
c           MAKE GRAPHICAL PICTOGRAM OF LEVELS AND TRANSITIONS
c     --------------------------------------------------------------
      subroutine pictlevels()
      implicit none
c
#include "common_lines.h"
c
      integer i,k,nlinmax,istart,iend
      character*83 pict
c
      istart=14
      nlinmax = min(nlines,69)
      do k=nlevels,1,-1
          write(pict(1:10),403) enerlevcm(k)
  403     format(1F10.5)
          do i=11,istart
              pict(i:i) = ' '
          enddo
          iend = 1
          do i=1,nlinmax
              if((k.ge.lev_down(i)).and.(k.le.lev_up(i))) then
   14             pict(i+istart:i+istart) = '|'
                  iend = i+istart
              else
                  pict(i+istart:i+istart) = ' '
              endif
          enddo
          do i=nlinmax+istart+1,83
              pict(i:i) = ' ' 
          enddo
          write(*,*) pict (1:iend+1)
ccc  410     format(A83)
      enddo
c
      end





c     --------------------------------------------------------------
c                       THE VOIGT LINE PROFILE 
c     
c     This function evaluates a line profile including the effects
c     of natural broadening, and turbulence. 
c
c     ARGUMENTS:
c       nu00             = Center freq in Hz. If specified <0 then 
c                          the output nu will be nu-abs(nu00), so as to
c                          avoid numerical problems later.
c       gammaks          = Natural line width in km/s
c       aks              = Turbulent line width in km/s
c       passbks          = The width of the sample interval 
c                           (from line center to edge) in km/s
c       nfr              = The amout
c
c     OUTPUT:
c       nu               = The frequency points at which the line is given
c       phi              = The line profile 
c       wgt              = The integration weights
c
c     WARNING: The weights are chosen rather bluntly, so no 
c              mathematical justification can be given for them.
c              But the norm is exactly 1.
c
c     WARNING: Right now the Voigt is simply a Gaussian, since the
c              integral turned out not to be so simple and I have
c              no time to do it properly now (14-11-98). So only
c              the aks is the thing that gives you broadning: the
c              microturbulence.
c     --------------------------------------------------------------
      subroutine make_voigt_profile(nu00,gammaks,aks,passbks,
     %                              nfr,nu,phi,wgt)
      implicit none
      doubleprecision nu00,gammaks,aks,passbks
      doubleprecision phi(SZ_LINEPROFILE),nu(SZ_LINEPROFILE)
      doubleprecision wgt(SZ_LINEPROFILE)
      integer nfr
      doubleprecision voigtfunc
c
      integer i
      doubleprecision gamma,aa,passb,nu1,nu2,dnu,nu0
      doubleprecision norm,dum,avoigt,uvoigt
c
      nu0 = abs(nu00)
c
c     Check line freq
c
      if(nu0.eq.0.d0) then
          write(*,*) 'Problem in Voigt(): nu0=0 !'
          stop 13
      endif
c
c     If nfr is too large or too small, then stop
c
      if((nfr.gt.SZ_LINEPROFILE).or.(nfr.lt.1)) then
          write(*,*) 'Number of frequencies for this line is ',
     %          'out of range: ',nfr
          write(*,*) 'Perhaps recompile with larger SZ_LINEPROFILE'
          stop 13
      endif
c     
c     First convert the line widths from km/s to Hz
c
      gamma = 3.33567d-6 * nu0 * gammaks      ! in Hz
      aa    = 3.33567d-6 * nu0 * aks          ! in Hz
      passb = 3.33567d-6 * nu0 * passbks      ! in Hz
c
c     If only one point in the line, then the function is assumed to
c     be a block with width gamma+aa from center to edge, so the
c     function will be normalized to 0.5d0/(gamma+aa)
c
      if(nfr.eq.1) then
          nu(1)  = nu0
          phi(1) = 0.5d0/(gamma+aa)
          wgt(1) = 2.0d0*(gamma+aa)
      endif
c
c     Else, do the real thing. Make a loop over the line profile
c      (for the moment it is simply a Gaussian)
c
c                        1       -((nu-nu0)/aa)^2
c         phi(nu) = ----------- e
c                   aa sqrt(pi)
c
c
      dum    = 0.d0
      norm   = 0.56419583546d0 / aa
c      avoigt = 0.0795774715459d0 * gamma / aa
      nu1    = 0.d0 - passb       ! nu0 already subtracted
      nu2    = 0.d0 + passb       ! nu0 already subtracted
      dnu    = 2*passb / (nfr-1.d0)
      do i=1,nfr
          nu(i)  = nu1 + (i-1)*dnu
          uvoigt = ( nu(i) - 0.d0 ) / aa
c          phi(i) = norm * voigtfunc(avoigt,uvoigt)
          phi(i) = norm * exp(-uvoigt**2)
          dum    = dum + phi(i)
      enddo
c
c     Determine the weights very crudely for now: equal weights
c     for each point.
c
      dum = 1.d0 / dum 
      do i=1,nfr
          wgt(i) = dum
      enddo      
c
      end


c     --------------------------------------------------------------
c                       THE VOIGT LINE PROFILE 
c     
c     This function evaluates a line profile including the effects
c     of natural broadening, and turbulence. 
c
c     ARGUMENTS:
c       nu00             = Center freq in Hz. If specified <0 then 
c                          the output nu will be nu-abs(nu00), so as to
c                          avoid numerical problems later.
c       gammaks          = Natural line width in km/s
c       aks              = Turbulent line width in km/s
c
c     OUTPUT:
c       phi              = The line profile 
c
c     WARNING: Right now the Voigt is simply a Gaussian, since the
c              integral turned out not to be so simple and I have
c              no time to do it properly now (14-11-98). So only
c              the aks is the thing that gives you broadning: the
c              microturbulence.
c     --------------------------------------------------------------
      function voigt_profile(nu00,gammaks,aks,dnu)
      implicit none
      doubleprecision voigt_profile
      doubleprecision nu00,gammaks,aks,dnu
c      doubleprecision voigtfunc
c
      doubleprecision gamma,aa,nu0,norm
      doubleprecision uvoigt
c
      nu0 = abs(nu00)
c
c     Check line freq
c
      if(nu0.eq.0.d0) then
          write(*,*) 'Problem in Voigt(): nu0=0 !'
          stop 13
      endif
c     
c     First convert the line widths from km/s to Hz
c
c      gamma = 3.33567d-6 * nu0 * gammaks      ! in Hz
      aa    = 3.33567d-6 * nu0 * aks          ! in Hz
      norm  = 0.56419583546d0 / aa
c
c     The line profile is:
c
c                        1       -((nu-nu0)/aa)^2
c         phi(nu) = ----------- e
c                   aa sqrt(pi)
c
      uvoigt = dnu / aa
      voigt_profile = norm * exp(-uvoigt**2)
c
      return
      end





c     --------------------------------------------------------------
c                  EVALUATION OF THE VOIGT FUNCTION
c     --------------------------------------------------------------
      function voigtfunc(a,u)
      implicit none
      doubleprecision a,u,voigtfunc
      common/voigtint/aa,uu
      doubleprecision aa,uu
      external voigtintgr
      aa=a
      uu=u
      call qgausn(voigtintgr,-2.d0,2.d0,20,voigtfunc)
      return
      end
c
      function voigtintgr(y)
      implicit none
      doubleprecision y,voigtintgr
      common/voigtint/aa,uu
      doubleprecision aa,uu
      voigtintgr = 0.318309886184d0 * aa * exp(-y*y) / 
     %             ( aa*aa + (uu-y)*(uu-y) )
      write(*,*) 'y = ',y,' intgr = ',voigtintgr
      return
      end


c     --------------------------------------------------------------
c                PLOT VOIGT LINE PROFILE (TEST PURPOSE)
c     --------------------------------------------------------------
      subroutine plotvoigt()
      implicit none
c
      doubleprecision nu0,gammaks,aks,passbks
      integer i,nfr
      doubleprecision phi(SZ_LINEPROFILE),nu(SZ_LINEPROFILE)
      doubleprecision wgt(SZ_LINEPROFILE)
c
      write(*,*) 'nu0'
      read(*,*) nu0
      write(*,*) 'nfr'
      read(*,*) nfr
      write(*,*) 'gamma (km/s)'
      read(*,*) gammaks
      write(*,*) 'a (km/s)'
      read(*,*) aks
      write(*,*) 'passband (km/s)'
      read(*,*) passbks
c
      call make_voigt_profile(nu0,gammaks,aks,passbks,nfr,nu,phi,wgt)
c
      do i=1,nfr
          write(*,*) nu(i),phi(i),wgt(i)
      enddo
c
      end

c     --------------------------------------------------------------
c                CONVERT JBARPLANCK INTO JBARLINE
c
c     The Jbarplanck is the \bar J in units of the planck function
c     at the frequency of the line. So we have to compute the 
c     line frequency from the energy levels, and then the planck
c     function there, and then compute the real bar J.
c     --------------------------------------------------------------
      subroutine convert_jbar(temp)
      implicit none
      doubleprecision temp
c
#include "common_lines.h"
c
      integer i
      doubleprecision nu,bplanck
c
      do i=1,nlines
c
c         Compute the line frequency
c          
      nu = 1.509160d26 * ( enerlevel(lev_up(i)) - 
     %                    enerlevel(lev_down(i)) ) 
      bplanck = 1.47455253991d-47 * nu * nu * nu / 
     %              ( exp(4.79915984645d-11*nu/temp) - 1.d0 )
      Jbarline(i) = Jbarplanck(i) * bplanck
      enddo
c
      end



c     --------------------------------------------------------------
c                     COPY THE STORED JBAR TO BACKUP
c
c     istore1 = -1 means from line_jbar_field
c     istore2 = -1 means to line_jbar_field
c     --------------------------------------------------------------
c$$$      subroutine line_jbar_copy(istore1,istore2)
c$$$      implicit none
c$$$      integer istore1,istore2
c$$$c
c$$$#include "common_grid.h"
c$$$#include "common_source.h"
c$$$#include "common_lines.h"
c$$$c
c$$$      integer it,ir,iline,itmax
c$$$c
c$$$      itmax   = irsi_frsizey
c$$$c
c$$$#     ifdef MIRROR_THETA
c$$$      itmax = (irsi_frsizey+1)/2
c$$$#     endif
c$$$#     ifdef RADGRID_ONEDIM
c$$$      itmax = 1
c$$$#     endif
c$$$      if((istore1.eq.-1).and.(istore2.ne.-1)) then
c$$$          do ir=1,irsi_frsizex
c$$$              do it=1,itmax
c$$$                  do iline=1,nlines
c$$$                      line_jbar_field_bk(iline,it,ir,istore2) = 
c$$$     %                          line_jbar_field(iline,it,ir)
c$$$                  enddo              
c$$$              enddo
c$$$          enddo
c$$$      elseif((istore1.ne.-1).and.(istore2.eq.-1)) then
c$$$          do ir=1,irsi_frsizex
c$$$              do it=1,itmax
c$$$                  do iline=1,nlines
c$$$                      line_jbar_field(iline,it,ir) = 
c$$$     %                      line_jbar_field_bk(iline,it,ir,istore1)
c$$$                  enddo              
c$$$              enddo
c$$$          enddo          
c$$$      elseif((istore1.ne.-1).and.(istore2.ne.-1)) then
c$$$          do ir=1,irsi_frsizex
c$$$              do it=1,itmax
c$$$                  do iline=1,nlines
c$$$                      line_jbar_field_bk(iline,it,ir,istore2) = 
c$$$     %                      line_jbar_field_bk(iline,it,ir,istore1)
c$$$                  enddo              
c$$$              enddo
c$$$          enddo                    
c$$$      else 
c$$$          write(*,*) 'Error in line_jbar_copy'
c$$$          stop 13
c$$$      endif
c$$$c
c$$$      end

#endif /* LINE_NG_ACCEL */





c     ==============================================================
c            ROUTINES FOR VELOCITY-DEPENDENT LINE PROFILES
c
c     So far, the routines in this file are only for media with
c     zero systematic velocity. The line profile is then isotropic.
c     Also the line profile was presumed to be a global one.
c     But when systematic velocities are present, the line profile
c     depends on angle, and the integral over I should be done in
c     both Omega and nu. The routines below allow for this systematic
c     velocity, and also for varying line profile. 
c
c     WARNING: For velocity dependent lines it can be very important
c              to chop the short characteristics into pieces so as
c              to make sure that, as the photon angles change along
c              a short characteristic, the correct line profile
c              value is taken. (priv.comm. v.Noort, Hubeny, Lanza).
c              But so far this 'chopping' is not yet implemented!
c              It is rather cumbersome to implement... 
c
c     NOTE: THIS MODULE IS NOT FINISHED YET, SO DONT ACTIVATE
c           Still to be done:
c             - Add the appropriate routines for copying intens to
c                intensity_array, and resetting the counter when
c                the next line is started
c             - Implement calls to these routines in the iteration.F
c               module.
c             - Add routines for reading and storing the velocity
c               into the arrays (well, this can be done in medium.F)
c             - Add switch in code so that if-not-set, the usual
c               line stuff is done, and if-set, the velocity-dependent
c               stuff is done.
c             - Make sure that if that switch is put on, but some
c               things are not initialized, the routine will 
c               automatically reproduce the results as if the usual
c               line transfer (without velocities) was done.
c           (02-04-99) 
c     ==============================================================
#ifdef LINE_VELOCITIES

c     --------------------------------------------------------------
c          INTEGRATE I OVER OMEGA AND LINE PROFILE TO GET JBAR
c         For systematic velocities: use this i.s.o. calc_jbar()
c     
c     This routine computes the $\bar J$ for a given line, from the
c     angle-dependent intensity I(\Omega,\nu). If _no_ systematic
c     velocities are present, then it is better to use the function
c     calc_jbar(), and undefine LINE_VELOCITIES in the Makefile. 
c     When there _are_ systematic velocities $\vec v$ in the medium,
c     one must define the constant LINE_VELOCITIES (make sure that
c     you have sufficient RAM memory!), and use calc_velo_jbar()
c     instead. The line profile $\phi$ will be given by:
c
c       \phi_ij(\Omega,\nu) = \phi^0_ij(\nu-\nu_ij\Omega\cdot\vec v)
c
c     The function calc_velo_jbar() computes the double integral:
c
c                     / /
c             bar J = | | I \phi_ij(\Omega,\nu) d\nu d\Omega
c                     / /
c
c     Note that (in contrast to calc_jbar) this function immediately
c     integrates the _intensity_ I, not the mean-intensity J. The
c     LINE_VELOCITIES constant (if set) will make sure that the
c     intensity I is being stored and not dropped after each 
c     frequency (as it is done usually, after having computed J).
c     But since the angular-dependence of I requires a large chunk 
c     of memory, the complete I(R,\Theta,\mu,\phi,\nu) is only 
c     stored for a limited set of frequencies: those frequencies 
c     that are in the same passband for one line. As soon as all 
c     the frequency-samples of the line have been done, this function
c     calc_velo_jbar() should be called, so that the arrays for 
c     the storage of I(R,\Theta,\mu,\phi,\nu) are re-usable for the
c     next line. For a typical molecule this saves about a factor of
c     10 in memory, which could mean the difference between 'go' or
c     'no-go' for the problem you wish to solve in 2-D. This memory
c     saving might be a bit confusing though.
c
c     There is another thing that this routine does, as opposed
c     to the calc_jbar(): it locally computes the line profile,
c     so that space-dependent line widths can be used. The other
c     routine uses a pre-existing line profile, set up by the line.F
c     module prior to starting the radiative transfer. Here we 
c     actually compute this on the fly. 
c
c     This routine will be called by calc_level_populations().
c
c     WARNING: Since we 'recycle' the intensity_array, it is of 
c              crucial importance to be sure that the iline argument
c              in fact corresponds to the line for which the 
c              radiation field is now stored in intensity_array.
c              Otherwise the integral is done on the wrong intensity.
c
c     NOTE:    The 'recycling' obviously excludes line blanketing
c     --------------------------------------------------------------
      subroutine calc_velo_jbar(iline)
      implicit none
      integer iline
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
#include "common_lines.h"
#include "common_setup.h"
c
      doubleprecision omega_dot_v,line_integr_phi_domega_dnu
c
      doubleprecision aks
      integer ir,it,ithetamax
c
c     Now the mirror things...
c
      ithetamax = irsi_frsizey
c
#     ifdef MIRROR_THETA
      ithetamax = (irsi_frsizey+1)/2
#     endif
#     ifdef RADGRID_ONEDIM
      ithetamax = 1
#     endif
c
c     For 2-D runs the intensity_array must be kept in check (for it
c     for to become too large). It has another frequency dimension.
c     Check here if the line_nfreq(iline) is not too big.
c
      if(line_nfreq(iline).gt.SZ_LINEPROFILE_SMALL) then
          write(*,*) 'ERROR: line.F/calc_velo_jbar():'
          write(*,*) '       iline = ',iline,'  nrfeq = ',
     %                    line_nfreq(iline)
          write(*,*) '       Too many frequencies in this line!'
          if(line_nfreq(iline).le.SZ_LINEPROFILE) then
              write(*,*) '       To save memory, the maximum ',
     %               'amount of frequencies for the '
              write(*,*) '       ESC/MESC algorithm is smaller ',
     %               'than the maximum for the spectrum computation'
              write(*,*) '       Maximum for ESC/MESC = ',
     %                 SZ_LINEPROFILE_SMALL
              write(*,*) '       Modify SZ_LINEPROFILE_SMALL in ',
     %                'the Makefile if necessary'
#             ifndef RADGRID_ONEDIM
              write(*,*) '       But beware: in 2-D this may cause ',
     %                  'memory problems!!'
#             endif
          endif
          stop 13
      endif
c
c     Loop over space
c
      do ir=1,irsi_frsizex
          do it=1,ithetamax
c
c             Call a function to compute the integral
c             \int I \phi(\Omega,\nu) d\Omega d\nu
c
              line_jbar_field(iline,it,ir) = 
     %            line_integr_phi_domega_dnu(iline,
     %                ir,it,medium_arr_velocity(1,it,ir),
     %                intensity_array(1,1,-FRSIZE_MU_HALF,it,ir))
c
          enddo
      enddo
c      
      end


c     --------------------------------------------------------------
c                 (SUBROUTINE USED BY CALC_VELO_JBAR)
c
c     This routine computes the inner product between the angle
c     and the velocity, divided by c. The velocity is in [cm/s] 
c     and has three components:
c
c       v(1)    =  radial velocity (outward: >0)
c       v(2)    =  velocity towards negative Theta (i.e. toward NP)
c       v(3)    =  velocity towards positive Phi (i.e. counter
c                    clockwise if viewed from Theta=0)
c     --------------------------------------------------------------
      function omega_dot_v(mu,phi,v)
      implicit none
      doubleprecision omega_dot_v
      doubleprecision mu,phi,v(3)
      doubleprecision dummy
c
      if(mu.gt.1.d0) stop 393
c
      omega_dot_v = 3.335668d-11 * ( mu * v(1) + sqrt(1.d0-mu*mu) * 
     %              ( v(2) * sin(phi) + v(3) * cos(phi) ) )
c
      return
      end



c     --------------------------------------------------------------
c              COMPUTE DOUBLE INTEGRAL OVER OMEGA,NU
c
c     This subroutine computes the double integral over Omega and nu
c     of any quantity, tabulated in a 3-dimensional array in
c     mu_i,phi_j,nu_k, times the line profile phi shifted according
c     to the Doppler shift. Make sure that the line profile has
c     been computed globally (i.e. line_phi() and line_wgt()), or
c     locally (if LINE_LOCAL_WIDTH is set), and that the splines
c     for these profiles have been computed. All of this is done
c     automatically by the routine line_velo_prep_profiles().
c
c     The present function computes:
c
c       1    / /
c      ----  | | func(\Omega,\nu) \phi(\Omega,\nu) d\Omega d\nu
c      4 pi  / /
c
c     over the line pass band.
c
c     NOTE: It might be confusing that there are two 'phi' variables.
c           the philine() is the line profile and the phi is the angular
c           coordinate.
c     NOTE: This is a very subtle function, please check it again and
c           again, because errors may not show up conspicuously.
c     NOTE: Remember that when MIRROR_PHI is set, the intensity is
c           stored in a compact way, which messes up the indexing
c           of the func in the iphi-index. This is the reason for
c           the weird iphioff term.
c     NOTE: if LINE_LOCAL_WIDTH is not set, then the arguments ir and
c           it are irrelevant and ignored. 
c     --------------------------------------------------------------
      function line_integr_phi_domega_dnu(iline,ir,it,velocity,func)
      implicit none
      doubleprecision line_integr_phi_domega_dnu
c
      doubleprecision velocity(3)
      doubleprecision func(SZ_LINEPROFILE_SMALL,1:FRSIZE_PHI_SMALL,
     %                     -FRSIZE_MU_HALF:FRSIZE_MU_HALF)
      integer iline,ir,it
c
#include "common_grid.h"
#include "common_lines.h"
c
      doubleprecision omega_dot_v,line_integr_phi_dnu
c
      doubleprecision mu,phi
      doubleprecision moment,phiweight,dummy
      integer imu,iphi,iangset,iphimin,iphimax,iphihalf,iphioff,k
c
c     Mirror things
c
      iangset   = 1
      iphimin   = 1
      iphimax   = nrphi(iangset)
      iphihalf  = nrphi(iangset) / 2
      iphioff   = 0
c
#     ifdef MIRROR_PHI
      iphimin   = nrphi(iangset) / 4 + 1
      iphimax   = 3 * nrphi(iangset) / 4 
      iphioff   = - nrphi(iangset) / 4
      if(abs(velocity(3)).gt.1.d-10) then
          write(*,*) 'ERROR: If velocity(3).ne.0 then ',
     %               'no MIRROR_PHI allowed'
          stop 13
      endif
#     endif
#     ifdef RADGRID_ONEDIM
      iphimin   = 1
      iphimax   = 1
      iphioff   = 0
#     endif
c
c     Integration weight in phi angle
c
      phiweight = 1.d0 / (1.d0*(iphimax-iphimin+1))
c
c     Reset all values
c         
      moment = 0.d0
      dummy  = 0.d0
c
c     If a mu=0 gridpoint exists in the mu-space, then 
c     start with this one.
c
      if(anggrid_mu_zero.eq.1) then
          imu = 0
          do iphi=iphimin,iphimax
c
c             Find the angles
c     
              mu    = rmu(imu,iangset)
              phi   = rphi(iphi,iangset)
c
c             Perform the frequency integration
c
              dummy = dummy + line_integr_phi_dnu(iline,ir,it,
     %                func(1,iphi+iphioff,imu),mu,phi,velocity)
c
          enddo
      endif
      dummy  = dummy * phiweight
      moment = moment + dummy * 
     %    (rmu_i(1,iangset)-rmu_i(-1,iangset))
c
c     Now the other mu-gridpoints
c
      do imu=1,nrmu(iangset)
          dummy = 0.d0
          do iphi=iphimin,iphimax
c
c             Find the angles
c     
              mu    = rmu(imu,iangset)
              phi   = rphi(iphi,iangset)
c
c             Perform the frequency integration
c
              dummy = dummy + line_integr_phi_dnu(iline,ir,it,
     %                func(1,iphi+iphioff,imu),mu,phi,velocity)
c
c             Find the angles
c     
              mu    = rmu(-imu,iangset)
              phi   = rphi(iphi,iangset)
c
c             Perform the frequency integration
c
              dummy = dummy + line_integr_phi_dnu(iline,ir,it,
     %                func(1,iphi+iphioff,-imu),mu,phi,velocity)
c
          enddo
          dummy  = dummy * phiweight
          moment = moment + dummy * 
     %        (rmu_i(imu+1,iangset)-rmu_i(imu,iangset))
      enddo
c
c     Return the integral
c
      line_integr_phi_domega_dnu = moment * 0.5d0
c
      return
      end




c     --------------------------------------------------------------
c            (SUBROUTINE USED BY LINE_INTEGR_PHI_DOMEGA_DNU)
c
c     This subroutine is used by line_integr_phi_dnu to perform
c     the d\nu part of the double integral of a function over the
c     velocity dependent line profile along d\Omega and d\nu.
c
c     --------------------------------------------------------------
      function line_integr_phi_dnu(iline,ir,it,func,mu,phi,velocity)
      implicit none
      doubleprecision line_integr_phi_dnu,line_velo_phi
      doubleprecision func(SZ_LINEPROFILE_SMALL),mu,phi,velocity(3)
      integer iline,ir,it
c
c#include "common_grid.h"
#include "common_lines.h"
c
      doubleprecision dummy,phiprof,dnu1,phtot,wgt
      integer i,nfr
c
c     Prepare...
c
      nfr   = line_nfreq(iline)
      dummy = 0.d0
      phtot = 0.d0
      wgt   = 0.d0
c
c     Loop over frequencies
c
      do i=1,nfr
c
c         Find the line profile function
c
          wgt     = line_wgt(i,iline)
          phiprof = line_velo_phi(iline,i,ir,it,mu,phi,velocity)
c
c         Now add a step to the integral
c
          dummy = dummy + wgt * phiprof * func(i)
c
c         And add to the \int \phi d\nu for strict normalization
c
          phtot = phtot + wgt * phiprof
c
      enddo
c
c     Now make sure that the dummy is normalized to phtot, to
c     counter numerical errors in the interpolation that could
c     ruin the normalization of the line profile function
c
      dummy = dummy / phtot
c
c     Now return the integral of func*phi(\Omega,\nu) over d\nu
c
      line_integr_phi_dnu = dummy 
c
      return
      end



c     --------------------------------------------------------------
c                      RETURN THE PHI(OMEGA,NU)
c
c     This routine returns the line profile phi(\Omega,\nu) at a
c     given gridpoint ir,it, and frequency point inu, in a given
c     direction mu,phi, for a given velocity. It does so by taking
c     the line profile as stored in the locprof_phi array, and its
c     spline locprof_phi2. The formula:
c
c        \phi(\Omega,\nu) = \phi^0( \nu - \nu_0 \Omega\cdot\vec v )
c
c
c     NOTE: If the line profile is global (i.e. LINE_LOCAL_WIDTH 
c           not set), then the ir and it are irrelevant arguments
c           and need not be specified: just pass e.g. 0 for them.
c
c     --------------------------------------------------------------
      function line_velo_phi(iline,ilfr,ir,it,mu,phi,velocity)
      implicit none
      doubleprecision line_velo_phi,omega_dot_v,voigt_profile
c
#include "common_grid.h"
#include "common_lines.h"
c
      integer iline,ilfr,ir,it,nfr
      doubleprecision mu,phi,velocity(3),phiprof,aks
c
      doubleprecision dnu1,dvmu
c
c     Copy some information, for convenience
c
      nfr = line_nfreq(iline)
c
c     First compute the inner product between Omega and v
c
      dvmu = omega_dot_v(mu,phi,velocity)
c
c     Then compute the Doppler shifted dnu1
c
      dnu1 = line_dnu(ilfr,iline) - line_nu0(iline) * dvmu
c
c     Now interpolate the line profile at this shifted frequency
c
#     ifdef LINE_LOCAL_WIDTH
#     ifdef LINE_PROFILE_ARRAY
      write(*,*) 'ERROR: ARRAY FOR LINE PROFILE IS REMOVED'
      stop 911
#     endif
#     else
      write(*,*) 'ERROR: GLOBAL LINE PROFILE IS REMOVED'
      stop 911
#     endif
c
c     Local, compute on-the-fly
c
      aks = locprof_linewidth(iline,it,ir)
c
c         MODIFIED 19.03.07 (to make sure that only gt counts; 
c                   modified also in line.F)
c      if(aks.lt.1.d-5) then
      if(aks.le.0.d0) then
          aks = line_fix_linewidth(iline)
      endif
      phiprof = voigt_profile(linefreq(iline),0.d0,aks,dnu1)
c
c     Return the answer
c
      line_velo_phi = phiprof
c
      return
      end




c     --------------------------------------------------------------
c       PREPARE ALL THE LOCAL/GLOBAL LINE PROFILES AND ABUNDANCES
c
c     This routine computes the local line profiles and stores
c     then. Also it will produce splines for all of them, and store
c     these splines. This will make it easy for the transfer 
c     algorithm to quickly find the appropriate line profile 
c     function. Also the abundances will be copied locally from the
c     global abun_molec, unless the flag locabun_done is set,
c     meaning that some part of the program has already specified
c     the local abundances. 
c
c     NOTE: if LINE_LOCAL_WIDTH is not set, then only the global 
c           line profile is computed and the spline produced.
c     NOTE: The make_voigt_profile() routine in fact modifies the dnu 
c           as well, but since we always pass the same line_rangewidth, 
c           this is not dangerous.
c     --------------------------------------------------------------
      subroutine line_velo_prep_profiles_abun(icheck)
      implicit none
c
      integer icheck
c
#include "common_grid.h"
#include "common_lines.h"
#include "common_setup.h"
c
      doubleprecision aks,velmax,veltot,rwidthmin
      doubleprecision dvkmsmax,samplerr
      integer ir,it,ithetamax
      integer iline
c
      aksmax = 0.d0
      aksmin = 1.d5
c
c     Now the mirror things...
c
      ithetamax = irsi_frsizey
c
#     ifdef MIRROR_THETA
      ithetamax = (irsi_frsizey+1)/2
#     endif
#     ifdef RADGRID_ONEDIM
      ithetamax = 1
#     endif
c
c     Find out if we want local or global profiles
c
#ifdef LINE_LOCAL_WIDTH
c
c     Local profiles. Loop over space and lines.
c
      do ir=1,irsi_frsizex
          do it=1,ithetamax
c
c             Find the absolute value of the local velocity.
c             This is needed for the verification of whether
c             the line profile fits in the passband and whether
c             it is sampled sufficiently well.
c
              veltot = sqrt( medium_arr_velocity(1,it,ir)**2 +
     %                       medium_arr_velocity(2,it,ir)**2 +
     %                       medium_arr_velocity(3,it,ir)**2 )
c
c             Now produce the line profiles for all lines
c
              do iline=1,nlines
c
c                 First compute the line profile width
c     
                  aks = locprof_linewidth(iline,it,ir)
c
c                 MODIFIED 19.03.07 (to make sure that only gt counts; 
c                   modified also in line.F)
c                  if(aks.lt.1.d-5) then
                  if(aks.le.0.d0) then
                      aks = line_fix_linewidth(iline)
                  endif
                  if(aks.gt.aksmax) then 
                      aksmax = aks
                  endif
                  if(aks.lt.aksmin) then 
                      aksmin = aks
                  endif
c
c                 Now produce the local line profile and store it 
c                 into an array (if the LINE_PROFILE_ARRAY is set;
c                 if not, then the profile is assumed to be computed
c                 on-the-fly).
c
#                 ifdef LINE_PROFILE_ARRAY
c
                  call make_voigt_profile(linefreq(iline),0.d0,aks,
     %                line_rangewidth(iline),
     %                line_nfreq(iline),line_dnu(1,iline),
     %                locprof_phi(1,iline,it,ir),
     %                line_wgt(1,iline))
c      
c                 Produce a spline through the line profile
c
                  call spline(line_dnu(1,iline),
     %                locprof_phi(1,iline,it,ir), 
     %                line_nfreq(iline),0.d0,0.d0,
     %                locprof_phi2(1,iline,it,ir))
c
#                 endif /* LINE_PROFILE_ARRAY */
c
c                 Check if line profile fits into the passband
c
                  if(icheck.ne.0) then
                      if(0.5d0*aks+veltot*1.d-5.gt.
     %                     line_rangewidth(iline)) then
                          write(*,*) 'FATAL ERROR: ',
     %                               'Passband too small!!!'
                          write(*,*) '  Position: ir=',ir,
     %                         ' it=',it,' iline=',iline
                          stop 13
                      endif
                      if((1.4d0*aks+veltot*1.d-5.gt.
     %                     line_rangewidth(iline)).and.
     %                     (line_warn_lineout.eq.0)) then
                          write(*,*) 'WARNING: Passband too small'
                          write(*,*) '  Position: ir=',ir,' it=',it,
     %                         ' iline=',iline
                          write(*,*) '  Suppressing further messages...'
                          line_warn_lineout = 1
                      endif
c
c                     Check if line profile is sampled sufficiently
c                     accurately
c
                      if(line_rangewidth(iline)/line_nfreq(iline)
     %                     .gt.2.d0*aks) then
                          write(*,*) 'FATAL ERROR: Line freq sampling ',
     %                         'spacing is larger than line width'
                          write(*,*) 'This is assured to result in ',
     %                         'erroneous answers. Aborting...'
                          write(*,*) '  Position: ir=',ir,' it=',it,
     %                         ' iline=',iline
                          stop 13
                      endif
                      if((line_rangewidth(iline)/line_nfreq(iline)
     %                     .gt.0.7d0*aks).and.
     %                     (line_warn_badsampling.eq.0)) then
                          write(*,*) 'WARNING: Line freq sampling ',
     %                         'spacing is larger than line width'
                          write(*,*) '  Position: ir=',ir,' it=',it,
     %                         ' iline=',iline
                          write(*,*) '  Suppressing further messages...'
                          line_warn_badsampling = 1
                      endif
                  endif
              enddo
          enddo
      enddo
#else
c
c     Global profiles
c      
c     First compute the line profile
c     
      do iline=1,nlines
          aks = line_fix_linewidth(iline)
          if(aks.gt.aksmax) then 
              aksmax = aks
          endif
          call make_voigt_profile(linefreq(iline),0.d0,aks,
     %        line_rangewidth(iline),
     %        line_nfreq(iline),line_dnu(1,iline),
     %        line_phi(1,iline),
     %        line_wgt(1,iline))
c     
c         Produce a spline through the line profile
c     
          call spline(line_dnu(1,iline),
     %        line_phi(1,iline), 
     %        line_nfreq(iline),0.d0,0.d0,
     %        line_phi2(1,iline))
c      
      enddo
#endif
c
c     Now the local/global abundances
c
#ifdef LINE_LOCAL_ABUNDANCE
      if(locabun_done.eq.0) then
          do ir=1,irsi_frsizex
              do it=1,ithetamax
                  locabun_abund_mol(it,ir)  = abund_molec
                  locabun_abund_coll(it,ir) = abund_colpartner
              enddo
          enddo
      endif
#endif
c
c     Now print some info about line profile, passband and
c     frequency(velocity) sampling.
c
      velmax = 0.d0
#     ifdef LINE_VELOCITIES
      do ir=1,irsi_frsizex
          do it=1,ithetamax
              veltot = sqrt( medium_arr_velocity(1,it,ir)**2 +
     %                       medium_arr_velocity(2,it,ir)**2 +
     %                       medium_arr_velocity(3,it,ir)**2 )
              if(veltot.gt.velmax) then
                  velmax = veltot
              endif
          enddo
      enddo
#     endif
      write(*,'(a23,f7.3,a5)') 
     %     'Maximum gas velocity = ',velmax*1.d-5,' km/s'
      write(*,'(a,f7.3,a)') 'Minimum intr. line width   = ',
     %     aksmin,' km/s'
      write(*,'(a,f7.3,a)') 'Maximum intr. line width   = ',
     %     aksmax,' km/s'
      rwidthmin=1.d8
      dvkmsmax =-1.d0
      samplerr =-1.d0
      do iline=1,nlines
          if(rwidthmin.gt.line_rangewidth(iline)) then
              rwidthmin = line_rangewidth(iline)
          endif
          if(dvkmsmax.lt.line_rangewidth(iline)/line_nfreq(iline)) then
              dvkmsmax = line_rangewidth(iline)/line_nfreq(iline)
          endif
#         ifndef LINE_LOCAL_WIDTH
          if(samplerr.lt.line_rangewidth(iline)/
     %       (line_nfreq(iline)*line_fix_linewidth(iline))) then
              samplerr = line_rangewidth(iline)/
     %            (line_nfreq(iline)*line_fix_linewidth(iline))
          endif
#         endif
      enddo
      write(*,'(a,f7.3,a)') 'Passband (smallest)  = ',rwidthmin,' km/s'
      write(*,'(a,f7.3,a)') 'Delta V  (biggest)   = ',dvkmsmax,' km/s'
c
c     If line profiles are global, then we still have not yet
c     checked the sampling and the passband. Do it here...
c
#ifndef LINE_LOCAL_WIDTH
      if(icheck.ne.0) then
          if(0.5d0*aksmax+velmax*1.d-5.gt.rwidthmin) then
              write(*,*) 'FATAL ERROR: Passband too small!!!'
              stop 13
          endif
          if((1.4d0*aksmax+velmax*1.d-5.gt.rwidthmin).and.
     %         (line_warn_lineout.eq.0)) then
              write(*,*) 'WARNING: Passband too small'
              line_warn_lineout = 1
          endif
          if(samplerr.gt.2.d0) then
              write(*,*) 'FATAL ERROR: Line freq sampling ',
     %             'spacing is larger than line width'
              write(*,*) 'This is assured to result in ',
     %             'erroneous answers. Aborting...'
              stop 13
          endif
          if((samplerr.gt.0.7d0).and.(line_warn_badsampling.eq.0)) then
              write(*,*) 'WARNING: Line freq sampling ',
     %             'spacing is larger than line width'
              line_warn_badsampling = 1
          endif
      endif
c     
#endif
c
      end



c     --------------------------------------------------------------
c             COMPUTE THE ANGLE DEPENDENT SRC AND ALPHA
c     --------------------------------------------------------------
      subroutine line_velo_src_alp(inu,ir,it,mu,phi,src,alp)
      implicit none
      integer inu,ir,it
      doubleprecision mu,phi,src,alp
c
      doubleprecision line_velo_phi
c
#include "common_grid.h"
#include "common_radfield.h"
#include "common_source.h"
#include "common_lines.h"
#include "common_setup.h"
c
      integer itmax,ilevel,ilfr,iline
      doubleprecision dummy,nup,nlow,nmol,dens,temp,phiprof,molabun
c
      itmax   = irsi_frsizey
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
      if(it.gt.itmax) then
          write(*,*) 'ERROR: internal inconsistency in line.F!'
          write(*,*) it,itmax
          stop 643
      endif
c
c     Find the density and temperature
c
      call get_phys_var(ir,it)
      temp = medium_el_temp
      dens = medium_rho
c
c     Get the molecular abundance
c
#     ifdef LINE_LOCAL_ABUNDANCE
c
      molabun = locabun_abund_mol(it,ir)
#     else
      molabun = abund_molec
#     endif
c     
c     Find out in which line and at which freq of line
c     we are here
c
      iline = line_iline(inu)
      ilfr  = line_ilifreq(inu)
c
c     Find the upper and lower level population
c
      nup  = line_level_popul(lev_up(iline),it,ir)
      nlow = line_level_popul(lev_down(iline),it,ir)
c
c     Find now the real nup and nlow in 1/cm^3
c
      nmol = molabun * dens / (umass_av*1.6726d-24) 
      nup  = nup * nmol
      nlow = nlow * nmol
c
c     Find the value of the line profile at this point
c     and for this velocity
c
      phiprof = line_velo_phi(iline,ilfr,ir,it,
     %               mu,phi,medium_arr_velocity(1,it,ir))
c
c     The j_nu is now:
c
c                   h nu_0
c            j_nu = ------ n_up A_ud phi(Omega,nu)
c                    4 pi
c
      src  = 5.27296241956d-28 * linefreq(iline) * nup * 
     %          Aud_field(iline,it,ir) * phiprof
c
c     The alpha_nu is now:
c
c                       h nu_0
c            alpha_nu = ------ ( n_down B_du - n_up B_ud ) phi(Omega,nu)
c                        4 pi
c
      alp  = 5.27296241956d-28 * linefreq(iline) * phiprof *
     %          ( nlow * Bdu_field(iline,it,ir) - 
     %             nup * Bud_field(iline,it,ir) )
c
c     Now check if the opacity remains positive. If not, then
c     we are dealig with a maser, and we should take measures 
c     to keep the stuff from blowing up violently.
c
#     ifndef LINE_NO_LIMIT_MASERS
      if(alp.lt.0.d0) then 
          alp = 1.e-37
          line_warn_opacity_neg = 1
          line_warn_neg_opacity = 1
      endif
#     endif
c
      return
      end



c     --------------------------------------------------------------
c              COPY ABUNDANCES TO THE GLOBAL ARRAYS
c     --------------------------------------------------------------
      subroutine copy_abundances(ir,it)
      implicit none
      integer ir,it
c
#include "common_grid.h"
#include "common_lines.h"
c
      abund_molec      = locabun_abund_mol(it,ir)  
      abund_colpartner = locabun_abund_coll(it,ir)      
c     
      end


c     --------------------------------------------------------------
c           CHECK IF THE VELOCITY STEPS ARE NOT TOO LARGE
c     --------------------------------------------------------------
#     ifdef LINE_VELOCITIES
      subroutine check_velo_step_linewidth()
      implicit none
c
#include "common_grid.h"
#include "common_lines.h"
#include "common_setup.h"
c
      doubleprecision dvr,dvtheta,dvphi,dvmax,aks
      integer iline,ir,it,ithetamax,k
c
c     Now the mirror things...
c
#ifdef MIRROR_THETA
      ithetamax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      ithetamax = 1
#endif
c      
      do ir=2,irsi_frsizex-1
          do it=1,ithetamax
              dvr     = 0.5d0 * ( medium_arr_velocity(1,it,ir+1) - 
     %                            medium_arr_velocity(1,it,ir-1) )
              dvtheta = 0.5d0 * ( medium_arr_velocity(2,it,ir+1) -
     %                            medium_arr_velocity(2,it,ir-1) )
              dvphi   = 0.5d0 * ( medium_arr_velocity(3,it,ir+1) -
     %                            medium_arr_velocity(3,it,ir-1) )
              dvmax   = max(dvr,dvtheta,dvphi)
              aksmin = 1.d33
              do iline=1,nlines
                  aks    = locprof_linewidth(iline,it,ir)
c
c         MODIFIED 19.03.07 (to make sure that only gt counts; 
c                   modified also in line.F)
c                  if(aks.lt.1.d-5) then
                  if(aks.le.0.d0) then
                      aks = line_fix_linewidth(iline)
                  endif
                  if(aks.lt.aksmin) then
                      aksmin = aks
                  endif
              enddo
c             KMP: No need for this warning in RADLite, since this is automatically taken care of
c              if((dvmax.gt.2.0d0*aksmin*1.d5).and.
c     %           (line_warn_velostep.ne.2)) then
c                  write(*,*) '-----------------------------------'
c                  write(*,*) 'HIGH ALERT: Velocity step is MUCH'
c                  write(*,*) '  larger than line width!!!!'
c                  write(*,*) '  Position: ir=',ir,' it=',it
c                  write(*,*) '  Factor = ',dvmax/(aksmin*1.d5)
c                  write(*,*) '  Suppressing further messages...'                  
c                  write(*,*) '-----------------------------------'
c                  line_warn_velostep = 2
c              endif
c              if((dvmax.gt.0.7d0*aksmin*1.d5).and.
c     %           (line_warn_velostep.eq.0)) then
c                  write(*,*) 'WARNING: Velocity step ',
c     %                          'is larger than line width'
c                  write(*,*) '  Position: ir=',ir,' it=',it
c                  write(*,*) '  Suppressing further messages...'                  
c                  line_warn_velostep = 1
c              endif
          enddo
      enddo
c
      end
#     endif


#ifdef LINE_SRC_INTERPOL_LEVPOP
c     --------------------------------------------------------------
c             GET THE LEVELPOPS AND TEMPERATURE AT POINT
c
c     This routine finds the N*nup,N*ndown and the temperature
c     at a given gridpoint.
c
c     REMARK:
c       This routine was meant for use by find_levpop_and_temp(), but
c       I already solved the original problem by another method. 
c       However, since in the future I should do the interpolation
c       for line transfer in this neat way anyway, I keep this routine
c       here, though inactivated. (13-May-1999)
c     --------------------------------------------------------------
      subroutine get_levpop_temp(iline,ir,it,nup,ndown,temp)
      implicit none
      integer iline,ir,it
      doubleprecision nup,ndown,temp
c
#include "common_grid.h"
#include "common_lines.h"
#include "common_setup.h"
c
      doubleprecision dummy,nmol,molabun,dens
      integer itmax
c
      itmax   = irsi_frsizey
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
      if(it.gt.itmax) then
          write(*,*) 'ERROR: internal inconsistency in line.F!'
          stop 645
      endif
c
c     Find the density and temperature
c
      call get_phys_var(ir,it)
      temp = medium_el_temp
      dens = medium_rho
c
c     Get the molecular abundance
c
#     ifdef LINE_LOCAL_ABUNDANCE
      molabun = locabun_abund_mol(it,ir)
#     else
      molabun = abund_molec
#     endif
c
c     Find the upper and lower level population
c
      nup   = line_level_popul(lev_up(iline),it,ir)
      ndown = line_level_popul(lev_down(iline),it,ir)
c
c     Find now the real nup and nlow in 1/cm^3
c
      nmol  = molabun * dens / (umass_av*1.6726d-24) 
      nup   = nup * nmol
      ndown = ndown * nmol
c
      return
      end
#endif




c     --------------------------------------------------------------
c     COMPUTE THE DUST THERMAL EMISSION/ABSORPTION FOR LINE TRANSFER
c
c     We compute this only for a complete line, not for each 
c     individual frequency bin within the line passband. 
c     Note that only the dust thermal emission/absorption is done
c     here, not the scattering. That is a much more difficult thing.
c
c     >>>> ADDED KEES/KLAUS MARCH 2007 <<<<
c
c     --------------------------------------------------------------
      subroutine line_dust_compute_src_alp(iline,ir,it,src,alp)
      implicit none
      integer iline,ir,it
      doubleprecision bplanck
#ifndef LINE_FREQDEP_DUST
      doubleprecision src,alp
#else /* LINE_FREQDEP_DUST */
      doubleprecision src(SZ_LINEPROFILE),alp(SZ_LINEPROFILE)
#endif /* LINE_FREQDEP_DUST */

c
#include "common_grid.h"
#include "common_radfield.h"
#include "common_source.h"
#include "common_dust.h"
#include "common_lines.h"
#include "common_setup.h"
c
      integer itmax,inudust,isize,ispec,k
      doubleprecision dummy,dens,temp,freq,kappawgt,wgt
c
      itmax   = irsi_frsizey
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
      if(it.gt.itmax) then
          write(*,*) 'ERROR: internal inconsistency in line.F!'
          stop 644
      endif
c     
c     Find the density and temperature
c
      call get_phys_var(ir,it)
      temp = medium_el_temp
      dens = medium_rho
c
c     NEW 27.03.06: Combine line and dust, a check
c
      if(freq_nr.ne.icont_freq_nr) then
          write(*,*) 'INTERNAL ERROR 8292:',
     %         ' Nr of freqs not OK...'
          write(*,*) freq_nr, icont_freq_nr
          stop
      endif
c
c     >>>> ADDED KEES/KLAUS MARCH 2007 <<<<
c
#ifndef LINE_FREQDEP_DUST
c
c     ORIGINAL METHOD: Dust continuum constant over each line bandwidth
c
c     Find out which frequency we are...
c
      freq = linefreq(iline)
c
c     ...and figure out where in the line-dust frequency array
c     this is...
c
      call hunt(cont_freq_nu,freq_nr,freq,inudust)
      if(inudust.eq.0) then
          src = 0.d0
          alp = 0.d0
      elseif(inudust.eq.freq_nr) then
          src = 0.d0
          alp = 0.d0
      else
          wgt = (freq-cont_freq_nu(inudust)) /
     %         (cont_freq_nu(inudust+1)-cont_freq_nu(inudust))
          src = 0.d0
          alp = 0.d0
c
c         Add the thermal emission / absorption
c
          do ispec=1,dust_nr_species
              do isize=1,dust_nr_size(ispec)
                  kappawgt = 
     %              wgt * dust_kappawgt_abs(inudust+1,1,isize,ispec) +
     %             (1.d0-wgt) * dust_kappawgt_abs(inudust,1,isize,ispec)
                  src = src + dust_rho(ispec,it,ir) * kappawgt * 
     %                   bplanck(dust_temp(isize,ispec,it,ir),freq)
                  alp = alp + dust_rho(ispec,it,ir) * kappawgt 
              enddo
          enddo
c
c         Add the scattering/PAH source / extinction
c
          src = src + wgt * scati_src(inudust+1,it,ir)
     %               + (1.d0-wgt) * scati_src(inudust,it,ir)
c          src = src + (wgt * scati_src(inudust+1,it,ir)
c     %               + (1.d0-wgt) * scati_src(inudust,it,ir))*0.5d0
c3333333
          do ispec=1,dust_nr_species
              do isize=1,dust_nr_size(ispec)
                  kappawgt = 
     %                 wgt * dust_kappawgt_scat(inudust+1,1,
     %                       isize,ispec) +
     %                 (1.d0-wgt) * dust_kappawgt_scat(inudust,1,
     %                       isize,ispec)
                  alp = alp + dust_rho(ispec,it,ir) * kappawgt 
              enddo
          enddo
      endif
c
#else /* ifndef LINE_FREQDEP_DUST */
c
c     NEW METHOD: Dust continuum variable over each line bandwidth
c     Note that this doesn't work right now because the line_nu array has
c     not yet been initialized with proper values. To do this, it is necessary
c     to redesign the code such that the global frequency is not overwritten! KMP January 2009
c
      do k=1,line_nfreq(iline)
c
c         Find out which frequency we are...
c
          freq = line_nu(k,iline)

c
c         ...and figure out where in the line-dust frequency array
c         this is...
c
          call hunt(freq_nu,freq_nr,freq,inudust)
          if(inudust.eq.0) then
              src(k) = 0.d0
              alp(k) = 0.d0
          elseif(inudust.eq.freq_nr) then
              src(k) = 0.d0
              alp(k) = 0.d0
          else
              wgt = (freq-freq_nu(inudust)) /
     %             (freq_nu(inudust+1)-freq_nu(inudust))
              src(k) = 0.d0
              alp(k) = 0.d0
c
c             Add the thermal emission / absorption
c
              do ispec=1,dust_nr_species
                  do isize=1,dust_nr_size(ispec)
                      kappawgt = 
     %                     wgt * dust_kappawgt_abs(inudust+1,
     %                                             1,isize,ispec) +
     %                     (1.d0-wgt) * dust_kappawgt_abs(inudust,
     %                                             1,isize,ispec)
                      src(k) = src(k) + dust_rho(ispec,it,ir)*kappawgt* 
     %                     bplanck(dust_temp(isize,ispec,it,ir),freq)
                      alp(k) = alp(k) + dust_rho(ispec,it,ir)*kappawgt 
                  enddo
              enddo
c
c             Add the scattering/PAH source / extinction
c
              src(k) = src(k) + wgt * scati_src(inudust+1,it,ir)
     %               + (1.d0-wgt) * scati_src(inudust,it,ir)
              do ispec=1,dust_nr_species
                  do isize=1,dust_nr_size(ispec)
                      kappawgt = 
     %                     wgt * dust_kappawgt_scat(inudust+1,1,
     %                     isize,ispec) +
     %                     (1.d0-wgt) * dust_kappawgt_scat(inudust,1,
     %                     isize,ispec)
                      alp(k) = alp(k) + dust_rho(ispec,it,ir)*kappawgt
                  enddo
              enddo
          endif
      enddo
c

#endif /* ifndef LINE_FREQDEP_DUST */

c
      return
      end



c     --------------------------------------------------------------
c             COMPUTE ALL SRC AND ALP FOR DUST IN LINES
c      
c      >>>> ADDED KEES/KLAUS MARCH 2007 <<<<
c
c     --------------------------------------------------------------
      subroutine global_prepare_line_dust()
      implicit none
c
#include "common_grid.h"
#include "common_lines.h"
#include "common_dust.h"
c
      integer ir,it,inu,itmax,ispec,isize,iline,k
      doubleprecision dummy,dens,temp,nu,bplanck
#ifndef LINE_FREQDEP_DUST
      doubleprecision src,alp
#else /* LINE_FREQDEP_DUST */
      doubleprecision src(SZ_LINEPROFILE),alp(SZ_LINEPROFILE)
#endif /* LINE_FREQDEP_DUST */
c
c     The usual spatial grid stuff
c
      itmax = irsi_frsizey
c
#ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#endif
#ifdef RADGRID_ONEDIM
      itmax = 1
#endif
c
c     Now loop over space and lines
c
#ifndef LINE_FREQDEP_DUST
       write(*,*) 'WARNING: Using grey dust for the line!'
       write(*,*) 'If you want frequency dependent dust'
       write(*,*) 'set LINE_FREQDEP_DUST in configure.h'
#endif /* LINE_FREQDEP_DUST */

      do ir=1,irsi_frsizex
          do it=1,itmax
              do iline=1,nlines
                  call line_dust_compute_src_alp(iline,ir,it,src,alp)
#ifndef LINE_FREQDEP_DUST
                  line_dust_src(1,iline,it,ir) = src
                  line_dust_alp(1,iline,it,ir) = alp
#else /* LINE_FREQDEP_DUST */
                  do k=1,line_nfreq(iline)
                      line_dust_src(k,iline,it,ir) = src(k)
                      line_dust_alp(k,iline,it,ir) = alp(k)
                  enddo
#endif /* LINE_FREQDEP_DUST */
              enddo
          enddo
      enddo
c
c     Now set the flag that this task is done
c
      line_warn_dust_done = 1
c
      return
      end



c     --------------------------------------------------------------
c         THE ROUTINE FOR QUICK RETRIEVAL OF SRC AND ALP OF DUST
c      
c      >>>> ADDED KEES/KLAUS MARCH 2007 <<<<
c
c     --------------------------------------------------------------
      subroutine line_dust_src_alp(inu,ir,it,src,alp)
      implicit none
      integer inu,ir,it
      doubleprecision src,alp
c
#include "common_grid.h"
#include "common_lines.h"
c
      doubleprecision dummy
      integer iline,ilfr
c
      if(line_warn_dust_done.eq.0) then
          write(*,*) 'ERROR: Dust sources for line transfer'
          write(*,*) '       not yet prepared...'
          stop 14
      endif
c
      iline = line_iline(inu)
#ifndef LINE_FREQDEP_DUST
      src   = line_dust_src(1,iline,it,ir)
      alp   = line_dust_alp(1,iline,it,ir)
#else /* LINE_FREQDEP_DUST */
      ilfr  = line_ilifreq(inu)
      src   = line_dust_src(ilfr,iline,it,ir)
      alp   = line_dust_alp(ilfr,iline,it,ir)
#endif /* LINE_FREQDEP_DUST */
c
      return
      end



#endif /* LINE_VELOCITIES */
c     ============== END OF VELO DEPENDENT LINE PART ===============



c     --------------------------------------------------------------
c               SET THE STAR BOUNDARY CONDITION FOR LINES
c                         FOR THE TELESCOPE
c
c     NEW: 30-03-06: Special for the telescope routine, and for
c                    a single line passband only.
c     --------------------------------------------------------------
      subroutine line_init_passband_starbc(nfr,itr_freq)
      implicit none
c
      integer nfr,itr_freq
c
#include "common_grid.h"
#include "common_boundary.h"
#include "common_lines.h"
c
      doubleprecision freq,wgt
      integer inu,inudust
c
c     First a check
c
      if(icont_freq_nr.eq.0) then
          write(*,*) 'ERROR: Cannot use line stellar BC without having'
          write(*,*) '       read the stellar spectrum.'
          stop
      endif
      if(itr_freq+nfr-1.gt.FRSIZE_FREQ) stop 83931
      if(itr_freq+nfr-1.gt.SZ_LINEPROFILE) stop 83932
c
c     Loop over wavelength
c
      do inu=itr_freq,itr_freq+nfr-1
c
c         Copy
c
          freq = freq_nu(inu)
c
c         ...figure out where in the line-dust frequency array
c         this is, and plug the stellar spectrum there...
c
          call hunt(cont_freq_nu,icont_freq_nr,freq,inudust)
          
          if(inudust.eq.0) then
              radbnd_starspec_line(inu) = 0.d0
          elseif(inudust.eq.icont_freq_nr) then
              radbnd_starspec_line(inu) = 0.d0
          else
              wgt = (freq-cont_freq_nu(inudust)) /
     %             (cont_freq_nu(inudust+1)-cont_freq_nu(inudust))
              radbnd_starspec_line(inu) = 
     %             (1.d0-wgt)*radbnd_cont_starspec(inudust) + 
     %             wgt*radbnd_cont_starspec(inudust+1)
          endif
      enddo
c
      end


c     --------------------------------------------------------------
c               SET THE OUTER BOUNDARY CONDITION FOR LINES
c                         FOR THE TELESCOPE
c
c     NEW: 28-03-06: Special for the telescope routine, and for
c                    a single line passband only.
c     --------------------------------------------------------------
      subroutine line_init_passband_outerbc(nfr,itr_freq)
      implicit none
c
      integer nfr,itr_freq
c
#include "common_grid.h"
#include "common_boundary.h"
#include "common_lines.h"
c
      doubleprecision freq,wgt
      integer inu,inudust
c
c     First a check
c
      if(icont_freq_nr.eq.0) then
          write(*,*) 'ERROR: Cannot use line outer BC without having'
          write(*,*) '       read the interstellar spectrum.'
          stop
      endif
      if(itr_freq+nfr-1.gt.FRSIZE_FREQ) stop 83921
      if(itr_freq+nfr-1.gt.SZ_LINEPROFILE) stop 83922
c
c     Now a loop over the wavelength
c
      do inu=itr_freq,itr_freq+nfr-1
c
c         Copy: Klaus: Is this correct that this is freq_nu and not
c                      cont_freq_nu?
c
          freq = freq_nu(inu)
c
c         ...figure out where in the line-dust frequency array
c         this is, and plug the stellar spectrum there...
c
          call hunt(cont_freq_nu,icont_freq_nr,freq,inudust)
          if(inudust.eq.0) then
              radbnd_interstellfield(inu) = 0.d0
          elseif(inudust.eq.freq_nr) then
              radbnd_interstellfield(inu) = 0.d0
          else
              wgt = (freq-cont_freq_nu(inudust)) /
     %             (cont_freq_nu(inudust+1)-cont_freq_nu(inudust))
              radbnd_interstellfield(inu) = 
     %             (1.d0-wgt)*radbnd_cont_interstellfield(inudust)+
     %             wgt*radbnd_cont_interstellfield(inudust+1)
          endif
      enddo
c
      end


c     ==============================================================
c           A FEW ROUTINES THAT MAY BE USED WHEN LINES ARE OFF
c     ==============================================================


c     --------------------------------------------------------------
c                      READ UMASS FROM LINE.INP
c
c     This routine may be used if for instance line processes are
c     switched off, but the dust module wants to be able to convert
c     n_H2 to rho. So it reads the line.inp to find the umass.
c     --------------------------------------------------------------
      subroutine read_line_umass()
      implicit none
c
#include "common_grid.h"
#include "common_lines.h"
c
      integer dum
c
      open(unit=3,file='line.inp',status='unknown',err=701)
      read(3,*) dum
      read(3,*) dum
      read(3,*) umass_av
      close(3)
      goto 710
  701 continue
      write(*,*) 'ERROR: Could not find line.inp'
      write(*,*) '       Need this file to find umass_av'
      stop 13
  710 continue
c
      end



c     ==============================================================
c             ROUTINES FOR THE PRECISION LINE INTEGRATION
c
c                 >>>> ADDED KEES/KLAUS MARCH 2007 <<<<
c
c     These are the routines necessary for the precision integration
c     for the telescope.F routines that allow very narrow intrinsic
c     line width without gridding problems. This method only works
c     for the formal integration of the RT equation, i.e. it does
c     not work for the non-LTE stuff, because for non-LTE we need
c     (in addition to the formal RT integration) also an accurate
c     integration of the intensity as a function of Omega and nu,
c     i.e. a precise integration to obtain $\bar J$ of the line.
c     This integration suffers from the same gridding problem as
c     the formal RT integration, but this time not in R- and Theta-
c     grid, but in mu- and phi-grid. The line may be narrower than
c     the mu- and phi- grid spacing, as well as the nu-spacing. 
c     
c     ==============================================================

c     --------------------------------------------------------------
c                 GET THE DUST SRC,ALP AND LINE PROJ VELO
c     --------------------------------------------------------------
      subroutine get_line_dust_values(icr,ir,it,dr,dt,mu,phi,iline,inu,
     %                 src_dust,alp_dust,dvmu,linewidth,Nup,Ndown)
      implicit none
c
      integer ir,it,iline,inu,icr,idr,idt,ir0,it0
      doubleprecision dr,dt,mu,phi,src_dust,alp_dust,dvmu,drr,dtr
      doubleprecision linewidth,Nup,Ndown
c
#include "common_lines.h"
#include "common_grid.h"
#include "common_setup.h"
c
      integer i,ii
      doubleprecision src(0:3),alp(0:3),lw(0:3),n(0:3)
      doubleprecision vel(0:3),velocity(3)
      doubleprecision rclip(0:1),thetaclip(0:1)
      integer indexr(0:1),indext(0:1)
      doubleprecision theta,r
      doubleprecision omega_dot_v
      doubleprecision molpg
c
      if((dr.lt.0.d0).or.(dr.gt.1.d0)) then
          write(*,*) 'ERROR: Erroneous dr found'
          stop 6024
      endif
      if((dt.lt.0.d0).or.(dt.gt.1.d0)) then
          write(*,*) 'ERROR: Erroneous dt found'
          stop 6023
      endif
c
c     Define molecules/gram
c

      molpg = 1.d0/(umass_av*1.6726d-24)
c
c     Find the radial stencil
c

c
c    CHECK: What are r and theta actually used for???
c
      if(dr.gt.0.d0) then
         indexr(0) = ir
         indexr(1) = ir + 1
         rclip(0)  = rsi_x_c(indexr(0),1)
         rclip(1)  = rsi_x_c(indexr(1),1)
         if(indexr(1).gt.irsi_frsizex) then
            indexr(1) = irsi_frsizex
         endif
         r = rclip(0) + dr * ( rclip(1) - rclip(0) )
      else
         indexr(0) = ir
         indexr(1) = ir - 1
         rclip(0)  = rsi_x_c(indexr(0),1)
         rclip(1)  = rsi_x_c(indexr(1),1)
         if(indexr(1).lt.1) then
            indexr(1) = 1
         endif
         r = rclip(0) + dr * ( rclip(1) - rclip(0) )
      endif
c
c     Now complete the 2-point Theta index stencil.
c
      if(dt.gt.0.d0) then
         indext(0) = it 
         indext(1) = it + 1
       else
         indext(0) = it 
         indext(1) = it - 1
      endif
      thetaclip(0) = rsi_x_c(indext(0),2)
      thetaclip(1) = rsi_x_c(indext(1),2)
      theta = thetaclip(0) + dt * ( thetaclip(1) - thetaclip(0) )
      indext(0)   = ridx_it(indext(0))
      indext(1)   = ridx_it(indext(1))
c
#ifndef LINE_FREQDEP_DUST
      ii     = 1
#else
      ii     = line_ilifreq(inu)
#endif
#ifndef LINE_LOCAL_WIDTH
      write(*,*) 'ERROR: The compilation with global line width'
      write(*,*) '       (i.e. NOT setting LINE_LOCAL_WIDTH)'
      write(*,*) '       is obsolete and disactivated.'
      write(*,*) '       >>>> ADDED KEES/KLAUS MARCH 2007 <<<<'
      stop 10245
#endif
c
      if(icr.eq.1) then
c
c         The point 'P' lies on an R=const gridline
c
          src(0)   = line_dust_src(ii,iline,indext(0),indexr(0))
          src(1)   = line_dust_src(ii,iline,indext(1),indexr(0))
          src_dust = (1.d0-dt)*src(0) + dt*src(1)
          alp(0)   = line_dust_alp(ii,iline,indext(0),indexr(0))
          alp(1)   = line_dust_alp(ii,iline,indext(1),indexr(0))
          alp_dust = (1.d0-dt)*alp(0) + dt*alp(1)
c
          lw(0)    = locprof_linewidth(iline,indext(0),indexr(0))
          lw(1)    = locprof_linewidth(iline,indext(1),indexr(0))
          linewidth= (1.d0-dt)*lw(0) + dt*lw(1)
c
          N(0)     = line_level_popul(lev_up(iline),
     %                              indext(0),indexr(0))
     %                   * locabun_abund_mol(indext(0),indexr(0))
     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
          N(1)     = line_level_popul(lev_up(iline),
     %                                indext(1),indexr(0))
     %                   * locabun_abund_mol(indext(1),indexr(0))
     %                   * medium_arr_rho(indext(1),indexr(0))*molpg
          Nup      = (1.d0-dt)*N(0) + dt*N(1)
c          write(*,*) N(0),N(1),Nup,icr
c         55555         
          N(0)     = line_level_popul(lev_down(iline),
     %                                indext(0),indexr(0))
     %                   * locabun_abund_mol(indext(0),indexr(0))
     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
          N(1)     = line_level_popul(lev_down(iline),
     %                                indext(1),indexr(0))
     %                   * locabun_abund_mol(indext(1),indexr(0))
     %                   * medium_arr_rho(indext(1),indexr(0))*molpg
          Ndown    = (1.d0-dt)*N(0) + dt*N(1)
c
          do i=1,3
              vel(0)      = medium_arr_velocity(i,indext(0),indexr(0))
              vel(1)      = medium_arr_velocity(i,indext(1),indexr(0))
              velocity(i) = (1.d0-dt)*vel(0) + dt*vel(1)
          enddo
c
      elseif(icr.eq.2) then
c
c         The point 'P' lies on an Theta=const gridline
c
          src(0)   = line_dust_src(ii,iline,indext(0),indexr(0))
          src(1)   = line_dust_src(ii,iline,indext(0),indexr(1))
          src_dust = (1.d0-dr)*src(0) + dr*src(1)
          alp(0)   = line_dust_alp(ii,iline,indext(0),indexr(0))
          alp(1)   = line_dust_alp(ii,iline,indext(0),indexr(1))
          alp_dust = (1.d0-dr)*alp(0) + dr*alp(1)
c
          lw(0)    = locprof_linewidth(iline,indext(0),indexr(0))
          lw(1)    = locprof_linewidth(iline,indext(0),indexr(1))
          linewidth= (1.d0-dr)*lw(0) + dr*lw(1)
c
          N(0)     = line_level_popul(lev_up(iline),
     %                                indext(0),indexr(0))
     %                   * locabun_abund_mol(indext(0),indexr(0))
     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
          N(1)     = line_level_popul(lev_up(iline),
     %                                indext(0),indexr(1))
     %                   * locabun_abund_mol(indext(0),indexr(1))
     %                   * medium_arr_rho(indext(0),indexr(1))*molpg
          Nup      = (1.d0-dr)*N(0) + dr*N(1)
c          write(*,*) N(0),N(1),icr
c         55555         
          N(0)     = line_level_popul(lev_down(iline),
     %                                indext(0),indexr(0))
     %                   * locabun_abund_mol(indext(0),indexr(0))
     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
          N(1)     = line_level_popul(lev_down(iline),
     %                                indext(0),indexr(1))
     %                   * locabun_abund_mol(indext(0),indexr(1))
     %                   * medium_arr_rho(indext(0),indexr(1))*molpg
          Ndown    = (1.d0-dr)*N(0) + dr*N(1)
c
          do i=1,3
              vel(0)      = medium_arr_velocity(i,indext(0),indexr(0))
              vel(1)      = medium_arr_velocity(i,indext(0),indexr(1))
              velocity(i) = (1.d0-dr)*vel(0) + dr*vel(1)
          enddo
c
      else
c
c         The point 'P' lies somewhere in between all gridlines
c
          src(0)   = line_dust_src(ii,iline,indext(0),indexr(0))
          src(1)   = line_dust_src(ii,iline,indext(1),indexr(0))
          src(2)   = line_dust_src(ii,iline,indext(0),indexr(1))
          src(3)   = line_dust_src(ii,iline,indext(1),indexr(1))
          src_dust = (1.d0-dr) * ((1.d0-dt)*src(0) + dt*src(1)) 
     %              +       dr * ((1.d0-dt)*src(2) + dt*src(3)) 
          alp(0)   = line_dust_alp(ii,iline,indext(0),indexr(0))
          alp(1)   = line_dust_alp(ii,iline,indext(1),indexr(0))
          alp(2)   = line_dust_alp(ii,iline,indext(0),indexr(1))
          alp(3)   = line_dust_alp(ii,iline,indext(1),indexr(1))
          alp_dust = (1.d0-dr) * ((1.d0-dt)*alp(0) + dt*alp(1)) 
     %              +       dr * ((1.d0-dt)*alp(2) + dt*alp(3)) 
c
          lw(0)    = locprof_linewidth(iline,indext(0),indexr(0))
          lw(1)    = locprof_linewidth(iline,indext(1),indexr(0))
          lw(2)    = locprof_linewidth(iline,indext(0),indexr(1))
          lw(3)    = locprof_linewidth(iline,indext(1),indexr(1))
          linewidth= (1.d0-dr) * ((1.d0-dt)*lw(0) + dt*lw(1)) 
     %              +       dr * ((1.d0-dt)*lw(2) + dt*lw(3)) 
c
          N(0)     = line_level_popul(lev_up(iline),
     %                                indext(0),indexr(0))
     %                   * locabun_abund_mol(indext(0),indexr(0))
     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
          N(1)     = line_level_popul(lev_up(iline),
     %                                indext(1),indexr(0))
     %                   * locabun_abund_mol(indext(1),indexr(0))
     %                   * medium_arr_rho(indext(1),indexr(0))*molpg
          N(2)     = line_level_popul(lev_up(iline),
     %                                indext(0),indexr(1))
     %                   * locabun_abund_mol(indext(0),indexr(1))
     %                   * medium_arr_rho(indext(0),indexr(1))*molpg
          N(3)     = line_level_popul(lev_up(iline),
     %                                indext(1),indexr(1))
     %                   * locabun_abund_mol(indext(1),indexr(1))
     %                   * medium_arr_rho(indext(1),indexr(1))*molpg
          Nup      = (1.d0-dr) * ((1.d0-dt)*N(0) + dt*N(1)) 
     %              +       dr * ((1.d0-dt)*N(2) + dt*N(3)) 
          N(0)     = line_level_popul(lev_down(iline),
     %                                indext(0),indexr(0))
     %                   * locabun_abund_mol(indext(0),indexr(0))
     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
          N(1)     = line_level_popul(lev_down(iline),
     %                                indext(1),indexr(0))
     %                   * locabun_abund_mol(indext(1),indexr(0))
     %                   * medium_arr_rho(indext(1),indexr(0))*molpg
          N(2)     = line_level_popul(lev_down(iline),
     %                                indext(0),indexr(1))
     %                   * locabun_abund_mol(indext(0),indexr(1))
     %                   * medium_arr_rho(indext(0),indexr(1))*molpg
          N(3)     = line_level_popul(lev_down(iline),
     %                                indext(1),indexr(1))
     %                   * locabun_abund_mol(indext(1),indexr(1))
     %                   * medium_arr_rho(indext(1),indexr(1))*molpg
          Ndown    = (1.d0-dr) * ((1.d0-dt)*N(0) + dt*N(1)) 
     %              +       dr * ((1.d0-dt)*N(2) + dt*N(3)) 
c
          do i=1,3
              vel(0)      = medium_arr_velocity(i,indext(0),indexr(0))
              vel(1)      = medium_arr_velocity(i,indext(1),indexr(0))
              vel(2)      = medium_arr_velocity(i,indext(0),indexr(1))
              vel(3)      = medium_arr_velocity(i,indext(1),indexr(1))
              velocity(i) = (1.d0-dr) * ((1.d0-dt)*vel(0)+dt*vel(1)) 
     %              +              dr * ((1.d0-dt)*vel(2)+dt*vel(3)) 
          enddo
c
      endif
c
c     First compute the inner product between Omega and v
c
      dvmu = omega_dot_v(mu,phi,velocity)
c
c     Return src_dust, alp_dust and dvmu
c
      return
      end

c
c Higher order interpolation version of the original get_line_dust_values
c
c$$$      subroutine get_line_dust_values_2nd(icr,ir,it,dr,dt,mu,phi,iline,
c$$$     %     inu,src_dust,alp_dust,dvmu,linewidth,Nup,Ndown)
c$$$      implicit none
c$$$c
c$$$      integer ir,it,iline,inu,icr
c$$$      doubleprecision dr,dt,mu,phi,src_dust,alp_dust,dvmu
c$$$      doubleprecision linewidth,Nup,Ndown
c$$$c
c$$$#include "common_lines.h"
c$$$#include "common_grid.h"
c$$$#include "common_setup.h"
c$$$c
c$$$      integer i,ii
c$$$      doubleprecision src(0:3),alp(0:3),lw(0:3),n(0:3)
c$$$      doubleprecision vel(0:3),velocity(3)
c$$$      doubleprecision rclip(0:3),thetaclip(0:3)
c$$$      integer indexr(0:3),indext(0:3)
c$$$      doubleprecision theta,r
c$$$      doubleprecision omega_dot_v
c$$$      doubleprecision molpg
c$$$c
c$$$      if((dr.lt.0.d0).or.(dr.gt.1.d0)) then
c$$$          write(*,*) 'ERROR: Erroneous dr found'
c$$$          stop 6024
c$$$      endif
c$$$      if((dt.lt.0.d0).or.(dt.gt.1.d0)) then
c$$$          write(*,*) 'ERROR: Erroneous dt found'
c$$$          stop 6023
c$$$      endif
c$$$c
c$$$c     Define molecules/gram
c$$$c
c$$$
c$$$      molpg = 1.d0/(umass_av*1.6726d-24)
c$$$c
c$$$c     Find the radial stencil
c$$$c
c$$$
c$$$c
c$$$c    CHECK: What are r and theta actually used for???
c$$$c
c$$$      indexr(0) = ir
c$$$      indexr(1) = ir + 1
c$$$      indexr(2) = ir + 2
c$$$      
c$$$      rclip(0)  = rsi_x_c(indexr(0),1)
c$$$      rclip(1)  = rsi_x_c(indexr(1),1)
c$$$      rclip(2)  = rsi_x_c(indexr(2),1)
c$$$      if(indexr(1).gt.irsi_frsizex) then
c$$$          indexr(1) = irsi_frsizex
c$$$      endif
c$$$      if(indexr(2).gt.irsi_frsizex) then
c$$$          indexr(2) = irsi_frsizex
c$$$      endif
c$$$      r = rclip(0) + dr * ( rclip(1) - rclip(0) )
c$$$c
c$$$c     Now the theta stencil
c$$$c
c$$$      indext(0) = it 
c$$$      indext(1) = it + 1
c$$$      thetaclip(0) = rsi_x_c(indext(0),2)
c$$$      thetaclip(1) = rsi_x_c(indext(1),2)
c$$$      theta = thetaclip(0) + dt * ( thetaclip(1) - thetaclip(0) )
c$$$      indext(0)   = ridx_it(indext(0))
c$$$      indext(1)   = ridx_it(indext(1))
c$$$c
c$$$#ifndef LINE_FREQDEP_DUST
c$$$      ii     = 1
c$$$#else
c$$$      ii     = line_ilifreq(inu)
c$$$#endif
c$$$#ifndef LINE_LOCAL_WIDTH
c$$$      write(*,*) 'ERROR: The compilation with global line width'
c$$$      write(*,*) '       (i.e. NOT setting LINE_LOCAL_WIDTH)'
c$$$      write(*,*) '       is obsolete and disactivated.'
c$$$      write(*,*) '       >>>> ADDED KEES/KLAUS MARCH 2007 <<<<'
c$$$      stop 10245
c$$$#endif
c$$$c
c$$$      if(icr.eq.1) then
c$$$c
c$$$c         The point 'P' lies on an R=const gridline
c$$$c
c$$$          src(0)   = line_dust_src(ii,iline,indext(0),indexr(0))
c$$$          src(1)   = line_dust_src(ii,iline,indext(1),indexr(0))
c$$$          src_dust = (1.d0-dt)*src(0) + dt*src(1)
c$$$          alp(0)   = line_dust_alp(ii,iline,indext(0),indexr(0))
c$$$          alp(1)   = line_dust_alp(ii,iline,indext(1),indexr(0))
c$$$          alp_dust = (1.d0-dt)*alp(0) + dt*alp(1)
c$$$c
c$$$          lw(0)    = locprof_linewidth(iline,indext(0),indexr(0))
c$$$          lw(1)    = locprof_linewidth(iline,indext(1),indexr(0))
c$$$          linewidth= (1.d0-dt)*lw(0) + dt*lw(1)
c$$$c
c$$$          N(0)     = line_level_popul(lev_up(iline),
c$$$     %                              indext(0),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(0))
c$$$     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
c$$$          N(1)     = line_level_popul(lev_up(iline),
c$$$     %                                indext(1),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(1),indexr(0))
c$$$     %                   * medium_arr_rho(indext(1),indexr(0))*molpg
c$$$          Nup      = (1.d0-dt)*N(0) + dt*N(1)
c$$$          N(0)     = line_level_popul(lev_down(iline),
c$$$     %                                indext(0),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(0))
c$$$     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
c$$$          N(1)     = line_level_popul(lev_down(iline),
c$$$     %                                indext(1),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(1),indexr(0))
c$$$     %                   * medium_arr_rho(indext(1),indexr(0))*molpg
c$$$          Ndown    = (1.d0-dt)*N(0) + dt*N(1)
c$$$c
c$$$          do i=1,3
c$$$              vel(0)      = medium_arr_velocity(i,indext(0),indexr(0))
c$$$              vel(1)      = medium_arr_velocity(i,indext(1),indexr(0))
c$$$              velocity(i) = (1.d0-dt)*vel(0) + dt*vel(1)
c$$$          enddo
c$$$c
c$$$      elseif(icr.eq.2) then
c$$$c
c$$$c         The point 'P' lies on an Theta=const gridline
c$$$c
c$$$          src(0)   = line_dust_src(ii,iline,indext(0),indexr(0))
c$$$          src(1)   = line_dust_src(ii,iline,indext(0),indexr(1))
c$$$          src_dust = (1.d0-dr)*src(0) + dr*src(1)
c$$$          alp(0)   = line_dust_alp(ii,iline,indext(0),indexr(0))
c$$$          alp(1)   = line_dust_alp(ii,iline,indext(0),indexr(1))
c$$$          alp_dust = (1.d0-dr)*alp(0) + dr*alp(1)
c$$$c
c$$$          lw(0)    = locprof_linewidth(iline,indext(0),indexr(0))
c$$$          lw(1)    = locprof_linewidth(iline,indext(0),indexr(1))
c$$$          linewidth= (1.d0-dr)*lw(0) + dr*lw(1)
c$$$c
c$$$          N(0)     = line_level_popul(lev_up(iline),
c$$$     %                                indext(0),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(0))
c$$$     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
c$$$          N(1)     = line_level_popul(lev_up(iline),
c$$$     %                                indext(0),indexr(1))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(1))
c$$$     %                   * medium_arr_rho(indext(0),indexr(1))*molpg
c$$$          Nup      = (1.d0-dr)*N(0) + dr*N(1)
c$$$          N(0)     = line_level_popul(lev_down(iline),
c$$$     %                                indext(0),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(0))
c$$$     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
c$$$          N(1)     = line_level_popul(lev_down(iline),
c$$$     %                                indext(0),indexr(1))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(1))
c$$$     %                   * medium_arr_rho(indext(0),indexr(1))*molpg
c$$$          Ndown    = (1.d0-dr)*N(0) + dr*N(1)
c$$$c
c$$$          do i=1,3
c$$$              vel(0)      = medium_arr_velocity(i,indext(0),indexr(0))
c$$$              vel(1)      = medium_arr_velocity(i,indext(0),indexr(1))
c$$$              velocity(i) = (1.d0-dr)*vel(0) + dr*vel(1)
c$$$          enddo
c$$$c
c$$$      else
c$$$c
c$$$c         The point 'P' lies somewhere in between all gridlines
c$$$c
c$$$          src(0)   = line_dust_src(ii,iline,indext(0),indexr(0))
c$$$          src(1)   = line_dust_src(ii,iline,indext(1),indexr(0))
c$$$          src(2)   = line_dust_src(ii,iline,indext(0),indexr(1))
c$$$          src(3)   = line_dust_src(ii,iline,indext(1),indexr(1))
c$$$          src_dust = (1.d0-dr) * ((1.d0-dt)*src(0) + dt*src(1)) 
c$$$     %              +       dr * ((1.d0-dt)*src(2) + dt*src(3)) 
c$$$          alp(0)   = line_dust_alp(ii,iline,indext(0),indexr(0))
c$$$          alp(1)   = line_dust_alp(ii,iline,indext(1),indexr(0))
c$$$          alp(2)   = line_dust_alp(ii,iline,indext(0),indexr(1))
c$$$          alp(3)   = line_dust_alp(ii,iline,indext(1),indexr(1))
c$$$          alp_dust = (1.d0-dr) * ((1.d0-dt)*alp(0) + dt*alp(1)) 
c$$$     %              +       dr * ((1.d0-dt)*alp(2) + dt*alp(3)) 
c$$$c
c$$$          lw(0)    = locprof_linewidth(iline,indext(0),indexr(0))
c$$$          lw(1)    = locprof_linewidth(iline,indext(1),indexr(0))
c$$$          lw(2)    = locprof_linewidth(iline,indext(0),indexr(1))
c$$$          lw(3)    = locprof_linewidth(iline,indext(1),indexr(1))
c$$$          linewidth= (1.d0-dr) * ((1.d0-dt)*lw(0) + dt*lw(1)) 
c$$$     %              +       dr * ((1.d0-dt)*lw(2) + dt*lw(3)) 
c$$$c
c$$$          N(0)     = line_level_popul(lev_up(iline),
c$$$     %                                indext(0),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(0))
c$$$     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
c$$$          N(1)     = line_level_popul(lev_up(iline),
c$$$     %                                indext(1),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(1),indexr(0))
c$$$     %                   * medium_arr_rho(indext(1),indexr(0))*molpg
c$$$          N(2)     = line_level_popul(lev_up(iline),
c$$$     %                                indext(0),indexr(1))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(1))
c$$$     %                   * medium_arr_rho(indext(0),indexr(1))*molpg
c$$$          N(3)     = line_level_popul(lev_up(iline),
c$$$     %                                indext(1),indexr(1))
c$$$     %                   * locabun_abund_mol(indext(1),indexr(1))
c$$$     %                   * medium_arr_rho(indext(1),indexr(1))*molpg
c$$$          Nup      = (1.d0-dr) * ((1.d0-dt)*N(0) + dt*N(1)) 
c$$$     %              +       dr * ((1.d0-dt)*N(2) + dt*N(3)) 
c$$$          N(0)     = line_level_popul(lev_down(iline),
c$$$     %                                indext(0),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(0))
c$$$     %                   * medium_arr_rho(indext(0),indexr(0))*molpg
c$$$          N(1)     = line_level_popul(lev_down(iline),
c$$$     %                                indext(1),indexr(0))
c$$$     %                   * locabun_abund_mol(indext(1),indexr(0))
c$$$     %                   * medium_arr_rho(indext(1),indexr(0))*molpg
c$$$          N(2)     = line_level_popul(lev_down(iline),
c$$$     %                                indext(0),indexr(1))
c$$$     %                   * locabun_abund_mol(indext(0),indexr(1))
c$$$     %                   * medium_arr_rho(indext(0),indexr(1))*molpg
c$$$          N(3)     = line_level_popul(lev_down(iline),
c$$$     %                                indext(1),indexr(1))
c$$$     %                   * locabun_abund_mol(indext(1),indexr(1))
c$$$     %                   * medium_arr_rho(indext(1),indexr(1))*molpg
c$$$          Ndown    = (1.d0-dr) * ((1.d0-dt)*N(0) + dt*N(1)) 
c$$$     %              +       dr * ((1.d0-dt)*N(2) + dt*N(3)) 
c$$$c
c$$$          do i=1,3
c$$$              vel(0)      = medium_arr_velocity(i,indext(0),indexr(0))
c$$$              vel(1)      = medium_arr_velocity(i,indext(1),indexr(0))
c$$$              vel(2)      = medium_arr_velocity(i,indext(0),indexr(1))
c$$$              vel(3)      = medium_arr_velocity(i,indext(1),indexr(1))
c$$$              velocity(i) = (1.d0-dr) * ((1.d0-dt)*vel(0)+dt*vel(1)) 
c$$$     %              +              dr * ((1.d0-dt)*vel(2)+dt*vel(3)) 
c$$$          enddo
c$$$c
c$$$      endif
c$$$c
c$$$c     First compute the inner product between Omega and v
c$$$c
c$$$      dvmu = omega_dot_v(mu,phi,velocity)
c$$$c
c$$$c     Return src_dust, alp_dust and dvmu
c$$$c
c$$$      return
c$$$      end





c     --------------------------------------------------------------
c                    BASIC LINE ELEMENT INTEGRATOR 
c
c     This routine does first order integration without sub-grid
c     model. 
c
c     ARGUMENTS: 
c       intensity  = Input AND output intensity for the formal RT 
c                    integration along this small ray segment
c       ds         = Path length of the ray segment
c       srcd0/1    = The dust emissitivity j_nu at points 0/1
c       alpd0/1    = The dust extinction alpha_nu at points 0/1
c       lw0/1      = The FWHM line width (in km/s) at points 0/1
c       dvmu0/1    = The Omega*v/c inner product at pionts 0/1
c                    where Omega is the direction of the radiation
c       nup0/1     = The upper level population [1/cm^3] at points 0/1
c       ndown0/1   = The lower level population [1/cm^3] at points 0/1
c       inu        = The frequency gridpoint index where we are now
c                    integrating the RT equation
c       iline      = The index of the line we are now integrating
c       char_tau   = (For analysis only) The total optical depth so far
c       char_emis  = (For analysis only) The total emitted radiation so far
c       init       = 1 --> Locally compute srcl1 AND srcl0 as well 
c                          as alpl1/0 etc. This is the setting you 
c                          need for the first ray element of the ray.
c                    0 --> Assume srcl0, alpl0, phiprof0 to be given
c                          (see arguments below). This is useful to
c                          speed up the code for all ray elements 
c                          except the first one, since srcl1 of the
c                          previous ray element becomes srcl0 of the
c                          next one.
c                    NOTE: This routine will automatically set init=0
c                          just before returning. So the resetting is 
c                          done automatic. You only must put init=1
c                          at the very start of the big ray.
c       phiprof0   = The phi line profile at point 0
c                    But this becomes the phiprof1 as a return-value
c       srcl0      = The j_nu_line at point 0
c                    But this becomes the j_nu_line at point 1 
c                    as a return-value
c       alpl0      = The alpha_nu_line at point 0
c                    But this becomes the alpha_nu_line at point 1 
c                    as a return-value
c
c     NOTES:
c       line_nu0   = Rest frame line center [Hz]
c       line_dnu   = Frequency we are looking at MINUS line_nu0 [Hz]
c       dnu0/1     = Comoving-frame frequency minus line_nu0
c     --------------------------------------------------------------
      subroutine integrate_element_linedust(intensity,ds,
     %             srcd0,srcd1,alpd0,alpd1,
     %             lw0,lw1,dvmu0,dvmu1,
     %             nup0,nup1,ndown0,ndown1,
     %             char_tau,char_emis,
     %             inu,iline,init,
     %             phiprof0,srcl0,alpl0,iray,
     %             r,theta,tr_s_is)
      implicit none
c
#include "common_lines.h"
c
      doubleprecision srcd0,srcd1,alpd0,alpd1,dvmu0,dvmu1,lw0,lw1
      doubleprecision nup0,nup1,ndown0,ndown1,intensity,ds
      doubleprecision char_tau,char_emis
      doubleprecision r,theta,tr_s_is      !r,theta included for debugging
      integer inu,iline,ilfr,init,iray,is  !iray,is included for debugging
c
      doubleprecision lwav,dnu0,dnu1,phiprof0,phiprof1
      doubleprecision voigt_profile,srcl0,srcl1,alpl0,alpl1
      doubleprecision src0,src1,alp0,alp1
      doubleprecision qdr_src_2
c
c     Make average of line width
c
      lwav  = 0.5 * ( lw0 + lw1 )
c
c     Get the frequency bin in the line passband from the 
c     global frequency bin index inu. In the current compilation
c     setup we will then get ilfr=inu
c
      ilfr  = line_ilifreq(inu)
c
c     Compute the line src and line alp at 0 and 1
c
c     ...First for point 0
c
c        Compute the Doppler shifted dnu
c
      dnu0 = line_dnu(ilfr,iline) - line_nu0(iline) * dvmu0
      dnu1 = line_dnu(ilfr,iline) - line_nu0(iline) * dvmu1
c
c     Local, compute on-the-fly
c
      if(init.eq.1) then
          phiprof0 = voigt_profile(linefreq(iline),0.d0,lwav,dnu0)
      endif
      phiprof1 = voigt_profile(linefreq(iline),0.d0,lwav,dnu1)
c
c     The j_nu is now:
c
c                   h nu_0
c            j_nu = ------ n_up A_ud phi(Omega,nu)
c                    4 pi
c
      if(init.eq.1) then
          srcl0 = 5.27296241956d-28 * linefreq(iline) * Nup0 * 
     %              Aud(iline) * phiprof0
      endif
      srcl1 = 5.27296241956d-28 * linefreq(iline) * Nup1 * 
     %          Aud(iline) * phiprof1
c     
c     The alpha_nu is now:
c
c                       h nu_0
c            alpha_nu = ------ ( n_down B_du - n_up B_ud ) phi(Omega,nu)
c                        4 pi
c
      if(init.eq.1) then
          alpl0 = 5.27296241956d-28 * linefreq(iline) * phiprof0 *
     %          ( Ndown0 * Bdu(iline) - Nup0 * Bud(iline) )
      endif
      alpl1 = 5.27296241956d-28 * linefreq(iline) * phiprof1 *
     %          ( Ndown1 * Bdu(iline) - Nup1 * Bud(iline) )
c
c     Now call the integrator
c
      src0  = srcd0 + srcl0
      src1  = srcd1 + srcl1
      alp0  = alpd0 + alpl0
      alp1  = alpd1 + alpl1

      intensity = qdr_src_2(intensity,src0,alp0,src1,alp1,ds)
      char_tau  = char_tau + 0.5d0 * ( alp0 + alp1 ) * ds
      char_emis = char_emis + 0.5d0 * ( src0 + src1 ) * ds
c
c     Very occasionally you may want to figure out how the 
c     integration along the characteristics has proceeded
c     For this, we can dump the intermediate result to a file
c
c#     ifdef MONITOR_CHARINT
c         write(4,'(11(E11.4,1X))') intensity,src0,src1,alp0,alp1,
c     %        char_tau,char_emis,r,theta,tr_s_is,ds
c#     endif

c
c     Swap srcl1->srcl0 etc
c
      phiprof0 = phiprof1
      srcl0    = srcl1
      alpl0    = alpl1
c
c     Put init to 0
c
      init = 0
c
c     Done...
c

      end


c     --------------------------------------------------------------
c                    CLEVER LINE ELEMENT INTEGRATOR 
c     
c     This routine is like integrate_element_linedust() (in fact,
c     it calls that routine) but it makes sure to do subgridding
c     in case the line passes by too quickly.
c
c     ARGUMENTS: See integrate_element_linedust() routine
c     --------------------------------------------------------------
      subroutine clever_integrate_element_linedust(intensity,ds,
     %             srcd0,srcd1,alpd0,alpd1,
     %             lw0,lw1,dvmu0,dvmu1,
     %             nup0,nup1,ndown0,ndown1,
     %             char_tau,char_emis,
     %             inu,iline,init,
     %             phiprof0,srcl0,alpl0,iray,
     %             r,theta,tr_s_is)
      implicit none
c
#include "common_lines.h"
c
      doubleprecision srcd0,srcd1,alpd0,alpd1,dvmu0,dvmu1,lw0,lw1
      doubleprecision srcl0,srcl1,alpl0,alpl1
      doubleprecision r,theta,tr_s_is
      doubleprecision nup0,nup1,ndown0,ndown1,intensity,ds
      doubleprecision char_tau,char_emis
      doubleprecision phiprof0
      integer inu,iline,init,i,iray   ! iray,is included for debugging
c
      doubleprecision crit_linecenter,crit_npt_fwhm
      parameter(crit_npt_fwhm=5.)     ! Nr of grid pts per FWHM of the line
      parameter(crit_linecenter=3.0)  ! Edge of refined grid in units of
c                                     !   FWHM from line center
      integer lg_nrmax,lg_nr
      parameter(lg_nrmax=2*crit_linecenter*crit_npt_fwhm+1)
      doubleprecision lg_s(lg_nrmax+2),lg_int(lg_nrmax+2)
      doubleprecision lg_dvmu(lg_nrmax+2),lg_phi(lg_nrmax+2)
      doubleprecision lg_nup(lg_nrmax+2),lg_ndown(lg_nrmax+2)
      doubleprecision lg_srcd(lg_nrmax+2),lg_alpd(lg_nrmax+2)
      doubleprecision lg_ds
      doubleprecision s,s_c,dls,lw
      doubleprecision eps,epsp
      doubleprecision ds_over_deltal_s
      doubleprecision sright, sleft
      integer ilfr
c
c     Get the frequency bin in the line passband from the 
c     global frequency bin index inu. In the current compilation
c     setup we will then get ilfr=inu
c
      ilfr  = line_ilifreq(inu)
#ifdef SUBGRID
c     If subgridding is turned on: 
c
c     Check whether we need subgridding or not
c
c     FIRST CONDITION: The `line width' in space (due to a gradient
c            in the doppler shift) must be smaller than a few times
c            the length of the integration path element.
c 
c     ...Calculate ds / Delta_l s (Delta_l s is the width of 
c        the line as seen along the integration path; ds is the
c        length of the integration element; note that this is 
c        infinite when the velocity does not change over the 
c        integration path)
c
c           ds       | U_1 - U_0 |  lambda_rest    | U_1 - U_0 | nu_rest
c        --------- = -------------------------- =  ---------------------
c        Delta_l s        Delta_l lambda                 Delta_l nu
c
c     where U is here now (v.Omega)/c. NOTE that Delta_l lambda is the
c     intrinsic line widths, and Delta_l nu is the same, but in Hz.
c     Therefore we can say that Delta_l lambda / lambda_rest = 
c     Delta_l nu / nu_rest (without minus sign, because Delta_l lambda
c     as well as Delta_l nu are positive definite).
c
c      ds_over_deltal_s = abs((dvmu1-dvmu0) * line_nu0(iline) /
c     %                      line_dnu(ilfr,iline))

      lw               = 0.5d0 * ( lw0 + lw1 ) ! NOTE: Approximation

      ds_over_deltal_s = abs((dvmu1-dvmu0) /
     %                      (lw/2.99792458d5))
c
c     Now check if we MIGHT need the subgridding
c
c      write(*,*) 2.d0*crit_linecenter*ds_over_deltal_s,
c     %     dvmu1-dvmu0,lw0/2.9979d5
      if(2.d0*crit_linecenter*ds_over_deltal_s.gt.1.d0) then
c
c          OK, we might need subgridding, UNLESS the line is not even
c          interacting with our wavelength point. So check for this.
c
c         ...Calculate s_c (which is the location along the integration
c            path where the line center crosses the wavelength of the ray)
c
c                (s_1-s_0)*(1-nu_rest/nu) + s0*U_1-s1*U_0 
c          s_c = ------------------------------------------------ 
c                                U_1 - U_0
c
c         where U is here now (v.Omega)/c, and we choose s_0=0 and s_1=ds. 
c         We have:
c
c          (1-nu_rest/nu) = (nu-nu_rest)/nu \simeq (nu-nu_rest)/nu_rest
c                         = line_dnu/line_nu0
c
          s_c = ds * ( (line_dnu(ilfr,iline)/line_nu0(iline) ) - 
     %                dvmu0) / ( dvmu1 - dvmu0 )
c
c         ...Also calculate the Delta_l s
c
          dls = ds / ds_over_deltal_s
c
c         Compute the left- and right- boundaries of the to-be-refined
c         region.
c
          sright = s_c + crit_linecenter*dls
          sleft  = s_c - crit_linecenter*dls
          if((sright.gt.0.d0).and.(sleft.lt.ds)) then

c
c             WE NEED SUBGRIDDING!
c
              lg_ds       = (sright-sleft)/(lg_nrmax-1.d0)
c
c
c             First set up the grid and the local quantities on the
c             grid
c
c             First set up the start of the (big) ray element as 
c             the first point of the integration
c             
              lg_nr       = 1
              lg_s(1)     = 0.d0  
              lg_nup(1)   = Nup0
              lg_ndown(1) = Ndown0
              lg_dvmu(1)  = dvmu0
              lg_srcd(1)  = srcd0
              lg_alpd(1)  = alpd0
c
c             Then add the refined grid points
c
              do i=1,lg_nrmax
c                  s  = sright + (i-1)*lg_ds
                  s = sleft + (i-1)*lg_ds
                  if((s.gt.0.d0).and.(s.lt.ds)) then
c
c                     This refined grid point lies inside of the
c                     current (big) ray element, so we must include
c                     it
c
                      eps             = s/ds
                      epsp            = 1.d0-eps
                      lg_nr           = lg_nr + 1
                      lg_s(lg_nr)     = s
                      lg_nup(lg_nr)   = epsp*Nup0+eps*Nup1
                      lg_ndown(lg_nr) = epsp*Ndown0+eps*Ndown1
                      lg_dvmu(lg_nr)  = epsp*dvmu0+eps*dvmu1
                      lg_srcd(lg_nr)  = epsp*srcd0+eps*srcd1
                      lg_alpd(lg_nr)  = epsp*alpd0+eps*alpd1
                  endif
              enddo
c
c             Finally add the final point, which is the end of the 
c             (big) ray element
c
              lg_nr           = lg_nr + 1
              lg_s(lg_nr)     = ds
              lg_nup(lg_nr)   = Nup1
              lg_ndown(lg_nr) = Ndown1
              lg_dvmu(lg_nr)  = dvmu1
              lg_srcd(lg_nr)  = srcd1
              lg_alpd(lg_nr)  = alpd1
c
c             Now do the integration
c
              lg_phi(1)   = phiprof0    ! NOTE: For the very first
c                                       !  ray-element this might be 
c                                       !  wrong, but it is for debug
c                                       !  purposes-only anyway
              do i=2,lg_nr
                  lg_ds = lg_s(i)-lg_s(i-1)
                  call integrate_element_linedust(intensity,lg_ds,
     %                 lg_srcd(i-1),lg_srcd(i),
     %                 lg_alpd(i-1),lg_alpd(i),
     %                 lw,lw,lg_dvmu(i-1),lg_dvmu(i),
     %                 lg_nup(i-1),lg_nup(i),
     %                 lg_ndown(i-1),lg_ndown(i),
     %                 char_tau,char_emis,
     %                 inu,iline,init,
     %                 phiprof0,srcl0,alpl0,iray,
     %                 r,theta,tr_s_is)  
c                  if(inu.eq.100) write(*,*) 
c     %            lg_nr,lg_s(i),alpl0
              enddo
c              if(lg_nr.eq.lg_nrmax) write(*,*) 
c     %        '-------'
c              write(*,*) lg_s
c              write(*,*) lg_nr
c              write(*,*) '---------'
c     
c             Done...
c
c              srcl0 = 0.d0
c              alpl0 = 0.d0
              return
          endif
      endif
#endif
c
c     DO NORMAL ONE-STEP INTEGRATION
c
      call integrate_element_linedust(intensity,ds,
     %             srcd0,srcd1,alpd0,alpd1,
     %             lw0,lw1,dvmu0,dvmu1,
     %             nup0,nup1,ndown0,ndown1,
     %             char_tau,char_emis,
     %             inu,iline,init,
     %             phiprof0,srcl0,alpl0,iray,
     %             r,theta,tr_s_is)    
c
      end
