c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c      ( MODULE: IMPLEMENT THE STELLAR LIGHT AS A CENTRAL BOUNDARY)
c
c                            Leiden, Sep 1998
c                             C.P. Dullemond   
c
c     =============================================================


#define STAR_F_VERSION "Trans_05oct98"
#include "main.h"


c     =============================================================
c            THE BOUNDARY CONDITIONS ON THE INTENSITY FIELD
c     =============================================================

#     ifndef NO_INTENSITY_STORAGE

c     --------------------------------------------------------------
c                    INNER BOUNDARY : STELLAR SURFACE
c
c     If in the center of the spherical coordinate system there is
c     a star of some kind, then here it can be installed. This
c     could be just a source of radiation, or it could be a complete
c     1-D radiative transfer algorithm of the stellar atosphere with
c     irradiation from the 2-D transfer, or it could be anything else.
c     The user may install any kind of thing here. Currently this 
c     routine simply copies the given stellar spectrum into the 
c     intensity array. This can be done in the innermost radial 
c     coordinate shell (ir=1), meaning that the stellar radius is 
c     precisely the same as the inner edge of the coordinate system. 
c     The other option (by setting ir=0) is to insert these values 
c     into the inner ghost gridpoints (lying inside the innermost
c     gridline), so that they represent the stellar surface of a star
c     with smaller radius than the inner edge of the coordinate 
c     system. The transfer routines will then take care of the 
c     proper geometrical factors to include this radiation into the
c     transfer. 
c
c     ARGUMENTS:
c      ir :   The radial index of the intensity array at which the
c             central source should be inserted. This can (reasonably)
c             be either 0 (the auxiliary inner gridpoints, used for 
c             the stellar surface) or 1 (the real inner gridpoints,
c             i.e. the inner edge of the coordinate system).
c      inu:   The index of the frequency bin of the spectrum of the
c             star that should be used as the inner boundary condition.
c     --------------------------------------------------------------
      subroutine central_source(ir,inu)
      implicit none
      integer ir,inu
c
#include "common_grid.h"
#include "common_radfield.h"
c#include "common_units.h"
#include "common_boundary.h"
c
      integer it,imu,iphi,iphimin,iphimax,iphioff
      integer ithetamax,imumax,iangset
      doubleprecision int
c
      iangset = 1
c
      imumax  = nrmu(iangset)
      iphimin = 1
      iphimax = nrphi(iangset)
      iphioff = 0
      ithetamax = irsi_frsizey
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
      ithetamax=1
      iphimin=1
      iphimax=1
#endif
c
      do it=1,ithetamax
          do imu=1,imumax
              do iphi=iphimin,iphimax
                  intens(iphi+iphioff,imu,it,ir) = 
     %                radbnd_starspec(inu)
              enddo
          enddo
      enddo
c
      end


c     --------------------------------------------------------------
c                    INNER BOUNDARY : MIRROR
c
c     This inner boundary is a reflecting mirror which copies the
c     intensity of the downward moving radiation and copies it into
c     the upward moving radiation.
c     --------------------------------------------------------------
      subroutine central_mirror(inu)
      implicit none
      integer inu
c
#include "common_grid.h"
#include "common_radfield.h"
c#include "common_units.h"
#include "common_boundary.h"
c
      integer it,imu,iphi,iphimin,iphimax,iphioff
      integer ithetamax,imumax,iangset
      doubleprecision int
c
      iangset = 1
c
      imumax  = nrmu(iangset)
      iphimin = 1
      iphimax = nrphi(iangset)
      iphioff = 0
      ithetamax = irsi_frsizey
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
      ithetamax=1
      iphimin=1
      iphimax=1
#endif
c
      do it=1,ithetamax
          do imu=1,imumax
              do iphi=iphimin,iphimax
                  intens(iphi+iphioff,imu,it,1) = 
     %               intens(iphi+iphioff,-imu,it,1)
              enddo
          enddo
      enddo
c
      end


c     --------------------------------------------------------------
c                    OUTER EQUILIBRIUM CONDITION
c
c     For 1-D calculations it is handy to have an outer boundary
c     condition that enables the solution to get to equilibrium 
c     there. That's what I do here. It's simply inserting the
c     value of j/alpah (at R=R_{imax-1})into the intensity at 
c     R=R_imax. 
c
c     REMARK: Be sure that the sources are set, in other words
c             that the make_lte_src() routine has been called
c             and that the solve_statequil_and_scat() 
c             has been called as well.
c     --------------------------------------------------------------
      subroutine outer_equilibrium(inu)
      implicit none
c
      doubleprecision calc_src_pnt,calc_alp_pnt
      integer inu
c
#include "common_grid.h"
#include "common_radfield.h"
c#include "common_units.h"
#include "common_boundary.h"
c
      integer ir,it,imu,iphi,iphimin,iphimax,iphioff
      integer ithetamax,imumax,iangset
      doubleprecision int
c
      ir      = irsi_frsizex
      iangset = 1
c
      imumax  = nrmu(iangset)
      iphimin = 1
      iphimax = nrphi(iangset)
      iphioff = 0
      ithetamax = irsi_frsizey
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
      ithetamax=1
      iphimin=1
      iphimax=1
#endif
c
      do it=1,ithetamax
          do imu=-imumax,0
              do iphi=iphimin,iphimax
                  intens(iphi+iphioff,imu,it,ir) = 
     %        calc_src_pnt(inu,ir,it,rmu(imu,iangset),
     %                     rphi(iphi,iangset)) /
     %        calc_alp_pnt(inu,ir,it,rmu(imu,iangset),
     %                     rphi(iphi,iangset))
              enddo
          enddo
      enddo      
c
      end


c     --------------------------------------------------------------
c              OUTER BC: MICROWAVE BACKGROUND RADIATION
c
c     This is just a hack, so don't take this subroutine too serious.
c     I use the common_lines.h file to extract the appropriate data.
c     --------------------------------------------------------------
      subroutine outer_microwave_backgr(inu)
      implicit none
c
      integer inu
c
#include "common_grid.h"
#include "common_radfield.h"
c#include "common_units.h"
#include "common_boundary.h"
#include "common_lines.h"
#include "common_source.h"
c
      integer ir,it,imu,iphi,iphimin,iphimax,iphioff
      integer ithetamax,imumax,iangset,iline
      doubleprecision int,freq,temp
c
c     Are we doing lines?
c
#ifdef INCLUDE_LINES
      if(iradproc_line.eq.0) then
          write(*,*) 'WARNING: Microwave BG only for lines'
          return
      endif
c
      ir      = irsi_frsizex
      iangset = 1
c
      imumax  = nrmu(iangset)
      iphimin = 1
      iphimax = nrphi(iangset)
      iphioff = 0
      ithetamax = irsi_frsizey
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
      ithetamax=1
      iphimin=1
      iphimax=1
#endif
c
c     First find out in which line we are
c
      iline = line_iline(inu)
      freq  = linefreq(iline)
      temp  = TEMPCMB               ! = 2.728d0
c
      do it=1,ithetamax
          do imu=-imumax,0
              do iphi=iphimin,iphimax
                  intens(iphi+iphioff,imu,it,ir) = 
     %                  1.47455253991d-47 * (freq**3) /
     %              ( exp(4.7991598d-11*freq/temp) - 1.d0 )
              enddo
          enddo
      enddo      
#else
      write(*,*) 'Cosmic Microwave BG not active'
      stop 12
#endif
c
      end


c     --------------------------------------------------------------
c              OUTER INTERSTELLAR RADIATION FIELD CONDITION
c
c     This routine makes it possible to fix the input radiation at
c     the outer edge to a given intensity as a function of freq.
c     This is useful for calculations of prestellar cloud cores
c     and class 0 objects in which the external radiation field is
c     dominant over the central source at the outer regions of
c     the cloud.
c     --------------------------------------------------------------
      subroutine outer_interstellar_radfield(inu)
      implicit none
c
      doubleprecision calc_src_pnt,calc_alp_pnt
      integer inu
c
#include "common_grid.h"
#include "common_radfield.h"
c#include "common_units.h"
#include "common_boundary.h"
c
      integer ir,it,imu,iphi,iphimin,iphimax,iphioff
      integer ithetamax,imumax,iangset
      doubleprecision int
c
      ir      = irsi_frsizex
      iangset = 1
c
      imumax  = nrmu(iangset)
      iphimin = 1
      iphimax = nrphi(iangset)
      iphioff = 0
      ithetamax = irsi_frsizey
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
      ithetamax=1
      iphimin=1
      iphimax=1
#endif
c
      do it=1,ithetamax
          do imu=-imumax,0
              do iphi=iphimin,iphimax
                  intens(iphi+iphioff,imu,it,ir) = 
     %                 radbnd_interstellfield(inu)
              enddo
          enddo
      enddo      
c
      end


c     --------------------------------------------------------------
c                        VACUUM OUTER SPACE BC
c     --------------------------------------------------------------
      subroutine outer_vacuum(inu)
      implicit none
c
      integer inu
c
#include "common_grid.h"
#include "common_radfield.h"
c#include "common_units.h"
#include "common_boundary.h"
#include "common_source.h"
c
      integer ir,it,imu,iphi,iphimin,iphimax,iphioff
      integer ithetamax,imumax,iangset
      doubleprecision int,freq,temp
c
c     Are we doing lines?
c
      ir      = irsi_frsizex
      iangset = 1
c
      imumax  = nrmu(iangset)
      iphimin = 1
      iphimax = nrphi(iangset)
      iphioff = 0
      ithetamax = irsi_frsizey
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
      ithetamax=1
      iphimin=1
      iphimax=1
#endif
      do it=1,ithetamax
          do imu=-imumax,0
              do iphi=iphimin,iphimax
                  intens(iphi+iphioff,imu,it,ir) = 0.d0
              enddo
          enddo
      enddo      
c
      end






#endif /* ifndef NO_INTENSITY_STORAGE */





c     =============================================================
c      GENERAL ROUTINES FOR BOUNDARY CONDITIONS AND STAR SPECTRUM
c     =============================================================

c     --------------------------------------------------------------
c                  READ THE STELLAR SPECTRUM FROM FILE
c
c     This routine reads the stellar spectrum into memory. This can
c     be done in the old way (intensity at stellar surface) or in 
c     the new way (compatible with spectrum.dat output: i.e. flux
c     at 1 pc.), depending on which file exists in the current 
c     directory. The starspec.inp contains the intensity at the
c     stellar surface. The starspectrum.inp (which is the new style
c     and is to be preferred) contains the flux at 1 pc. Note that
c     the frequencies must coincide with those specified in the
c     frequency.inp file.
c     --------------------------------------------------------------
      subroutine read_starspec()
      implicit none
c
#include "common_boundary.h"
#include "common_grid.h"
c
      integer i,iformat
      doubleprecision dummy,frq
      logical fexold,fexnew
c
c     First check if there is a starinfo.inp file to read the
c     stellar radius
c
      call read_starinfo()
c
c     Check old and new style
c
      inquire(file='starspec.inp',EXIST=fexold)
      inquire(file='starspectrum.inp',EXIST=fexnew)
      if(fexold.and.fexnew) then
          write(*,*) 'PROBLEM: The file starspec.inp is the old'
          write(*,*) '         style star file. The file ',
     %         'starspectrum.inp'
          write(*,*) '         is the new style. Both files are'
          write(*,*) '         present in this directory. Please'
          write(*,*) '         make sure only one exists.'
          stop 13
      endif
c
      if(fexold) then
c
c         Do old way: I_nu(R=R*)
c         
          open(unit=3,file='starspec.inp',status='old',err=701)
          read(3,*) i
          if(i.ne.icont_freq_nr) then
              write(*,*) 'Sorry, stellar spectrum in starspec.inp'
              write(*,*) 'has not the same number of freqs as ',
     %             'other input files.'
              stop 13
          endif
          do i=1,icont_freq_nr
              read(3,*) dummy 
              radbnd_cont_starspec(i) = dummy
          enddo
          close(3)
      elseif(fexnew) then
c         
c         Do modern way: F_nu(R=1pc)
c         
          open(unit=3,file='starspectrum.inp',status='old',err=702)
          read(3,*) i
          if(i.ne.icont_freq_nr) then
              write(*,*) 'Sorry, stellar spectrum in ',
     %                   'starspectrum.inp'
              write(*,*) 'does not have the same number of ',
     %             'freqs as other input files.'
              stop 13
          endif
          if(i.eq.1) then
              read(3,*) frq,dummy
              radbnd_cont_starspec(1) = 3.0308410d36 * 
     %             dummy / radbnd_rstar**2
          else
              do i=1,icont_freq_nr
                  read(3,*) frq,dummy 
                  if(abs(frq-cont_freq_nu(i))/
     %                 (frq+cont_freq_nu(i)).gt.1.d-3) then
                      write(*,*) 'PROBLEM: Frequency grid of ',
     %                     'stellar spectrum unequal to ',
     %                     'frequency.inp'
                      write(*,*) frq,cont_freq_nu(i)
                      stop 13
                  endif
c                 NOTE: pc^2/pi = 3.0308410d36
                  radbnd_cont_starspec(i) = 3.0308410d36 * 
     %                 dummy / radbnd_rstar**2
              enddo
              close(3)
          endif
      else
          write(*,*) 'PROBLEM: Trying to read either starspec.inp'
          write(*,*) '         or starspectrum.inp. But neither of '
          write(*,*) '         these files is present in this ',
     %         'directory.'
          stop 13
      endif
c     
c     Compute total luminosity of the star
c     
      if(freq_nr.gt.1) then
          radbnd_lstar = 0.d0
          do i=2,icont_freq_nr
              radbnd_lstar = radbnd_lstar + 
     %             0.5d0 * 3.14159265359d0 *
     %             ( radbnd_cont_starspec(i)+
     %               radbnd_cont_starspec(i-1) ) *
     %             abs( cont_freq_nu(i) - cont_freq_nu(i-1) ) 
          enddo
          radbnd_lstar = radbnd_lstar * 
     %         1.25663706144d1 * radbnd_rstar**2 
      else
          radbnd_lstar = 39.4784176045 * radbnd_rstar**2 * 
     %         radbnd_cont_starspec(1)
      endif
c
c     Flag that spectrum was read
c
      iradbnd_read_starspec = 1
c
c     Final stuff
c
      goto 710
 701  continue
      write(*,*) 'Could not open file starspec.inp'
      stop 13
 702  continue
      write(*,*) 'Could not open file starspectrum.inp'
      stop 13
 710  continue
      return
      end


c     --------------------------------------------------------------
c                          READ STAR INFO
c     --------------------------------------------------------------
      subroutine read_starinfo()
      implicit none
c
#include "common_boundary.h"
#include "common_grid.h"
c
      integer iformat
      logical fexstar
c
      if(iradbnd_read_rstar.ne.0) then
          return
      endif
      inquire(file='starinfo.inp',EXIST=fexstar)
      if(fexstar) then
          if(radbnd_rstar.ne.0.d0) then
              write(*,*) '---------------------------------------------'
              write(*,*) 'WARNING: Found starinfo.inp file, while'
              write(*,*) '         Rstar is set elsewhere as well.'
              write(*,*) '         I now take the value in starinfo.inp'
              write(*,*) '---------------------------------------------'
          endif
          write(*,*) 'Reading starinfo.inp'
          open(unit=3,file='starinfo.inp',status='old')
          read(3,*) iformat
          read(3,*) radbnd_rstar
          close(3)
      else
          if(radbnd_rstar.eq.0.d0) then
              if(rsi_x_c(1,1).eq.0.d0) then
                  write(*,*) 'ERROR: While reading stellar spectrum'
                  write(*,*) '       from starspectrum.inp, I find'
                  write(*,*) '       that rstar is 0 (=unset).'
                  write(*,*) '       This cannot be correct. '
                  write(*,*) '       Warn author!'
                  stop 55
              endif
              radbnd_rstar = rsi_x_c(1,1)*1.d-3 ! Small enough
          endif
      endif
      iradbnd_read_rstar = 1
c
      end





c     --------------------------------------------------------------
c                    SET UP THE STELLAR SPECTRUM
c     --------------------------------------------------------------
      subroutine setup_starspec(iline)
      implicit none
c
#include "common_source.h"
#include "common_boundary.h"
#include "common_grid.h"
c
      integer inu,nf,iline
c
c     NEW: 28.03.06: If the frequency grid is not yet read in, then
c     do it now. This is useful in case of line transfer with pumping
c     by the star.
c
      if(icont_freq_nr.eq.0) then
          if(iline.ne.0) then
              call create_freqgrid(1d0,1d2,nf,-1,1)
          else
              stop 79623
          endif
      endif
c
c     If the inner boundary condition does not require a spectrum
c     to be read in, then put rstar to zero and forget about reading
c     spectrum
c
      if((iradbnd_in_itype.ne.1).and.(iradbnd_in_itype.ne.2)) then
          radbnd_rstar = 0.d0
          return
      endif
c
c     Else read spectrum
c
      call read_starspec()
c
c     NEW: 28.03.06: Copy, in case line overwrites dust

      if(iline.eq.0) then
          do inu=1,freq_nr 
              radbnd_starspec(inu) = radbnd_cont_starspec(inu) 
          enddo
      endif
c
      end



c     --------------------------------------------------------------
c                 READ THE INTERSTELLAR RADIATION FIELD
c     --------------------------------------------------------------
      subroutine read_interstellfield(iline)
      implicit none
c
#include "common_grid.h"
#include "common_boundary.h"
c
      integer inu,nf,iline
      logical fex
c
c     NEW: 28-03-06: If frequency.inp not yet read, then read it now.
c
      if(icont_freq_nr.eq.0) then
          if(iline.ne.0) then
              call create_freqgrid(1d0,1d2,nf,-1,1)
          else
              stop 79624
          endif
      endif
c
c     Now read the interstellar radiation field
c
      inquire(file='interstellfield.inp',exist=fex)
      if(.not.fex) then
          write(*,*) 'ERROR: You have chosen outer BC type 3.'
          write(*,*) '    But I cannot find the file '
          write(*,*) '    interstellfield.inp...'
          stop 13
      endif
      open(unit=1,file='interstellfield.inp',status='old')
      read(1,*) nf
      if(nf.ne.icont_freq_nr) then
          write(*,*) 'ERROR: interstellfield.inp inconsistent nr'
          write(*,*) '    of frequencies compared to frequency.inp' 
          stop 13
      endif
      do inu=1,nf
          read(1,*) radbnd_cont_interstellfield(inu)
c
c         NEW: 28.03.06: Copy, in case line overwrites dust
c
          if(iline.eq.0) then
              radbnd_interstellfield(inu) = 
     %              radbnd_cont_interstellfield(inu)
          endif
      enddo
      close(1)
c
      end
