c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c                    ( MODULE: DUST OPACITY PACKAGE ) 
c
c                           Leiden, March 1999
c                    C.P. Dullemond / Alex de Koter
c
c     =============================================================
c
c     The general idea of this module:
c     
c     Rather than computing the emission and opacities from the
c     medium_rho and medium_temp, here we have a collection of
c     different dust species, each consisting of grain size 
c     distribution which is modeled by sampling the distribution.
c     The density of these species must be specified in the array
c     dust_rho(), specifying the density of each of the species
c     of dust particles. Note that the density is specified for
c     each species separately, but NOT for each grainsize separately.
c     The grainsizes are integrated to yield an average grain 
c     weight, so that the average grain density can be computed,
c     which is supposed to be done by the routines in the module
c     called "dustopacities.F", which at the time of writing this
c     is still a separate program that produces an input file
c     called "dustopac.inp", which the dust.F module can read in.
c     The dust temperatures, on the other hand, _are_ specified for
c     each species AND grain size. 
c
c     This module can compute the scattering and absorbtion for 
c     setups involving varying dust species abundances over space. 
c     Also it can solve the temperatures of the dust particles.
c
c     NOTE: The solving of the temperatures is not yet compatible
c           with ALI.
c
c       Leiden, 11 april 1999
c     =============================================================

#define DUST_F_VERSION "26apr00"
#include "main.h"

#ifdef INCLUDE_DUST  /* Compile only if necessary */

c#define DUST_PRIVATE_EXTERNALS

c     =============================================================
c
c                 THE DUST DATA AND DUST OPACITY PART
c
c     This part of the dust module one can use also in other 
c     programs. It requires externally:
c
c       - A frequency array freq_nu(1..freq_nr) in which the 
c         frequencies of the dust opacity files are given. This
c         must be located in the common_grid.h file.
c       - A subroutine that reads these frequencies (if you dont
c         have such a routine, uncomment the line with the 
c         #define DUST_PRIVATE_EXTERNALS, and the routine
c         read_frequencies() will be provided to you, see below). 
c       - Define the following constants:
c           DUST_SPECIES_MAX   Maximum nr of dust species
c           DUST_SIZE_MAX      Maximum nr of grain sizes per species
c           DUST_TRANGE_MAX    Maximum nr of temperature ranges 
c                              for temperature-dependent opacities
c           INCLUDE_DUST       Just define this to activate module
c
c     Then, to use the module, do these things:
c
c       1) Read the frequency grid [call read_frequencies()]
c       2) Read the dust opacity data [call read_dustdata()]
c       3) To get the dust kappa for a certain species ispec,
c          a grain size isize, a temperature temp, call the
c          function find_dust_kappa().
c
c     =============================================================



c     --------------------------------------------------------------
c                       READ THE DUST DATA
c
c     This is the main routine for the input of all relevant dust
c     data, and the preprocessing of this data, if necessary.
c
c     REMARK: dust_done_read is a flag signalling to the program
c             that we have (fortunately) not forgotten to read in
c             the dust data. This is only for internal consistency
c             checking.
c     --------------------------------------------------------------
      subroutine read_dustdata()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
c
      integer ispec,inu,isize,iformat,idum,idustfile,idum2
      doubleprecision dummy,temp0,temp1,dtemp0,dtemp1,dum3,mp
      character*80 comstring,filename,base,ext
      logical fex_fwgt
      parameter(mp=1.6726d-24)
c
c     Check if the transfer routine has bothered at all about the
c     amount of expected dust species to model.
c
#     ifdef DUST_RADICAL      
      if(dust_setup_nrspecies.eq.0) then
          write(*,*) 'POSSIBLE PROBLEM in read_dustdata():'
          write(*,*) '  It seems that the dust_setup_nrspecies'
          write(*,*) '  has not been set. It seems to me that this'
          write(*,*) '  implies that the calc_phys_var() routine'
          write(*,*) '  has not bothered to think about how many'
          write(*,*) '  dust species should be included. For '
          write(*,*) '  safety I quit here. If you are sure that'
          write(*,*) '  everything is okay, simply set the '
          write(*,*) '  dust_setup_nrspecies to the number of'
          write(*,*) '  expected dust species, or put it to -1'
          write(*,*) '  to leave it up to the read_dustdata() routine.'
          stop 99
      endif
#     endif
c
c     Now open the master file for the dust
c
      open(unit=3,file='dustopac.inp',status='old',err=701)
      read(3,*) iformat
      read(3,*) dust_nr_species
      read(3,*) comstring      
      if(dust_nr_species.gt.DUST_SPECIES_MAX) then
          write(*,*) 'ERROR in file dustopac.inp:',
     %            ' Too many species specified.'
          write(*,*) 'Aborting...'
          stop
      endif
c
c     Now a loop over the dust species
c
      do ispec=1,dust_nr_species
          read(3,*) idum
c
c         NEWSTUFF 12-12-04
c
          if(iformat.ge.2) then
             read(3,*) idum2
c
c            NEWSTUFF 29-12-04
c            MODIFIED 11-04-05
c
             if(idum2.ne.0) then 
                 dust_quantum(ispec) = 1
                 if(idum2.eq.2) then
                     read(3,*) dum3
                     dust_mgrain(ispec) = dum3*12*mp
                 elseif(idum2.eq.3) then
                     read(3,*) dum3
                     dust_mgrain(ispec) = dum3*12*mp
                     read(3,*) dum3
                     dust_tmax(ispec) = dum3
                 else
                     dust_mgrain(ispec) = 0
                 endif
             else
                 dust_quantum(ispec)=0
             endif
          else
             dust_quantum(ispec)=0
          endif
          read(3,*) idustfile
          base='dustopac_'
          ext ='.inp'
          call make_indexed_filename(base,idustfile,ext,filename)
c
          if(idum.eq.-1) then
c
c             Simple temperature-independent opacities
c
c              write(*,*) filename
              call read_dustopac_file(ispec,filename)
          elseif(idum.eq.-2) then
c
c             Include the minimal and maximal temperature for
c             this species.
c

              read(3,*) temp0
              read(3,*) temp1
              dust_tmin(ispec)  = temp0
              dust_tmax(ispec)  = temp1
              dust_dtmin(ispec) = 0.d0
              dust_dtmax(ispec) = 0.d0
c              if((temp0.ge.0.d0).or.(temp1.gt.0.d0)) then
c                  iradproc_dust_tminmax = 1
c              endif
              call read_dustopac_file(ispec,filename)
          elseif(idum.eq.-3) then
c
c             Not only tmin/tmax, but also smooth switch-off
c

              read(3,*) temp0
              read(3,*) dtemp0
              read(3,*) temp1
              read(3,*) dtemp1
              dust_tmin(ispec)  = temp0
              dust_tmax(ispec)  = temp1
              dust_dtmin(ispec) = abs(dtemp0)
              dust_dtmax(ispec) = abs(dtemp1)
c              if((temp0.ge.0.d0).or.(temp1.gt.0.d0)) then
c                  iradproc_dust_tminmax = 1
c              endif
              call read_dustopac_file(ispec,filename)
          elseif(idum.eq.-4) then
c
c             This is a near-future mode in which the dust opacities
c             can be really temperature-dependent
c
              write(*,*) 'Sorry, dust input type ',idum,' not yet',
     %             'finished.'
              stop 55
          else
              write(*,*) 'While reading dustopac.inp: other ',
     %         'input modes than -1 not yet implemented'
              stop 13
          endif
          read(3,*) comstring                
      enddo
c
      close(3)
c
c     Check if the nr of dust species found here agrees with
c     the number that the transfer routines expect. If the
c     dust_setup_nrspecies=-1 then this is a signal of the
c     transfer routine that the number of dust species is in
c     fact determined by the dustopac.inp file itself.
c
      if((dust_setup_nrspecies.gt.0).and.
     %    (dust_setup_nrspecies.ne.dust_nr_species)) then
          write(*,*) 'INCONSISTENCY, read_dustdata():'
          write(*,*) '  The dust file dustopac.inp specifies ',
     %                 dust_nr_species,' dust species.'
          write(*,*) '  The radiative transfer setup expects ',
     %                 dust_setup_nrspecies,' dust species.'
          write(*,*) '  (=dust_setup_nrspecies).'
          write(*,*) '  So: simulation and dust files do not '
          write(*,*) '  agree on the number of dust species.'
          write(*,*) '  modify either of the two. Now I quit.'
          stop 99
      endif
      if(dust_setup_nrspecies.lt.0) then
          dust_setup_nrspecies = dust_nr_species
      endif
c
      goto 710
  701 continue
      write(*,*) 'Could not open file dustopac.inp'
      stop 13
  710 continue
c
c     Now read the dust opacity frequency weights, if they have 
c     not yet been specified (which can be seed from dust_frwgt_read)
c
      inquire(file='dustopac_fwgt.inp',EXIST=fex_fwgt)
      if(fex_fwgt.and.(dust_frwgt_read.eq.0)) then
          open(unit=1,file='dustopac_fwgt.inp',status='old',err=701)
          read(1,*) idum
          if(idum.ne.freq_nr) then
              write(*,*) 'ERROR: file dustopac_fwgt.inp ',
     %            'is incompatible with frequency.dat'
              stop 13
          endif
          do inu=1,idum
              read(1,*) dummy
              dust_freq_wgt(inu) = dummy
          enddo
          close(1)
          dust_frwgt_read = 1
      endif
c
c     If the frequency weights have still not been specified (for 
c     whatever reason), then we generate them here
c
#     ifndef EXCLUDE_DUST_ITERATION
      if(dust_frwgt_read.eq.0) then
          call make_dust_freq_weights()
          dust_frwgt_read = 1
      endif
#     endif
c
      dust_done_read = 1
      return
      end



c     --------------------------------------------------------------
c                     READ THE DUST OPACITY FILES 
c
c     This routines read the ready-to-use data from the file
c     dustopac.inp. This file contains the dust opacities and
c     weights for the frequencies which are used in this simulation.
c     There are two possible ways these tables can be given. One
c     way (the old way) is to give a table that is INdependent off
c     the temperature. The other is to specify the opacities for
c     a series of temperature ranges (in each range the opacity 
c     remains independent off temperature, but may differ from range
c     to range). The latter option is automatically selected by
c     starting the file with the number -1. See below how the format
c     is then.
c
c     NOTE:   The dust opacities are assumed to have been computed 
c             elsewhere and written readily into this file. 
c
c     REMARK: This routine will be called by read_dustdata(), which
c             will keep unit=3 open.
c     --------------------------------------------------------------
      subroutine read_dustopac_file(ispec,filename)
      implicit none
      character*80 filename
      integer ispec
c
#include "common_dust.h"
#include "common_grid.h"
c
      integer ifr,isize,nsize,itemp,ntemp,idum
      doubleprecision dummy,a,b,bold
      logical readtrange
c
      ntemp = 1
      readtrange = .false.
c
#     ifdef DUST_RADICAL
      if(freq_grid_type.ne.-1) then
          write(*,*) 'ERROR: In radlite.inp the frequency grid type'
          write(*,*) '       is specified as ',freq_grid_type
          write(*,*) '       But when dust opacity files are used,'
          write(*,*) '       the frequency grid MUST be read from '
          write(*,*) '       the file "frequency.inp", since the '
          write(*,*) '       opacity files dustopac_*.inp contain '
          write(*,*) '       the opacities at those frequencies. '
          write(*,*) '          In the future perhaps we will allow'
          write(*,*) '       other frequency grids, and use some kind'
          write(*,*) '       of interpolation to map the opactities '
          write(*,*) '       from the frequency.inp grid onto another '
          write(*,*) '       frequency grid. But for now: specify '
          write(*,*) '       frequency grid type -1 in radlite.inp'
          stop 13
      endif
#     endif
c
      open(unit=1,file=filename,status='old',err=701)
      read(1,*) ifr,nsize
      if(ifr.eq.-1) then
c
c         ifr=-1 is a signal to say that the dust opacity file consists
c         of opacities at a series of temperature points. 
c
          readtrange = .true.
          ifr=nsize
          read(1,*) nsize,idum,ntemp
c
c         Check also if the DUST_OPAC_TEMPDEP is set.
c
#         ifndef DUST_OPAC_TEMPDEP
          write(*,*) 'PROBLEM: The code is compiled without the'
          write(*,*) '         DUST_OPAC_TEMPDEP option. But the '
          write(*,*) '         dust opacity file',filename,
     %          ' in fact has'
          write(*,*) '         multiple temperature ranges.'
          write(*,*) '         Recompile with DUST_OPAC_TEMPDEP in'
          write(*,*) '         the configure.h file!'
          stop 13
#         endif
      else
c
c         Else, there is only one (fixed) opacity table for the dust,
c         independent on temperature (with the exception that there is
c         a Tmin and Tmax within which this species can survive). 
c
          ntemp = 1
      endif
      dust_nr_size(ispec) = nsize
      dust_nr_temp(ispec) = ntemp
      if(ifr.ne.freq_nr) then
          write(*,*) 'Number of frequencies in ',filename
          write(*,*) 'not equal to the number in frequency.inp'
          write(*,*) ifr,freq_nr
          stop 13
      endif
      if(nsize.gt.DUST_SIZE_MAX) then
          write(*,*) 'Number of dust sizes for species ',ispec
          write(*,*) 'exceeds maximum as specified in Makefile'
          stop 13
      endif
      if(ispec.gt.DUST_SPECIES_MAX) then
          write(*,*) 'Dust species identifier is larger than '
          write(*,*) 'maximum as specified in Makefile'
          stop 13
      endif
      if(ntemp.gt.DUST_TRANGE_MAX) then
          write(*,*) 'Sorry, the nr of dust opacity temperature'
          write(*,*) 'points is larger than DUST_TRANGE_MAX'
          stop 13
      endif
c
c     Now read the ready-to-use kappa*wgt values. These values
c     include already the relative abundances of the various
c     dust sizes.
c
c     NOTE: the temperature ranges should (for the moment) be 
c           exactly matching: the upper temp of the previous range
c           should equal the lower temp of the next. In the future
c           this will no longer be necessary. Then there may be a
c           gap in between: linear interpolation will then be used
c           to bridge the gap. By taking the lower and upper temp
c           of the same temp range to be equal, one naturally 
c           reproduces the other interesting possibility: gradually
c           varying opacities as a function of temperature, with
c           linear interpolation. 
c
c     Make sure that, if the dust opacity table is simply a constant
c     opacity (as a function of temperature), that the first temperature
c     range is from 0 to infinity (=1d33). 
c
      dust_temprange_low(1,ispec)      = 0.d0
      dust_temprange_high(1,ispec)     = 1.d33
c
c     Now read tables in
c
      bold = 0.d0
      do itemp=1,ntemp
          if(readtrange) then
              read(1,*) a,b
              dust_temprange_low(itemp,ispec)  = a
              dust_temprange_high(itemp,ispec) = b
              dust_opacity_tempdep = 1
              if(a.lt.bold) then
                  write(*,*) 'ERROR while reading ',filename
                  write(*,*) '    for dust species ',ispec
                  write(*,*) '    The minimum temperature for '
                  write(*,*) '    opacity temp range ',itemp,' is '
                  write(*,*) '    smaller than the maximum for '
                  write(*,*) '    opacity temp range ',itemp-1
                  stop 13
              endif
              if(abs((a-bold)/(a+bold)).gt.1.d-3) then
                  write(*,*) 'ERROR in dust.F'
                  write(*,*) '    Sorry: I can not yet handle'
                  write(*,*) '    temperature-dependent dust opacity'
                  write(*,*) '    tables with gaps in temperature.'
                  write(*,*) '    So Tmax of one temp-range must equal'
                  write(*,*) '    Tmin of the next temp-range, in the'
                  write(*,*) '    file ',filename
                  stop 13
              endif
              if(a.lt.0.d0) then
                  write(*,*) 'PROBLEM: while reading ',filename
                  write(*,*) '    found negative minimum temperature'
                  stop 13
              endif
              if(b.lt.0.d0) then
                  write(*,*) 'PROBLEM: while reading ',filename
                  write(*,*) '    found negative maximum temperature'
                  stop 13
              endif
              if(b.lt.a) then
                  write(*,*) 'PROBLEM: while reading ',filename
                  write(*,*) '    minimum temp is larger than max'
                  stop 13
              endif
              bold = b
          else
              if(itemp.gt.1) then
                  write(*,*) 'BUG: Inconsistency in dust.F/',
     %                            'read_dustopac_file()'
                  stop 13
              endif
          endif
          do ifr=1,freq_nr
              do isize=1,nsize
                  read(1,*) dummy
                  dust_kappawgt_abs(ifr,itemp,isize,ispec) = dummy
              enddo
          enddo
          do ifr=1,freq_nr
              do isize=1,nsize
                  read(1,*) dummy
                  dust_kappawgt_scat(ifr,itemp,isize,ispec) = dummy
              enddo
          enddo
      enddo
c
c     Now, if there was a minimum or maximum temperature of the
c     dust selected, put this into the t-range database
c
      if(dust_tmax(ispec).gt.1.d-2) then
          dust_temprange_high(ntemp,ispec) = min(dust_tmax(ispec),
     %                          dust_temprange_high(ntemp,ispec))
          dust_opacity_tempdep = 1
      endif
      if(dust_tmin(ispec).gt.0.d0) then
c       >>> bugfix 06-08-00: .ge. --> .gt. <<<
          dust_temprange_low(1,ispec) = max(dust_tmin(ispec),
     %                          dust_temprange_low(1,ispec))
          dust_opacity_tempdep = 1
      endif
c
c     Done...
c     
      close(1)
c
      goto 710
  701 continue
      write(*,*) 'Could not open file ',filename
      stop 13
  710 continue
      return
      end


      
c     --------------------------------------------------------------
c                         FIND DUST OPACITY
c     
c     This function returns the dust opacity of a given species,
c     size, and at a given frequency-index. It takes consistently
c     into account the temperature range in which the dust grains
c     can exist.
c
c     --------------------------------------------------------------
      function find_dust_kappa(inu,isize,ispec,temp,iabs,iscat)
      implicit none
      doubleprecision temp,find_dust_kappa
      integer inu,isize,ispec,iabs,iscat
c
#include "common_grid.h"
#include "common_dust.h"
c
      doubleprecision tmin,tmax,dtmin,dtmax,plindex
      doubleprecision condense,fact,omfact
      integer ntemp,itlo,ithi
c
      condense = 1.d0
      ntemp = dust_nr_temp(ispec)
      tmin  = dust_temprange_low(1,ispec)
      tmax  = dust_temprange_high(ntemp,ispec)
      dtmin = dust_dtmin(ispec)
      dtmax = dust_dtmax(ispec)
      itlo  = 1
      ithi  = 1
      fact  = 1.d0
c
c     Now take care of temperature regimes
c
      if(temp.lt.tmin) then
          if(dtmin.eq.0.d0) then
c             
c             Hard switch-off at low temperature
c         
              condense = 0.d0
          else
c         
c             Smooth switch-off at low temperature
c         
              plindex  = log(1.d0-(dtmin/tmin))
              condense = (tmin/temp)**plindex
              if(condense.gt.1.d0) then
                  stop 9932
              endif
          endif
      elseif(temp.ge.tmax) then
          if((dtmax.eq.0.d0).and.(tmax.gt.0.d0)) then
c         
c             Hard switch-off at high temperature
c         
              condense = 0.d0
          else
c         
c             Smooth switch-off at high temperature
c         
              plindex  = log(1.d0+(dtmax/tmax))
              condense = (tmin/temp)**plindex
              if(condense.gt.1.d0) then
                  stop 9932
              endif              
          endif
      else
          if(ntemp.gt.1) then
              call hunt(dust_temprange_low(1,ispec),ntemp,temp,itlo)
              if(temp.gt.dust_temprange_high(itlo,ispec)) then
c
c                 Consistency check
c
                  if(itlo.ge.ntemp) then
                      write(*,*) 'BUG: Inconsistency in dust opacity'
                      write(*,*) '     temperature range....'
                      stop 13
                  endif
c
c                 We are in the linear interpolation regime
c
                  ithi = itlo + 1
                  fact = ( temp - dust_temprange_high(itlo,ispec) ) 
     %                   / ( dust_temprange_low(ithi,ispec) -
     %                       dust_temprange_high(itlo,ispec) )
              else
c
c                 We are within a range of constant opacity
c
                  ithi = itlo
                  fact = 1.d0
              endif
          endif
      endif
      omfact = 1.d0 - fact
c
c     Now get the opacities
c
      find_dust_kappa = 0.d0
      if(iabs.ne.0) then
          find_dust_kappa = find_dust_kappa +
     %           omfact * dust_kappawgt_abs(inu,itlo,isize,ispec) +
     %           fact * dust_kappawgt_abs(inu,ithi,isize,ispec)
      endif
      if(iscat.ne.0) then
          find_dust_kappa = find_dust_kappa +
     %           omfact * dust_kappawgt_scat(inu,itlo,isize,ispec) +
     %           fact * dust_kappawgt_scat(inu,ithi,isize,ispec)
      endif
c
c     If either iscat or iabs .lt.0 then do not incorporate the
c     evaporation of the grain. This can be handy to ensure that
c     at every temperature (even for the unphysically large ones)
c     a dust temperature can be computed. Simply for consistency.
c
      if((iabs.ge.0).and.(iscat.ge.0)) then
          find_dust_kappa = find_dust_kappa * condense
      endif
c
      return
      end




c     --------------------------------------------------------------
c                COMPUTE THE DUST FREQUENCY WEIGHTS
c
c     This routine computes the integration weights for the integrals
c     over the frequency. In principle this is a very general thing
c     and should not be here in this specialized module. But for now
c     this is the way we do it. Be sure to have generated the frequency
c     grid before calling this routine.
c
c     The frequency weight is such that
c                        __
c         /              \
c         | f(nu) dnu  =  >  f(nu_i) wgt_i
c         /              /_
c
c     We use the simple trapezium rule for the integral. 
c
c     --------------------------------------------------------------
      subroutine make_dust_freq_weights()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
c
      integer ifreq
c
      if(freq_nr.lt.1) then
          write(*,*) 'Zero or negative Number of freqs impossible'
          stop 99
      elseif(freq_nr.eq.1) then
          dust_freq_wgt(1) = 1.d0
cc          write(*,*) 'PROBLEM in dust module:'
cc          write(*,*) '  Attempt to construct frequency integration'
cc          write(*,*) '  weights for case of 1-frequency computation.'
cc          write(*,*) '  This is inconsistent. Are you sure you have '
cc          write(*,*) '  not accidently forgotten to switch off the '
cc          write(*,*) '  dust-temperature solver? This can be done '
cc          write(*,*) '  by setting iradproc_dust to -1 instead of 1.'
cc          stop 13
          return
      elseif((freq_nr.ge.2).and.(freq_nr.le.5)) then
          write(*,*) 'WARNING for dust temperature solver: '
          write(*,*) '  Too few frequency bins for consistently solving'
          write(*,*) '  dust temperatures! Continuuing while hoping '
          write(*,*) '  for the best....'
          dust_warn_few_freqs = 1
      else
          dust_freq_wgt(1) = 0.5d0 * abs( freq_nu(2) - freq_nu(1) )
          dust_freq_wgt(freq_nr) = 0.5d0 * abs( freq_nu(freq_nr) 
     %                              - freq_nu(freq_nr-1) )
          do ifreq=2,freq_nr-1
              dust_freq_wgt(ifreq) = 0.5d0 * 
     %             abs( freq_nu(ifreq+1) - freq_nu(ifreq-1) )
          enddo
          dust_frwgt_read = 1
      endif
c     
      end





#ifdef DUST_PRIVATE_EXTERNALS

c     --------------------------------------------------------------
c         READ THE FREQUENCY ARRAY (FOR USE OUTSIDE OF RADICAL)
c
c     If this module is used outside of RADICAL, then this module
c     itself will provide a routine for reading the frequency grid.
c     --------------------------------------------------------------
      subroutine read_frequencies()
      implicit none
c
#include "common_grid.h"
c
      integer inu
      doubleprecision frnu,dfrnu
c
      open(unit=1,file='frequency.inp',status='old',err=701)
      read(1,*) freq_nr
      if(freq_nr.gt.FRSIZE_FREQ) then
          write(*,*) 'ERROR: frequency.inp has more gridpoints'
          write(*,*) '       than FRSIZE_FREQ. '
          write(*,*) '       Recompile with larger FRSIZE_FREQ'
          stop 13
      endif          
      do inu=1,freq_nr
          read(1,*) freq_nu(inu)
      enddo
      close(1)
      goto 710
  701 continue
      write(*,*) 'Could not open file frequency.inp'
      stop 13
  710 continue
      return
      end


#endif





c     ==============================================================
c
c                        RADICAL-SPECIFIC ROUTINES
c
c     These routines are specific to the application as a module 
c     within RADICAL, the 2-D transfer program. If you wish to use
c     the dust opacity module in other programs, then remove the
c     stuff below.
c
c     ==============================================================

#     ifdef DUST_RADICAL

c     --------------------------------------------------------------
c                      THE DUST CONDENSATION SWITCH
c
c     This routine takes care of the minimum and maximum temperatures
c     at which the particular dust species can exist, and thereby
c     provide opacity. One way is to simply abruptly switch off the
c     opacity beyond tmax and below tmin. Another way is to multiply
c     by a powerlaw reduction factor, so that beyond tmax and below
c     tmin, the opacity drops off as a powerlaw as a function of 
c     temperature. The e-folding delta-temp with which this opacity
c     should then drop is determined by the dtmax and dtmin.
c     --------------------------------------------------------------
      function dust_condense(isize,ispec,it,ir)
      implicit none
      doubleprecision dust_condense
c
#include "common_grid.h"
#include "common_dust.h"
c
      doubleprecision temp,tmin,tmax,dtmin,dtmax,plindex
      integer it,ir,ispec,isize,ntemp,itemp
c
      dust_condense = 1.d0
      ntemp = dust_nr_temp(ispec)
      tmin  = dust_temprange_low(1,ispec)
      tmax  = dust_temprange_high(ntemp,ispec)
      dtmin = dust_dtmin(ispec)
      dtmax = dust_dtmax(ispec)
      temp  = dust_temp(isize,ispec,it,ir)
      if(dtmin.eq.0.d0) then
c
c         Hard switch-off at low temperature
c
          if(temp.lt.tmin) then
              dust_condense = 0.d0
          endif
      else
c
c         Smooth switch-off at low temperature
c
          plindex       = log(1.d0-(dtmin/tmin))
          dust_condense = (tmin/temp)**plindex
          if(dust_condense.gt.1.d0) then
              stop 9932
          endif
      endif
      if((dtmax.eq.0.d0).and.(tmax.gt.0.d0)) then
c     
c         Hard switch-off at high temperature
c     
          if(temp.gt.tmax) then
              dust_condense = 0.d0
          endif
      else
c     
c         Smooth switch-off at high temperature
c     
          plindex       = log(1.d0+(dtmax/tmax))
          dust_condense = (tmin/temp)**plindex
          if(dust_condense.gt.1.d0) then
              stop 9932
          endif              
      endif
c
c     If ntemp > 1, then we need to take care of the different
c     temperature regimes. 
c
      if(ntemp.gt.1) then
          write(*,*) 'ERROR: This should not have occured'
          write(*,*) '       Please warn me (Dullemond).'
          stop 54345
      endif
c
      return
      end




#ifndef EXCLUDE_DUST_ITERATION


c     --------------------------------------------------------------
c                 COMPUTE THE EMISSION FROM THE DUST
c
c     For given dust temperatures we can compute the total emissivity
c     from the dust grains. We do this for all the dust grains and
c     add them up.
c
c     NOTE: The philosophy of the RADICAL program is that only the
c           static (non-changing) sources are put into slte_src and
c           slte_alp. All other things (including the dust emissivities
c           since the dust temperatures are changing) are put into the
c           scat_src and scat_alp. The "scat" should therefore not be
c           taken as literally "scattering", but more like all things
c           that are varying per iteration.
c     --------------------------------------------------------------
      subroutine calc_emis_dust()
      implicit none
c
#include "common_grid.h"
#include "common_radfield.h"
#include "common_source.h"
#include "common_dust.h"
#include "common_setup.h"
c
      integer ir,it,inu,itmax,ispec,isize
      doubleprecision dummy,dens,nu,bplanck,dust_condense
      doubleprecision kappa,temp,find_dust_kappa
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
c     If grey opacity mode (freq-integrated), then it's simple
c
      if(freq_nr.eq.1) then
          do ir=1,irsi_frsizex
              do it=1,itmax
                  do ispec=1,dust_nr_species
                      do isize=1,dust_nr_size(ispec)
c
c                         Dust temperature
c
                          temp = dust_temp(isize,ispec,it,ir)
c
c                         Opacity
c
#                         ifndef DUST_OPAC_TEMPDEP
                          kappa = dust_condense(isize,ispec,it,ir) *
     %                         dust_kappawgt_abs(1,1,isize,ispec)
#                         else
                          kappa = find_dust_kappa(1,isize,ispec,
     %                                 temp,1,0)
#                         endif
c
c                         The emission: alpha * (sigma/pi) * T^4
c
                          scati_src(1,it,ir) = scati_src(1,it,ir) +
     %                      dust_rho(ispec,it,ir) * kappa *
     %                      1.80491254763d-5 * temp**4
                      enddo
                  enddo
              enddo
          enddo
c
c         Return 
c
          return
      endif
c
c     Full frequency dependent stuff
c
      do ir=1,irsi_frsizex
          do it=1,itmax
c
c             Now loop over all frequencies
c
              do inu=1,freq_nr
c
                  nu = freq_nu(inu)
c
c                 Now a loop over dust species / sizes
c
                  do ispec=1,dust_nr_species
                      do isize=1,dust_nr_size(ispec)
c
c                         Dust temperature
c
                          temp = dust_temp(isize,ispec,it,ir)
c
c                         Opacity
c
#                         ifndef DUST_OPAC_TEMPDEP
                          kappa = dust_condense(isize,ispec,it,ir) *
     %                         dust_kappawgt_abs(inu,1,isize,ispec)
#                         else
                          kappa = find_dust_kappa(inu,isize,ispec,
     %                                 temp,1,0)
#                         endif
c
c                         The emission
c
                          scati_src(inu,it,ir) = scati_src(inu,it,ir) +
     %                      dust_rho(ispec,it,ir) * kappa *
     %                      bplanck(temp,nu)
                      enddo
                  enddo
              enddo
          enddo
      enddo
c
      return
      end



c     --------------------------------------------------------------
c                    INTEGRAL OF EMISSION-ABSORPTION
c     --------------------------------------------------------------
      function dust_integral_emisabs(temp)
      implicit none
c
      common/dustintegral/di_ir,di_it,di_isize,di_ispec
      integer di_ir,di_it,di_isize,di_ispec
      doubleprecision temp,dust_integral_emisabs,bplanck,kappa
      doubleprecision find_dust_kappa
c
#include "common_grid.h"
#include "common_radfield.h"
#include "common_dust.h"
#include "common_main.h"
c     
      integer inu
      doubleprecision dummy,corr
c
c     If the flux conservation trick is active, then put corr to 
c     the appropriate value, else to 1.0
c
      if(iter_fluxcons.le.0.d0) then
          corr=1.d0
      else
          corr=1.d0-fl_corrfact(di_it,di_ir)
c          write(*,*) 'corr=',corr
      endif
c
c     Now evaluate the integral
c
      dummy = 0.d0
      do inu=1,freq_nr
#         ifndef DUST_OPAC_TEMPDEP
             kappa = dust_kappawgt_abs(inu,1,di_isize,di_ispec) 
#         else
             kappa = find_dust_kappa(inu,di_isize,di_ispec,temp,-1,0)
#         endif
          dummy = dummy + kappa * dust_freq_wgt(inu)
     %      * ( bplanck(temp,freq_nu(inu)) -
     %             corr*intmom_0(inu,di_it,di_ir) )
      enddo
ccc
cc      if(iter_fluxcons.gt.0.d0) then
cc          dum2 = 0.d0
cc          do inu=1,freq_nr
cc#             ifndef DUST_OPAC_TEMPDEP
cc              kappa = dust_kappawgt_abs(inu,1,di_isize,di_ispec) 
cc#             else
cc              kappa = find_dust_kappa(inu,di_isize,di_ispec,temp,-1,0)
cc#             endif
cc              dum2 = dum2 + kappa * dust_freq_wgt(inu)
cc     %             * intmom_0(inu,di_it,di_ir)
cc          enddo
cc          dum2 = dum2 * fl_corrfact(di_it,di_ir)
cc      endif
c
      dust_integral_emisabs = dummy
      return
      end



c     --------------------------------------------------------------
c          SOLVE ONE SINGLE DUST SPECIES AT ONE SINGLE POSITION
c     --------------------------------------------------------------
      function solve_one_dust_temp(isize,ispec,ir,it)
      implicit none
      doubleprecision solve_one_dust_temp,dust_integral_emisabs,zbrent
      integer isize,ispec,ir,it
      doubleprecision dens,temp
      external dust_integral_emisabs
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
#include "common_dust.h"
#include "common_main.h"
#include "common_setup.h"
#include "common_check.h"
c
      common/dustintegral/di_ir,di_it,di_isize,di_ispec
      integer di_ir,di_it,di_isize,di_ispec
c
      doubleprecision dum1,dum2
      integer i
c
c     Put the ir, it, ispec, isize into the common
c
      di_ir      = ir
      di_it      = it
      di_isize   = isize
      di_ispec   = ispec
c
c     Find the density
c
      call get_phys_var(ir,it)
c
c     If in grey mode (freq-integrated) then the answer is simple
c
      if(freq_nr.eq.1) then
          if(iter_fluxcons.le.0.d0) then
              solve_one_dust_temp = 1.53421406189d1 * 
     %             intmom_0(1,it,ir)**0.25
          else
cc              write(*,*) fl_corrfact(it,ir)
              solve_one_dust_temp = 1.53421406189d1 * 
     %             (intmom_0(1,it,ir)*(1.d0-fl_corrfact(it,ir)))**0.25
          endif
          return
      endif
c
c     Now check if there is indeed sufficient radiation for a
c     solution.
c
      dum1 = dust_integral_emisabs(DUST_TEMP_MIN)
      dum2 = dust_integral_emisabs(DUST_TEMP_MAX)
      if(dum1.gt.0.d0) then
ccc
ccc    (removed this on 15-02-01 to avoid crash of CG97 models)
ccc
ccc          write(*,*) 'Error: dust temperature < 0.01 K'
ccc          write(*,*) '       at ir=',ir,' it=',it
ccc          write(*,*) 'The mean intensities are:'
ccc          write(*,*) (intmom_0(i,di_it,di_ir),i=1,freq_nr)
ccc#         ifdef SAFETY_CHECKS_ACTIVE
ccc          if(selfch_taumax.gt.SAFETY_CHECK_TAUMAX) then
ccc              write(*,*) '--------------------------------------------'
ccc              write(*,*) 'Possible cause:'
ccc              write(*,*) '  Found rather high one-cell optical depths'
ccc              write(*,*) '  Largest optical depth found is:',
ccc     %                     selfch_taumax
ccc              write(*,*) '--------------------------------------------'
ccc          endif
ccc#         endif
ccc          stop 1
c
c         Dust temperature too low. Simply put it to zero and return.
c
          dust_warn_zero_temp = dust_warn_zero_temp + 1
          solve_one_dust_temp = 0.d0
          return
      endif
      if(dum2.lt.0.d0) then
          write(*,*) 'Error: dust temperature > ',DUST_TEMP_MAX,' K'
          write(*,*) '       at ir=',ir,' it=',it
          stop 1
      endif
c
c     Now solve for the dust temperature
c
cccc      solve_one_dust_temp = zbrent(dust_integral_emisabs,
cccc     %                              2.d-2,3.d3,1.d-1) 
      solve_one_dust_temp = zbrent(dust_integral_emisabs,
     %          DUST_TEMP_MIN,DUST_TEMP_MAX,DUST_TEMP_ACCUR)
c
      return
      end



c     --------------------------------------------------------------
c                         SOLVE DUST TEMPERATURES
c
c     This routine solves the temperatures of all the dust species
c     and sizes, for a given mean intensity int_moment_0(). 
c     --------------------------------------------------------------
      subroutine solve_dust_temp()
      implicit none
c
      doubleprecision solve_one_dust_temp
c
#include "common_grid.h"
#include "common_source.h"
#include "common_radfield.h"
#include "common_dust.h"
#include "common_setup.h"

#ifdef PARALLEL
!parallel headerfile
#include "mpif.h"

      integer ierr        !error code
      integer MyRank      !rank of this process
      integer NrProc      !number of processes
      integer MyNrR       !number of rows for this process
      integer BeginR      !start of rows
      integer EndR        !end of rows
      integer SendNrR     !nr of rows that are send to this process
      integer SendBeginR  !start of rows of sending process
      integer SendEndR    !end of rows domain of sending process
      integer NrDustTemp  !number of dust_temp items that will be send or recv.
      integer i,j,k,l,p,m !counters

      !buffer for sending and receiving data
      doubleprecision SendRecvArray(1:DUST_SIZE_MAX*DUST_SPECIES_MAX*
     %                                FRSIZE_Y)
#endif
c
      integer ir,it,inu,itmax,ispec,isize
      doubleprecision dummy,dens,temp,nu
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

#ifdef PARALLEL
      !get rank of this process and the number of processes
      call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, NrProc, ierr)

      !decompose R domain
      MyNrR = irsi_frsizex/NrProc
      if((MyRank+1).LE.(Mod(irsi_frsizex, NrProc))) then
         MyNrR = MyNrR + 1
         BeginR = MyNrR*MyRank + 1
         EndR = BeginR + MyNrR - 1
      else
         EndR = irsi_frsizex - ((NrProc-MyRank-1)*MyNrR)
         BeginR = EndR - MyNrR + 1
      end if
#endif
c
c     Loop over space
c
#ifdef PARALLEL
      do ir=BeginR,EndR
#else
      do ir=1,irsi_frsizex
#endif
          do it=1,itmax
c
c             Loop over dust species / sizes
c
              do ispec=1,dust_nr_species
                  do isize=1,dust_nr_size(ispec)
                      dust_temp(isize,ispec,it,ir) = 
     %                 solve_one_dust_temp(isize,ispec,ir,it)
                  enddo
              enddo
          enddo
      enddo

#ifdef PARALLEL
      !determine the number of dust_temp items that will be sent or received
      NrDustTemp = 0
      do i=1, dust_nr_species
         NrDustTemp = NrDustTemp + dust_nr_size(i)
      end do
      NrDustTemp = NrDustTemp*itmax

      !broadcast and receive data
      do p=0, NrProc-1
         if(p.EQ.MyRank) then
           do i=BeginR, EndR
             !put dust temperature data in send-receive array
             m=0
             do j=1, itmax
                do k=1, dust_nr_species
                   do l=1, dust_nr_size(k)
                      m = m+1 
                      SendRecvArray(m) = dust_temp(l,k,j,i)
                   end do
                end do
             end do

             !broadcast data
             call MPI_Bcast(SendRecvArray(1), NrDustTemp, 
     %               MPI_DOUBLE_PRECISION, p, MPI_COMM_WORLD, ierr)
           end do

         else
            !determine how much and which rows will be sent
            SendNrR = irsi_frsizex/NrProc
            if((p+1).LE.(Mod(irsi_frsizex, NrProc))) then
               SendNrR = SendNrR + 1
               SendBeginR = SendNrR*p + 1
               SendEndR = SendBeginR + SendNrR - 1
            else
               SendEndR = irsi_frsizex - ((NrProc-p-1)*SendNrR)
               SendBeginR = SendEndR - SendNrR + 1
            end if

            !receive data
            do i=SendBeginR, SendEndR
               call MPI_Bcast(SendRecvArray(1),NrDustTemp,
     %                 MPI_DOUBLE_PRECISION, p, MPI_COMM_WORLD, ierr)

               !put received data into dust temperature array
               m=0
               do j=1, itmax
                  do k=1, dust_nr_species
                     do l=1, dust_nr_size(k)
                        m = m+1
                        dust_temp(l,k,j,i) = SendRecvArray(m)
                     end do
                  end do
               end do
            end do
 
        end if

      end do
#endif

      return
      end


c     --------------------------------------------------------------
c                    WRITE DUST TEMPERATURES TO FILES
c     --------------------------------------------------------------
      subroutine write_dust_temp(iter,flnr)
      implicit none
c
      doubleprecision solve_one_dust_temp
      integer flnr,iter
c
#include "common_grid.h"
#include "common_isave.h"
#include "common_source.h"
#include "common_radfield.h"
#include "common_dust.h"
#include "common_setup.h"
c
      integer ir,it,inu,itmax,ispec,isize,imirt
      integer filenr,nrdustgrain
      doubleprecision dummy,dens,temp,nu
      character*80 filename,base,ext
c
c     The usual spatial grid stuff
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
c     Make filename
c
      filenr = flnr
      if(filenr.eq.-1) then
          write(*,*) 'No append for dusttemp.dat'
          stop 13
      elseif(filenr.eq.-2) then
          filename = 'dusttemp_final.dat'
      elseif(filenr.eq.-4) then
          filename = 'dusttemp_nocons.dat'
      elseif(filenr.eq.-5) then
          filename = 'dusttemp_cons.dat'
      else
          base='dusttemp_'
          ext ='.dat'
          call make_indexed_filename(base,filenr,ext,filename)
      endif
c
c     Open file
c
      open(unit=1,file=filename,status='unknown',err=701)
      write(1,*) dust_nr_species,irsi_frsizex,itmax,imirt
      write(1,*) ' '
c
c     Loop over dust species / sizes
c
      nrdustgrain = 0
      do ispec=1,dust_nr_species
          write(1,*) dust_nr_size(ispec)
          do isize=1,dust_nr_size(ispec)
              do ir=1,irsi_frsizex
                  do it=1,itmax
                      write(1,*) dust_temp(isize,ispec,it,ir)
                  enddo
              enddo
              write(1,*) ' '
              nrdustgrain = nrdustgrain + 1
          enddo
          write(1,*) ' '
      enddo
      close(1)
c
c     If this save was succesful, and if it was not the final save,
c     then put filenr into lastsave
c
      if(filenr.gt.0) then
          lastsave_dusttemp = filenr
      endif
c
c     Write dusttemp.info file to tell which files are written
c     so far and which is the current one
c
      open(unit=1,file='dusttemp.info',status='unknown',err=702)
c
c     The number of the latest file (-2 means: 'intmean_final.dat'),
c     some info on nr of dust species/size, and the iteration number.
c
      write(1,*) filenr
      write(1,*) nrdustgrain
      write(1,*) DUST_SIZE_MAX
      write(1,*) lastsave_dusttemp
      write(1,*) iter
      close(1)
c
c
      goto 710
  701 continue
      write(*,*) 'Could not open file ',filename
      stop 13
  702 continue
      write(*,*) 'Could not open file dusttemp.info'
      stop 13
  710 continue
      return
      end




c     --------------------------------------------------------------
c                   WRITE WHERE DUST EXISTS OR NOT
c     --------------------------------------------------------------
      subroutine write_dust_condens()
      implicit none
c
#include "common_grid.h"
c#include "common_source.h"
#include "common_dust.h"
c#include "common_setup.h"
c
      integer ir,it,itmax,ispec,isize,imirt,icond
      integer nrdustgrain
      doubleprecision dust_condense
      character*80 ch,filename
c
c     The usual spatial grid stuff
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
      filename = 'dustcondens.dat'
c
c     Open file
c
      open(unit=1,file=filename,status='unknown')
      write(1,*) dust_nr_species,irsi_frsizex,itmax,imirt
      write(1,*) ' '
c
c     Loop over dust species / sizes
c
      nrdustgrain = 0
      do ispec=1,dust_nr_species
          write(1,*) dust_nr_size(ispec)
          do isize=1,dust_nr_size(ispec)
              do ir=1,irsi_frsizex
                  do it=1,itmax
                      icond = dust_condense(isize,ispec,it,ir)
                      write(1,*) icond
                  enddo
              enddo
              write(1,*) ' '
              nrdustgrain = nrdustgrain + 1
          enddo
          write(1,*) ' '
      enddo
      close(1)
c
      end



c     --------------------------------------------------------------
c                DO THE DUST SCATTERING AND TEMPERATURES
c
c     This routine performs all the necessary steps that have to be
c     taken in between two iterations of RT for the dust. 
c
c     ARGUMENT:
c       iprecond       If 1, then use dust ALI for solving the 
c                      dust temperature and scattering.
c
c
c     REMARK: This routine is incompatible with all other processes.
c             This is because it erases the scati_src() arrays. This
c             should be changed in the future, but then all other
c             processes should also be changed.
c     --------------------------------------------------------------
c$$$      subroutine do_dust_temp_scat(iprecond)
c$$$      implicit none
c$$$      integer iprecond
c$$$c
c$$$#include "common_source.h"
c$$$#include "common_grid.h"
c$$$#include "common_dust.h"
c$$$c
c$$$      if(iradproc_dust.ne.0) then
c$$$c
c$$$c         Check if ALI and temp dep opacities
c$$$c
c$$$          if((dust_opacity_tempdep.gt.0).and.(iprecond.ne.0)) then
c$$$              write(*,*) 'Sorry, but for now the dust ALI is not'
c$$$              write(*,*) 'yet compatible with temperature-dependent'
c$$$              write(*,*) 'opacities.'
c$$$              stop 13
c$$$          endif
c$$$c
c$$$c         Compute the dust temperatures from the given 
c$$$c         mean intensities
c$$$c
c$$$          if(iradproc_dust.eq.1) then
c$$$              call solve_dust_temp()
c$$$          endif
c$$$c
c$$$c         Now the ALI stuff
c$$$c
c$$$#         ifdef DUST_ALI
c$$$          if(iprecond.ne.0) then
c$$$              call dust_ali_accelerate()
c$$$          endif
c$$$#         endif
c$$$c
c$$$c         If the tmin and tmax of the dust grains are active, then
c$$$c         we must recompute the opacitities
c$$$c
c$$$          if(dust_opacity_tempdep.gt.0) then
c$$$              write(*,*) '   (Recomputing dust opacities)'
c$$$              call make_scat_alpha(1,freq_nr)
c$$$          endif
c$$$c
c$$$c         Clear the scattering emission. This is necessary because
c$$$c         the various non-lte and scattering processes simply add their
c$$$c         contribution to the source, so we must first reset it.
c$$$c
c$$$          call scat_src_clear()
c$$$c
c$$$c         Compute the isotropic scattering contribution
c$$$c
c$$$          call do_dust_iso_scat()
c$$$c
c$$$c         Add the thermal emission of the dust
c$$$c
c$$$          call calc_emis_dust()
c$$$c
c$$$      endif
c$$$c
c$$$      end
c$$$





c     --------------------------------------------------------------
c                   DO THE DUST ISOTROPIC SCATTERING
c
c     This routine uses the mean intensity J (as can be found in the
c     array intmom_0()) to compute the scattering source term for
c     isotropic dust scattering, which is simply:
c 
c           j_{dustscat}(\nu) = \alpha_{dustscat}(\nu) J(\nu)
c
c     --------------------------------------------------------------
c$$$      subroutine do_dust_iso_scat()
c$$$      implicit none
c$$$c
c$$$#include "common_grid.h"
c$$$#include "common_radfield.h"
c$$$#include "common_source.h"
c$$$#include "common_dust.h"
c$$$c
c$$$      doubleprecision dustkappa,dust_condense,temp,find_dust_kappa
c$$$      integer ir,it,inu,itmax,ispec,isize
c$$$c
c$$$      itmax   = irsi_frsizey
c$$$c
c$$$#     ifdef MIRROR_THETA
c$$$      itmax = (irsi_frsizey+1)/2
c$$$#     endif
c$$$#     ifdef RADGRID_ONEDIM
c$$$      itmax = 1
c$$$#     endif
c$$$c
c$$$      do ir=1,irsi_frsizex
c$$$          do it=1,itmax
c$$$              do inu=1,freq_nr
c$$$                  do ispec=1,dust_nr_species
c$$$                      dustkappa = 0.d0
c$$$#                     ifndef DUST_OPAC_TEMPDEP
c$$$c
c$$$c                     The old way...
c$$$c
c$$$                      do isize=1,dust_nr_size(ispec)
c$$$                          dustkappa = dustkappa +
c$$$     %                       dust_kappawgt_scat(inu,1,isize,ispec) *
c$$$     %                       dust_condense(isize,ispec,it,ir)
c$$$                      enddo
c$$$#                     else
c$$$c
c$$$c                     The new way...
c$$$c
c$$$                      do isize=1,dust_nr_size(ispec)
c$$$                          temp = dust_temp(isize,ispec,it,ir)
c$$$                          dustkappa = dustkappa +
c$$$     %                       find_dust_kappa(inu,isize,ispec,temp,0,1)
c$$$                      enddo                      
c$$$#                     endif
c$$$                      scati_src(inu,it,ir) = scati_src(inu,it,ir) 
c$$$     %                   + dust_rho(ispec,it,ir) * 
c$$$     %                     dustkappa * intmom_0(inu,it,ir)
c$$$                  enddo
c$$$              enddo
c$$$          enddo
c$$$      enddo
c$$$c
c$$$      return
c$$$      end
c$$$      


#     ifndef NO_DUST_NG
c     --------------------------------------------------------------
c             COPY THE DUST TEMPERATURES IN BACKUP ARRAYS
c
c     For Ng acceleration the dust temperatures must be stored 
c     somewhere. 
c
c     istore1 = -1 means from dust_temp
c     istore2 = -1 means to dust_temp
c     --------------------------------------------------------------
      subroutine dusttemp_copy(istore1,istore2)
      implicit none
      integer istore1,istore2
c
#include "common_grid.h"
#include "common_dust.h"
c
      integer ir,it,inu,itmax,ispec,isize
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
      if((istore1.eq.-1).and.(istore2.ne.-1)) then
          do ir=1,irsi_frsizex
              do it=1,itmax
                  do ispec=1,dust_nr_species
                      do isize=1,dust_nr_size(ispec)
                          ngbk_dust_temp(isize,ispec,it,ir,istore2) = 
     %                         dust_temp(isize,ispec,it,ir)
                      enddo
                  enddo
              enddo
          enddo
      elseif((istore1.ne.-1).and.(istore2.eq.-1)) then
          do ir=1,irsi_frsizex
              do it=1,itmax
                  do ispec=1,dust_nr_species
                      do isize=1,dust_nr_size(ispec)
                          dust_temp(isize,ispec,it,ir) = 
     %                         ngbk_dust_temp(isize,ispec,it,ir,istore1)
                      enddo
                  enddo
              enddo
          enddo
      elseif((istore1.ne.-1).and.(istore2.ne.-1)) then
          do ir=1,irsi_frsizex
              do it=1,itmax
                  do ispec=1,dust_nr_species
                      do isize=1,dust_nr_size(ispec)
                          ngbk_dust_temp(isize,ispec,it,ir,istore2) = 
     %                         ngbk_dust_temp(isize,ispec,it,ir,istore1)
                      enddo
                  enddo
              enddo
          enddo
      else 
          write(*,*) 'Error in dusttemp_copy'
          stop 13
      endif
c
      end


c     --------------------------------------------------------------
c                     NG ACCELLERATION FOR DUST
c
c     This version of the Ng acceleration does the acceleration on
c     the dust temperatures.
c     --------------------------------------------------------------
      subroutine dust_ng_accellerate()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
c
      doubleprecision a1,a2,b1,b2,c1,c2,a,b,c,x0,x1,x2,x3,w
      doubleprecision dum1,dum2,dum3,dummy
      integer ir,it,inu,itmax,ispec,isize,number_invalid
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
      do ir=1,irsi_frsizex
          do it=1,itmax
              do ispec=1,dust_nr_species
                  do isize=1,dust_nr_size(ispec)
                      x0   = dust_temp(isize,ispec,it,ir)**4
                      if(x0.gt.0.d0) then
                          x1   = ngbk_dust_temp(isize,ispec,it,ir,1)**4
                          x2   = ngbk_dust_temp(isize,ispec,it,ir,2)**4
                          x3   = ngbk_dust_temp(isize,ispec,it,ir,3)**4
                          w    = 1.d0 / dust_temp(isize,ispec,it,ir)**4
                          dum1 = x0 - 2.d0*x1 + x2
                          dum2 = x0 - x1 - x2 + x3
                          dum3 = x0 - x1
                          a1   = a1 + dum1 * dum1 * w
                          b1   = b1 + dum1 * dum2 * w
                          c1   = c1 + dum1 * dum3 * w
                          a2   = b1
                          b2   = b2 + dum2 * dum2 * w
                          c2   = c2 + dum2 * dum3 * w
                      endif
                  enddo
              enddo
          enddo
      enddo
      a = ( c1*b2 - c2*b1 ) / ( a1*b2 - a2*b1 ) 
      b = ( c2*a1 - c1*a2 ) / ( a1*b2 - a2*b1 ) 
      c = 1.d0 - a - b
c
c     Check the numbers
c
      if((number_invalid(a).gt.0).or.(number_invalid(b).gt.0)) then
          write(*,*) 'ERROR in dust-Ng acceleration: invalid numbers'
          write(*,*) '    a = ',a,'  b = ',b
          write(*,*) '   a1 = ',a1,' a2 = ',a2
          write(*,*) '   b1 = ',b1,' b2 = ',b2
          write(*,*) '   c1 = ',c1,' c2 = ',c2
          stop 13
      endif
c
c     Write a message
c
      write(*,*) '       Ng : a = ',a,'  b = ',b
c
c     Now perform the accelleration step
c
      do ir=1,irsi_frsizex
          do it=1,itmax
              do ispec=1,dust_nr_species
                  do isize=1,dust_nr_size(ispec)
                      dummy = c * dust_temp(isize,ispec,it,ir)**4
     %                     + a * ngbk_dust_temp(isize,ispec,it,ir,1)**4
     %                     + b * ngbk_dust_temp(isize,ispec,it,ir,2)**4
                      if(dummy.lt.0.d0) then
                          write(*,*) 'Ng ERROR: negative T^4 ',
     %                         ir,it,ispec,isize
                          write(*,*) '   a,b,c = ',a,b,c
                          write(*,*) 'Skipping this Ng step'
                          goto 110
                      else
                          ngbk_dust_temp(isize,ispec,it,ir,3) =
     %                         sqrt(sqrt(dummy))
                      endif
                  enddo
              enddo
          enddo
      enddo
      do ir=1,irsi_frsizex
          do it=1,itmax
              do ispec=1,dust_nr_species
                  do isize=1,dust_nr_size(ispec)
                      dust_temp(isize,ispec,it,ir) = 
     %                     ngbk_dust_temp(isize,ispec,it,ir,3)
                  enddo
              enddo
          enddo
      enddo
c
 110  continue
      end

#endif /* IFNDEF NO_DUST_NG */



#endif /* IFNDEF EXCLUDE_DUST_ITERATION */

c     --------------------------------------------------------------
c              CHECK WHETHER THE DUST DENSITIES HAVE BEEN SET
c
c     Sometimes the user may forget to implement the proper setting
c     of the dust_rho() arrays. To avoid confusion why RADICAL does
c     not produce any result, it can be checked whether the dust_rho
c     array has any non-zero value.
c     --------------------------------------------------------------
      subroutine dust_check_grid()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
c
      integer ir,it,ispec,itmax,imirt
      logical founddens
c
      call set_mirror_limits_theta(itmax,imirt)
      founddens=.false.
      do ispec=1,DUST_SPECIES_MAX
          do ir=1,irsi_frsizex
              do it=1,itmax
                  if(dust_rho(ispec,it,ir).ne.0.d0) then
                      founddens = .true.
                  endif
              enddo
          enddo
      enddo
      if(.not.founddens) then
          write(*,*) 'ERROR: I presume that you wish to have '
          write(*,*) '       some dust emission/absorption/scatteirng'
          write(*,*) '       going on. But you have not specified the'
          write(*,*) '       dust densities. The entire dust density'
          write(*,*) '       array called dust_rho(ispec,it,ir) is'
          write(*,*) '       empty. Please check your setup in the '
          write(*,*) '       setup.F module to see if you set the'
          write(*,*) '       rho_density(ispec) array. '
          write(*,*) ''
          write(*,*) '       NOTE: Do NOT directly specify the '
          write(*,*) '       dust_rho(eispec,it,ir) array, since this'
          write(*,*) '       will be filled with copies of the '
          write(*,*) '       rho_density(ispec) array at each it,ir!'
          write(*,*) '       So, set rho_density(ispec) in the '
          write(*,*) '       setup.F/calc_phys_var() routine.'
          write(*,*) '       PS: The read_dust_density() routines in'
          write(*,*) '       fact DO directly specify dust_rho(), but'
          write(*,*) '       they set the flag setup_dustdens_read to 1'
          write(*,*) '       which disables the copying of the '
          write(*,*) '       rho_density into dust_rho. You could do so'
          write(*,*) '       as well, if you choose.'
          stop 13
      endif
c
      end


c     --------------------------------------------------------------
c                     COMPUTE THE TOTAL DUST MASS
c
c     Computes:
c
c              2 pi /rout        /+1
c          M = ---- |     d(r^3) | d(cos(theta)) rho_dust 
c               3   /rin         /-1
c
c     ARGUMENTS:
c       ispec             >0  The mass of this species only
c                         <0  The total mass of all species
c       icell             =0  Regard each gridpoint as cell-center
c                         =1  Regard gridpoints as cell-corner
c     --------------------------------------------------------------
      function total_dust_mass(ispec,icell)
      implicit none
      doubleprecision total_dust_mass,mass
      integer ispec,iispec,icell
c
#include "common_grid.h"
#include "common_dust.h"
c
      integer ir,it,itmax,imirt
      logical founddens
c
      call set_mirror_limits_theta(itmax,imirt)
      mass = 0.d0
      if(ispec.lt.0) then
c
c         All species summed up:
c     
          do iispec=1,dust_setup_nrspecies
              if(icell.eq.0) then 
                  do ir=1,irsi_frsizex
                      do it=1,itmax
                          mass = mass + 
     %                         2.09439510239d0 * 
     %                         dust_rho(iispec,it,ir) *
     %                         (rsi_x_i(ir+1,1)**3-rsi_x_i(ir,1)**3) *
     %                         abs(cos(rsi_x_i(it+1,2))-
     %                         cos(rsi_x_i(it,2)))
                      enddo
                  enddo
              else
                  do ir=1,irsi_frsizex-1
                      do it=1,itmax-1
                          mass = mass + 
     %                         2.09439510239d0 * 0.25 *
     %                         ( dust_rho(iispec,it,ir) +
     %                           dust_rho(iispec,it,ir+1) +
     %                           dust_rho(iispec,it+1,ir) +
     %                           dust_rho(iispec,it+1,ir+1) ) *
     %                         (rsi_x_c(ir+1,1)**3-rsi_x_c(ir,1)**3) *
     %                         abs(cos(rsi_x_c(it+1,2))-
     %                         cos(rsi_x_c(it,2)))
                      enddo
                  enddo
              endif
          enddo
          if(imirt.ne.0) then
              mass = mass * 2.d0
          endif
      else
c
c         Only this species
c     
          if(icell.eq.0) then 
              do ir=1,irsi_frsizex
                  do it=1,itmax
                      mass = mass + 
     %                     2.09439510239d0 * dust_rho(ispec,it,ir) *
     %                     (rsi_x_i(ir+1,1)**3-rsi_x_i(ir,1)**3) *
     %                     abs(cos(rsi_x_i(it+1,2))-cos(rsi_x_i(it,2)))
                  enddo
              enddo
          else
              do ir=1,irsi_frsizex-1
                  do it=1,itmax-1
                      mass = mass + 
     %                     2.09439510239d0 * 0.25 * 
     %                     ( dust_rho(ispec,it,ir) +
     %                       dust_rho(ispec,it,ir+1) +
     %                       dust_rho(ispec,it+1,ir) +
     %                       dust_rho(ispec,it+1,ir+1) ) *
     %                     (rsi_x_c(ir+1,1)**3-rsi_x_c(ir,1)**3) *
     %                     abs(cos(rsi_x_c(it+1,2))-cos(rsi_x_c(it,2)))
                  enddo
              enddo
          endif
          if(imirt.ne.0) then
              mass = mass * 2.d0
          endif
      endif
c
      total_dust_mass = mass
      return
      end








c     ==============================================================
c
c                         DUST-LOCAL-ALI ROUTINES
c
c     These routines are meant for the local dust ALI method. This
c     method is not as powerful as non-local ALI, but it at least
c     enables us to do high optical depths without waiting forever.
c     In the future, after gaining some experience with the local
c     ALI method I would like to try out the tridiagonal method
c     as well. 
c
c     ==============================================================

#ifdef DUST_ALI


c     --------------------------------------------------------------
c                     MAIN ROUTINE FOR DUST ALI
c
c     This routine is the main routine for dust ali. It will loop 
c     over space and solve the dust ali at every location, by 
c     calling the routines below.
c
c     INPUT:
C       dust_tempold      This is the temperature at the start 
c                         of the Lambda Iteration step
c       dust_meanintold   The mean intensity at the start of the 
c                         Lambda Iteration step
c
c     WARNING: For the moment we don't separate direct starlight
c              from diffuse radiation. So we put Q=0.
c
c     --------------------------------------------------------------
      subroutine dust_ali_accelerate()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
#include "common_source.h"
#include "common_radfield.h"
c
      doubleprecision kappaabs(FRSIZE_FREQ),kappascat(FRSIZE_FREQ)
      doubleprecision epsilon(FRSIZE_FREQ),temp,q,find_dust_kappa
      integer inu,ir,it,itmax
c
c     Check
c
      if((dust_nr_species.gt.1).or.(dust_nr_size(1).gt.1)) then
          write(*,*) 'Sorry, dust ALI works at the moment only'
          write(*,*) '   for one dust species and size'
          stop 13
      endif
c
c     The usual stuff
c
      itmax   = irsi_frsizey
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#     endif
#     ifdef RADGRID_ONEDIM
      itmax = 1
#     endif
c
c     Loop over all space
c     
      do ir=1,irsi_frsizex
          do it=1,itmax
              do inu=1,freq_nr
                  temp = dust_temp(1,1,it,ir)
                  kappaabs(inu) = find_dust_kappa(inu,1,1,temp,1,0)
                  kappascat(inu) = find_dust_kappa(inu,1,1,temp,0,1)
                  epsilon(inu) = kappaabs(inu)/
     %                         (kappaabs(inu)+kappascat(inu))
              enddo
              q=0.d0   ! No separation of direct/indirect
cc####################################
c              if(ali_approx_lambda(1,it,ir).gt.0.99) then
c                  write(*,*) 'LAMBDA LARGE ',ir,it,
c     %                  ali_approx_lambda(1,it,ir)
c              endif
c              ali_approx_lambda(1,it,ir) = 0.d0
cc####################################
              call solve_local_dust_ali(freq_nr,freq_nu,dust_freq_wgt,
     %             ali_approx_lambda(1,it,ir),kappaabs,epsilon,q,
     %             dust_tempold(it,ir),dust_meanintold(1,it,ir),
     %             dust_temp(1,1,it,ir),intmom_0(1,it,ir),
     %             dust_warn_zero_temp)
          enddo
      enddo
c
      end



c     --------------------------------------------------------------
c                  SOLVE THE LOCAL DUST-ALI EQUATION
c
c     This is the core routine for the local non-linearized dust
c     ALI accelerator. It works both in the usual mode and in the
c     frequency-integrated mode. This routine is completely self-
c     contained (no common blocks).
c
c     This routine can be called once the complete Lambda Iteration
c     procedure for a single iteration has been completed. At the
c     end of this iteration, the solve_local_dust_ali() can then
c     be solved 
c
c     ARGUMENTS:
c       nf            Nr of frequencies
c       lambdastar    The diagonal approx operator (see calc_ali_accell)
c       kappaabs      The absorption opacity
c       epsilon       The ratio kappa_abs/(kappa_scat+kappa_abs)
c       q             The external source (re-emitted photons)
c       tempold       The temperature before applying LI
c       meanintold    The mean intensity before applying LI
c       temp          The temperature from Lambda Iteration
c       meanint       The mean intensity from Lambda Iteration
c
c     RESULTS:
c       temp          The temperature from Accelerated Lambda Iteration
c       meanint       The mean intensity from Accelerated Lambda Iteration
c       iwarn         The return of dust_warn_zero_temp, i.e. the 
c                     counter of how many error: zero dust temp
c
c     --------------------------------------------------------------
      subroutine solve_local_dust_ali(nf,nu,dnu,lambdastar,kappaabs,
     %             epsilon,q,tempold,meanintold,temp,meanint,iwarn)
      implicit none
      integer nf,iwarn
      doubleprecision lambdastar(nf),epsilon(nf),temp,q
      doubleprecision meanint(nf),kappaabs(nf)
      doubleprecision meanintold(nf),tempold
      doubleprecision nu(nf),dnu(nf)
c
      common/dustalisolve/dae_meanint,dae_lambdastar,dae_epsilon,
     %       dae_kappaabs,dae_rhs,dae_q,dae_nu,dae_dnu
      doubleprecision dae_meanint(FRSIZE_FREQ)
      doubleprecision dae_lambdastar(FRSIZE_FREQ)
      doubleprecision dae_epsilon(FRSIZE_FREQ)
      doubleprecision dae_kappaabs(FRSIZE_FREQ)
      doubleprecision dae_rhs(FRSIZE_FREQ)
      doubleprecision dae_nu(FRSIZE_FREQ)
      doubleprecision dae_dnu(FRSIZE_FREQ)
      doubleprecision dae_q
      common/idustalisolve/dae_nf
      integer dae_nf
c
      doubleprecision sspi,zbrent,dust_ali_equation
      doubleprecision bplanck,dum1,dum2
      parameter(sspi=1.80491254763d-5)     ! ss / pi
      external dust_ali_equation
c
      integer inu
c
c     Copy stuff to common blocks
c
      dae_nf=nf
      do inu=1,nf
          dae_lambdastar(inu) = lambdastar(inu)
          dae_kappaabs(inu) = kappaabs(inu)
          dae_epsilon(inu) = epsilon(inu)
          dae_meanint(inu) = meanint(inu)
          dae_nu(inu)      = nu(inu)
          dae_dnu(inu)     = dnu(inu)
      enddo
c
      if(nf.eq.1) then
c
c         Frequency-integrated mode
c
          dae_rhs(1) = meanint(1)-lambdastar(1)*
     %                  (epsilon(1)*sspi*tempold**4+
     %                   (1.d0-epsilon(1))*meanintold(1))
          dae_q      = q + kappaabs(1)*dae_rhs(1)
      else
c
c         Frequency-dependent mode
c
          dae_q      = q
          do inu=1,nf
              dae_rhs(inu) = meanint(inu)-lambdastar(inu)*
     %                  (epsilon(inu)*bplanck(temp,nu(inu))+
     %                   (1.d0-epsilon(inu))*meanintold(inu))
              dae_q      = dae_q + dnu(inu)*kappaabs(inu)*
     %                                dae_rhs(inu)
c################################
c   BUG HERE????  14-12-01
c################################
              write(*,*) 'ALI: BUG?? temp--->tempold!!'
          enddo
      endif
c
c     Check the behavior of the equation
c
      dum1 = dust_ali_equation(DUST_TEMP_MIN)
      dum2 = dust_ali_equation(DUST_TEMP_MAX)
      if(dum1.gt.0.d0) then
          iwarn = iwarn + 1
          temp = 0.d0
          return
      endif
      if(dum2.lt.0.d0) then
          write(*,*) 'Error in dust-ALI: '
          write(*,*) '  dust temperature > ',DUST_TEMP_MAX,' K'
          stop 1
      endif
c
c     Now use zbrent to solve the system
c
      temp     = zbrent(dust_ali_equation,DUST_TEMP_MIN,
     %                      DUST_TEMP_MAX,1.d-4*DUST_TEMP_MIN)
c
c     Return also the newly computed mean intensity
c
      do inu=1,nf
          meanint(inu) = dae_meanint(inu)
      enddo
c
      end



c     --------------------------------------------------------------
c                HELPER ROUTINE FOR SOLVING LOCAL ALI
c
c     This routine is called by the zbrent routine and returns the
c     thermal equilibrium equation for dust, using the ALI recipe
c     with a local operator. It includes (and solves for) scattering 
c     as well, and this is returned in the common block /dustali/. 
c     --------------------------------------------------------------
      function dust_ali_equation(temp)
      implicit none
      doubleprecision dust_ali_equation,temp
c
      common/dustalisolve/dae_meanint,dae_lambdastar,dae_epsilon,
     %       dae_kappaabs,dae_rhs,dae_q,dae_nu,dae_dnu
      doubleprecision dae_meanint(FRSIZE_FREQ)
      doubleprecision dae_lambdastar(FRSIZE_FREQ)
      doubleprecision dae_epsilon(FRSIZE_FREQ)
      doubleprecision dae_kappaabs(FRSIZE_FREQ)
      doubleprecision dae_rhs(FRSIZE_FREQ)
      doubleprecision dae_nu(FRSIZE_FREQ)
      doubleprecision dae_dnu(FRSIZE_FREQ)
      doubleprecision dae_q
      common/idustalisolve/dae_nf
      integer dae_nf
c
      doubleprecision bplanck,alifactor,result,sspi
      parameter(sspi=1.80491254763d-5)     ! ss / pi
      integer inu
c
      if(dae_nf.eq.1) then
c
c         Frequency-integrated mode
c
          alifactor      = 1.d0/(1.d0-dae_lambdastar(1)*
     %                           (1.d0-dae_epsilon(1)))
          dae_meanint(1) = alifactor * ( dae_rhs(1) +
     %                     dae_lambdastar(1)*dae_epsilon(1)*
     %                     sspi*temp**4 )
          result         = dae_q - dae_kappaabs(1) * ( 
     %                     (1.d0-dae_lambdastar(1)*
     %                     dae_epsilon(1))*sspi*temp**4 -
     %                     dae_lambdastar(1)*
     %                     (1.d0-dae_epsilon(1))*dae_meanint(1) )
      else
c
c         Frequency-dependent mode
c
          result = -dae_q
          do inu=1,dae_nf
              alifactor        = 1.d0/(1.d0-dae_lambdastar(inu)*
     %                                (1.d0-dae_epsilon(inu)))
              dae_meanint(inu) = alifactor * ( dae_rhs(inu) +
     %                           dae_lambdastar(inu)*dae_epsilon(inu)*
     %                           bplanck(temp,dae_nu(inu)) )
              result           = result + dae_dnu(inu)*
     %                           dae_kappaabs(inu)*((1.d0-
     %                           dae_lambdastar(inu)*dae_epsilon(inu))*
     %                           bplanck(temp,dae_nu(inu)) -
     %                           dae_lambdastar(inu)*(1.d0-
     %                           dae_epsilon(inu))*dae_meanint(inu))
          enddo
      endif
c
      dust_ali_equation = - result
      return
      end




c     --------------------------------------------------------------
c                   COPY THE DUST TEMP AND THE MEAN INT
c
c     This subroutine saves the dust temperature and the mean 
c     intensity of the radiation field to a special array. This is
c     necessary for the dust ALI algorithm. The ALI algorithm 
c     requires namely memory of the variables as they were before 
c     the LI step.
c     --------------------------------------------------------------
      subroutine copy_dust_temp_and_meanint()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
#include "common_source.h"
#include "common_radfield.h"
c
      integer itmax,it,ir,inu
c
c     The usual stuff
c
      itmax   = irsi_frsizey
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
#     endif
#     ifdef RADGRID_ONEDIM
      itmax = 1
#     endif
c
c     Loop over all space
c     
      do ir=1,irsi_frsizex
          do it=1,itmax
              do inu=1,freq_nr
                  dust_meanintold(inu,it,ir) = intmom_0(inu,it,ir)
              enddo
              dust_tempold(it,ir) = dust_temp(1,1,it,ir)
          enddo
      enddo
c     
      end




#endif /* IFDEF DUST_ALI */

#endif /* IFDEF DUST_RADICAL */



#else /* IFDEF INCLUDE_DUST */
      subroutine dust_dummy()
      write(*,*) 'hello'
      end
#endif /* IFDEF INCLUDE_DUST */






