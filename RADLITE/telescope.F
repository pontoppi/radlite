

c     =============================================================
c                      2-D FORMAL TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c           (MODULE: RENDER IMAGE OF OBJECT FOR EARTH OBSERVER)
c
c                            Leiden, Aug 1998
c                             C.P. Dullemond   
c
c     This module contains all the routines necessary to compute 
c     the image of the object as it is observerd on earth by a
c     telescope. One can compute the image itself, and one can also 
c     compute the spectrum only, by integrating the image for each
c     frequency. This module has its own set of (long) characteristics
c     and performs the transfer using the source terms as given by
c     the routines in the module source.F . For non-LTE systems, one
c     should first perform the transfer within the system, using the
c     module transfer.F and the routines therein. For LTE systems, 
c     the present module suffices to compute the image and spectrum
c     of the source.
c
c     The routines of interest to the user are:
c      - setup_rays_circular():
c           Set up the rays for the transfer. Must be called 
c           before use of the other routines.
c      - make_image_circular():
c           Produce the image as seen in a telescope.
c      - calc_freq_flux_observer(): 
c           Produce a spectrum of the source, in erg / s cm^2 Hz .
c      - write_image_circular():
c           Write the image to a file.
c      
c     =============================================================
c                           IMPORTANT NOTE:
c
c     The parameterization used here stems from a previous version
c     of the code, which used a weird parameterization. Instead of
c     completely rewriting the routine for the long characteristics,
c     I decided to leave it.
c
c       - This module:     z0,x0,thetainf,s
c       - Rest of program: z00,b,thetainf,s00
c
c     Relations:
c
c          b^2 = x0^2 + z0^2*(sin(thetainf))^2
c          s00 = s + z0*cos(thetainf)
c          z00 = z0*(sin(thetainf))^2
c
c     =============================================================

#define TELESCOPE_F_VERSION "Trans_16mar1999"
#include "main.h"

#ifdef LONGCHAR_EXTRA_POINTS
#define LONGCHAR_EX_PTS_RADIUS
#define LONGCHAR_EX_PTS_THETA
#endif


c     =============================================================
c              ROUTINES FOR SPECTRA AND FOR CIRCULAR IMAGES 
c
c     The spectrum is computed by first computing the image of 
c     the object at each frequency, and integrating over it. In 
c     order to be sure to treat the emission from small radii 
c     equally well as emission from large radii, we compute the
c     images _not_ with a rectangular grid (like a CCD), but instead
c     with a circular grid. The image consists of pixels arranged
c     in concentric circles around the center. There are more circles
c     in the center than in the outer parts, so as to guarantee
c     sufficient resolution to resolve the innermost parts of the
c     object as well as the outer parts.
c
c     The circular images can be used for image plotting as well,
c     but perhaps the rectangular images (routines with the extension
c     _rectang, below) may be more useful.
c
c     =============================================================


c     --------------------------------------------------------
c               READ THE INFO FILE FOR THE LINE SPECTRUM
c     --------------------------------------------------------
      subroutine read_telesc_linespec()
      implicit none
c
#include "common_telescope.h"
#include "common_lines.h"
#include "common_main.h"
c
      integer iformat,iformtel,iformobj,isizespecifier,gotit,i
      character*80 string
      logical fexists
c
      inquire(file='linespectrum.inp',EXIST=fexists)
      if(.not.fexists) then
          write(*,*) 'ERROR: Could not find linespectrum.inp'
          write(*,*) '       I need this file to know how to'
          write(*,*) '       produce the line spectrum!'
          stop 13
      endif            
      open(unit=1,file='linespectrum.inp',status='old')
      read(1,*) iformat
      linespec_style = 0
      if(iformat.eq.1) then
          read(1,*) linespec_style
      endif
      if((iformat.lt.0).or.(iformat.gt.1)) then
          write(*,*) 'While reading linespectrum.inp:'
          write(*,*) '    Unknown format version : ',iformat
          stop 13
      endif
      read(1,*) string
      read(1,*) iformtel

      if(iformtel.eq.2) then
          read(1,*) string
c          read(1,*) linespec_beamshape
c          read(1,*) linespec_dishdiameter    ! [meter] 
          read(1,*) linespec_vmax            ! [km/s]
          read(1,*) linespec_dv              ! [km/s]
          read(1,*) string
          linespec_dishdiameter = linespec_dishdiameter * 100.d0
      else
          write(*,*) 'While reading linespectrum.inp:'
          write(*,*) '    Unknown telescope format : ',iformtel
          stop 13
      endif      

c     Read information about the object itself
c
      read(1,*) iformobj

      if(iformobj.eq.1) then
c
c         Format type 1 (8 May 1999) : Only spectrum
c
          read(1,*) linespec_distance
          read(1,*) linespec_incang
          linespec_incang = linespec_incang * 0.0174532925199d0
          read(1,*) linespec_radvelo
          read(1,*) linespec_nlines
          read(1,*) linespec_ilinestart
          linespec_command = 0
      elseif(iformobj.eq.2) then
c
c         Format type 2 
c
c     Read first the input filename for the molecular data.
c     I moved it here from line.inp, since I at some point want
c     to use only one input file.

         read(1,*) string
         gotit=0
         do i=1,80-3
            if(string(i:i+3).eq.'.dat') then
               molec_filename=string(1:i+3)
               molec_name=string(1:i-1)
               molec_name_len=i-1
               gotit = 1
            endif
         enddo
         if(gotit.eq.0) then
            write(*,*) 'Could not make sense of string ',string
            write(*,*) 'while reading the line.dat file'
            write(*,*) 'File name must end with .dat'
            stop 13
         endif
c
c        Read molecular data file
c
         call read_species_lambda(molec_filename,nlev_orig)

c        if i_top_level is set to 0 then take all the levels into account

c         Find out what to compute:
c
c            0 = Line spectrum
c            1 = Frequency-integrated image
c            2 = 3-D position-velocity data cube
c     

         read(1,*) linespec_command
         read(1,*) linespec_distance
         read(1,*) linespec_incang
         linespec_incang = linespec_incang * 0.0174532925199d0
         read(1,*) linespec_radvelo
         read(1,*) linespec_nlines
         read(1,*) linespec_ilinestart
         if(linespec_command.eq.0) then
            telesc_command=5
c            if(linespec_nlines.gt.1) then
c               telesc_command = 5
c            endif
c           linespec_command=0
         elseif(linespec_command.eq.1) then
            write(*,*) 'While reading linespectrum.inp:'
            write(*,*) 'Freq-integr image not yet implemented'
            stop 1
         elseif(linespec_command.eq.2) then
            telesc_command=6
            read(1,*) imr_nx
            read(1,*) imr_ny
            read(1,*) isizespecifier ! How are the sizes specified?
            read(1,*) imr_szimx
            read(1,*) imr_szimy
            imr_szimx = imr_szimx * 0.5d0 ! don't ask questions ;-)
            imr_szimy = imr_szimy * 0.5d0
            imr_spx   = imr_szimx / (imr_nx*0.5d0)
            imr_spy   = imr_szimy / (imr_ny*0.5d0)
            read(1,*) imr_phioff
            read(1,*) imr_xoff
            read(1,*) imr_yoff
            read(1,*) imrec_addstar
            if(isizespecifier.eq.0) then
               write(*,*) 'While reading linespectrum.inp:'
               write(*,*) 'Size specifier in arcsec not yet ',
     %              'implemented'
            elseif(isizespecifier.eq.1) then
               isizespecifier=1
            else
               write(*,*) 'While reading linespectrum.inp:'
               write(*,*) 'Dont recognize size spec ',
     %              isizespecifier
               stop 13                  
            endif
         else
            write(*,*) 'While reading linespectrum.inp:'
            write(*,*) 'Dont recognize command no. ',linespec_command
            stop 13
         endif
      else
         write(*,*) 'While reading linespectrum.inp:'
         write(*,*) '    Unknown object format : ',iformobj
         stop 13
      endif
c     
c     Now check some stuff
c     
      if((linespec_ilinestart+linespec_nlines-1.gt.nlines).or.
     %     (linespec_ilinestart.lt.1)) then
         write(*,*) 'ERROR while reading linespectrum.inp:'
         write(*,*) '      Attempting to make spectrum of '
         write(*,*) '      non-existing line!'
         stop 13
      endif
c     
c     Now postprocess some stuff if necessary
c     
      if(linespec_vmax.le.0.d0) then
         linespec_vmax = line_rangewidth(linespec_ilinestart)
      endif
      close(1)
c     
      end
      


c     --------------------------------------------------------
c             DO TRANSFER ALONG ALL PHOTON TRAJECTORIES
c                 FOR ONE FREQUENCY AND MAKE IMAGE
c 
c                     (CIRCULAR IMAGE VERSION)
c
c     This routine performs transfer along the set of photon
c     trajectories (long characteristics) as specified by the 
c     setup_rays_circular() routine. The results are put into an
c     array in the common block /imageinf/. This array can be
c     regarded as an image of the object as seen by an observer
c     at a far distance. 
c     
c     Since the transfer program is meant to be able to deal
c     with large ranges in radius, the image should also reflect
c     this large dynamic range. This is done by making the image
c     not on a rectangular grid in arcsec x arcsec, but instead
c     choosing the impact parameter of the rays (b, in cm) and
c     the angle on a circle around the object's origin (phi, 
c     in rad) as coordinates. One obtains the intensity I as
c     a function of R and phi. The IDL software (with this package)
c     will then transform this into a rectangular image. This
c     will then appear as if the image consists of concentric
c     circles, which is in fact precisely what it is. This way
c     we can be sure to have a high image resolution close to
c     the origin, so as to resolve all spatial scales. This
c     is crucial for the computation of the spectrum, if at
c     different frequencies the photosphere is at vastly 
c     different radii.
c
c     In any case, this image can subsequently be integrated to 
c     yield a luminosity in erg / s cm^2 Hz. In fact, the function
c     calc_freq_flux_observer() does just that. It first calls this
c     routine make_freq_image_circular(), and then integrates to
c     obtain the luminosity. This routine make_freq_image_circular()
c     can therefore be regarded as a part of calc_freq_flux_observer(),
c     but it can also be used independent of it, if one wants
c     to obtain actual images of the object as seen by an
c     observer. 
c
c     HOW TO USE:
c      - First set up the long characteristics necessary for the
c        computation of the spectrum/image. This is done by calling
c        the routine setup_rays_circular(). See the header of that routine
c        for more information.
c      - Be sure to have all the sources ready. The 
c        functions calc_src_pnt() and calc_alp_pnt() in the
c        module source.F should be properly set. See the 
c        comments in these function headers for further 
c        information. See also the header of setup_rays_circular().
c      - The inu refers to the index of the frequency array. 
c        so specify for which frequency you like to find the image.
c      
c     RESULT:
c     This routine produces an array called imcir_int(), that
c     contains the image as seen for a far away observer.
c
c     WARNING:
c     For non-LTE problems, you must first use the module
c     transfer.F to iterate to convergence, before using this
c     routine to compute the observed image at infinity.
c
c     SMALL_MEMORY MODE:
c     If the SMALL_MEMORY mode is on, the array imcir_int will
c     not contain all the frequencies separately, but rather 1
c     single frequency bin: ifr=1. The inu may still be any 
c     value, but the resulting image will be stored in the
c     array at the ifr=1 position (since there are no other
c     valid values).
c     
c     BUGFIX: There used to be a bug here (until 19-04-99) 
c     that caused the central star to be incorrectly incorporated.
c
c     BUGFIX 08.04.05 / 24.03.05: The new possibility of 
c     flexible extra rays has a minor problem. It sometimes
c     causes the simulation to abort even when it does not
c     have to. This is now fixed.
c     --------------------------------------------------------
      subroutine make_image_circular()
      implicit none
c
      integer inu
      doubleprecision charintline,dummy
c
#include "common_grid.h"
#include "common_telescope.h"
c     (Include common_boundary.h: Related to BUGFIX 08.04.05 / 24.03.05)
#include "common_boundary.h"
#include "common_lines.h"
c
      integer iy,iray,ir,iphi,iins,ifr
      doubleprecision dr,bb,velo,aksmax_c
      dimension velo(1:FRSIZE_FREQ)
c
c     Dependent on the SMALL_MEMORY macro, we use ifr=inu
c     or ifr=1
c
c#ifdef SMALL_MEMORY
c      ifr = 1
c#else
c      ifr = inu
c#endif
c
c     Check if the ray parameters have indeed been set up
c
      if(rays_ready.ne.321) then
          write(*,*) 'ERROR, make_image_circular(): ',
     %               'Ray paramters not yet set'
          stop 13
      endif
c
c     Check if rays_nrphi is not too large
c
      if(rays_nrphi.gt.FRSIZE_INF_PHI) then
          write(*,*) 'Sorry, rays_nrphi is too large 
     %                (in charintline())'
          stop 13
      endif
c
c     Then copy some information
c
      imcir_nrphi = rays_nrphi
      imcir_nrr   = rays_nrr

c
c     Initialize the max and min velocity arrays (max possible velocity is c=1)
c
      do iray=0,rays_amount
         minvel(iray) = 1
         maxvel(iray) = -1
      enddo
      do inu=1,freq_nr
         velo(inu) = 
     %       line_dnu(inu,ispec_currentline)/line_nu0(ispec_currentline)         
      enddo
      do ir=1,imcir_nrr
         do iphi=1,rays_nrphi
            imcir_cont(iphi,ir) = 0d0
         enddo
      enddo
      aksmax_c = aksmax/2.99792458d5

c
c     If we want to monitor the creation of the image, then
c     we must open file.
c

#     ifdef MONITOR_CHARINT
         write(*,*) '!!!!! WRITING THE BIG FILE !!!!!!'
         open(unit=4,file='charmonitor_radlite.dat',
     %        form='unformatted',
     %        access='stream')
         write(4) TRMAXFIL,imcir_nrphi*imcir_nrr
         if(imcir_nrphi*imcir_nrr.gt.500000) then
            write(*,*) 'safety stop in telescope.F'
            stop 1
         endif
#     endif

c
c     Put the impact parameters into an array. I presume here that they
c     are arranged as planned: for each radius of the grid one impact
c     parameter, plus possible some extra rays in the central hole. 
c     If other impact parameter choices are made, then the
c     stuff below is inapplicable. This array is not very important
c     anyway, but it might be handy to know how the family of rays
c     are arranged (for which phi and impact parameter b).
c
c     First check if the total amount of impact parameters is not larger
c     than FRSIZE_INF_R
c
c     [BUGFIX 24.02.09] Added rp_nrref to the next line
      if(rp_dbdr*(irsi_frsizex-1)+rp_nrrayextra+rp_nrref+1.gt.
     %           FRSIZE_INF_R ) then
         write(*,*) 'The total amount of impact parameters is ',
     %    'too large. Recompile with bigger FRSIZE_INF_R'
         stop 13
      endif
c
c     Now make an array of the impact parameters. These always
c     represent the center of each pixel. 
c
      imcir_r(0) = 0.d0
      do ir=1,rays_nrr
         imcir_r(ir) = rays_r(ir)
      enddo
c
c     Now the edges of the pixels, so that we can figure out
c     the solid angle that each pixel subtends.
c
      imcir_ri(0) = 0.d0
c     [BUGFIX 24.02.09] Added rp_nrref to the next line
      do ir=1,rp_dbdr*(irsi_frsizex-1)+rp_nrrayextra+rp_nrref
          imcir_ri(ir) = 0.5*(imcir_r(ir)+imcir_r(ir-1))
      enddo
c     [BUGFIX 24.02.09] Added rp_nrref to the next line
      imcir_ri(rp_dbdr*(irsi_frsizex-1)+1+rp_nrrayextra+rp_nrref) = 
     %    rsi_x_c(irsi_frsizex,1)
c
c     [Added 28-07-02]
c     Shift imcir_ri boundaries if next boundary lies much closer
c     to the imcir_r point. This is necessary for instance at the
c     inner edge where grid refinement can be very strong while the
c     'extra' circles lie far apart. Without this fix the contribution
c     of the circle at the inner edge is grossly overdone. This was
c     probably the cause of the funny problems when taking the 
c     spectrum of a disk around a Herbig Ae/Be star at incl=0.
c
c     [BUGFIX 24.02.09] Added rp_nrref to the next line
      do ir=1,rp_dbdr*(irsi_frsizex-1)+rp_nrrayextra+rp_nrref-1
          if(imcir_r(ir)-imcir_ri(ir).gt.
     %          2*(imcir_ri(ir+1)-imcir_r(ir))) then
              imcir_ri(ir) = imcir_r(ir)-2*(imcir_ri(ir+1)-imcir_r(ir))
          endif
      enddo
c
c     [BUGFIX 08.04.05 / 24.03.05]
c     Even if all the rays now lie nicely in between the star and
c     the inner boundary, the cell-interface of the inner extra ray
c     may still be smaller than R*. Correct for this case.
c
      if(imcir_ri(1).lt.radbnd_rstar) then
          if(imcir_r(1).lt.radbnd_rstar) then
              write(*,*) 'INTERNAL ERROR IN RAY-SETUP...'
              stop
          endif
          imcir_ri(1)=radbnd_rstar
      endif
c
c     Loop over all rays. First the single central ray / beam.
c
c     The central beam is larger than the star, so that part of the
c     beam will pass the star, while part of it will be on the star.
c     This means that the charintline should do transfer from R=infty
c     to the camera, but at the closest approach it should include
c     the intersection of the star and the stellar light
c
      iray = 1
#     ifndef TELESCOPE_TRAJARR
      call make_trajectory_c(rp_x0(iray),
     %     rp_z0(iray),rp_theta0(iray),rp_s0(iray),0)
#     endif
      
      do inu=1,freq_nr
		  
#        ifdef MONITOR_CHARINT
            if (inu.eq.2.or.inu.eq.freq_nr/2) then
               write(4) inu
		       write(4) iphi
		       write(4) imcir_ri(ir)
            endif
#        endif
		  
		  
         imcir_int(inu,1,0) = 
     %        charintline(ispec_currentline,inu,iray,imcir_ri(1))
c
c     Write tau of image-center to special variable
c
         char_tau_center = char_tau
c
c     Copy this for all phi, because they are all the same
c
         do iphi=2,rays_nrphi
            imcir_int(inu,iphi,0) = imcir_int(inu,1,0)
         enddo
      enddo
c
c     Then the rest of the rays
c
      iray = 2
      do ir=1,imcir_nrr
          do iphi=1,rays_nrphi


             
#            ifndef TELESCOPE_TRAJARR
             call make_trajectory_c(rp_x0(iray),
     %            rp_z0(iray),rp_theta0(iray),rp_s0(iray),0)
#            endif

c            IMPORTANT: at the first frequency, we do not yet know the min and max velocities - hence, here
c            we have to integrate ALL the rays.
             inu = 1 
             imcir_int(inu,iphi,ir) = 
     %            charintline(ispec_currentline,inu,iray,-1.d0)
c            Mark the ray as 'integrated'
             imcir_cmask(inu,iphi,ir) = 1
c            and if the line does not cross the ray, then also mark it as continuum.
             if(velo(inu).gt.maxvel(iray)+2.*aksmax_c.or.
     %            velo(inu).lt.minvel(iray)-2.*aksmax_c) then
                imcir_cont(iphi,ir) = 
     %               imcir_int(inu,iphi,ir)
             endif
                
c            Now we can do the rest of the frequencies without wasting precious CPU cycles
             do inu=2,freq_nr
				 
#ifdef NONREDUNDANT
c               Now check: is the line velocity within the velocity boundaries of the ray? 
                if(velo(inu).le.maxvel(iray)+2.*aksmax_c.and.
     %               velo(inu).ge.minvel(iray)-2.*aksmax_c) then
	 
#                  ifdef MONITOR_CHARINT
	                  if (inu.eq.2.or.inu.eq.freq_nr/2) then
						  write(4) inu
	                      write(4) iphi
	 	                  write(4) imcir_ri(ir)
	 		          endif
#                  endif
	 
                   imcir_int(inu,iphi,ir) = 
     %                  charintline(ispec_currentline,inu,iray,-1.d0)
c                  set the mask to show that the ray has been integrated (for diagnostic reasons).     
                   imcir_cmask(inu,iphi,ir) = 1

c               If the line does not cross the ray in a high density region, 
c               then use the continuum (if it has already been calculated). 
c               Note that this currently assumes that the continuum is the same across the line.
c               If the user wishes this to be done better - turn off the NONREDUNDANT switch and accept
c               the much longer processing times. 
                else
                   if(imcir_cont(iphi,ir).ne.0) then
c                  The continuum has already been calculated for this ray, so just copy it
                      imcir_int(inu,iphi,ir) = imcir_cont(iphi,ir) 
c                     Although we do need it if we are saving the full integration
#                     ifdef MONITOR_CHARINT
					     if (inu.eq.2.or.inu.eq.freq_nr/2) then
					  		write(4) inu
					  		write(4) iphi
					  		write(4) imcir_ri(ir)
					  	    dummy = charintline(ispec_currentline,inu,iray,-1.d0)
					  	 endif
#                     endif
					  
                   else
c                     The continuum has not yet been calculated for the ray, so do it now
#                     ifdef MONITOR_CHARINT
					     if (inu.eq.2.or.inu.eq.freq_nr/2) then
					   	    write(4) inu
					   		write(4) iphi
					   		write(4) imcir_ri(ir)
					   	 endif
#                     endif

                      imcir_int(inu,iphi,ir) = 
     %                     charintline(ispec_currentline,inu,iray,-1.d0)
c                     And save it for later use
                      imcir_cont(iphi,ir) = 
     %                     imcir_int(inu,iphi,ir)                      
                   endif				   
                endif
                
#else
                imcir_int(inu,iphi,ir) = 
     %               charintline(ispec_currentline,inu,iray,-1.d0)
#endif
             enddo
             iray = iray + 1
          enddo
       enddo
       
c
c     Close the debugging monitor file
c
#     ifdef MONITOR_CHARINT
         close(4)
#     endif
c
c
      end



c     --------------------------------------------------------------
c            GENERATE THE COMPLETE SET OF PHOTON TRAJECTORIES
c                    FOR THE TELESCOPE AT INFINITY
c 
c                      (CIRCULAR IMAGE VERSION)
c
c     This routine sets up the set of long characteristics needed
c     for making the image of the object as seen by an observer
c     at infinity at a certain inclination angle. Here, the 
c     image pixels will be arranged in concrentric circles on 
c     the sky, centered around the origin of the object. This
c     allows the image to have high resolution close to the
c     origin, so as to guarantee that emission from all radii
c     is treated equally. 
c
c     If the appropriate macro is set, this routine calls 
c     setup_all_trajectories() to create the rays. If they are 
c     not created, then only the parameters of the rays are set 
c     up. The creation should then be done as soon as each ray 
c     is needed.
c
c     This function should always be called at the beginning of the
c     program, because this information is static (independent on
c     any variable that might change in the course of the program)
c     and without it, the formal radiative transfer for obtaining
c     the image and/or spectrum at a far-away observer cannot 
c     be performed.
c
c     HOW TO USE:
c      - irmin is the index of the minimal impact parameter. The
c        index refers to bimpact_min = rsi_x_c(irmin,1). Advice:
c        take is irmin=1
c      - irmax is similar to irmin, but then for the maximum. 
c        Advice: take it irmax=irsi_frsizex . 
c      - istep is the step in the rsi_x_c array it should take.
c        Advice: take it 1
c      - anginf is the angle of the observer at infinity with 
c        respect to the polar axis. It is measured in steradian.
c      - nrphiinf is the amount of bins in the observer phi-
c        angle it should take. The observer phi-angle is explained 
c        in the image that the observer would see with a strong
c        enough telescope: if this image is coordinized by an theta
c        in sterrad from the center of the object and an angle around
c        that center point, then this latter angle is phi_inf. 
c        Taking a reasonable resolution in phi is good, even if in 
c        the end one wishes only to find the luminosity of the object,
c        because the phi-resolution tells how many sample rays are
c        being traced through the source at each impact parameter 
c        from the source. Low resolution in phi would mean 
c        undersampling and therefore a bad estimate of the luminosity.
c
c     COMMENT: In the telescope routine I use an old way of specifying
c              the parameters of the rays: z0,x0,Theta_inf, instead
c              of the modern way: R_0,z_0,Theta_inf. 
c     COMMENT: There is still an old option: if you put anginf<-30.
c              then the routine computes a more-or-less complete 
c              set of rays for all directions, and all radii.
c              This is not useful in the present context, so it is
c              advised not to use this option.
c     COMMENT: This is one of the hackiest routines I've ever built....
c              Yuckie... There is a parameter irstep which enables
c              you to skip radial points, and there is this parametr
c              telesc_dbdr to ADD extra radial points in comparison
c              with the radial grid. Bluhh.
c     COMMENT: Added the possibility for 1--D spherical symmetry
c              here. (30-12-98)
c     COMMENT: irstep is now inactivated. (03-05-99)
c
c
c     Relation between (x_0,z_0,theta_inf,s) and (R_0,zz_0,theta_inf,ss)
c     is:
c
c       b^2  = x_0^2 + z_0^2 ( 1 - cos^2(theta_inf) )
c       zz_0 = z_0 sin^2(theta_inf)
c       ss   = s + z_0 cos(theta_inf)
c
c     23/02/2009: KMP updated the setup_rays_circular to that developed by
c     Kees for raytrace. The last version of RADLite using the old RADICAL 
c     ray setup is version 20/02/2009. 
c     --------------------------------------------------------------
      subroutine setup_rays_circular(irmin,irmax,irstep,anginf,
     %                      nrphiinf,nrext,telesc_dbdr,rstar,
     %                      imethod,nrref)
      implicit none
      doubleprecision anginf,rstar
      integer irmin,irmax,irstep,nrphiinf,nrext,telesc_dbdr
      integer imethod,nrref
c
#include "common_grid.h"
#include "common_telescope.h"
c
      integer ix,iy,it0,iy_amount,iys,ir,cnt_inspnt,iins
      integer iradius,nrrextra
      dimension iy_amount(0:FRSIZE_INF_PHI)
      doubleprecision theta_c,r_c,theta0,sinth0,costh0
      doubleprecision zhat0,xhat0,zh02,dum,dphi,phi,dr
      doubleprecision epsxyz,epsrrr,refdum1,refdum2,refdum3
c      logical dynamic
      dimension zhat0(0:FRSIZE_INF_PHI,0:1)       ! it0 = 1 anyway!
      dimension xhat0(0:FRSIZE_INF_PHI,0:1)       ! it0 = 1 anyway!
      parameter(epsxyz=1.0d2*TELESC_EPS)
      parameter(epsrrr=1.0d3*TELESC_EPS)
c
c     Test
c     
      if(imethod.lt.0) then
          write(*,*) 'Negative imethod not allowed.'
          stop
      endif
c
c     Default
c     [BUGFIX 05.07.07]
c
      nrrextra = nrext
      if(nrrextra.lt.0) nrrextra = -nrrextra
c
c     Decide on which way the radial rays are arranged
c     
      if(imethod.eq.0) then
c
c         The original method from RADICAL
c
          if(nrext.gt.0) then
              nrrextra = nrext
              imethod = -1
          elseif(nrext.lt.0) then
              nrrextra = -nrext
              imethod = -2
          else
              write(*,*) 'ERROR telescope.F: ',
     %                   'Must have non-zero nrrextra'
              stop
          endif
      endif
c
c     Do some checks
c
      if(irstep.ne.1) then
          write(*,*) 'ERROR: setup_rays_circular():'
          write(*,*) '       Argument irstep is obsolete'
          write(*,*) '       It MUST be 1. Other value not accepted'
          stop 13
      endif
      if(nrphiinf.gt.FRSIZE_INF_PHI) then
          write(*,*) 'ERROR Telescope: nrphiinf larger than ',
     %                'FRSIZE_IN_PHI'
          write(*,*) 'Recompile with larger FRSIZE_IN_PHI'
          stop 13
      endif
      if(irmax.le.irmin) then
          write(*,*) 'ERROR Telescope: irmax.le.irmin'
          stop 13
      endif
c     [BUGFIX 05.07.07] Added nrref to the next line
      if((irmax-irmin+1)*telesc_dbdr+nrrextra+nrref.gt.
     %   FRSIZE_INF_R) then
          write(*,*) 'ERROR while making circular image / spectrum'
          write(*,*) '       Too many radial points in circular '
          write(*,*) '       image requested. Either decrease the '
          write(*,*) '       value of dbdr or nrextra in the input'
          write(*,*) '       file, or recompile with larger value'
          write(*,*) '       of FRSIZE_INF_R'
          stop 13
      endif
#ifdef RADGRID_ONEDIM
      if(nrphiinf.ne.1) then
          write(*,*) 'ERROR: telescope.F/setup_rays_circular():'
          write(*,*) '       In RADYGRID_ONEDIM mode, the amount ',
     %     'of phi-angles of the image at infinity should,'
          write(*,*) '       by virtue of spherical symmetry, be 1.',
     %     'Not ',nrphiinf 
          write(*,*) '       Internal error in program! Aborting...'
          stop 13
      endif
#endif
c
      cnt_inspnt = 0 ! What is this for ???????
c
c     If the anginf is too close to the z-axis, some formulae
c     will diverge. Therefore it is required to be not smaller 
c     than 0.1. In the future this must be made smaller, but right
c     now I take it on the safe side.
c
      if(ABS(anginf).lt.1.d-1) then
          write(*,*) 'WARNING: In telescope.F/setup_rays_circular():'
          write(*,*) '         For an observer at inclination angle ',
     %           'smaller than 0.1 radians,'
          write(*,*) '         the current routine may make ',
     %     'round-off errors.'
          write(*,*) '         Therefore I have put the inclination ',
     %     'angle at Theta=0.1 rads. Hope you agree...'
          anginf = 0.1*ABS(anginf)/anginf
      endif
c
      if((nrrextra.gt.0).and.(irmin.ne.1)) then
         write(*,*) 'ERROR: In telescope.F/setup_rays_circular():'
         write(*,*) '       In You request extra rays at the ',
     %           'inner edge while you dont have irmin=1.'
         write(*,*) '       This is contradictory. Aborting...'
         stop 13
      endif
c
c
c     Make the rays only for one single impact parameter.
c     Later use this information for the loop over the impact
c     parameter: the r_inf in the (r_inf,phi_inf) image.
c
#ifdef RADGRID_ONEDIM
c
      stop
ccc
ccc     In case of 1--D spherically symmetric transfer, we take
ccc     z_0=0 and Theta0=pi/2
ccc
cc      it0        = 1
cc      theta0     = 1.5707963268d0
cc      sinth0     = 1.d0
cc      costh0     = 0.d0
cc      rays_nrphi = 1
cc      dphi       = 6.28318530718d0
cc      phi        = 0.d0
c
#else      
c
c     In case of 2--D transfer use the Theta0 as given by
c     the user
c     
      it0        = 1
      theta0     = anginf + 1.d-4
      sinth0     = sin(theta0)
      costh0     = cos(theta0)
      rays_nrphi = nrphiinf
      dphi       = 6.28318530718d0 / (1.d0*nrphiinf)
      phi        = 0.5 * dphi
#endif
c     
c     Now let's really make the set of rays for the image at
c     infinity... This consists of a series of concentric rings
c     of pixels around the center of the object. Each pixel of 
c     the image corresponds to one ray defined here. This
c     circular image has coordinates r (in [cm], is the radius
c     of each ring in the circular image) and phi (in [rad], is
c     the angle along each concentric ring). 
c
c     First find out how the z0 and x0 change as a function of 
c     the phi-angle of the circular image. z0 and x0 are also
c     proportional to the radius in the image, but here we
c     only determine the dependence on phi. Hence we call
c     these zhat0 and xhat0, to indicate that the radius is
c     divided out.
c     
      iys = 1
      do iys=1,nrphiinf
          zhat0(iys,it0) = - sin(phi) / sinth0
          zh02           = zhat0(iys,it0)*zhat0(iys,it0)
          dum            = 1.d0 - zh02*sinth0*sinth0
          dum            = dum + epsxyz
          if(dum.ge.0.d0) then
              if(cos(phi).gt.0.d0) then
                  xhat0(iys,it0) = sqrt( dum )
              else
                  xhat0(iys,it0) = -sqrt( dum )
              endif
          else
              write(*,*) 'ERROR in setup_rays_circular'
              stop 13
          endif
          phi = phi + dphi
      enddo
      iy_amount(it0) = nrphiinf
c     
c     Now determine the real x0 and z0, this time as a function
c     of both phi and r. The image has concentric rings, AND one
c     central circular pixel, corresponding to the line of sight
c     that goes right through the central star (so to speak).
c
c     First make this single radial ray right through the center itself
c     
      ir             = 1
      iradius        = 0
      rp_x0(ir)      = 0.d0
      rp_z0(ir)      = 0.d0
      rp_theta0(ir)  = theta0
      rp_iserth0(ir) = it0
      rp_iserrad(ir) = 0
      rp_s0(ir)      = 1.d30
      rays_r(iradius)= 0.d0
c     
      ir      = 2
      iradius = 1
c     
c     Next make the extra rays, which are the rays with impact parametes
c     smaller than the radius of the inner edge of the (R,Theta) spatial
c     grid of the object itself. 
c     
      if(imethod.lt.0) then
c
c         These are the original RADICAL-style ways to arrange these
c         rays. But be careful that there are drawbacks for near-face-on
c         spectra for very optically thick disks with vertical inner walls.
c
          do ix=1,nrrextra 
              if(iradius.gt.FRSIZE_INF_R) then
                  write(*,*) 'BUG in telescope.F/setup_rays_circular()'
                  write(*,*) '    iradius exceeds FRSIZE_INF_R'
                  write(*,*) '    (at extra inner rays)'
                  stop 155
              endif
              if(imethod.eq.-2) then
c                 
c                 The newer way, which is important for simulations of
c                 e.g. T Tauri stars and Brown Dwarfs [13.01.05]
c                 
                  if(rstar.gt.rsi_x_c(1,1)) stop 83991
                  r_c = ( ix * ( rsi_x_c(1,1) - rstar ) / 
     %                 (nrrextra+1.d0) ) + rstar

              elseif(imethod.eq.-1) then
c                 
c                 The old way
c                 
                  r_c = ix * rsi_x_c(1,1) / (nrrextra+1.d0)
              else
                  write(*,*) 'ERROR: Do not know imethod =',imethod
                  write(*,*) '   as a method for arraning rays...'
                  stop
              endif
              rays_r(iradius) = r_c
c     
c             Now a loop over theta_c (cell centers)
c             but only the ones that have a trajectory
c             at minimum radius here
c         
              if(iy_amount(it0).ge.1) then
                  do iys=1,iy_amount(it0)
                      rp_x0(ir)      = r_c * xhat0(iys,it0)
                      rp_z0(ir)      = r_c * zhat0(iys,it0)
                      rp_theta0(ir)  = theta0
                      rp_iserth0(ir) = it0
                      rp_iserrad(ir) = ix
                      rp_s0(ir)      = 1.d30
                      ir = ir + 1
                      if(ir.ge.MAXRAYNR) then
                          write(*,*) 'Exceeded maximum number ',
     %                         'of rays!!'
                          stop 13
                      endif
                  enddo
              endif
              iradius = iradius + 1
          enddo
      elseif(imethod.gt.0) then
c
c         The newer (only RAYTRACE) methods for arranging the rays.
c         The idea here is to make sure that the inner edge is better
c         sampled. Also, at some point, I want to introduce here the
c         method that allows the star to be sampled in a better way,
c         in case the inner edge of the disk/envelope is not so far 
c         away from the stellar surface.
c
          if(nrref.le.0) then
              write(*,*) 'ERROR: If new method for ray-arrangement'
              write(*,*) '  is chosen, then nrref must be set>0, too'
              stop
          endif
          refdum1 = 0.5d0
          refdum2 = 0.0d0
c         max scale height of inner rim
          refdum3 = 0.3d0 
          do ix=1,nrrextra+nrref 
              if(iradius.gt.FRSIZE_INF_R) then
                  write(*,*) 'BUG in telescope.F/setup_rays_',
     %                 'circular()'
                  write(*,*) '    iradius exceeds FRSIZE_INF_R'
                  write(*,*) '    (at extra inner rays)'
                  stop 155
              endif
              if(imethod.eq.1) then
c
c                 A method by which a minimal grid fineness is guaranteed,
c                 but also refinement near the inner edge is done.
c             
                  if(rstar.gt.rsi_x_c(1,1)) stop 91991
                  if(ix.le.nrrextra) then
                      r_c = ( (ix-1) * ( rsi_x_c(1,1) - rstar ) / 
     %                     nrrextra ) + rstar
                  else
                      refdum2 = refdum2 + refdum1
                      refdum1 = refdum1/2
                      r_c = ( ( nrrextra + refdum2 ) 
     %                     * ( rsi_x_c(1,1) - rstar ) / 
     %                     (nrrextra+1.d0) ) + rstar
                  endif
              else
                  write(*,*) 'ERROR: Do not know imethod =',imethod
                  write(*,*) '   as a method for arranging rays...'
                  stop
              endif
              rays_r(iradius) = r_c
c     
c             Now a loop over theta_c (cell centers)
c             but only the ones that have a trajectory
c             at minimum radius here
c             
              if(iy_amount(it0).ge.1) then
                  do iys=1,iy_amount(it0)
                      rp_x0(ir)      = r_c * xhat0(iys,it0)
                      rp_z0(ir)      = r_c * zhat0(iys,it0)
                      rp_theta0(ir)  = theta0
                      rp_iserth0(ir) = it0
                      rp_iserrad(ir) = ix
                      rp_s0(ir)      = 1.d30
                      ir = ir + 1
                      if(ir.ge.MAXRAYNR) then
                          write(*,*) 'Exceeded maximum number ',
     %                         'of rays!!'
                          stop 13
                      endif
                  enddo
              endif
              iradius = iradius + 1
          enddo
      else
          stop 33720
      endif
c     
c     Now make the usual rays, which are the rays with impact parameters
c     larger than the inner edge, but smaller than the outer edge of the
c     spatial object. These rays will contain most of the information of 
c     the image, once the telescope imaging routine has done it's job.
c     
      irmax = max(irmax,irsi_frsizex)
      do ix=irmin,irmax-1,irstep
          if(iradius.gt.FRSIZE_INF_R) then
              write(*,*) 'BUG in telescope.F/setup_rays_circular()'
              write(*,*) '    iradius exceeds FRSIZE_INF_R'
              write(*,*) '    (at usual rays)'
              stop 155
          endif
ccc          r_c = rsi_x_c(ix,1) * (1.d0 + 1.d-3)  ! (bug: 15-11-00)
          r_c = rsi_x_c(ix,1) * (1.d0 + epsrrr )
          rays_r(iradius) = r_c
c     
c         Now a loop over theta_c (cell centers)
c         but only the ones that have a trajectory
c         at minimum radius here
c         
          if(iy_amount(it0).ge.1) then
              do iys=1,iy_amount(it0)
                  rp_x0(ir)      = r_c * xhat0(iys,it0)
                  rp_z0(ir)      = r_c * zhat0(iys,it0)
                  rp_theta0(ir)  = theta0
                  rp_iserth0(ir) = it0
                  rp_iserrad(ir) = ix
                  rp_s0(ir)      = 1.d30
                  ir = ir + 1
                  if(ir.ge.MAXRAYNR) then
                      write(*,*) 'Exceeded maximum number ',
     %                    'of rays!!'
                      stop 13
                  endif
              enddo
          endif
          iradius = iradius + 1
c
c         If I am not at the last point, and if the telesc_dbdr
c         gt 1 then I will add a few points in between this
c         R_i and R_i+1.... And I will simply count the amount
c         of added points....
c
          if((ix.lt.irmax).and.(telesc_dbdr.gt.1)) then
             dr = ( rsi_x_c(ix+1,1) - rsi_x_c(ix,1) ) / 
     %                 ( 1.d0*telesc_dbdr )
             do iins = 1,telesc_dbdr-1
                if(iradius.gt.FRSIZE_INF_R) then
                     write(*,*) 'BUG in telescope.F/setup_',
     %                          'rays_circular()'
                     write(*,*) '    iradius exceeds FRSIZE_INF_R'
                     write(*,*) '    (at secondary rays)'
                     stop 155
                endif
                r_c = rsi_x_c(ix,1) + iins*dr
                rays_r(iradius) = r_c
                if(iy_amount(it0).ge.1) then
                   do iys=1,iy_amount(it0)
                      rp_x0(ir)      = r_c * xhat0(iys,it0)
                      rp_z0(ir)      = r_c * zhat0(iys,it0)
                      rp_theta0(ir)  = theta0
                      rp_iserth0(ir) = it0
                      rp_iserrad(ir) = ix
                      rp_s0(ir)      = 1.d30
                      ir = ir + 1
                      if(ir.ge.MAXRAYNR) then
                         write(*,*) 'Exceeded maximum ',
     %                        'number of rays!!'
                         stop 13
                      endif
                   enddo
                endif
                iradius = iradius + 1
             enddo
             cnt_inspnt = cnt_inspnt + (telesc_dbdr-1) ! ??????
          endif
      enddo
      rays_amount = ir - 1
      rays_nrr    = iradius - 1
c
c     The telescope imaging routine sometimes likes to know which
c     pixels belong where (although the routine for the transfer 
c     along the rays, charint(), does not care about this info, but
c     the routine for producing spectra from the images DOES need
c     to know this!). So store some information about the setup of 
c     the rays into the common block
c
      rp_nrrayextra = nrrextra
c     [BUGFIX 05.07.07] Added rp_nrref = nrref to the next line
      rp_nrref      = nrref
      rp_dbdr       = telesc_dbdr
      rp_irmin      = irmin
      rp_irmax      = irmax
      rp_irstep     = irstep
      rp_thetainf   = theta0    ! used to be: anginf
c
c     Put flag up to say that the ray parameters have been set up
c     This is an internal consistency check. 
c
      rays_ready = 321
c
c     Now render the rays... This takes some CPU time...
c     Only do this if the TELESCOPE_TRAJARRAY option is set, meaning
c     that we wish to pre-compute all trajectories. For high-resolution
c     images this seems not so good, because it is too memory expensive.
c     So per default for the telescope it is better to switch off this
c     option (i.e. do not define TELESCOPE_TRAJARRAY).
c
#ifdef TELESCOPE_TRAJARRAY
      write(*,*) 'Rendering trajectories'
      call setup_all_trajectories(1)
#endif
c
c     Call the reindexation routine, necessary for the interpolations
c     
c     WARNING: This is actually a hack, because iangset=1 has nothing to
c              do with the camera functioning. But the reindexation
c              is also necessary for the theta mirroring, so that's
c              why we need to do it anyway.
c
      call make_index(1)
c
c     Reset imethod if negative
c
      if(imethod.lt.0) then 
          imethod=0
      endif

      return
      end


c     --------------------------------------------------------
c           WRITE THE IMAGE AT INFINITY TO A FILE
c
c     Write the image to a file. This is done for all frequencies,
c     if they are all in memory. 
c
c     SMALL_MEMORY MODE:
c       If the SMALL_MEMORY mode is on, the array imcir_int will
c       not contain all the frequencies separately, but rather 1
c       single frequency bin: ifr=1. The inu may still be any 
c       value, but the resulting image will be stored in the
c       array at the ifr=1 position (since there are no other
c       valid values).
c     --------------------------------------------------------
      subroutine write_image_circular(ifreq)
      implicit none
c
      integer ifreq
c
#include "common_grid.h"
#include "common_telescope.h"
c
      integer ir,ip,inu,inu0,inu1
c
c     Dependent on the macro SMALL_MEMORY we save all frequencies
c     or just a single snapshot.
c     
#     ifdef SMALL_MEMORY
      inu0 = 1
      inu1 = 1
#     else
      if(ifreq.gt.0) then
          inu0 = ifreq
          inu1 = ifreq
      else
          inu0 = 1
          inu1 = freq_nr
      endif
#     endif      
c
c     First write a file containing the set of impact parameters
c     of the telescope image. Normally this is equal the the R_i
c     of the computational grid, but sometimes you want extra
c     resoltution or extra rays towards the center.
c
      open(unit=1,file='iminf_b.dat',status='unknown')
      write(1,*) imcir_nrr
      do ir=1,imcir_nrr
          write(1,*) imcir_r(ir)
      enddo
      write(1,*) ' '
      do ir=1,imcir_nrr+1
          write(1,*) imcir_ri(ir)
      enddo
      close(1)
c
c     Open the file and write a header. The '0' is for the 
c     possible mirror symmetry in phi, which is not yet 
c     implemented.
c
      open(unit=1,file='imageinf.dat',status='unknown')
      write(1,*) imcir_nrr,imcir_nrphi,0,inu1+1-inu0
c
c     Write the main stuff
c
      write(1,*) ' '
      do inu=inu0,inu1
          do ir=0,imcir_nrr
              write(1,*) (imcir_int(inu,ip,ir),ip=1,imcir_nrphi)
          enddo
          write(1,*) ' '
      enddo      
c
      close(1)
      end



c     --------------------------------------------------------
c             INTEGRATE THE CIRCULAR IMAGE AT INFINITY 
c
c     This function computes the actual flux as seen
c     by an observer at a certain distance from the object.
c     For line transfer it can also convolve the circular image
c     with a beam pattern, before integration, and thus yield
c     the beam--averaged intensity.
c
c     This function automatically calls the routine for producing
c     the circular image: make_freq_image_circular(inu).
c     Then it integrates this image over the solid 
c     angle on the sky to obtain the luminosity in ergs/s.cm^2 
c     at the location of the observer. Note that the image
c     produced by make_freq_image_circular() is an image as
c     a function of (R,phi), i.e. circular coordinates on the
c     sky, centered around the origin of the object. This is
c     to allow accurate treatment of radiation emerging from
c     all radii, small or large alike. 
c
c     ARGUMENTS:
c       inu        = Index of frequency bin
c       distance   = Distance of source to earth in cm (not
c                    used when convolution is performed).
c       iprecalc   = 0 : Image at this freq NOT yet computed, 
c                        so calc_freq_flux_observer() will call it.
c                    1 : Image at this freq already computed and
c                        in memory. So need not recompute this.
c       iconvolve  = 0 : Produce the integral of the image, i.e. flux
c                    1 : Convolve the intensity with a function given
c                        in the array imcir_convfunc(), and return the
c                        average intensity.
c
c     WARNING:
c       For non-LTE problems, you must first use the module
c       transfer.F to iterate to convergence, before using this
c       routine to compute the observed spectrum at infinity.
c       Otherwise the source function and opacity are not yet
c       at their local statistical equilibrium value.
c
c     BUGFIX: 
c       The wrong solid angle was used for the flux...
c       And besides, for all-except-one bin, the solid angle
c       was not even multiplied!!
c       This is now fixed (19-04-99)
c
c     ADDED: conv is now returned. (06.04.06)
c     --------------------------------------------------------
      subroutine calc_freq_flux_observer(distance,iprecalc,
     %                                 iconvolve,rapert,conv)
      implicit none
      doubleprecision distance,rapert
      doubleprecision conv,velo
      integer inu,iprecalc,iconvolve
c
c#include "common_units.h"
#include "common_telescope.h"
#include "common_grid.h"
#include "common_lines.h"
c
      doubleprecision slum,surf,dslum
      integer ir,ip,ifr
c
c     First perform transfer along all the rays (if the circular
c     image has not yet been made before this call)
c
      if(iprecalc.ne.1) then
          call make_image_circular
      endif

c
c
c    In case one would like to save the circular image

      do inu=1,freq_nr
         
#ifdef SAVE_IMCIR
         if(inu.eq.1)then
            write(9,*) line_nu0(ispec_currentline)
            write(9,*) imcir_nrphi, imcir_nrr
            do ir=1,imcir_nrr+1
               write(9,'(E12.6)') imcir_ri(ir)
            enddo
            do ir=1,imcir_nrr+1
               write(9,'(E12.6)') imcir_r(ir)
            enddo
         endif
         write(9,*) '       '
         velo = 
     %      line_dnu(inu,ispec_currentline)/line_nu0(ispec_currentline)
         write(9,*) velo*2.99792458d5
         write(9,*) '       '   
         write(9,*) imcir_int(inu,1,0)
         write(9,*) '       '   
         do ip=1,imcir_nrphi
            do ir=1,imcir_nrr
               write(9,'(E10.4,2X,i1)') imcir_int(inu,ip,ir), 
     %              imcir_cmask(inu,ip,ir)
            enddo
         enddo
         
#endif

c
c     Next add all the contributions up. This is done by first
c     order integration. I integrate over the impact parameter B:
c
c                    1   /   /
c              L = ----  | B | I dphi dB      [ erg / s cm^2 Hz]
c                  D^2   /   /                                    
c
c     where D=distance. 
c
c     So, start the loop in frequency. Reset the flux and the convolution
c     function integral. 
c
         slum  = 0.d0
         conv  = 0.d0
c
c     First the purely radial ray
c
         surf  = 3.14159265359d0 * (imcir_ri(1)**2)
         dslum = surf * imcir_int(inu,1,0)
         slum  = slum + dslum
c
c     Then the rest of the rays, with non-nil impact parameter
c      
         do ir=1,imcir_nrr
c
c         Check the aperture
c
            if(imcir_ri(ir).lt.rapert) then
c
c             Calculate surface between circle of radius imcir_ri(ir+1) 
c             and circle of radius imcir_ri(ir)
c            
               if(imcir_ri(ir+1).lt.rapert) then
                  surf = 3.14159265359d0 *  
     %                 (imcir_ri(ir+1)**2-imcir_ri(ir)**2) 
               else
                  surf = 3.14159265359d0 *  
     %                 (rapert**2-imcir_ri(ir)**2) 
               endif
c
c             Now compute the average of I over phi
c
               dslum = 0.d0
               do ip=1,imcir_nrphi
                  dslum = dslum + imcir_int(inu,ip,ir)
               enddo
               dslum = dslum / (1.d0*imcir_nrphi)
               dslum = dslum * surf

c
c             Next add to slum
c
               slum = slum + dslum
c
            endif
         enddo
         
         spec_flux_observer(inu) = slum / (distance**2)

      enddo
      return
      end


c     =============================================================
c                      ROUTINES FOR LINE SPECTRA
c     
c     These routines are specially for line transfer. 
c     =============================================================

c     --------------------------------------------------------
c            HIGH-RESOLUTION SPECTRUM FOR A SINGLE LINE
c
c     This routine will compute the temperature profile of a
c     line, for a given beam centered around the origin of the
c     object under consideration. It can also compute absolute
c     flux.
c
c     ARGUMENTS:
c       iline      = Index of the line
c       distance   = Distance of source to earth in cm
c       beamsize   < 0 : Compute total flux instead of temperature
c                  > 0 = Beam size in " (arcsec) over which to convolve.
c       passband   = Width of passband in km/s from line center
c       nfr        = Number of frequency points in passband (left-right)
c
c
c
c     --------------------------------------------------------
      subroutine calc_line_spectrum(iline,distance,passband,nfr)
      implicit none
c
      doubleprecision distance,passband
      doubleprecision conv
      doubleprecision telesc_airy_beam,telesc_beam_size
      integer iline,nfr,itr_freq
      common/calclinespec/idone_print_rminmax
      integer idone_print_rminmax
      character*80 filename
      character*12 ch
c
c#include "common_units.h"
#include "common_grid.h"
#include "common_telescope.h"
#include "common_source.h"
#include "common_lines.h"
#include "common_boundary.h"
c
      doubleprecision bbeam,beamsize
      integer inu,ir,iconvolve
c
      call telescope_check_safety_numbers()
c
c     Consistency check
c
#     ifndef SMALL_MEMORY
      write(*,*) 'ERROR: telescope.F/calc_line_spectrum()'
      write(*,*) '       SMALL_MEMORY must be set!'
      stop 13
#     endif
#     ifndef INCLUDE_LINES
      write(*,*) 'ERROR: telescope.F/calc_line_spectrum()'
      write(*,*) '       INCLUDE_LINES must be set!'
      stop 13
#     endif
#     ifndef LINE_VELOCITIES
      write(*,*) 'ERROR: telescope.F/calc_line_spectrum()'
      write(*,*) '       LINE_VELOCITIES must be set!'
      stop 13
#     endif
#     ifdef LINE_PROFILE_ARRAY
      write(*,*) 'ERROR: telescope.F/calc_line_spectrum()'
      write(*,*) '       LINE_PROFILE_ARRAY not allowed!'
      stop 13      
#     endif
      if(iradproc_line.eq.0) then
          write(*,*) 'ERROR: telescope.F/calc_line_spectrum()'
          write(*,*) '       iradproc_line must be set!'
          stop 13
      endif
      if(nfr.gt.SZ_LINEPROFILE) then
          write(*,*) 'ERROR: telescope.F/calc_line_spectrum()'
          write(*,*) '       SZ_LINEPROFILE too small'
          write(*,*) '       Requested ',nfr,' frequencies per',
     %      ' line, while'
          write(*,*) '       only ',SZ_LINEPROFILE,' are allowed.'
          write(*,*) '       Recompile with SZ_LINEPROFILE at ',
     %      'least ',nfr
          stop 13
      endif
      if(nfr.gt.FRSIZE_FREQ) then
          write(*,*) 'ERROR: telescope.F/calc_line_spectrum()'
          write(*,*) '       FRSIZE_FREQ too small'
          write(*,*) '       Requested ',nfr,' frequencies per',
     %      ' line, while'
          write(*,*) '       only ',FRSIZE_FREQ,' are allowed.'
          write(*,*) '       Recompile with FRSIZE_FREQ at ',
     %      'least ',nfr
          stop 13
      endif
c
c     Initialize maser detector
c
      maserflag = 0
c
c     Tell the rest of the telescope module that we're working
c     on line nr. iline, and all the other info
c
      ispec_currentline = iline
c
c
c     Set up the passband
c
      itr_freq = 1
      call line_setup_passband(iline,passband,nfr,itr_freq)
c
c     Prepare the star spectrum for this passband
c
c     NOTE: This routine is not used for the Non-LTE stuff
c           (for that, see routine line_init_starbc())
c
      call line_init_passband_starbc(nfr,1)
c
c     Prepare the outer boundary condition for this passband
c
c     NOTE: This routine is not used for the Non-LTE stuff
c           (for that, see routine line_init_outerbc())
c
      if(iradbnd_out_itype.eq.3) then
          call line_init_passband_outerbc(nfr,1)
      else
          do inu=1,nfr
              radbnd_interstellfield(inu) = 0.d0
          enddo
      endif
c
c     Prepare the transitions...
c
      call global_prepare_transitions()
c
c     Copy the frequencies into the arrays belonging to the spectrum
c     common block
c
      ispec_nrfreq = freq_nr
      do inu=1,freq_nr
          spec_freq(inu) = freq_nu(inu)
      enddo
c
c     Always switch off convolution (KMP 2008/11/22)
c
      iconvolve = 0
      linespec_flag1 = 1
c
c     Now compute the luminosity for each frequency and copy into the
c     spectrum array
c

c
c     If the user wants to output the circular image, prepare the file here. 
#ifdef SAVE_IMCIR 
      write(ch,'(i2)') iline
#ifdef COMPILER_INTEL
      if(iline.lt.10) then
         filename = 'lineposvelcirc_'//
     %        molec_name(1:molec_name_len)//
     %        '_'//ch(2:2)//'.dat'
      elseif(iline.lt.100) then
         filename = 'lineposvelcirc_'//
     %        molec_name(1:molec_name_len)//
     %        '_'//ch(1:2)//'.dat'
      endif
#else
      if(iline.lt.10) then
         filename = 'lineposvelcirc_'//
     %        molec_name(1:molec_name_len)//
     %        '_'//ch(2:2)//'.dat'
      elseif(iline.lt.100) then
         filename = 'lineposvelcirc_'//
     %        molec_name(1:molec_name_len)//
     %        '_'//ch(1:2)//'.dat'
      endif
#endif
      open(unit=9,file=filename,
     %     status='unknown')
      write(9,*) nfr
#endif

c
c         Generate the flux
c
      call calc_freq_flux_observer(distance,0,iconvolve,
     %     1d99,conv)
c
c         ADDED 06.04.06
c
c      linespec_beamradian(inu) = conv
c
c         If close to line center, then save the center-beam 
c         optical depth
c
c          if(inu.gt.1) then
c              if(line_dnu(inu,iline)*line_dnu(inu-1,iline)
c     %                            .le.0.d0) then
c                  char_tau_linecenter(iline) = char_tau_center
c              endif
c          endif
c

#ifdef SAVE_IMCIR 
      close(9)
#endif

c
      IF (maserflag.eq.1) THEN 
         write(*,'(A35,I10)') 'WARNING: Masing detected!'
         write(*,*) 'Will only warn once for this transition.'
      ENDIF

      end


c     --------------------------------------------------------
c           WRITE THE HEADER OF THE LINE SPECTRUM FILE
c     
c     For use in combination with write_line_spectrum() below.
c     --------------------------------------------------------
      subroutine header_line_spectrum(nlinespec,nfrmax,distance,
     %                        radvelo,anginf,iformout)
      implicit none
c
      integer nlinespec,nfrmax,iformout
      doubleprecision distance,radvelo,anginf
c
#include "common_telescope.h"
#include "common_lines.h"
c
      linespec_filename = 'linespectrum_'//
     %               molec_name(1:molec_name_len)//'.dat'
      open(unit=1,file=linespec_filename,status='unknown')
      write(1,'(i2)') 1                                 ! File format version
      write(1,'(i2)') iformout                          ! 0: spectrum=temperatures
      write(1,'(a80)') molec_name                       ! Name of the molecule
      write(1,'(a80)') molec_filename                   ! Name of the molec data file
      write(1,'(i10)') nlinespec,nfrmax 
      write(1,'(e12.4,e12.4,f7.3)') distance/3.08572d18,! Distance in [pc]
     %           radvelo,                               ! Radial velocity [km/s]
     %           anginf*57.2957795132                   ! Inclination angle in [degrees]
      close(1)
      end




c     --------------------------------------------------------
c                      WRITE LINE SPECTRUM
c
c     This routine writes the arrays spec_flux_observer and
c     spec_freq to a file, in such a way that it is easy to 
c     interpret as a line spectrum. 
c
c       iformout = 0    Output spec is `temperature'
c       iformout = 1    Output spec is F_nu
c
c      NOTE: This routine APPENDS the spectrum of a line to 
c            the file linespectrum_<molecname>.dat. Use 
c            header_line_spectrum() to write the header.
c     --------------------------------------------------------
      subroutine write_line_spectrum(iformout)
      implicit none
c
      integer iformout
c
c#include "common_units.h"
#include "common_grid.h"
#include "common_telescope.h"
#include "common_lines.h"
c
      doubleprecision velo,temp
      integer inu
c
c     Check if the spectrum has been computed
c
      if(ispec_nrfreq.eq.0) then
          write(*,*) 'ERROR: telescope.F/write_spectrum()'
          write(*,*) '       Attempt to write spectrum, but ',
     %          'spectrum not yet computed!'
          stop 13
      endif
c
c     Write file for frequencies and luminosities
c
      open(unit=1,file=linespec_filename,status='old',
     %                OPENACC='append')
      write(1,*) 
      write(1,'(i5,i5)') lev_up(ispec_currentline),
     %     lev_down(ispec_currentline)
      write(1,'(e14.9)') linefreq(ispec_currentline)
      write(1,'(e10.5)') linespec_beamsize
      write(1,'(i5)') ispec_nrfreq
      write(1,*) ' '
c
c     Write the spectrum
c
      if(iformout.eq.0) then
c
c         The output is: velocity , temperature. The temperature is
c         computed using Rayleigh-Jeans Law
c
          if(linespec_flag1.eq.2) then
              do inu=ispec_nrfreq,1,-1
                  velo = spec_freq(inu) - linefreq(ispec_currentline)
                  velo = - 2.99792458d5 * velo / 
     %                 linefreq(ispec_currentline)
                  velo = velo + linespec_radvelo
                  temp = 3.25465503368d36 / 
     %                   linefreq(ispec_currentline)**2
                  temp = temp * spec_flux_observer(inu)
                  if(temp.lt.1.d-50) then
                      temp = 0.d0
                  endif
                  write(1,10) velo,temp
              enddo
          else
              write(*,*) 'Outputting temperature instead of Fnu only'
              write(*,*) '   possible in convolution mode (put beam'
              write(*,*) '   size or dish size >0 in linespectrum.inp'
              stop 52987
          endif
      elseif(iformout.eq.1) then
c
c         Output is normal flux
c
          if(linespec_flag1.eq.1) then
c
c             The spec_flux_observer is the actual F_nu, so just dump
c             this.
c
              do inu=ispec_nrfreq,1,-1
                  velo = spec_freq(inu) - linefreq(ispec_currentline)
                  velo = - 2.99792458d5 * velo / linefreq(ispec_currentline)
                  velo = velo + linespec_radvelo
                  write(1,10) velo,spec_flux_observer(inu)
              enddo
          elseif(linespec_flag1.eq.2) then
c
c             The spec_flux_observer is an average intensity over the
c             beam. So to get flux output multiply by the beam surface.
c
c             ADDED: 06.04.06
c
              do inu=ispec_nrfreq,1,-1
                  velo = spec_freq(inu) - linefreq(ispec_currentline)
                  velo = - 2.9979d5 * velo / linefreq(ispec_currentline)
                  velo = velo + linespec_radvelo
                  write(1,10) velo,spec_flux_observer(inu)*
     %               linespec_beamradian(inu)
              enddo
          else
              write(*,*) 'ERROR: Dont know linespec_flat1 = ',
     %                    linespec_flag1
              stop 52986
          endif
      else
          stop 18273
      endif
c   10 format(E21.13,1X,E21.13)
   10 format(E13.6,1X,E13.6)
      close(1)
      end




c     --------------------------------------------------------
c            3-D RECTANGULAR POSITION-VELOCITY DATA CUBE
c
c     This routine is specially for line transfer. It will 
c     first set up everything necessary for line transfer.
c     Then it will compute images of the object at all the
c     frequency points in the line, and append each of them
c     to the file lineposvel_*.dat (*=iline). The output
c     will be in Kelvin (brightness temperature). The
c     optical depth through the object will be the second
c     column of the file.
c
c     ARGUMENTS:
c       iline      = Index of the line
c       distance   = Distance of source to earth in cm
c       passband   = Width of passband in km/s from line center
c       nfr        = Number of frequency points in passband (left-right)
c
c     --------------------------------------------------------
      subroutine calc_write_line_posvel(iline,distance,
     %                              passband,nfr)
      implicit none
c
      doubleprecision distance,passband,velo(SZ_LINEPROFILE)
      integer iline,nfr,itr_freq,ifreq,ix,iy
      common/calclinespec/idone_print_rminmax
      integer idone_print_rminmax
      character*80 filename
      character*12 ch
      doubleprecision temp
c
#ifndef RADGRID_ONEDIM
c
c#include "common_units.h"
#include "common_grid.h"
#include "common_telescope.h"
#include "common_source.h"
#include "common_lines.h"
c
      integer inu,ir,iconvolve
c
c     Consistency checks
c
      call telescope_check_safety_numbers()
c
#     ifdef TELESCOPE_NOIMTAU
      write(*,*) 'ERROR: telescope.F/calc_write_line_posvel()'
      write(*,*) '       Cannot define TELESCOPE_NOIMTAU!'
      stop 13
#     endif
#     ifndef LONGCHAR_EXTRA_POINTS
      write(*,*) 'WARNING: While making pos-vel diagram:'
      write(*,*) '         LONGCHAR_EXTRA_POINTS must be defined'
      write(*,*) '         to avoid grainy images'
      write(*,*) '         STOP disabled, Klaus 6/4/07'
c      stop 1
#     endif
#     ifndef SMALL_MEMORY
      write(*,*) 'ERROR: telescope.F/calc_write_line_posvel()'
      write(*,*) '       SMALL_MEMORY must be set!'
      stop 13
#     endif
#     ifndef INCLUDE_LINES
      write(*,*) 'ERROR: telescope.F/calc_write_line_posvel()'
      write(*,*) '       INCLUDE_LINES must be set!'
      stop 13
#     endif
#     ifndef LINE_VELOCITIES
      write(*,*) 'ERROR: telescope.F/calc_write_line_posvel()'
      write(*,*) '       LINE_VELOCITIES must be set!'
      stop 13
#     endif
#     ifdef LINE_PROFILE_ARRAY
      write(*,*) 'ERROR: telescope.F/calc_write_line_posvel()'
      write(*,*) '       LINE_PROFILE_ARRAY not allowed!'
      stop 13      
#     endif
      if(iradproc_line.eq.0) then
          write(*,*) 'ERROR: telescope.F/calc_write_line_posvel()'
          write(*,*) '       iradproc_line must be set!'
          stop 13
      endif
      if(nfr.gt.SZ_LINEPROFILE) then
          write(*,*) 'ERROR: telescope.F/calc_write_line_posvel()'
          write(*,*) '       SZ_LINEPROFILE too small'
          write(*,*) '       Requested ',nfr,' frequencies per',
     %      ' line, while'
          write(*,*) '       only ',SZ_LINEPROFILE,' are allowed.'
          write(*,*) '       Recompile with SZ_LINEPROFILE at ',
     %      'least ',nfr
          stop 13
      endif
      if(nfr.gt.FRSIZE_FREQ) then
          write(*,*) 'ERROR: telescope.F/calc_write_line_posvel()'
          write(*,*) '       FRSIZE_FREQ too small'
          write(*,*) '       Requested ',nfr,' frequencies per',
     %      ' line, while'
          write(*,*) '       only ',FRSIZE_FREQ,' are allowed.'
          write(*,*) '       Recompile with FRSIZE_FREQ at ',
     %      'least ',nfr
          stop 13
      endif
c
c     Initalize maser detector
c
      maserflag = 0
c
c     Tell the rest of the telescope module that we're working
c     on line nr. iline, and all the other info
c
      ispec_currentline = iline
c
c     Set up the passband
c
      itr_freq = 1
      call line_setup_passband(iline,passband,nfr,itr_freq)
c
c     Prepare the transitions...
c
      call global_prepare_transitions()
c
c     Copy the frequencies into the arrays belonging to the spectrum
c     common block
c
      ispec_nrfreq = freq_nr
      do inu=1,freq_nr
          spec_freq(inu) = freq_nu(inu)
      enddo
c
c     Open file and write header
c
      write(ch,'(i2)') iline
#ifdef COMPILER_INTEL
      if(iline.lt.10) then
          filename = 'lineposvel_'//
     %               molec_name(1:molec_name_len)//
     %               '_'//ch(2:2)//'.dat'
      elseif(iline.lt.100) then
          filename = 'lineposvel_'//
     %               molec_name(1:molec_name_len)//
     %               '_'//ch(1:2)//'.dat'
      else
          stop 177
      endif
#else
      if(iline.lt.10) then
          filename = 'lineposvel_'//
     %               molec_name(1:molec_name_len)//
     %               '_'//ch(2:2)//'.dat'
      elseif(iline.lt.100) then
          filename = 'lineposvel_'//
     %               molec_name(1:molec_name_len)//
     %               '_'//ch(2:3)//'.dat'
      else
          stop 177
      endif
#endif
      open(unit=1,file=filename,status='unknown')
      write(1,*) 
      write(1,'(I2)') 1                        ! File format version
      write(1,'(A80)') molec_name              ! Name of the molecule
      write(1,'(A80)') molec_filename          ! Name of the molec data file
      write(1,'(ES12.4,ES12.4,F7.3)') 
     %     linespec_distance/3.08572d18,       ! Distance in [pc]
     %           linespec_radvelo,             ! Radial velocity [km/s]
     %           linespec_incang*57.2957795132 ! Inclination angle in [degrees]
      write(1,'(I5,I5)') lev_up(ispec_currentline),    ! Up-level of this line
     %           lev_down(ispec_currentline)   ! Down-level of this line
      write(1,'(ES12.4)') linefreq(ispec_currentline)   ! Frequency at line center
      write(1,'(I5)') ispec_nrfreq                  ! Nr of frequencies in this line
      write(1,'(I5,I5,ES12.4,ES12.4,F7.3,ES12.4,ES12.4)') 
     %           imr_nx,                       ! Nr of pixels in x-direction
     %           imr_ny,                       ! Nr of pixels in y-direction
     %           imr_spx,                      ! Pixel size in x-direction
     %           imr_spy,                      ! Pixel size in y-direction
     %           imr_phioff,                   ! Rotation angle of image
     %           imr_xoff,                     ! X-offset of image
     %           imr_yoff                      ! Y-offset of image
      write(1,*) ' '      
c
c     Write the velocities
c
      do inu=1,freq_nr
          velo(inu) = spec_freq(inu) - linefreq(iline)
          velo(inu) = - 2.99792458d5 * velo(inu) / linefreq(iline)
          velo(inu) = velo(inu) + linespec_radvelo
      enddo
      write(1,'(E10.3,1X)') (velo(inu),inu=1,freq_nr)
      write(1,*) ' '      
c
c     Now make the image
c     
      call make_freq_image_rectang
c     
c     And write the image and the line-of-sight 
c     optical depth
c     
      do inu=1,freq_nr
         do iy=1,imrec_nry
            do ix=1,imrec_nrx
               temp = 3.25465503368d36 / linefreq(iline)**2
               temp = temp * imrec_int(inu,ix,iy)
               write(1,'(ES12.5,1X,ES12.5,1X)') 
     %              temp,imrec_tau(inu,ix,iy)
            enddo
         enddo
         write(1,*) ' '
      enddo
c     
      close(1)
      
      IF (maserflag.eq.1) THEN 
         write(*,*) 'WARNING: Masing detected!'
         write(*,*) 'Will only warn once for this transition.'
      ENDIF

c
#else
      write(*,*) 'ERROR: Pos-vel only for 2-D'
      stop 13
#endif /* ifndef RADGRID_ONEDIM */
      return
      end




c     =============================================================
c               ROUTINES FOR THE RECTANGULAR IMAGES
c
c     While the circular images (above) are very robust, in that
c     they always resolve all detail, they may be less handy for
c     actual representation of an image of the object. A rectangular
c     grid is more useful. Here are the routines for doing so.
c     =============================================================


#     ifndef RADGRID_ONEDIM
c     --------------------------------------------------------
c             DO TRANSFER ALONG ALL PHOTON TRAJECTORIES
c                 FOR ONE FREQUENCY AND MAKE IMAGE
c 
c                   (RECTANGULAR IMAGE VERSION)
c
c     This routine is similar to the routine
c              make_freq_image_circular(inu)
c     but this time the image is rectangular.
c
c     REMARK: If the rectangular image has a resolution that
c             resolves the central star, then the routine is
c             still not capable of dealing with that!!!
c     --------------------------------------------------------
      subroutine make_freq_image_rectang
      implicit none
c
      integer inu
      doubleprecision charintline
c
#include "common_grid.h"
#include "common_telescope.h"
#include "common_radfield.h"
#include "common_boundary.h"
#include "common_lines.h"
c
      integer ix,iy,iray,ir,iphi,iins,ifr,ii,istar
      doubleprecision dr,bb,rmax,srat,srat1,dummy
c
c     De istar option
c
c      write(*,*) '---------------------'
c      write(*,*) 'WARNING: istar.ne.0!!'
c      write(*,*) '---------------------'
c      istar = 1
      istar = 0 
c
c     Check if the ray parameters have indeed been set up
c
      if(rays_ready.ne.321) then
          write(*,*) 'ERROR, make_freq_image_rectang(): ',
     %               'Ray paramters not yet set'
          stop 13
      endif
c
c     Reset error counter
c
      transerr_overshoot = 0
c
c     Check if rays_nrphi is not too large
c
      if(rays_nrx.gt.FRSIZE_INF_X) then
          write(*,*) 'Sorry, rays_nrx is too large'
          stop 13
      endif
      if(rays_nry.gt.FRSIZE_INF_Y) then
          write(*,*) 'Sorry, rays_nry is too large'
          stop 13
      endif
c
c     Then copy some information
c
      imrec_nrx       = rays_nrx
      imrec_nry       = rays_nry
      imrec_sizepix_x = rays_sizepix_x
      imrec_sizepix_y = rays_sizepix_y
c
c     Loop over all rays of the rectangular image
c
      iray = 1
      rmax = rsi_x_c(irsi_frsizex,1)
      do ix=1,imrec_nrx
         do iy=1,imrec_nry
            if(rp_b(iray).lt.0.999d0*rmax) then
#              ifndef TELESCOPE_TRAJARR
               call make_trajectory_c(rp_x0(iray),
     %              rp_z0(iray),rp_theta0(iray),rp_s0(iray),istar)
#              endif
               do inu=1,freq_nr
                  imrec_int(inu,ix,iy) = 
     %                 charintline(ispec_currentline,inu,iray,0.d0)
#                 ifndef TELESCOPE_NOIMTAU
                  imrec_tau(inu,ix,iy) = char_tau
#                 endif
               enddo
            else
               do inu=1,freq_nr
                  if(iradbnd_out_itype.eq.3) then
                     imrec_int(inu,ix,iy) = radbnd_interstellfield(inu)
                  else
                     imrec_int(inu,ix,iy) = 0.d0
                  endif
#                 ifndef TELESCOPE_NOIMTAU
                  imrec_tau(inu,ix,iy) = 0.d0
#                 endif
               enddo
            endif
            iray = iray + 1
         enddo
      enddo
c
c     Artificially add the central star to the image, if/when it is 
c     not resolved. (This feature was added 20 June 2000, on request of 
c     Jeroen Bouwman, Alex de Koter and Roy van Boekel, University of 
c     Amsterdam)
c
      if((imrec_addstar.gt.0).and.(imrec_starunres.gt.0)) then
c
c         First do the transfer along this radial ray
c
#         ifndef TELESCOPE_TRAJARR
          call make_trajectory_c(rp_x0(iray),
     %      rp_z0(iray),rp_theta0(iray),rp_s0(iray),0)
#         endif
c
c         Compute the ratio of the angular surface of the
c         star over the angular surface of four pixels (since the
c         stellar flux has to be averaged over four pixels).
c
          srat  = ( 3.14159265d0 * radbnd_rstar**2 ) / 
     %            ( 4.d0 * rays_sizepix_x * rays_sizepix_y )
          srat1 = 1.d0 - srat  
c          
c         Now loop over line frequencies
c
          do inu=1,freq_nr
c     
c            Find the intensity at the location of the star
c     
             dummy = charintline(ispec_currentline,inu,iray,0.d0)
c     
c            Scale the intensity to the smeared-out intensity over the
c            four central pixels
c
             dummy = dummy * srat
c
c            Add the central star to the center of the image. This is
c            only going to work for images that are star-centered, i.e.
c            no offsets. This has been checked beforehand by the
c            setup_rays_rectangular() subroutine. 
c     
             ix = imrec_nrx / 2
             iy = imrec_nry / 2
             imrec_int(inu,ix,iy)     = dummy + 
     %            srat1 * imrec_int(inu,ix,iy)     
             imrec_int(inu,ix+1,iy)   = dummy + 
     %            srat1 * imrec_int(inu,ix+1,iy)   
             imrec_int(inu,ix,iy+1)   = dummy + 
     %            srat1 * imrec_int(inu,ix,iy+1)   
             imrec_int(inu,ix+1,iy+1) = dummy + 
     %            srat1 * imrec_int(inu,ix+1,iy+1)
          enddo
          iray = iray + 1
       endif
c
c
c
      if(transerr_overshoot.gt.0) then
          write(*,*) ' ' 
          write(*,*) '-----------------------------------------------' 
          write(*,*) '|  WARNING:                                   |'
          write(*,*) '|  While making the image there were serious  |'
          write(*,*) '|  overshoots in the quadrature.              |'
          write(*,*) '|  These have been corrected automatically,   |'
          write(*,*) '|  but they indicate poor grid resolution!    |'
          write(*,*) '-----------------------------------------------' 
      endif
c
      end
#     endif


#     ifndef RADGRID_ONEDIM
c     --------------------------------------------------------------
c            GENERATE THE COMPLETE SET OF PHOTON TRAJECTORIES
c                    FOR THE TELESCOPE AT INFINITY
c 
c                    (RECTANGULAR IMAGE VERSION)
c
c     Similar to setup_rays_circular(), but this time for the
c     rectangular image.
c     --------------------------------------------------------------
      subroutine setup_rays_rectang(nrx,nry,sizepix_x,sizepix_y,
     %                       anginf,phioffset,xoffset,yoffset)
      implicit none
      doubleprecision anginf,sizepix_x,sizepix_y
      doubleprecision phioffset,xoffset,yoffset
      integer nrx,nry
c
#include "common_grid.h"
#include "common_telescope.h"
#include "common_boundary.h"
c
      integer ix,iy,it0,iys,ir,cnt_inspnt,iins,nrxhalf,nryhalf
      doubleprecision x_c,y_c,p_c,r_c,theta0,sinth0,costh0
      doubleprecision xh0,zh0,zh02,dum
c
c     Do some checks
c
#     ifdef RADGRID_ONEDIM
      write(*,*) 'ERROR: in 1-D only circular image allowed, not ',
     %           'rectangular image'
      stop 13
#     endif
      if(nrx.gt.FRSIZE_INF_X) then
          write(*,*) 'ERROR Telescope: nrx larger than ',
     %                'FRSIZE_INF_X'
          write(*,*) 'Recompile with larger FRSIZE_INF_X'
          stop 13
      endif
      if(nry.gt.FRSIZE_INF_Y) then
          write(*,*) 'ERROR Telescope: nry larger than ',
     %                'FRSIZE_INF_Y'
          write(*,*) 'Recompile with larger FRSIZE_INF_Y'
          stop 13
      endif
      if((nrx+1)/2.ne.nrx/2) then
          write(*,*) 'ERROR Telescope: nrx must be even'
          stop 13
      endif
      if((nry+1)/2.ne.nry/2) then
          write(*,*) 'ERROR Telescope: nry must be even'
          stop 13
      endif
      if(sizepix_x.le.0.d0) then
          write(*,*) 'ERROR: telescope.F/setup_rays_rectang():'
          write(*,*) 'sizepix_x.le.0'
          stop 13
      endif
      if(sizepix_y.le.0.d0) then
          write(*,*) 'ERROR: telescope.F/setup_rays_rectang():'
          write(*,*) 'sizepix_y.le.0'
          stop 13
      endif
c
c     If the anginf is too close to the z-axis, some formulae
c     will diverge. Therefore it is required to be not smaller 
c     than 0.1. In the future this must be made smaller, but right
c     now I take it on the safe side.
c
      if(anginf.lt.1.d-1) then
          write(*,*) 'WARNING: In telescope.F/setup_rays_circular():'
          write(*,*) '         For an observer at inclination angle ',
     %           'smaller than 0.1,'
          write(*,*) '         the current routine may make ',
     %     'round-off errors.'
          write(*,*) '         Therefore I have put the inclination ',
     %     'angle at Theta=0.1. Hope you agree...'
          anginf = 0.1
      endif
c
c     Use the Theta0 as given by the user
c     
      it0            = 1
      theta0         = anginf + 1.d-4
      sinth0         = sin(theta0)
      costh0         = cos(theta0)
      rays_nrx       = nrx
      rays_nry       = nry
      rays_sizepix_x = sizepix_x
      rays_sizepix_y = sizepix_y
c     
      ir       = 1
      nrxhalf  = nrx / 2
      nryhalf  = nry / 2
c
      do ix=1,nrx
          do iy=1,nry
              x_c = ( ix - nrxhalf - 0.5d0 ) * sizepix_x - xoffset
              y_c = ( iy - nryhalf - 0.5d0 ) * sizepix_y - yoffset
              r_c = sqrt( x_c**2 + y_c**2 )
              if(x_c.eq.0.d0) then
                  x_c = x_c + 0.001d0*sizepix_x
              endif
              p_c = atan( y_c / x_c ) - phioffset
              if(x_c.lt.0.d0) then
                  p_c = p_c + 3.14159265359d0
              endif
   17         continue
              if(p_c.lt.0.d0) then
                  p_c = p_c + 6.28318530718d0
                  goto 17
              endif
   18         continue
              if(p_c.ge.6.28318530718d0) then
                  p_c = p_c - 6.28318530718d0
                  goto 18
              endif
              zh0  = - sin(p_c) / sinth0
              zh02 = zh0*zh0
              dum  = 1.d0 - zh02*sinth0*sinth0
              dum  = dum + 1d-4
              if(dum.lt.0.d0) then
                  write(*,*) 'ERROR in setup_rays_rectang'
                  stop 13
              endif
ccc              xh0  = sqrt( dum )
              if(cos(p_c).gt.0.d0) then
                  xh0 = sqrt( dum )
              else
                  xh0 = -sqrt( dum )
              endif
c
c             Now store the results in the arrays
c
              rp_x0(ir)      = r_c * xh0
              rp_z0(ir)      = r_c * zh0
              rp_theta0(ir)  = theta0
              rp_s0(ir)      = 1.d30
              rp_b(ir)       = r_c
c
c             Increase pointer
c
              ir = ir + 1
              if(ir.ge.MAXRAYNR) then
                  write(*,*) 'Exceeded maximum number ',
     %                'of rays!!'
                  stop 13
              endif
c
          enddo
      enddo
c
c     Add a central ray, to be able to artificially add the central
c     star to the image, even though it is not resolved. (This feature 
c     was added 20 June 2000, on request of Jeroen Bouwman,
c     Alex de Koter and Roy van Boekel, University of Amsterdam)
c
      if(imrec_addstar.gt.0) then
c
c         For the moment we implement this feature only for the
c         case of star-centred images (i.e. xoffset=0 and yoffset=0)
c
          if((xoffset.ne.0.d0).or.(yoffset.ne.0.d0)) then
              write(*,*) 'PROBLEM: The artificial addition of the ',
     %                    'unresolved central star in the'
              write(*,*) '         rectangular images is an added ',
     %                    'special feature that has been '
              write(*,*) '         implemented only for the special ',
     %                    'case when the image is centered on the'
              write(*,*) '         central star: xoffset=yoffset=0.'
              write(*,*) '         So please put these offsets --> 0'
              stop 13
          endif
c
c         Add the single central ray
c
          rp_x0(ir)      = 0.d0
          rp_z0(ir)      = 0.d0
          rp_theta0(ir)  = theta0
          rp_s0(ir)      = 1.d30
c
c         Increase pointer
c
          ir = ir + 1
          if(ir.ge.MAXRAYNR) then
              write(*,*) 'Exceeded maximum number ',
     %             'of rays!!'
              stop 13
          endif
c
c         Even though we always (for safety) add the central beam, we 
c         may not need it. So check whether or not the artificial 
c         addition of the star is necessary (i.e. if the star is 
c         resolved or not).
c
          if(sizepix_x**2+sizepix_y**2.gt.radbnd_rstar**2) then
c
c             Okay, star is indeed unresolved, so signal to the
c             imager that the central beam must be used and added
c             to the image.
c
              imrec_starunres = 1
          elseif((sizepix_x.gt.0.4*radbnd_rstar).or.
     %           (sizepix_y.gt.0.4*radbnd_rstar)) then
c
c             Dangerous marginal resolving of central star...
c
              write(*,*) 'WARNING: Central star is only marginally'
              write(*,*) '         resolved by the rectangular image'
              write(*,*) '         pixels. So: the star flux is '
              write(*,*) '         not very reliable!'
              imrec_starunres = 0
          else
              imrec_starunres = 0
          endif
      endif
c
c     Count the rays
c
      rays_amount = ir - 1
c
c     Put flag up to say that the ray parameters have been set up
c     This is an internal consistency check. 
c
      rays_ready = 321
c
c     Now render the rays... This takes some CPU time...
c     Only do this if the TELESCOPE_TRAJARRAY option is set, meaning
c     that we wish to pre-compute all trajectories. For high-resolution
c     images this seems not so good, because it is too memory expensive.
c     So per default for the telescope it is better to switch off this
c     option (i.e. do not define TELESCOPE_TRAJARRAY).
c
#ifdef TELESCOPE_TRAJARRAY
      write(*,*) 'Rendering trajectories'
c      write(*,*) 'THIS OPTION IS OUT OF ORDER FOR RECTANGULAR IMAGE'
c      stop 13
      call setup_all_trajectories(1)
#endif
c
c     Call the reindexation routine, necessary for the interpolations
c     
c     WARNING: This is actually a hack, because iangset=1 has nothing to
c              do with the camera functioning. But the reindexation
c              is also necessary for the theta mirroring, so that's
c              why we need to do it anyway.
c
      call make_index(1)
c
      return
      end
#     endif





#     ifndef RADGRID_ONEDIM
c     --------------------------------------------------------
c               WRITE THE IMAGE AT INFINITY TO A FILE
c 
c                    (RECTANGULAR IMAGE VERSION)
c
c     Write the image to a file. Similar to circular version
c     --------------------------------------------------------
      subroutine write_image_rectang(ifreq)
      implicit none
c
      integer ifreq
c
#include "common_grid.h"
#include "common_telescope.h"
c
      integer ix,iy,inu,inu0,inu1
c
c     Dependent on the macro SMALL_MEMORY we save all frequencies
c     or just a single snapshot.
c     
#     ifdef SMALL_MEMORY
      inu0 = 1
      inu1 = 1
#     else
      if(ifreq.gt.0) then
          inu0 = ifreq
          inu1 = ifreq
      else
          inu0 = 1
          inu1 = freq_nr
      endif
#     endif      
c
c     Open the file and write a header. The '0' is for the 
c     possible mirror symmetry in phi, which is not yet 
c     implemented.
c
      open(unit=1,file='image.dat',status='unknown')
      write(1,*) imrec_nrx,imrec_nry,inu1+1-inu0
      write(1,*) imrec_sizepix_x,imrec_sizepix_y
c
c     Write the main stuff
c
      write(1,*) ' '
      do inu=inu0,inu1
          do iy=1,imrec_nry
              do ix=1,imrec_nrx
                  if(imrec_int(inu,ix,iy).gt.1d-97) then
                      write(1,*) imrec_int(inu,ix,iy)
                  else
                      write(1,*) 0.d0
                  endif
              enddo
          enddo
          write(1,*) ' '
      enddo      
      close(1)
c
c     Write the optical depth map to file
c
#     ifndef TELESCOPE_NOIMTAU
      open(unit=1,file='imtau.dat',status='unknown')
      write(1,*) imrec_nrx,imrec_nry,inu1+1-inu0
      write(1,*) imrec_sizepix_x,imrec_sizepix_y
c
c     Write the main stuff
c
      write(1,*) ' '
      do inu=inu0,inu1
          do iy=1,imrec_nry
              do ix=1,imrec_nrx
                  if(imrec_tau(inu,ix,iy).gt.1d-97) then
                      write(1,*) imrec_tau(inu,ix,iy)
                  else
                      write(1,*) 0.d0
                  endif
              enddo
          enddo
          write(1,*) ' '
      enddo      
      close(1)
#     endif
c
      end
#     endif




c     =============================================================
c                  ROUTINES FOR INTERNAL USE ONLY
c     =============================================================


c     =============================================================
c
c                     EXPLANATION OF GRID ARRANGEMENT
c
c     I have irsi_frsizex times irsi_frsizey cells on my grid. But on 
c     the border I have put extra cells in order to be able to properly 
c     impose boundary conditions. So the x index of the array runs from 
c     0 to irsi_frsizex+1 where ix=0 is the left extra boundary cell and 
c     ix=irsi_frsizex+1 is the right extra boundary cell. I number the 
c     cell interfaces from 0 to FrSize+2, where interface 0 is the left 
c     boundary of cell 0 and interface irsi_frsizex+2 is the right 
c     interface of cel irsi_frsizex+1 (for the Y direction similar).
c
c                  Example: x-grid for irsi_frsizex = 3
c                                                                            
c             bndry cell     cell        cell        cell      bndry cell    
c            |           |           |           |           |           |   
c     ix=    |     0     |     1     |     2     |     3     |     4     |   
c          Intrf       Intrf       Intrf       Intrf       Intrf       Intrf 
c            0           1           2           3           4           5   
c
c     In 2-D one has 2-D cells, 1-D cell boundaries (which are line-elements)
c     and 0-D cell nodes (which are the begin and end points of the
c     line elements). The cells/cell-centers, and the cell boundaries are 
c     numbered as:
c
c            .---(0,0)---.---(1,0)---.---(2,0)---.---(3,0)---.---(4,0)---.
c            |           |           |           |           |           |
c            |           |           |           |           |           |
c          (0,0) (0,0) (1,0) (1,0) (2,0) (2,0) (3,0) (3,0) (4,0) (4,0) (5,0)
c            |           |           |           |           |           |
c            |           |           |           |           |           |
c            .---(0,1)---.---(1,1)---.---(2,1)---.---(3,1)---.---(4,1)---.
c            |           |           |           |           |           |
c            |           |           |           |           |           |
c          (0,1) (0,1) (1,1) (1,1) (2,1) (2,1) (3,1) (3,1) (4,1) (4,1) (5,1)
c            |           |           |           |           |           |
c            |           |           |           |           |           |
c            .---(0,2)---.---(1,2)---.---(2,2)---.---(3,2)---.---(4,2)---.
c            |           |           |           |           |           |
c            |           |           |           |           |           |
c          (0,2) (0,2) (1,2) (1,2) (2,2) (2,2) (3,2) (3,2) (4,2) (4,2) (5,2)
c            |           |           |           |           |           |
c            |           |           |           |           |           |
c            .---(0,3)---.---(1,3)---.---(2,3)---.---(3,3)---.---(4,3)---.
c
c     The numbering of the cell nodes is not depicted here, but it is 
c     basically similar, starting with (0,0) in the far upper left.
c
c     In the transfer routine we also use a mesh that is spaced precisely
c     in between the above mesh lines. This is used for collecting all the
c     intensities and computing the scatterin integral
c
c                  |           |           |           |           |
c                  |           |           |           |           |
c                  |           |           |           |           |
c            ----(0,0)-(1,0)-(1,0)-(2,0)-(2,0)-(3,0)-(3,0)-(4,0)-(4,0)----
c                  |           |           |           |           |
c                  |           |           |           |           |
c                (0,1)       (1,1)       (2,1)       (3,1)       (4,1)       
c                  |           |           |           |           |
c                  |           |           |           |           |
c            ----(0,1)-(1,1)-(1,1)-(1,1)-(2,1)-(1,1)-(3,1)-(1,1)-(4,1)----
c                  |           |           |           |           |
c                  |           |           |           |           |
c                (0,2)       (1,2)       (2,2)       (3,2)       (4,2)       
c                  |           |           |           |           |
c                  |           |           |           |           |
c            ----(0,2)-(1,2)-(1,2)-(1,2)-(2,2)-(1,2)-(3,2)-(1,2)-(4,2)----
c                  |           |           |           |           |
c                  |           |           |           |           |
c                  |           |           |           |           |
c
c     Example of numbering of theta in the actual spherical geometry.
c     This is actually the same as described above, but then depicted
c     in the correct geometry. Take irsi_frsizey.eq.4. Then we have 
c     the following cells & cell boundaries:
c
c                 1                                         
c        0        |        2       The cells on the left of the z-axis
c         .    0  |  1   /         are the ghost cells in the theta 
c        XX .     |    /           direction. They have no physical
c       XXXXX .   |  /   2         meaning, and that's why the cell
c       XXXXXXX . |/               boundary is dotted .....
c       XXXXXXXXX .--------- 3
c       XXXXXXX . |\               The XXX means that there is not only
c       XXXXX .   |  \   3         no physical space there, but also
c        XX .     |    \           no ghost cell. Only the cells 1, 2,
c         .    5  |  4   \         3 and 4 are physical. 
c        6        |        4                                      
c                 5                                            
c
c     =============================================================



c     -------------------------------------------------------------
c                  MAKE PHOTON TRAJECTORY IN I-GRID
c
c     This routine finds the points where a photon trajectory in
c     the 3-D space crosses the (R,theta) grid surfaces. It computes
c     the length of each of the trajecory elements in between the
c     crossings. In the (R,theta) plane each photon moves as a
c     hyperbolic curve. So we wish to find the crossing between 
c     hyperbolic curves in the (R,theta) plane and the lines 
c     R=const and theta=const. 
c
c     A more detailed explanation... Start with Cartesian coordinate
c     system (x,y,z). The z-axis is the axis of axial symmetry. Any
c     photon trajectory is a function x(s),y(s),z(s), where s is the
c     path length along the line. By applying axial symmetry, and by 
c     gauging s in a particular way, one can parameterize the line as
c
c       x = x_0
c       y = sin(theta_0) * s
c       z = z_0 + cos(theta_0) * s
c
c     The theta_0 is  the angle at s --> + infinity.
c     The corresponding (R,theta,phi) are
c
c       R     = sqrt( x_0^2 + z_0^2 + s^2 + 2 * z_0 * cos(theta_0) * s )
c       theta = taninv ( sqrt( x_0^2 + sin^2(theta_0) * s^2) / 
c                         ( z_0 + cos(theta_0) * s ) )
c       phi   = taninv ( sin(theta_0) * s / x_0 )
c
c     Since we have axial symmetry, the phi is irrelevant, and we'll
c     ignore is from now on. 
c
c     The curves in the (R,theta) plane are hyperbolic curves. Define
c     D = sqrt( x^2 + y^2 ), then the hyperbolic curve is simply
c
c        D^2 - tan^2(theta_0) * ( z - z_0 )^2 = x_0^2
c
c     which is clearly a hyperbolic curve with locus at the z-axis.
c     The minimal distance to the origin is
c 
c       R_min = sqrt( x_0^2 + z_0^2 (1-cos^2(theta_0)) )
c
c     So the hyperbolic curve crosses the circles of constant radius
c     R for all R obeying
c
c        R > R_min = sqrt( x_0^2 + z_0^2 (1-cos^2(theta_0)) )
c
c     The hyperbolic curve crosses lines of constant theta, as long as 
c     theta obeys 
c
c        - | theta_0 | < theta < | theta_0 |
c
c     Within this domain we wish to find the crossing points of the
c     hyperbolic photon trajectory with the circles of constant radius
c     and the (radially pointing) lines of constant theta. The discrete 
c     set of circles and radial lines I am talking about are the 
c     grid lines of the computational domain.
c
c     The grid that I use here is not necessarily the same grid as the
c     hydrodynamic grid. The reason is that this will be computationally
c     too expensive. The radiative transfer will be done on a much 
c     coarser grid. 
c
c
c                       WHAT IS BEING CALCULATED:
c     
c     The result of this routine is a 1-D array of elements. Each element
c     describes a crossing with either a R=const circle or a theta=const
c     line, dependent on which comes first along the photon trajectory.
c     Whether it is a R=const crossing or theta=const crossing, is stored
c     in th_icross(is). Suppose it is a crossing with an R=const circle.
c     Then the th_irad(is) gives the number i of the grid circle R_i, and
c     th_itheta(is) gives the number i of the theta line element that is
c     crossed. The line elements are numbered the same as the cell centers
c     are (see above for a pictogram of the grid numbering). So line element
c     2 is in between node 2 and node 3. 
c
c       tr_icross(is)          = 1 --> crossing of radius=const
c                              = 2 --> crossing of theta=const
c
c     REMARK:  The grid for radiation always spans theta=[0,pi]. The
c              routines for radiative transfer do not care about a 
c              possible equatorial mirror symmetry, only when the actual
c              source terms are computed. If mirror symmetry is present
c              then one simply chooses the set of trajectories only for
c              upwards moving photons for example. Of course, when 
c              scattering and emission/absorption are done, then 
c              one should take into account whether or not the mirror
c              symmetry is present. 
c
c     WARNING: This routine expects that there exists exactly one 
c              theta=Pi/2 cell boundary. 
c
c     WARNING: The amount of grid cells in theta should be even,
c              since we include both quadrants, and do not allow
c              a cell to be half-up half-down the equator.
c
c     REMARK:  Added the possibility for 1--D spherical symmetry
c              here. (30-12-98)
c
c     REMARK:  Added (mu,phi) here, in particular for the line transfer
c              (01-05-99)
c
c     REMARK:  From now on, the responsibility of including the central
c              star is solely on the imager. In the make_trajectory
c              routine the istar must always be 0!
c
c     REMARK:  Bugfix: rmint is now rsi_x_c(1,1)*(1.0+1.d-5), instead
c              of rsi_x_c(1,1)*(1.0-1.d-5). Also I added rmaxt and
c              rmaxr, in place of rmax. 
c              (27-06-00)
c
c     REMARK:  All PI changed to pi, with parameter(pi=PICONST)
c              And all 3.14 --> pi and 6.28 --> 2.d0*pi
c              (12-09-00)
c
c     REMARK:  All safety numbers (1.d-6 and 1.d-5) changed to eps
c              and epsplus respectively. 
c              (12-09-00)
c     
c     -------------------------------------------------------------
      subroutine make_trajectory_c(x0,z0,theta0,send,istar)
      implicit none
      doubleprecision x0,z0,theta0,send
      integer istar
c
#include "common_grid.h"
#include "common_boundary.h"
#include "common_telescope.h"
c
      doubleprecision th_radius,th_theta,th_s
      dimension th_radius(0:RAYSIZE)
      dimension th_theta(0:RAYSIZE)
      dimension th_s(0:RAYSIZE)
      integer th_ir,th_itheta
      dimension th_ir(0:RAYSIZE)
      dimension th_itheta(0:RAYSIZE)
      doubleprecision r_radius,r_theta,r_s
      dimension r_radius(0:RAYSIZE)
      dimension r_theta(0:RAYSIZE)
      dimension r_s(0:RAYSIZE)
      integer r_ir,r_itheta
      dimension r_ir(0:RAYSIZE)
      dimension r_itheta(0:RAYSIZE)
      doubleprecision ex_radius,ex_theta,ex_s,ex_wksp
      dimension ex_radius(RAYEXPT)
      dimension ex_theta(RAYEXPT)
      dimension ex_wksp(RAYEXPT)
      dimension ex_s(RAYEXPT)
      integer ex_ir,ex_itheta,ex_iwksp,ex_iiwksp
      dimension ex_ir(RAYEXPT)
      dimension ex_itheta(RAYEXPT)
      dimension ex_iwksp(RAYEXPT)
      dimension ex_iiwksp(RAYEXPT)
c
      integer ix,iy,iyeq,iyend,is,iss,ist,isr,iad,irng
      integer isex,nrex,isnr,isdblnr,isrt
      integer ir_min,ith_amount,ir_amount
      doubleprecision pitheta0,r,theta,pitheta,tanth2,sdiscr
      doubleprecision costh0,sinth0,costh02,sinth02,a,b,c
      doubleprecision s0,s1,s2,s3,bimpact,sar1,sar2,dum1,dum2,ds
      doubleprecision rmaxr,rmaxt,rminr,rmint,sbeg,sprev,rr
      dimension sar1(0:RAYSIZE)
      dimension sar2(0:RAYSIZE)
      integer iyar,iup
      dimension iyar(0:RAYSIZE)
      doubleprecision rrgrid(FRSIZE_MAX+2),ttgrid(FRSIZE_MAX+2)
      doubleprecision snew,znew,bnew,sinphi,dummy
c
      doubleprecision pi
      parameter(pi=PICONST)
c
      doubleprecision eps,epsplus
      parameter(eps=TELESC_EPS)
      parameter(epsplus=1.d1*TELESC_EPS)
c
      if(istar.ne.0) then
          write(*,*) 'SAFETY STOP: make_trajectory_c called() with'
          write(*,*) '             argument istar.ne.0.'
          write(*,*) '             NOT ALLOWED ANYMORE.'
          write(*,*) '             Consult the author C.P. Dullemond'
          stop 13
      endif
c
      do ix=1,irsi_frsizex
          rrgrid(ix) = rsi_x_c(ix,1)
      enddo
      do iy=1,irsi_frsizey
          ttgrid(iy) = rsi_x_c(iy,2)
      enddo
c
c     Compute 0.5*pi-theta0 for convenience
c
      pitheta0 = 0.5d0*pi - theta0
c
c     Compute cos(theta0) and sin(theta0) for later convenience
c     
      costh0   = cos(theta0)
      sinth0   = sin(theta0)
      costh02  = costh0**2
      sinth02  = sinth0**2
c
c     If in 1--D mode (RADGRID_ONEDIM defined) then we must
c     skip the search for Theta crossings, since there will
c     be none, and searching for it might crash the system.
c
#ifdef RADGRID_ONEDIM
c
c     If 1--D mode: count the amount of theta--crossings as 0
c     This will produce the right tr_amount, and also automatically
c     assure that in the sorting routine, the Theta crossings
c     will be skipped altogether.
c
      ith_amount = 0
c
#else
c
c     For 2--D mode: find the Theta crossings
c
c     Compute the iy for which rsi_x_c(iy,2) = Pi/2, and
c     the iyend for which iyend - iy is opposite to the 
c     equator
c
      iyeq  = irsi_frsizey / 2
      iyend = irsi_frsizey  
c
c     First I am going to do the crossings with the theta=const
c     cell boundaries. Treat the cases z_0<0 and z_0>0 separately.
c     Take as an example z_0>0. There are 3 ranges of theta. There 
c     exists a theta_crit such that for theta<theta_crit there is
c     no intersection between the hyperbola and the theta=const
c     line. This is depicted here:
c
c        |           /       .......
c        |          /    ....              
c        |         /  ...        
c        |        / ..         
c        | theta /..                   
c        | crit /.              The meaning of theta_crit
c        |__   /.               -------------------------
c        |  \ /.                The / is the theta=const line.
c        |   /  .               The ... is the photon trajectory.
c        |  /    .              The theta of the theta=const line 
c        | /      ..            is marginally too small to 
c        |/         ..          intersect. This is the theta_crit
c        .------------...-----         
c        |               ....          
c        |                          
c
c     For theta_crit < theta < |theta0|, one has two intersections 
c     as depected below. 
c
c        |                          / 
c        |                        /.......
c        |                  ..../..
c        |              ....  / ^          
c        |           ...    /   |      
c        |         ..     /   crossing 2        
c        |       ..     /             
c        |      .     /               Two intersections (same quadrant)  
c        |     .    /                 ---------------------------------
c        |    .   /                   The / is the theta=const line.
c        |     ./ <-- crossing 1      The ... is the photon trajectory.
c        |    / .                     There are two intersections in the
c        |  /    ..                   same quadrant. The photon moves
c        |/        ..                 upwards in this diagram.
c        .-----------...-----         
c        |\             ....          
c        |  \              ....         
c
c     For |theta0| < theta < pi-|theta0| there are again two 
c
c        |          ...           
c        |        ..                  Two intersection (other quadrant)
c        |      ..          ____/     ---------------------------------
c        |     .       ____/          The / is the theta=const line. Two
c        |    .   ____/               lines are shown: for theta and for
c        |   .___/                    pi-theta. Each has an intersection.
c        |__/.                        The ... is the photon trajectory.
c        .__-.---------------         
c        |  \_.__  
c        |     . \____ 
c        |      ..    \____ 
c        |        ..       \____
c        |          ..          \   
c
c
c     Okay, now let's go... First we make a loop over theta, 
c     from 0 to 0.5*pi. We do not yet discriminate between 
c     theta and pi-theta. We treat them on the same footing
c     in this first loop. In this loop we neglect the equator.
c
      is = 1
      do iy=1,iyeq
          theta   = rsi_x_c(iy,2)
          pitheta = 0.5d0 * pi - theta
          tanth2  = ( tan(theta) )**2
          a       = tanth2 * costh02 - sinth02
          b       = 2.d0 * tanth2 * costh0 * z0
          c       = tanth2 * z0 * z0 - x0 * x0
          sdiscr  = b*b-4.d0*a*c
c
c         Check if there are solutions, either for
c         theta or pi-theta.
c          
          if(sdiscr.gt.0.d0) then
c
c             Square root of discriminant
c
              sdiscr = sqrt(sdiscr)
c
c             Find the s position along the trajectory of the
c             two possible intersections
c             
              iyar(is) = iy
              sar1(is) = ( -b - sdiscr ) / (2.d0*a)
              sar2(is) = ( -b + sdiscr ) / (2.d0*a)
c
c             Order them properly
c
              if(sar1(is).gt.sar2(is)) then
                  dum1     = sar1(is)
                  sar1(is) = sar2(is)
                  sar2(is) = dum1
              endif
c
c             Increase is
c
              is = is + 1
c
          endif
c
      enddo
c
      isnr = is - 1
c
c     Count the number of times this condition is satisfied:
c             |0.5*pi-theta| > |0.5*pi-theta0|
c     This condition means that the two intersections 
c     are for the same theta (not one for theta and one
c     for pi-theta, but both for the same theta).
c
      isdblnr = 0
      do is=1,isnr
          pitheta = 0.5d0*pi - rsi_x_c(iyar(is),2)
          if(abs(pitheta).gt.abs(pitheta0)) then
              isdblnr = isdblnr + 1
          endif
      enddo
c
c     Now collect the solutions in the correct order along
c     the photon trajectory.
c
      if(pitheta0.gt.0.d0) then
          iup = 1
      else
          iup = 0
      endif
      if(z0*pitheta0.gt.0.d0) then
c
c         Upward moving photon & double solutions are in 
c         upper quadrant OR downward moving photon & 
c         double solution are in lower quadrant. 
c
c         In order to communicate here, let's do as if 
c         photon is moving upwards. But note that in all
c         comments here one can also simultaneously replace 
c         upwards with downwards and upper quadrant with 
c         lower quadrant.
c
c         We start in the lower quadrant, moving upwards /
c         (or upper quandrant moving downwards)
c             
c         Loop over the lower quadrant. Skip the 
c         double solutions here. 
c
          is = 1 + isdblnr
          do iss=1,isnr-isdblnr
              th_radius(iss) = sqrt( x0*x0 + 
     %            z0*z0 + sar1(is)*sar1(is) + 
     %            2.d0*z0*costh0*sar1(is) )
              th_s(iss)      = sar1(is)
              if(iup.eq.1) then 
                  th_itheta(iss) = iyend + 1 - iyar(is) ! bugfix aug 98
              else
                  th_itheta(iss) = iyar(is)
              endif
              th_theta(iss)  = rsi_x_c(th_itheta(iss),2)
              is = is + 1
          enddo
cc
cc         Now the equator
cc     
c          th_radius(isnr-isdblnr+1) = sqrt( x0*x0 + 
c     %        z0*z0 * ( (1.d0/costh02) - 1.d0 ) )
c          th_s(isnr-isdblnr+1)      = (- z0) / costh0
c          th_itheta(isnr-isdblnr+1) = iyeq
c          th_theta(isnr-isdblnr+1)  = 0.5*pi
c
c         Now the upper quadrant. First the single solutions
c
          is = isnr
          do iss=isnr-isdblnr+1,2*isnr-2*isdblnr
              th_radius(iss) = sqrt( x0*x0 + 
     %            z0*z0 + sar2(is)*sar2(is) + 
     %            2.d0*z0*costh0*sar2(is) )
              th_s(iss)      = sar2(is)
              if(iup.eq.1) then 
                  th_itheta(iss) = iyar(is)
              else
                  th_itheta(iss) = iyend + 1 - iyar(is) ! bugfix aug 98
              endif
              th_theta(iss)  = rsi_x_c(th_itheta(iss),2)
              is = is - 1
          enddo
c
c         Now do the double solutions, if any
c
          if(isdblnr.gt.0) then
              do iss=2*isnr+1-2*isdblnr,2*isnr-isdblnr   
                  th_radius(iss) = sqrt( x0*x0 + 
     %                z0*z0 + sar1(is)*sar1(is) + 
     %                2.d0*z0*costh0*sar1(is) )
                  th_s(iss)      = sar1(is)
                  if(iup.eq.1) then 
                      th_itheta(iss) = iyar(is)
                  else
                      th_itheta(iss) = iyend + 1 - iyar(is) ! bugfix aug 98
                  endif
                  th_theta(iss)  = rsi_x_c(th_itheta(iss),2)
                  is = is - 1
              enddo
              is = is + 1
              do iss=2*isnr+1-isdblnr,2*isnr
                  th_radius(iss) = sqrt( x0*x0 + 
     %                z0*z0 + sar2(is)*sar2(is) + 
     %                2.d0*z0*costh0*sar2(is) )
                  th_s(iss)      = sar2(is)
                  if(iup.eq.1) then 
                      th_itheta(iss) = iyar(is)
                  else
                      th_itheta(iss) = iyend + 1 - iyar(is) ! bugfix aug 98
                  endif
                  th_theta(iss)  = rsi_x_c(th_itheta(iss),2)
                  is = is + 1
              enddo
          endif
      else
c
c     The double solutions are in lower quadrant
c     and we start in the lower quadrant, moving 
c     upwards. 
c     
c     Loop over the lower quadrant. Start with the 
c     double solutions here. 
c
          is = isdblnr
          if(isdblnr.gt.0) then
              do iss=1,isdblnr
                  th_radius(iss) = sqrt( x0*x0 + 
     %                z0*z0 + sar1(is)*sar1(is) + 
     %                2.d0*z0*costh0*sar1(is) )
                  th_s(iss)      = sar1(is)
                  if(iup.eq.1) then 
                      th_itheta(iss) = iyend + 1 - iyar(is) ! bugfix aug 98
                  else
                      th_itheta(iss) = iyar(is)
                  endif
                  th_theta(iss)  = rsi_x_c(th_itheta(iss),2)
                  is = is - 1
              enddo
              is = is + 1
              do iss=isdblnr+1,2*isdblnr
                  th_radius(iss) = sqrt( x0*x0 + 
     %                z0*z0 + sar2(is)*sar2(is) + 
     %                2.d0*z0*costh0*sar2(is) )
                  th_s(iss)      = sar2(is)
                  if(iup.eq.1) then 
                      th_itheta(iss) = iyend + 1 - iyar(is) ! bugfix aug 98
                  else
                      th_itheta(iss) = iyar(is)
                  endif
                  th_theta(iss)  = rsi_x_c(th_itheta(iss),2)
                  is = is + 1
              enddo
          endif
c     
c         Next the single solutions
c
          is = isdblnr + 1
          do iss=2*isdblnr+1,isnr+isdblnr
              th_radius(iss) =  sqrt( x0*x0 + 
     %            z0*z0 + sar1(is)*sar1(is) + 
     %            2.d0*z0*costh0*sar1(is) )
              th_s(iss)      = sar1(is)
              if(iup.eq.1) then 
                  th_itheta(iss) = iyend + 1 - iyar(is) ! bugfix aug 98
              else
                  th_itheta(iss) = iyar(is)
              endif
              th_theta(iss)  = rsi_x_c(th_itheta(iss),2)
              is = is + 1
          enddo
cc     
cc         Now the equator
cc
c          th_radius(isnr+isdblnr+1) = sqrt( x0*x0 + 
c     %        z0*z0 * ( (1.d0/costh02) - 1.d0 ) )
c          th_s(isnr+isdblnr+1)      = (- z0) / costh0
c          th_itheta(isnr+isdblnr+1) = iyeq
c          th_theta(isnr+isdblnr+1)  = 0.5*pi
c
c         Now the upper quadrant. Skip the double solutions.
c
          is = isnr 
          do iss=isnr+isdblnr+1,2*isnr
              th_radius(iss) = sqrt( x0*x0 + 
     %            z0*z0 + sar2(is)*sar2(is) + 
     %            2.d0*z0*costh0*sar2(is) )
              th_s(iss)      = sar2(is)
              if(iup.eq.1) then 
                  th_itheta(iss) = iyar(is)
              else
                  th_itheta(iss) = iyend + 1 - iyar(is)  ! bugfix aug 98
              endif
              th_theta(iss)  = rsi_x_c(th_itheta(iss),2)
              is = is - 1
          enddo
      endif
c
c     The number of theta crossings
c
      ith_amount = 2 * isnr
c
c     Now find the indexes of the line elements that are crossed,
c     or in other words: which index ix belongs to each of the 
c     th_radius?
c
      do is=1,ith_amount
          r         = th_radius(is)
          if(r.lt.rsi_x_c(1,1)) then
              th_ir(is) = 0
          else if (r.gt.rsi_x_c(irsi_frsizex,1)) then
              th_ir(is) = irsi_frsizex
          else
              call hunt(rrgrid,irsi_frsizex,r,ix)
              th_ir(is) = ix
          endif
      enddo
c
c     Now I've got all the theta=const crossings with the cell
c     boundaries. 
c
#endif /* The switch between 1--D and 2--D */
c
c     Now let's do the same for the radius
c     
c        |              ....         This shows the upper               
c        |           ...             quadrant in which the          
c        |         ..                hyperbolic photon trajectory   
c    R_i |----__ ..                  (dotted line) crosses with     
c        |      \__<-- crossing 1    R=R_i (for some i). This is 
c        |     .   \_                the example with 4 r-gridcells
c        |    .      \               only. The two crossings both 
c        |     .      \              have a theta. These are 
c        |      .      \             r_theta(1) and r_theta(2).
c        |       ..    |             Normally there will be many
c        |         ..  | crossing 2  more crossings. The order 
c        .-----------...-----        (crossing 1 and 2) depends  
c        |             | ....        on whether the photon moves
c        |             |             up or down. Here: down.
c
c     Find the ir_min. This is the smallest radius which 
c     intersects with the hyperbola
c
      bimpact = sqrt( x0*x0 + z0*z0 * ( 1.d0 - costh02 ) )
      do ix=1,irsi_frsizex
          if(rsi_x_c(ix,1).gt.bimpact) then 
              goto 30
          endif
      enddo
      write(*,*) 'Cannot find minimal approach radius!'
      write(*,*) 'Maybe the trajectory lies outside of domain!'
      write(*,*) 'Aborting...'
      write(*,*) ' x0 = ',x0,' z0 = ',z0,' th0 = ',theta0
      write(*,*) ' b  = ',bimpact,' Rmin = ',rsi_x_c(1,1),
     %                         ' Rmax = ',rsi_x_c(irsi_frsizex,1)
      stop 13
   30 continue
      ir_min = ix
      if(bimpact.le.radbnd_rstar) then
cc################################
c          SHOULD THIS NOT BE ELIMINATED? (20 JUNE 2000)
c          write(*,*) 'Bimpact.le.rstar: bimpact=',bimpact
cc################################
          ir_min = 0
      endif
c
c     Save the impact parameter and the theta0 to the tr_ common.
c     This is used in the imager, to include the central star 
c     without really having to include it as an r-crossing.
c
      tr_b    = bimpact
      tr_th0  = theta0
c
c     Compute the amount of crossings expected. Note that we 
c     assume that the nearest approach (bimpact) is never EXACTLY
c     equal to an R=const cell boundary. So there will always
c     be an even number of crossings.
c     
      ir_amount = 2 * ( irsi_frsizex + 1 - ir_min )
c
c     Now make a loop over all R_i within the domain
c     Each of these R_i should have an intersection 
c     with the photon trajectory. We'll compute at which
c     theta this occurs.
c
      is = 1
      do ix=irsi_frsizex,ir_min,-1
          if(ix.eq.0) then
              r   = radbnd_rstar
          else
              r   = rsi_x_c(ix,1)
          endif
          a       = 1.d0
          b       = 2.d0 * costh0 * z0
          c       = x0*x0 + z0*z0 -r*r
          sdiscr  = b*b-4.d0*a*c+eps*b*b
          if(sdiscr.lt.0.d0) then
              write(*,*) 'BUG IN CODE make_trajectory_i() 2 !'
              write(*,*) 'Input to make_trajectory = ',x0,z0,
     %                     theta0,send,istar
              write(*,*) 'Current position in loop = ',ix,r,a,b,c,sdiscr
              stop 13
          endif
          sdiscr = sqrt(sdiscr)
c
c         Find the s position along the trajectory of the
c         two possible intersections
c
          s1 = ( -b - sdiscr ) / (2.d0*a)
          s2 = ( -b + sdiscr ) / (2.d0*a)
c
c         Now compute the corresponding theta, and store
c         them
c
          r_theta(is)               = atan(sqrt(x0*x0+sinth02*s1*s1)
     %                                 /(z0+costh0*s1))
          if(r_theta(is).lt.0.d0) then
              r_theta(is) = r_theta(is) + pi
          endif
          r_s(is)                   = s1
          r_ir(is)                  = ix
          r_radius(is)              = r
c
          r_theta(ir_amount+1-is)   = atan(sqrt(x0*x0+sinth02*s2*s2)
     %                                 /(z0+costh0*s2))
          if(r_theta(ir_amount+1-is).lt.0.d0) then
              r_theta(ir_amount+1-is) = 
     %             r_theta(ir_amount+1-is) + pi
          endif
          r_s(ir_amount+1-is)       = s2
          r_ir(ir_amount+1-is)      = ix
          r_radius(ir_amount+1-is)  = r
c
c         Increase is
c
          is = is + 1
c
      enddo
c
c     Now find the number of the line elements that are crossed,
c     or in other words: which iy belongs to the r_theta?
c
      do is=1,ir_amount
          theta    = r_theta(is)
          call hunt(ttgrid,irsi_frsizey,theta,iy)
          r_itheta(is) = iy
      enddo
c
c     Now we've got both the theta crossings and the radius
c     crossings. Next we should order them, and create an 
c     array of crossings and line elements.
c
c     First the sorting and array of crossings
c
      tr_amount = ir_amount + ith_amount
      ist = 1
      isr = 1
      th_s(ith_amount+1) = 1.d30
      r_s(ir_amount+1)   = 1.d30
      rmaxt = rsi_x_c(irsi_frsizex,1) * (1.d0-epsplus)
      rmaxr = rsi_x_c(irsi_frsizex,1) * (1.d0+epsplus)
      rmint = rsi_x_c(1,1) * (1.d0+epsplus)
      rminr = rsi_x_c(1,1) * (1.d0-epsplus)
c      if(istar.ne.0) then
c          rminr = 0.d0
c      endif
      is = 1
      sbeg  = -1.d30
      sprev = -1.d30
c
c     Check if surface of star is reached, and if so, include
c
      if((bimpact.le.radbnd_rstar).and.(istar.ne.0)) then
c
c         Start of SC is at stellar surface
c
c       (was ok: 15-05-01)   sbeg = sqrt( radbnd_rstar**2 - bimpact**2 ) 
c       (was ok: 15-05-01)   sbeg = sbeg - z0 * costh0
c       (bug: 15-05-01)     sbeg = 0.d0 - z0 * costh0
c
          dum2 = sqrt( radbnd_rstar**2 - bimpact**2 ) - z0 * costh0
c
c         But only if stellar surface nearer than send, are we
c         going to reach the stellar surface. Else (and this happens 
c         only when used for long characteristics, see longchar.F), 
c         we should do as if nothing happened.
c     
          if(dum2.lt.send) then
              sbeg  = dum2
              rminr = 0.5*radbnd_rstar    ! 0.5 is safety nr
          endif
      endif
c
c     Before we proceed, let's check if we do not exceed the
c     array bounds
c
      if(tr_amount.gt.RAYSMN) then
          write(*,*) 'ERROR: telescope.F/make_trajectory_c()'
          write(*,*) '       Exceeded RAYSIZE limit'
          stop 13
      endif
c
c     At the point where the trajectory reaches an R-extremum or
c     a Theta-extremum, we must make sure that the extremum is in
c     fact represented with a point. 
c
c        \       \   . 
c         \      .#
c          |   .   |
c          |  #    |
c          |   .   |
c         /      .#
c                /   .
c
c        R_i   R_i+1
c
c     If not, then we may end up getting fringes in the image. 
c     Here we prepare here these 1 or 2 extra points.
c
#     ifdef LONGCHAR_EXTRA_POINTS
      if(RAYEXPT.lt.2+4*RAYADPT*RAYRNPT) then
          write(*,*) 'ERROR: RAYEXPT too small'
          stop 199
      endif
      do isex=1,RAYEXPT
         ex_s(isex) = 1.d30
      enddo
c      ex_s(1)      = 1.d30
c      ex_s(2)      = 1.d30
c      ex_s(3)      = 1.d30
      isex = 1
     
c
c     Is there a radial extremum? If so, include...
c
#     ifdef LONGCHAR_EX_PTS_RADIUS
      if(bimpact.gt.rsi_x_c(1,1)) then
c
c        The extremum itself:
c
         s2              = 0.d0 - z0 * costh0
         ex_s(isex)      = s2
         ex_radius(isex) = bimpact
         ex_theta(isex)  = atan(sqrt(x0*x0+sinth02*s2*s2)
     %                                 /(z0+costh0*s2))
         if(ex_theta(isex).lt.0.d0) then
             ex_theta(isex) = ex_theta(isex) + pi
         endif
         ex_s(isex)      = s2
         ex_ir(isex)     = ir_min - 1
         call hunt(ttgrid,irsi_frsizey,ex_theta(isex),iy)
         ex_itheta(isex) = iy
         isex = isex + 1
c
c        Some additional points around it:
c
cc#ifdef UNDEF
         s3 = s2
         call hunt0(r_s,ir_amount,s2,isrt)
         do irng=1,min(2*RAYRNPT,4)
             if(irng.eq.1) then
                 s1 = s3
                 s2 = r_s(isrt+1)
             elseif(irng.eq.2) then
                 s1 = r_s(isrt)
                 s2 = s3
             elseif(irng.eq.3) then
                 s1 = r_s(isrt+1)
                 s2 = r_s(isrt+2)
             elseif(irng.eq.4) then
                 s1 = r_s(isrt-1)
                 s2 = r_s(isrt)
             endif
             ds = ( s2 - s1 ) / (1.d0+1.d0*RAYADPT)
             do iad=1,RAYADPT
                 s0              = iad*ds + s1
                 ex_s(isex)      = s0
                 rr              = sqrt( x0*x0 + z0*z0 + 
     %                              s0*s0 + 2.d0*z0*costh0*s0 )
                 ex_radius(isex) = rr
                 ex_theta(isex)  = atan(sqrt(x0*x0+sinth02*s0*s0)
     %                              /(z0+costh0*s0))
                 if(ex_theta(isex).lt.0.d0) then
                     ex_theta(isex) = ex_theta(isex) + pi
                 endif
                 call hunt(rrgrid,irsi_frsizex,ex_radius(isex),ix)
                 ex_ir(isex)     = ix 
                 call hunt(ttgrid,irsi_frsizey,ex_theta(isex),iy)
                 ex_itheta(isex) = iy
                 isex = isex + 1
             enddo
         enddo
cc#endif
      endif
#endif
c
c     Is there a theta extremum? If so, include...
c
#     ifdef LONGCHAR_EX_PTS_THETA
#     ifndef RADGRID_ONEDIM
      s2      = x0**2 * costh0 / ( z0 * sinth02 )
      rr      = sqrt( x0*x0 + z0*z0 + s2*s2 + 2.d0*z0*costh0*s2 )
      if((rr.gt.rsi_x_c(1,1)).and.(rr.lt.rsi_x_c(irsi_frsizex,1))) then
c
c        The extremum itself:
c
         ex_s(isex)      = s2
         ex_radius(isex) = rr
         ex_theta(isex)  = atan(sqrt(x0*x0+sinth02*s2*s2)
     %                                 /(z0+costh0*s2))
         if(ex_theta(isex).lt.0.d0) then
             ex_theta(isex) = ex_theta(isex) + pi
         endif
cc         if((ex_theta(isex).le.rsi_x_c(1,2)).or.
cc     %      (ex_theta(isex).ge.rsi_x_c(irsi_frsizey,2))) then
cc             ex_s(isex) = 1e33
cc             goto 78
cc         endif
         ex_s(isex)      = s2
         call hunt(rrgrid,irsi_frsizex,ex_radius(isex),ix)
         if((ix.eq.0).or.(ix.eq.irsi_frsizex)) then
             write(*,*) 'telescope.F/make_traject_t(): hunt failed'
             stop 192
         endif
         ex_ir(isex)     = ix 
         call hunt(ttgrid,irsi_frsizey,ex_theta(isex),iy)
         ex_itheta(isex) = iy
         isex = isex + 1
c
c        Some additional points around it:
c
cc#ifdef UNDEF
         s3 = s2
         call hunt0(th_s,ith_amount,s3,isrt)
         if((isrt-RAYRNPT+1.lt.1).or.
     %      (isrt+RAYRNPT.gt.ith_amount)) then 
             goto 78
         endif
         do irng=1,min(2*RAYRNPT,4)
             if(irng.eq.1) then
                 s1 = s3
                 s2 = th_s(isrt+1)
             elseif(irng.eq.2) then
                 s1 = th_s(isrt)
                 s2 = s3
             elseif(irng.eq.3) then
                 s1 = th_s(isrt+1)
                 s2 = th_s(isrt+2)
             elseif(irng.eq.4) then
                 s1 = th_s(isrt-1)
                 s2 = th_s(isrt)
             else
                 stop 176
             endif
             if(s2.eq.s1) then 
                 write(*,*) 's2=s1'
                 stop 987
             endif
             if(s2.lt.s1) then 
                 write(*,*) 's2<s1'
                 write(*,*) irng,s1,s2
                 write(*,*) ith_amount,s3,isrt
                 write(*,*) th_s(isrt-1),th_s(isrt),
     %                       th_s(isrt+1),th_s(isrt+2)
                 stop 988
             endif
             ds = ( s2 - s1 ) / (1.d0+1.d0*RAYADPT)
             do iad=1,RAYADPT
                 s0              = iad*ds + s1
                 ex_s(isex)      = s0
                 rr              = sqrt( x0*x0 + z0*z0 + 
     %                              s0*s0 + 2.d0*z0*costh0*s0 )
                 ex_radius(isex) = rr
                 ex_theta(isex)  = atan(sqrt(x0*x0+sinth02*s0*s0)
     %                              /(z0+costh0*s0))
                 if(ex_theta(isex).lt.0.d0) then
                     ex_theta(isex) = ex_theta(isex) + pi
                 endif
                 call hunt(rrgrid,irsi_frsizex,ex_radius(isex),ix)
                 if((ix.eq.0).or.(ix.eq.irsi_frsizex)) then
                     goto 79
                 endif
                 ex_ir(isex)     = ix 
                 call hunt(ttgrid,irsi_frsizey,ex_theta(isex),iy)
                 ex_itheta(isex) = iy
                 isex = isex + 1
   79            continue
             enddo
         enddo
cc#endif
      endif
   78 continue
#     endif /* ifndef RADGRID_ONEDIM */
#     endif
c
c
      nrex = isex - 1
      isex = 1
c
c     Sort these in ascending order in ex_s
c
      call ray_sort(nrex,ex_s,ex_radius,ex_theta,ex_ir,ex_itheta, 
     %              ex_wksp,ex_iwksp,ex_iiwksp)
c

#     endif /* ifdef LONGCHAR_EXTRA_POINTS */
c
c     Now make a loop and sort the points
c
      do iss=1,tr_amount 
#         ifdef LONGCHAR_EXTRA_POINTS
  310     continue
          if(ex_s(isex).lt.min(th_s(ist),r_s(isr))) then
             if((ex_radius(isex).le.rmaxt).and.
     %          (ex_radius(isex).ge.rmint).and.
     %          ((ex_s(isex)-sprev).gt.
     %               eps*ex_radius(isex)).and.
     %          (ex_s(isex).ge.sbeg).and.
     %          (ex_s(isex).le.send)) then
                 tr_icross(is)  = 3
                 tr_radius(is)  = ex_radius(isex)
                 tr_theta(is)   = ex_theta(isex)
                 tr_iradius(is) = ex_ir(isex)
                 tr_itheta(is)  = ex_itheta(isex)
                 tr_s(is)       = ex_s(isex)
                 sprev          = tr_s(is)
                 is  = is + 1
             endif
             isex = isex + 1
             goto 310
         endif
         
c     Sanity check for array boundary crossing - should not happen.
         if(isex.gt.RAYEXPT) then
            write(*,*) 'Internal exception error',isex,RAYEXPT
            STOP 83459
         endif
#       endif
         
         if(th_s(ist).lt.r_s(isr)) then
c
c             A crossing with theta=const
c
            if((th_radius(ist).le.rmaxt).and.
     %           (th_radius(ist).ge.rmint).and.
     %           ((th_s(ist)-sprev).gt.eps*th_radius(ist)).and.
     %           (th_s(ist).ge.sbeg).and.
     %           (th_s(ist).le.send)) then
               tr_icross(is)  = 2
               tr_radius(is)  = th_radius(ist)
               tr_theta(is)   = th_theta(ist)
               tr_iradius(is) = th_ir(ist)
               tr_itheta(is)  = th_itheta(ist)
               tr_s(is)       = th_s(ist)
               sprev          = tr_s(is)
               is  = is + 1
            endif
            ist = ist + 1
         else
c     
c             A crossing with R=const
c
            if((r_radius(isr).le.rmaxr).and.
     %           (r_radius(isr).ge.rminr).and.
     %           ((r_s(isr)-sprev).gt.eps*r_radius(isr)).and.
     %           (r_s(isr).ge.sbeg).and.
     %           (r_s(isr).le.send)) then
               tr_icross(is) = 1
               tr_radius(is)  = r_radius(isr)
               tr_theta(is)   = r_theta(isr)
               tr_iradius(is) = r_ir(isr)
               tr_itheta(is)  = r_itheta(isr)
               tr_s(is)       = r_s(isr)
               sprev          = tr_s(is)
               is  = is + 1
            endif
            isr = isr + 1
         endif
      enddo
      tr_amount = is - 1
c     
c     Now that we have all the points of the trajectory, we are
c     going to add the local photon-angles (mu,phi). 
c
c
c     First convert z0,x0,s to znew,bnew,snew (i.e. the coordinates
c     used in the paper Dullemond & Turolla 1999 and by the short
c     characteristics module). 
c
      znew = z0 * sinth02
      bnew = sqrt(x0**2 + z0**2 * sinth02)
c
c     Now do all points along the characteristic
c
      do is=1,tr_amount 
          snew       = tr_s(is) + z0 * costh0
          tr_mu(is)  = snew / sqrt(bnew**2+snew**2)
          dummy      = bnew * sqrt(bnew**2+snew**2-
     %                    (znew+snew*costh0)**2)
          if(dummy.gt.0.d0) then
              sinphi = ( bnew**2 * costh0 - znew*snew ) / dummy
          else
              sinphi = epsplus
          endif
          if(x0.lt.0.d0) then
              sinphi = asin(sinphi)
          else
              sinphi = pi - asin(sinphi) 
          endif
  103     continue
          if(sinphi.lt.0.d0) then
              sinphi = sinphi + 2.d0*pi
              goto 103
          endif
  104     continue
          if(sinphi.ge.2.d0*pi) then
              sinphi = sinphi - 2.d0*pi
              goto 104
          endif
          tr_phi(is) = sinphi
      enddo
c
      end



#ifdef TELESCOPE_TRAJARRAY
c     --------------------------------------------------------------
c                          CREATE ALL THE RAYS
c     --------------------------------------------------------------
      subroutine setup_all_trajectories(istar)
      implicit none
      integer istar
c
#include "common_grid.h"
#include "common_telescope.h"
c
      integer ir,is
c
c     Check if the ray parameters have indeed been set up
c
      if(rays_ready.ne.321) then
          write(*,*) 'ERROR, setup_all_trajectories(): ',
     %               'Ray paramters not yet set'
          stop 13
      endif
c
c     Loop over rays
c
      do ir=1,rays_amount
c
c         Make ray number ir
c
cc#########################
c          write(*,*) 'x0 = ',rp_x0(ir),'  z0 = ',rp_z0(ir),
c     %               '  t0 = ',rp_theta0(ir)
cc#########################
          call make_trajectory_c(rp_x0(ir),
     %              rp_z0(ir),rp_theta0(ir),rp_s0(ir),istar)
c
c         Copy this ray into the big array 
c         
          atr_amount(ir)         = tr_amount
          do is=1,tr_amount
              atr_radius(is,ir)  = tr_radius(is)
              atr_theta(is,ir)   = tr_theta(is)
              atr_iradius(is,ir) = tr_iradius(is)
              atr_itheta(is,ir)  = tr_itheta(is)
              atr_s(is,ir)       = tr_s(is)
              atr_mu(is,ir)      = tr_mu(is)
              atr_phi(is,ir)     = tr_phi(is)
              atr_icross(is,ir)  = tr_icross(is)
          enddo
          ael_amount(ir)         = el_amount
          do is=1,el_amount
c              ael_ix(is,ir)      = el_ix(is)
c              ael_iy(is,ir)      = el_iy(is)
              ael_ds(is,ir)      = el_ds(is)
          enddo
      enddo      
c
c     Set flag that the rays have been rendered
c
      atr_ready = 321
c
      end
#endif

#ifdef TELESCOPE_TRAJARRAY
c     --------------------------------------------------------------
c                      WRITE TRAJECTORIES TO FILE
c     --------------------------------------------------------------
      subroutine write_trajects()
      implicit none
c
#include "common_grid.h"
#include "common_telescope.h"
c
      integer tram_max,is,ir
c
c     Check if rays have been rendered
c
      if(atr_ready.ne.321) then
          write(*,*) 'ERROR, write_trajects: rays not yet rendered'
          stop 13
      endif
c
      open(unit=1,file='alltraject.dat',status='unknown')
      tram_max = 1
      do ir=1,rays_amount
          if(atr_amount(ir).gt.tram_max) then
              tram_max = atr_amount(ir)
          endif
      enddo
      write(1,*) rays_amount,tram_max
      write(1,*) ' '
      do ir=1,rays_amount
          write(1,*) atr_amount(ir)
          write(1,*) rp_iserrad(ir),rp_iserth0(ir)
          write(1,*) rp_x0(ir),rp_z0(ir),rp_theta0(ir)
          write(1,*) ' '
          do is=1,atr_amount(ir)
              write(1,1) atr_radius(is,ir),atr_theta(is,ir),
     %            atr_s(is,ir),atr_mu(is,ir),atr_phi(is,ir)
          enddo
    1     format(3(E13.6,1X))
          write(1,*) ' '
      enddo
      close(1)
c     
c      open(unit=1,file='allelement',status='unknown')
c      write(1,*) rays_amount
c      write(1,*) ' '
c      do ir=1,rays_amount
c          write(1,*) ael_amount(ir)
c          do is=1,ael_amount(ir)
c              write(1,2) ael_ix(is,ir),ael_iy(is,ir),
c     %            ael_ds(is,ir)
c          enddo
c    2     format(2(I4,1X),1(E13.6,1X))
c          write(1,*) ' '
c      enddo
c      close(1)
c      
      end
#else
c     --------------------------------------------------------------
c                      WRITE TRAJECTORIES TO FILE
c                            <<< DUMMY >>>
c     --------------------------------------------------------------
      subroutine write_trajects()
      implicit none
      write(*,*) 'NO TRAJECTORIES WRITTEN: WRITING DISABLED'
      end
#endif


c     --------------------------------------------------------
c            CHARACTERISTIC INTEGRATION OF (S-I)DTAU
c
c     This performs an integral along one of the trajectories
c     given in the common block /alltraject/. This is usually
c     a routine that is part of a more general routine for
c     computing the luminosity of the object under consideration.
c     Therefore I refrain from elaborate comments here.
c
c     Before calling this function, the trajectories should 
c     have been rendered, by calling setup_rays_circular().
c
c     The index inu is the index referring to which frequency
c     bin this transfer should be done. The only way this is
c     used here is that it uses this index to find the emission
c     and opacity terms, and the stellar spectrum. 
c
c     BUGFIX: For central beam, be sure to add starspec only
c             once!!!
c
c     ARGUMENTS:
c        inu     = index of the frequency bin
c        iray    = index of the pencil beam (a complete set
c                  of rays is set up beforehand, and this is
c                  the index of the present ray)
c        rbeam0  :  -1.d0    = Central star is ignored
c                    0.d0    = Central star is included as usual
c                  > 0.d0    = The radius of the central beam
c                              (if iray=1), necessary for the
c                              starspectrum inclusion
c
c     NOTE: STRONGLY MODIFIED, SPECIALLY FOR LINE TRANSFER!!
c
c     --------------------------------------------------------
      function charintline(iline,inu,iray,rbeam0)
      implicit none
      doubleprecision charintline,rbeam0,starfract
      integer inu,iray,iline
c
      doubleprecision find_src_t,find_alp_t,qdr_src_2
      doubleprecision bplanck,ssd_temp
c
#include "common_grid.h"
#include "common_telescope.h"
#include "common_boundary.h"
#include "common_lines.h"
#include "common_source.h"
#include "common_setup.h"
c
      doubleprecision sj,as,ds,r,theta,dr,dt,mu,phi
      doubleprecision srcd0,srcd1,alpd0,alpd1,charintdum,s0,s1
      doubleprecision temp,r_old,freq
      integer is,ir,it,icr,ir_old,it_old,icr_old,amount
      integer it_eq,istar_done,init
      doubleprecision phiprof0,srcl0,alpl0
      doubleprecision lw0,lw1,dvmu0,dvmu1,Nup0,Nup1,Ndown0,Ndown1
      integer ilfr
      doubleprecision dnu1,bmin,bmax
c
c     Check consistency for line transfer
c
c     NOTE: This check can be removed once the code is working fine
c           (KEES & KLAUS 22.03.07)
c
#     ifdef INCLUDE_LINES
      if((iradproc_line.ne.0).and.
     %     (line_warn_preptrans_done.eq.0)) then
          write(*,*) 'INTERNAL CONSISTENCY FAILURE: '
          write(*,*) '  telescope.F/charintline():'
          write(*,*) '  line transfer is active, but the transitions'
          write(*,*) '  have not been prepared yet.'
          write(*,*) '  Problem is easily solved by calling '
          write(*,*) '  global_prepare_transitions() before '
          write(*,*) '  using the telescope.F routines.'
          stop 13
      endif
#     endif
c
c     Some initializations
c
      charintline = 0.d0
      char_tau    = 0.d0
      istar_done  = 0
      init        = 1

c
c     Make preparations in case there might be an equatorial
c     infinitesimally thin optically thick disk present.
c
      it_eq  = (irsi_frsizey+1)/2
c
c     Find the starting point of the characteristic
c     
      icr    = tr_icross(1)
      r      = tr_radius(1)
      theta  = tr_theta(1)
      mu     = tr_mu(1)
      phi    = tr_phi(1)
      ir     = tr_iradius(1)
      it     = tr_itheta(1)
      dr     = ( r - rsi_x_c(ir,1) ) / 
     %    ( rsi_x_c(ir+1,1) - rsi_x_c(ir,1) )
      dt     = ( theta - rsi_x_c(it,2) ) / 
     %    ( rsi_x_c(it+1,2) - rsi_x_c(it,2) )
c     
c     Boundary conditions: both the inner and outer ones.
c
      if(ir.eq.0) then
c
c         If ir = 0, then we start at the stellar surface, so we must
c         use the star light as starting value of charintline.
c
c         NOTE: WE SHOULD CHECK IF THIS WORKS FINE ALSO FOR THE LINE
c               TRANSFER (KEES & KLAUS 22.03.07)
c
          charintline = radbnd_starspec_line(inu)
c
c         And put the src0 and alp0 to zero for safety (should not
c         be necessary because there is no integration necessary
c         in vacuum)
c
          
          srcd0 = 0.d0
          alpd0 = 0.d0
          lw0   = 0.d0
          dvmu0 = 0.d0
          Nup0  = 0.d0
          Ndown0= 0.d0
c
      else         !!! if(ir.eq.irsi_frsizex) then
c
c         If ir = irsi_frsizex then we are at the outer edge of the
c         domain. 
c
          if(iradbnd_out_itype.eq.2) then
c
c             Include microwave background radiation outer BC
c             (only for molecular lines)
c
              if (iradproc_line.ne.0) then
#             ifdef INCLUDE_LINES
                  iline = line_iline(inu)
                  freq  = linefreq(iline)
                  temp  = TEMPCMB    ! = 2.728d0
                  charintline = 1.47455253991d-47 * (freq**3) /
     %              ( exp(4.7991598d-11*freq/temp) - 1.d0 )
#             else 
                  stop 199
#             endif
              else
                  charintline = 0.d0
              endif
          elseif(iradbnd_out_itype.eq.0) then
c
c             We assume interstellar space to be perfectly black.
c
              charintline = 0.d0
c
          elseif(iradbnd_out_itype.eq.1) then
              write(*,*) 'Outer BC type 1 not allowed for telescope'
              stop 13              
          elseif(iradbnd_out_itype.eq.3) then
c
c             We include the interstellar radiation field here.
c
c             NOTE: CHECK IF THIS ALSO STILL WORKS FOR LINE TRANSFER
c                   (I DOUBT IT). KEES & KLAUS 22.03.07
c
              charintline = radbnd_interstellfield(inu)
c
          else
              write(*,*) 'Telecope: dont know this type of outer BC'
              stop 13              
          endif 
c
c         The starting opacity and source
c
          call get_line_dust_values(icr,ir,it,dr,dt,mu,phi,iline,inu,
     %                 srcd0,alpd0,dvmu0,lw0,Nup0,Ndown0)
c     
      endif
c
c     The loop along the ray
c
      amount = tr_amount
c
#     ifdef MONITOR_CHARINT
         if (inu.eq.2.or.inu.eq.freq_nr/2) then
   	        write(4) amount
         endif
#     endif
      do is=2,amount
c
c         Compute the length of the line element
c
          s1   = tr_s(is)
          s0   = tr_s(is-1)
          ds   = ( s1 - s0 ) 
c
c         Stupidity check
c
          if(ds.lt.0.d0) then
              stop 749
          endif
c
c         A backup of ir,it and icr
c
          ir_old  = ir
          it_old  = it
          icr_old = icr
          r_old   = r
c
c         Find the information for the next node
c
          icr    = tr_icross(is)
          r      = tr_radius(is)
          theta  = tr_theta(is)
          mu     = tr_mu(is)
          phi    = tr_phi(is)
          ir     = tr_iradius(is)
          it     = tr_itheta(is)
c
c         Stupidity check
c
c          write(*,*) is,icr,r,theta,mu,phi,ir,it,amount
          if(ir.eq.0) then
              write(*,*) 'INTERNAL ERROR: '
              write(*,*) '   The stellar surface is not done correctly'
              stop 7454
          endif
c
c         Find the relative location within the grid
c
          dr     = ( r - rsi_x_c(ir,1) ) / 
     %        ( rsi_x_c(ir+1,1) - rsi_x_c(ir,1) )
          dt     = ( theta - rsi_x_c(it,2) ) / 
     %        ( rsi_x_c(it+1,2) - rsi_x_c(it,2) )
c
c         If within inner radius, then ds=0, because no interaction
c
          if((ir.eq.1).and.(ir_old.eq.1).and.
     %       (icr.eq.1).and.(icr_old.eq.1)) then
              ds = 0.d0
              init = 1
          endif
          if((ir.eq.1).and.(ir_old.eq.0).and.
     %       (icr.eq.1).and.(icr_old.eq.1)) then
              ds = 0.d0
              init = 1
          endif
          if((ir.eq.0).and.(ir_old.eq.1)) then
              write(*,*) 'Huh?? Not possible... (charintline())'
              write(*,*) 'is=',is,' nrs=',amount,' ir=',ir,' it=',it,
     %                   ' icr=',icr,' r=',r,' theta=',theta,
     %                   ' mu=',mu,' phi=',phi
              stop 137
          endif
c
c         Implement the inner boundary condition here. But if we are
c         obviously not at the inner edge, then skip this...
c
          if(ir.gt.1) then
              goto 77
          endif
c
c         Ik snap niet die ir.eq.-1 mogelijkheid
c
          if(ir.eq.-1) then
              write(*,*) 'SAFETY STOP: charintline(): ir.eq.-1'
              stop 13
          endif
c
c         Switch inner boundary condition type
c
          if(iradbnd_in_itype.eq.0) then
              write(*,*) 'STOPPED: Inner BC type 0 is temporarily ',
     %               'disabled'
              write(*,*) '         Use type 2 instead, and choose R* ',
     %               'tiny.'
              stop 13
          elseif(iradbnd_in_itype.eq.1) then
c               >>>>>> CHANGED THIS ON 08-11-00 <<<<<<
c
c             BOUNDARY CONDITION OF BLACK BALL FILLING CENTRAL SPACE
c
              if((ir_old.eq.1).and.(icr_old.eq.1)) then
                  charintline = 0.d0
                  init = 1
              endif
          elseif(iradbnd_in_itype.eq.2) then 
c
c             BOUNDARY CONDITION OF VACUUM WITH CENTRAL STAR
c
              if(((ir.eq.1).or.(ir.eq.-1)).and.
     %             (istar_done.eq.0).and.(rbeam0.gt.0.d0)) then 
c
c                 If reemerging from the vacuum between star and inner 
c                 grid boundary then when rbeam0.gt.0 then we mix the 
c                 stellar flux with the present integrated intensity. 
c                 What actually happens is that for a small star, all 
c                 the rays will miss the star, except the single radial
c                 ray. In order to correctly include the star, the radial 
c                 ray must have a beam size, so that the relative 
c                 contribution of the star to the intensity can be 
c                 estimated. This is done with the rbeam0. 
c                 Of course, this should only be rbeam0.gt.0 for iray.eq.1. 
c                 Rbeam0 is the `impact parameter' of the edge of the beam.
c             
c                 REMARK: What should we do when the beam resolves the star?
c                         This seems not to be implemented yet...
c         
                  if(rbeam0.lt.radbnd_rstar) then
                      write(*,*) 'PROBLEM: The central beam of the',
     %                   'circular CCD has a size smaller than the'
                      write(*,*) '         stellar radius. This ',
     %                   'probably means that you specified the '
                      write(*,*) '         inner radius of the grid ',
     %                   'of the same (general) order of the stellar'
                      write(*,*) '         radius. '
                      write(*,*) '         While this is in principle ',
     %                   'okay, the present version of the telescope'
                      write(*,*) '         routine cannot deal with ',
     %                   'a situation in which the star is'
                      write(*,*) '         resolved. '
                      write(*,*) '         TIP: Choose in radical.inp ',
     %                   'the number of extra b_i '
                      write(*,*) '         inside the inner radius to ',
     %                   'a smaller value (worst case: choose 0).'
                      write(*,*) '         R_beam = ',rbeam0
                      write(*,*) '         R_*    = ',radbnd_rstar
                      stop 124
                  endif
                  starfract  = (radbnd_rstar/rbeam0)**2
                  charintline    = ( 1.d0 - starfract ) * charintline
     %                + starfract * radbnd_starspec_line(inu)
                  istar_done = 1
                  init = 1

              endif
              if((ir.eq.1).and.(ir_old.eq.1).and.(rbeam0.eq.0.d0)
     %             .and.(tr_b.le.radbnd_rstar)) then 

c
c                 Include the central star in the normal way, i.e.
c                 if the ray emerges from the central hole, and if
c                 it has an impact parameter smaller than the stellar
c                 radius, then it must have hit the star, and we put
c                 its intensity to the stellar intensity.
c                 (added 20 June 2000)
c
                  charintline    = radbnd_starspec_line(inu)
                  init = 1
c
              endif
          else
              write(*,*) 'Dont know this type of inner bc'
              stop 13
          endif
   77     continue
c
c         Consistency check
c
          if(ir.eq.0) then
              write(*,*) 'IR=0!! rbeam0=',rbeam0
              stop 13
          endif
c
c         Find the new src and alpha
c         
#ifdef INTERPOL_SRC_3 
          call get_line_dust_values_2nd(icr,ir,it,dr,dt,mu,phi,iline,
     %                 inu,srcd1,alpd1,dvmu1,lw1,Nup1,Ndown1)
#endif
#ifdef INTERPOL_SRC_1
          call get_line_dust_values(icr,ir,it,dr,dt,mu,phi,iline,inu,
     %                 srcd1,alpd1,dvmu1,lw1,Nup1,Ndown1)
#endif
c
c         Check whether the line crosses the ray for any frequency within the
c         current frequency interval.
c
c$$$          ilfr  = line_ilifreq(inu)
c$$$          dnu1  = line_dnu(ilfr,iline) - line_nu0(iline) * dvmu1
c$$$          if(ilfr.eq.1) then
c$$$             bmin  = (line_dnu(ilfr,iline)+line_dnu(ilfr,iline))/2.d0
c$$$          else
c$$$             bmin  = (line_dnu(ilfr,iline)+line_dnu(ilfr-1,iline))/2.d0
c$$$          endif
c$$$          if(ilfr.eq.freq_nr) then
c$$$             bmax  = (line_dnu(ilfr,iline)+line_dnu(ilfr,iline))/2.d0
c$$$          else
c$$$             bmax  = (line_dnu(ilfr+1,iline)+line_dnu(ilfr,iline))/2.d0
c$$$          endif
c$$$          if(dnu1.gt.bmin.and.dnu1.lt.bmax.and.iray.eq.1000) 
c$$$     %         write(*,*) bmin,dnu1,bmax,ilfr,dnu1*3d5/line_nu0(iline)
c$$$c     MARKER: ARF
c     
c         Now we can actually perform the formal
c         integral over this small piece of local
c         characteristic (r,theta and is are passed for debugging purposes)
c         
          call clever_integrate_element_linedust(charintline,ds,
     %         srcd0,srcd1,alpd0,alpd1,
     %         lw0,lw1,dvmu0,dvmu1,
     %         nup0,nup1,ndown0,ndown1,
     %         char_tau,char_emis,
     %         inu,iline,init,
     %         phiprof0,srcl0,alpl0,iray,
     %         r,theta,tr_s(is)) 

          if (dvmu0.lt.minvel(iray).and.nup0+ndown0.gt.LEVTHRES) then
             minvel(iray) = dvmu0             
          endif
          if (dvmu0.gt.maxvel(iray).and.nup0+ndown0.gt.LEVTHRES) then
             maxvel(iray) = dvmu0
          endif
          
          
c     
c         New is old for the dust, for the next step.
c         
          srcd0 = srcd1
          alpd0 = alpd1
c
c         Now swap also for the line stuff
c
          lw0    = lw1
          dvmu0  = dvmu1
          nup0   = nup1
          ndown0 = ndown1 
c
c         Note: for the line stuff the swapping for
c               phiprof0,srcl0,alpl0
c               is already done in the above subroutine.
c         

c
c     Check for masing (or lasing action) - this should not happen
c     - and if it does, the user should check the input level populations
c    
          IF (alpl0*ds.lt.-0.01) THEN
             maserflag = 1
          ENDIF

#         ifdef MONITOR_CHARINT
             if (inu.eq.2.or.inu.eq.freq_nr/2) then
		        write(4) r, theta, charintline
			 endif
#         endif

      enddo
c
c     If the density never gets above the threshold, make sure that
c     the ray is always counted as continuum.     
c 

      return
      end



c     --------------------------------------------------------
c          CHECK THE GRID SPACING VERSUS SAFETY NUMBERS
c
c     If the grid is too finely spaced in R, then the safety
c     numbers may have to be reduced. This routine checks if
c     the safety numbers are small enough to be safe.
c     --------------------------------------------------------
      subroutine telescope_check_safety_numbers()
      implicit none
c
#include "common_grid.h"
c
      integer ir
c
      do ir=1,irsi_frsizex-1
          if(rsi_x_c(ir+1,1)/rsi_x_c(ir,1)-1.d0.lt.1.d4*TELESC_EPS) then
              write(*,*) 'PROBLEM in telescope:'
              write(*,*) 'The radial grid resolution is too fine for'
              write(*,*) 'the safety number TELESC_EPS=',TELESC_EPS
              write(*,*) 'Tip: reduce TELESC_EPS, but not too small'
              write(*,*) 'since then machine precision errors can'
              write(*,*) 'occur.'
              stop 13
          endif
      enddo
c
      end


 
      
      

