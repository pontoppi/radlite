c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c                         ( MODULE: CREATE GRID )
c
c                          Muenchen, July 2000
c                             C.P. Dullemond   
c
c     =============================================================


#define GRID_F_VERSION "Trans_31Jul00"
#include "main.h"

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


c     --------------------------------------------------------------
c                          CREATE THE ANGULAR GRID 
c
c     This subroutine creates a grid in the angles mu and phi, at
c     each radius/position. The mu-gridding can be done in two 
c     different ways (see type below). If type=0, then the mu gridding
c     is regular. If type=1, then the mu-gridding is adapted to the
c     gridding in radius R (for spherical coordinates). 
c
c     The amount of gridpoints is determined by anggrid_frsizemu and
c     anggrid_frsizephi. 
c
c     IMPORTANT: If anggrid_frsizemu is an even number, then there will
c                be no gridpoint at mu=0, but if anggrid_frsizemu is an 
c                odd number, then there exists one gridpoint at mu=0.d0.
c                If there is no mu=0 gridpoint (anggrid_frsizemu even)
c                then the array element rmu(0) has no meaning. It will
c                be filled with the number 0.0, but that is meaningless.
c                In the intensity array the value of intens(*,0,*,*) will
c                therefore also be meaningless and should be skipped.
c
c     REQUIRES:
c       - The spatial gridding must be made before this call
c
c     USES: 
c       anggrid_drr       = The Delta R / R of the radial grid
c       anggrid_frsizemu  = The total amount of mu-gridpoints from mu=-1 to +1
c                           If this number is odd, then the middle mu(0)=0.d0
c                           while if it is even, then there exists no middle
c                           mu, and the mu(0)=0.d0 is ignored by the transfer
c                           routines.
c       anggrid_frsizephi = The total amount of phi-gridpoints between 0,2*pi
c
c     ARGUMENTS:
c       type            = The type of discretization of mu used:
c                         0 : Regular grid in mu
c                         1 : Grid in mu is adapted to grid in radius
c                         2 : Same as 1, but now the amount of mu angles
c                             (nrmu) will be such that the expected error
c                             (due to mu-binning at mu=+1,-1) will be kept
c                             below anggrid_muerr_max.
c       iangset         = The index of the set characteristics belonging
c                         to this value of Delta R / R
c 
c     OUTPUT: 
c       nrmu(iangset)   = The amount of actual mu gridpoints for mu>0. 
c                         This means that including mu=0 we have 
c                         2*nrmu+1 angles in mu. If mu=0 is not a gridpoint
c                         (i.e. if anggrid_frsizemu is even number) then
c                         of course this is 2*nrmu instead of 2*nrmu+1.
c       rmu(i,iangset)  = Array of mu griddings at each R_ix. The i index
c                         runs from -nrmu(ix) to nrmu(ix). If anggrid_frsizemu
c                         is an even number, then the i=0 index represents
c                         not a gridpoint. This array element will then be
c                         skipped by all routines.
c       rphi(i,iangset) = Array of phi griddings at each R_ix. 
c       anggrid_mu_zero = This is a flag that tells the rest of the program
c                         whether or not the amount of mu's is odd, so that
c                         (if odd) there exists a mu(0) = 0.d0. If it is
c                         even then the mu(0) is non-existent (ignored by
c                         the program).
c
c     --------------------------------------------------------------
      subroutine create_anggrid(type,iangset)
      implicit none
      integer iangset,type
c
#include "common_grid.h"
c
      doubleprecision sinth,bi,b,costh,phi,dphi,drr,mu1
      doubleprecision drr0,dummy,error,muibk
      integer imu,iphi,idumm
c
      drr0 = anggrid_drr
      if(anggrid_frsizemu.gt.FRSIZE_MU) then
          write(*,*) 'Too many angles in mu requested!'
          stop 13
      endif
      if(anggrid_frsizephi.gt.FRSIZE_PHI) then
          write(*,*) 'Too many angles in phi requested!'
          stop 13
      endif
c
      nrmu(iangset)  = anggrid_frsizemu / 2
      nrphi(iangset) = anggrid_frsizephi 
      if(anggrid_frsizemu/2.ne.(anggrid_frsizemu+1)/2) then
          anggrid_mu_zero = 1
      else
          anggrid_mu_zero = 0
      endif
c
      if(iextrmu(iangset).gt.nrmu(iangset)-1) then
          write(*,*) 'Too many extra mu requested around mu=0!'
          stop 13
      endif
c
c     Now select the type of mu gridding, and make the mu grid.
c
      if((type.eq.1).or.(type.eq.2)) then
c
c        The mu-grid is adapted to the R-grid
c
#ifdef COORD_CARTESIAN
          write(*,*) 'This code is compiled for Cartesian coordinates'
          write(*,*) 'Mu-grid type 1 is only usable for spherical',
     %      'coordinates. Use type 0 instead.'
          stop 13
#endif
c
c         For now, we demand the anggrid_frsizemu to be odd (i.e. we
c         demand there to be a mu(0)=0.d0 gridpoint).
c
          if(anggrid_mu_zero.eq.0) then
              write(*,*) 'For spherical coordinates we demand the ',
     %              'anggrid_frsizemu to be odd.'
              write(*,*) 'In other words: we demand there to be ',
     %              'a mu(0)=0.0 gridpoint.'
              write(*,*) 'Aborting... Sorry for the inconvenience.'
              write(*,*) 'Rerun with odd value in the input file'
              stop 13
          endif
c
c         Check if the mu-resolution not too high. If so, this would
c         cause (mild) troubles with the mu-binning, since this does
c         not expect to find more than 1 mu bin in range for the flux
c         bin. See transfer.F.
c
          if(dmuvdr.lt.1.d0) then
              write(*,*) 'PROBLEM: At the moment, the algorithm is not'
              write(*,*) '         yet suited for higher mu-resolution'
              write(*,*) '         than the tangent-ray configuration.'
              write(*,*) '         Put the Dmu wrt dR in radlite.inp'
              write(*,*) '         to 1 (tangent-ray) or higher (less'
              write(*,*) '         mu-resolution).'
              stop 13
          endif
c
c         Write a warning if coarse mu-resolution
c
          if(dmuvdr.gt.1.d0) then
              write(*,*) '    Choosing coarse mu-resolution'
              write(*,*) '    NOTE: This may (but not always does)'
              write(*,*) '          cause flux cons problems.'
          endif
c
c         Find the drr of the mu-grid. If dmuvdr.eq.1, then we have the
c         perfect tangent-ray method. If dmuvdr.gt.1, then we have a
c         coarser mu-resolution. If dmuvdr is 2,3,4 etc, then the mu
c         grid is tangent to every 2nd, 3rd, 4th etc grid circle. 
c
          drr = ((1.d0+drr0)**dmuvdr)-1.d0
          rmu(0,iangset) = 0.d0
c
c         First the extra mu's around mu=0
c    
          sinth = 1.d0 / ( 1.d0 + drr )
          costh = sqrt(1.d0-sinth*sinth)
          mu1   = costh
          do imu=1,iextrmu(iangset)
              rmu(imu,iangset)  =  imu*mu1/(iextrmu(iangset)+1)
              rmu(-imu,iangset) = -imu*mu1/(iextrmu(iangset)+1)
          enddo
c    
c         Now the regular ones (the 'tangent rays')
c    
          b           = 1.d0
          bi          = 1.d0 * sqrt( 1.d0 + drr)
c
ccccccccc          do imu=iextrmu(iangset)+1,nrmu(iangset)-1
c
c         Do a loop over imu. This is done in an oldfashioned way
c         because of the type.eq.2 mode possibility (see below).
c         
          imu         = iextrmu(iangset)+1
  611     continue
c         
c         Find the impact parameter b (for R=1)
c         
          b  = b / ( 1.d0 + drr )
          bi = bi / ( 1.d0 + drr )
c         
c         Store the mu angles
c          
          sinth  = b
          costh  = sqrt(1.d0-sinth*sinth)
          rmu(imu,iangset)  =  costh
          rmu(-imu,iangset) = -costh
          sinth  = bi
          costh  = sqrt(1.d0-sinth*sinth)
          rmu_i(imu,iangset)  =  costh
          rmu_i(-imu,iangset) = -costh
c         
c         Check if the present rmu is not too close to 1.0 .
c         If so, then ignore this mu, and consider previous mu 
c         as the last. Also reset the bi, for below. 
c         (added on 05-01-01)
c
          if(1.d0-rmu(imu,iangset).lt.1.d-3) then
              write(*,*) 'Too many mu-angles for accurate RT (mu too '
              write(*,*) '   close to 1). Using now nmu=',imu
              nrmu(iangset) = imu
              anggrid_frsizemu = 2*nrmu(iangset) + 1
              b  = b * ( 1.d0 + drr )
              bi = bi * ( 1.d0 + drr )
              goto 612
          endif
c
c         Continue loop or exit
c
          imu = imu + 1
          if(imu.le.nrmu(iangset)-1) then 
             goto 611
          endif
c
c         If the mu grid type is 2, then the amount of mu gridpoints
c         must be determined automatically here. This is done by
c         estimating the error made by choosing a certain amount of
c         mu-gridpoints, and demanding that the error is small enough.
c         If type.eq.2 then we actually wish to go on adding
c         mu-points (and correspondingly increase nrmu(iangset))
c
          if(type.eq.2) then
              dummy = rmu(nrmu(iangset)-1,iangset)
              error = abs(2.d0*abs(1.d0-dummy)/(1.d0-(dummy**2))-1.d0) 
              if(error.gt.anggrid_muerr_max) then
                  if(nrmu(iangset).lt.FRSIZE_MU_HALF) then
                      write(*,*) '    Adding mu to angle-grid. nrmu = ',
     %                    nrmu(iangset)+1
                      nrmu(iangset) = nrmu(iangset) + 1
                      anggrid_frsizemu = 2*nrmu(iangset) + 1
                      goto 611
                  else
                      write(*,*) '    WARNING: Could not add ',
     %                  'sufficient mu-angles to get error in check...'
                     write(*,*) '    Tip: Perhaps you have a very ',
     %                  'high spatial resolution Dr/r. The code will '
                     write(*,*) '         then need many mu angles ',
     %                  'as well... Solution: increase the value of'
                     write(*,*) '         FRSIZE_MU in the Makefile ',
     %                  'and recompile the code.'
                     stop 13
                  endif
              endif
          endif
c
 612      continue
c
c         Message
c
          dummy = rmu(nrmu(iangset)-1,iangset)
          error = abs(2.d0*abs(1.d0-dummy)/(1.d0-(dummy**2))-1.d0) 
          write(*,'(a,f8.4,a,e9.4)') '    Max |mu| = ',
     %         rmu(nrmu(iangset)-1,iangset),
     %         '  --> Error = ',error
          write(*,'(a,f8.4)') '    Min |mu| = ',rmu(1,iangset)
c
c         Now we finish the job by adding the mu=+1 and mu=-1 
c         gridpoint and all the interface mu's. 
c     
c         Middle mu_i=0
c    
          rmu_i(0,iangset)  = 0.d0
c    
c         The fore-last mu_i
c         
          bi = bi / ( 1.d0 + drr )
          sinth  = bi
          costh  = sqrt(1.d0-sinth*sinth)
          rmu_i(nrmu(iangset),iangset)  =  costh
          rmu_i(-nrmu(iangset),iangset) = -costh      
c         
c         The last angle is mu=+-1
c         
          rmu(nrmu(iangset),iangset)  =  1.d0
          rmu(-nrmu(iangset),iangset) = -1.d0
          rmu_i(nrmu(iangset)+1,iangset)  =  1.d0
          rmu_i(-nrmu(iangset)-1,iangset) = -1.d0
c    
c         Now the intermediate mu_i angles
c         
c          (FIX 08-01-01: special treatment of rmu_i(iextra+1).
c           First store muibk, then check if special treatment
c           possible, if so then do)
c
          muibk = rmu_i(iextrmu(iangset)+1,iangset)
          do imu=1,iextrmu(iangset)+1
              rmu_i(imu,iangset)  = 0.5 * 
     %          ( rmu(imu,iangset)  + rmu(imu-1,iangset) )
              rmu_i(-imu,iangset) = 0.5 * 
     %          ( rmu(-imu,iangset) + rmu(-imu+1,iangset) )
          enddo      
c
c         08-01-01: Special treatment
c
          imu = iextrmu(iangset) + 1
          if(muibk.gt.rmu(iextrmu(iangset),iangset)) then
              write(*,*) 'Note: special treatment of mu_i angles...'
              rmu_i(imu,iangset)  = muibk
              rmu_i(-imu,iangset) = -muibk
          endif
c
      elseif(type.eq.0) then
c
c         Regular mu-grid
c
          rmu(0,iangset) = 0.d0
          do imu=1,nrmu(iangset)
              rmu(imu,iangset)  =  imu/(1.0*nrmu(iangset))
              rmu(-imu,iangset) = -imu/(1.0*nrmu(iangset))
          enddo
c    
c         Middle mu_i=0
c    
          rmu_i(0,iangset)  = 0.d0
c    
c         Now the intermediate mu_i angles
c         
          do imu=1,nrmu(iangset)
              rmu_i(imu,iangset)  = 0.5 * 
     %          ( rmu(imu,iangset)  + rmu(imu-1,iangset) )
              rmu_i(-imu,iangset) = 0.5 * 
     %          ( rmu(-imu,iangset) + rmu(-imu+1,iangset) )
          enddo          
          rmu_i(nrmu(iangset)+1,iangset) = 1.0
          rmu_i(-nrmu(iangset)-1,iangset) = -1.0
c
c         Make correction if there exists nu mu(0)
c
          if(anggrid_mu_zero.eq.0) then
              rmu_i(1,iangset) = 0.d0
              rmu_i(-1,iangset) = 0.d0
          endif
      else
          write(*,*) 'Do not know mu-grid type ',type
          stop 13
      endif
c     
c     Next the phi angles (they have equal spacing)
c     Note that we add two ghostpoints on each side, just for 
c     convenience in the 3-rd order interpolation algorithms
c
#ifdef RADGRID_ONEDIM
      rphi(1,iangset) = 0.d0
      rphi(2,iangset) = 0.d0                ! Just in case...
      rphi_i(1,iangset) = 3.14159265359d0 
      rphi_i(0,iangset) = 0.d0              ! (Added 16-04-02)
      if(anggrid_frsizephi.ne.1) then
          write(*,*) '( 1-D mode: putting frsizephi=1 )'
          anggrid_frsizephi = 1
          nrphi(iangset) = 1
      endif
#else
      if(nrphi(iangset).eq.1) then
          write(*,*) 'ERROR: In 2-D mode I cannot accept '
          write(*,*) '       only 1 gripoint in phi.'
          stop 13
      endif
      if((nrphi(iangset)+3)/4.ne.(nrphi(iangset)/4)) then
          write(*,*) 'PROBLEM: In 2-D mode I would strongly '
          write(*,*) '       advise to have a nr of phi-points'
          write(*,*) '       that is divisible by 4.'
          write(*,*) '       Until possible later versions I will'
          write(*,*) '       now enforce this for safety..... Sorry.'
          stop 13
      endif
      dphi = 6.28318530718d0 / nrphi(iangset)
      phi  = -1.5d0 * dphi 
      do iphi=-1,nrphi(iangset)+2
          rphi(iphi,iangset) = phi
          phi = phi + dphi
      enddo
c     
c     Next the phi_i interfaces (Used only in the vetmethod.F module)
c     (BUGFIX: 16-04-02)
c     
      do iphi=1,nrphi(iangset)-1
          rphi_i(iphi,iangset) = 0.5 * 
     %            ( rphi(iphi,iangset) + rphi(iphi+1,iangset) )
      enddo
      rphi_i(0,iangset) = 0.5 * 
     %            ( rphi(nrphi(iangset),iangset) + rphi(1,iangset) 
     %              - 6.28318530718 )
      rphi_i(nrphi(iangset),iangset) = 0.5 * 
     %            ( rphi(nrphi(iangset),iangset) + rphi(1,iangset) 
     %              + 6.28318530718 )
#endif
c
c     Now a safety check...
c
      dummy = rmu(nrmu(iangset)-1,iangset)
      error = abs(2.d0*abs(1.d0-dummy)/(1.d0-(dummy**2))-1.d0) 
      if(anggrid_muerr_max.gt.0.d0) then 
          dummy = anggrid_muerr_max
      else
          dummy = 0.20d0
      endif
      if(error.gt.dummy) then
          write(*,*) '--------------------------------------'
          write(*,*) 'WARNING: Insufficient nr of mu-points.'
          write(*,*) 'The expected error is ',error,' which is'
          write(*,*) 'larger than the allowed',dummy
          write(*,*) 'Tip: Do you have a small dr/r '
          write(*,*) '     radial grid spacing? --> Take'
          write(*,*) '     more mu-points. The finer the'
          write(*,*) '     spacing in R, the larger the '
          write(*,*) '     nr of required mu points.'
          write(*,*) '--------------------------------------'
          iang_warn_mureslow = 1
      endif
cc      if(error.gt.0.1) then
cc          write(*,*) ' !!! PROBLEM: Too coarse mu-coverage !!!'
cc          iang_warn_mureslow = 1
cc      else
cc          iang_warn_mureslow = 0
cc      endif
#ifdef VERBOSE
      dummy = rmu(nrmu(iangset)-1,iangset)
      error = abs(2.d0*abs(1.d0-dummy)/(1.d0-(dummy**2))-1.d0) 
      write(*,*) 'Angle-Coverage:'
      write(*,*) '    Number of mu points: ',nrmu(iangset) 
      write(*,*) '    Highest mu  = ',dummy,
     %      ' ---> error = ',error
      dummy  = 0.d0
      idumm  = 0
      do imu=-nrmu(iangset),nrmu(iangset)-1
          if(abs(rmu(imu+1,iangset)-rmu(imu,iangset)).ge.dummy) then
              dummy  = abs(rmu(imu+1,iangset)-rmu(imu,iangset))
              idumm  = imu
          endif
      enddo
      write(*,*) '    Largest dmu = ',dummy, ' between mu = ',
     %     rmu(idumm,iangset),' and ',rmu(idumm+1,iangset)
#endif
c
      end


c     --------------------------------------------------------------
c                  COPY ANGULAR GRID TO OTHER IANGSET
c     --------------------------------------------------------------
      subroutine copy_anggrid(iangset1,iangset2)
      implicit none
c
      integer iangset1,iangset2
c
#include "common_grid.h"
c
      integer imu,iphi
c
      nrmu(iangset2)    = nrmu(iangset1)
      nrphi(iangset2)   = nrphi(iangset1)
      iextrmu(iangset2) = iextrmu(iangset1)
c
      do imu=-nrmu(iangset1),nrmu(iangset1)
          rmu(imu,iangset2)  = rmu(imu,iangset1)
      enddo
      do imu=-nrmu(iangset1)-1,nrmu(iangset1)+1
          rmu_i(nrmu(iangset1)+1,iangset2) = 
     %            rmu_i(nrmu(iangset1)+1,iangset1)
      enddo
      do iphi=-1,nrphi(iangset1)+2
          rphi(iphi,iangset2)   = rphi(iphi,iangset1)
      enddo
      do iphi=1,nrphi(iangset1)
          rphi_i(iphi,iangset2) = rphi_i(iphi,iangset1)
      enddo
c
      end


c     --------------------------------------------------------------
c                        CREATE THE GRIDPOINTS
c
c                [ THIS IS A COPY FROM THE ROE.F CODE, BUT 
c                  ADAPTED SPECIALLY FOR THE RADIATION ]
c
c     This routine sets up the arrays rsi_x_c,rsi_dx_c,rsi_x_i,rsi_dx_i.
c     These are used by the roe solver because they represent the
c     positions of the cell centers (_c) and cell interfaces (_i) 
c     in radius and in theta (resp (..,1) and (..,2)). Without
c     these arrays the data F0(ix,iy) have no meaning because the
c     ix and iy don't have meaning. 
c
c     Note that ix=0 and ix=irsi_frsizex+1 are also set by this 
c     routine and similar for iy=0 and iy=irsi_frsizey+1. These are
c     the `ghost cells' used for easy implementation of boundary 
c     conditions. For the y-direction this means that iy=0 represents 
c     a cell on the opposite side of the z-axis, and therefore does 
c     not represent any physical location.
c
c       x_c   = Position of cell centers
c       x_i   = Position of cell interfaces
c       dx_c  = Distance between two cell centers (defined at cell interface)
c       dx_i  = Distance between two cell interfaces (defined at cell center)
c
c     --------------------------------------------------------------
      subroutine create_spacegrid()
      implicit none
c
#include "common_grid.h"
c
      common/gridlinref/grsol_x0,grsol_x1,grsol_xr,grsol_a
      doubleprecision grsol_x0,grsol_x1,grsol_xr,grsol_a
      common/igridlinref/grsol_N,grsol_M
      integer grsol_N,grsol_M
c
      doubleprecision pi
      parameter(pi=PICONST)
c
      integer ix,iy
      integer itmax,imirt
      character*1024 string 
      doubleprecision dth,xpr,xprh,rr,thcur,dummy,thmaxerror
      logical yesouter
      external grid_linear_refine
      doubleprecision grid_linear_refine
c
c
c     First do the radial direction
c
      if((spacegrid_type_x.eq.'L').or.
     %   (spacegrid_type_x.eq.'l')) then
          write(*,*) 'ERROR: L and l type grid not allowed in RADICAL'
          stop 13
      else if(spacegrid_type_x.eq.'O') then
          xpr  = ( spacegrid_xi + spacegrid_dxi ) / spacegrid_xi
          xprh = sqrt(xpr)
c
c         Check whether spacegrid_xo is set. If so, then the outer edge
c         is roughly spacegrid_xo, and the irsi_frsizex is determined
c         from this. 
c
          if(spacegrid_xo.gt.spacegrid_xi) then
              irsi_frsizex = FRSIZE_X
              yesouter = .true.
          else
              yesouter = .false.
          endif
c
c         The cell centers
c
          rr  = spacegrid_xi / ( xpr * xprh )
          do ix=-1,irsi_frsizex+2
              rsi_x_c(ix,1)  = rr
              rr             = rr * xpr
              if(yesouter.and.(ix.gt.2)) then
                  if(rsi_x_c(ix-2,1).gt.spacegrid_xo) goto 555
                  if(ix.eq.irsi_frsizex+2) goto 555
              endif
          enddo
 555      continue
          if(yesouter) then
              irsi_frsizex = ix-2
              if(rsi_x_c(irsi_frsizex,1).lt.spacegrid_xo) then
                  write(*,*) 'ERROR: Trying to make spatial grid '
                  write(*,*) '       in between ',spacegrid_xi,
     %                       ' and ',spacegrid_xo
                  write(*,*) '       with dr/r= ',spacegrid_dxixi
                  write(*,*) '       Problem: need more gridpoints '
                  write(*,*) '       than allowed by array dimensions.'
                  write(*,*) '       Increase FRSIZE_X in configure.h'
                  stop 13
              endif
          endif
c
c         The cell interfaces
c
          rr  = spacegrid_xi / xpr
          do ix=0,irsi_frsizex+2
              rsi_x_i(ix,1)  = rr
              rr             = rr * xpr
          enddo
c
c         Now the dx_c (distance between cell centers)
c          
          do ix=0,irsi_frsizex+2
              rsi_dx_c(ix,1) = rsi_x_c(ix,1) - rsi_x_c(ix-1,1)
          enddo
          rsi_dx_c(irsi_frsizex+2,1) = 0.d0
c
c         and the dx_i (distance between cell interfaces)
c          
          do ix=0,irsi_frsizex+1
              rsi_dx_i(ix,1) = rsi_x_i(ix+1,1) - rsi_x_i(ix,1)
          enddo
c
c         Now the dx_cir (distance between cell center and its right cell interface)
c          
          do ix=0,irsi_frsizex+1
              rsi_dx_cir(ix,1) = rsi_x_i(ix+1,1) - rsi_x_c(ix,1)
          enddo
          rsi_dx_cir(irsi_frsizex+2,1) = 0.d0
c
c         Now the dx_cil (distance between cell center and its left cell interface)
c          
          do ix=0,irsi_frsizex+1
              rsi_dx_cil(ix,1) = rsi_x_c(ix,1) - rsi_x_i(ix,1)
          enddo
          rsi_dx_cil(irsi_frsizex+2,1) = 0.d0
c
      elseif(spacegrid_type_x.eq.'-') then
c
c         Read radial grid from file radius.inp, and make sure that
c         it is indeed logaritmically spaced.
c
          open(unit=3,file='radius.inp',status='old')
          read(3,*) irsi_frsizex
          if(irsi_frsizex.gt.FRSIZE_X) then
              write(*,*) 'ERROR: radius.inp has more radial gridpoints'
              write(*,*) '       than FRSIZE_X. '
              write(*,*) '       Recompile with larger FRSIZE_X'
              stop 13
          endif
          do ix=1,irsi_frsizex
              read(3,*) rsi_x_c(ix,1)
          enddo
          close(3)
          spacegrid_radius_read = 1
c
c         Check if log spaced 
c
#         ifndef NO_SHORT_CHARS
#         ifndef SHORTCHAR_NEWSTYLE
          xpr = rsi_x_c(2,1) / rsi_x_c(1,1)
          do ix=3,irsi_frsizex
              if(abs(rsi_x_c(ix,1)/rsi_x_c(ix-1,1)-xpr).gt.1.d-3) then
                  write(*,*) 'ERROR: Currently RADICAL is switched in'
                  write(*,*) '       old style Short Characteristic '
                  write(*,*) '       mode. When the SCs are switched on'
                  write(*,*) '       this old mode does not allow any'
                  write(*,*) '       other than equally spaced log R'
                  write(*,*) '       grid in radius. '
                  write(*,*) '       The radii provided in radius.inp'
                  write(*,*) '       are not exactly logarthmically'
                  write(*,*) '       spaced. '
                  write(*,*) '       When you switch on the define '
                  write(*,*) '       SHORTCHAR_NEWSTYLE, you have more'
                  write(*,*) '       freedom: then also factors 1/2^n'
                  write(*,*) '       in grid refinement in log(R) are'
                  write(*,*) '       allowed. '
                  write(*,*) '       If you use RADICAL only for ray'
                  write(*,*) '       tracing, then switch on the '
                  write(*,*) '       ONLY_RAY_TRACING macro, and then'
                  write(*,*) '       any R-grid is allowed!'
                  stop 13
              endif
          enddo
#         endif
#         endif
c
c         Set spacegrid_xi etc.
c
cc        (BUGFIX 04-09-00):
cc          spacegrid_xi  = rsi_x_c(1,1)
cc          spacegrid_dxi = rsi_x_c(2,1) - rsi_x_c(1,1)
cc          xpr  = ( spacegrid_xi + spacegrid_dxi ) / spacegrid_xi
cc          xprh = sqrt(xpr)
          xpr  = rsi_x_c(2,1) / rsi_x_c(1,1)
          xprh = sqrt(xpr)
          spacegrid_xi  = rsi_x_c(1,1) / xprh
          spacegrid_dxi = rsi_x_c(1,1) * xprh - rsi_x_c(1,1) / xprh
c     
c         The boundaries
c
          rsi_x_c(0,1)  = rsi_x_c(1,1)**2 / rsi_x_c(2,1)
          rsi_x_c(-1,1) = rsi_x_c(0,1)**2 / rsi_x_c(1,1)
          rsi_x_c(irsi_frsizex+1,1) = rsi_x_c(irsi_frsizex,1)**2 
     %                                   / rsi_x_c(irsi_frsizex-1,1)
          rsi_x_c(irsi_frsizex+2,1) = rsi_x_c(irsi_frsizex+1,1)**2 
     %                                   / rsi_x_c(irsi_frsizex,1)
c
c         The cell interfaces
c
cc        (MODIFIED 04-09-00)
cc          rr  = spacegrid_xi / xpr
cc          do ix=0,irsi_frsizex+2
cc              rsi_x_i(ix,1)  = rr
cc              rr             = rr * xpr
cc          enddo
          do ix=0,irsi_frsizex+2
              rsi_x_i(ix,1)  = sqrt(rsi_x_c(ix,1)*rsi_x_c(ix-1,1))
          enddo
c
c         Now the dx_c (distance between cell centers)
c          
          do ix=0,irsi_frsizex+2
              rsi_dx_c(ix,1) = rsi_x_c(ix,1) - rsi_x_c(ix-1,1)
          enddo
          rsi_dx_c(irsi_frsizex+2,1) = 0.d0
c
c         and the dx_i (distance between cell interfaces)
c          
          do ix=0,irsi_frsizex+1
              rsi_dx_i(ix,1) = rsi_x_i(ix+1,1) - rsi_x_i(ix,1)
          enddo
c
c         Now the dx_cir (distance between cell center and its right cell interface)
c          
          do ix=0,irsi_frsizex+1
              rsi_dx_cir(ix,1) = rsi_x_i(ix+1,1) - rsi_x_c(ix,1)
          enddo
          rsi_dx_cir(irsi_frsizex+2,1) = 0.d0
c
c         Now the dx_cil (distance between cell center and its left cell interface)
c          
          do ix=0,irsi_frsizex+1
              rsi_dx_cil(ix,1) = rsi_x_c(ix,1) - rsi_x_i(ix,1)
          enddo
          rsi_dx_cil(irsi_frsizex+2,1) = 0.d0
      else
          write(*,*) 'Other radial grid types than O, L or l ',
     %     'not yet implemented'
          write(*,*) 'ABORTING'
          stop 13
      endif
c
c     Next the theta direction
c
      if(spacegrid_type_y.eq.'l') then
          write(string,*) 'Theta is done using old grid spacing method'
c          call writelog(string)
          write(*,*) 'Theta is done using old grid spacing method'
c
c         Cell centers
c          
          do iy=-1,irsi_frsizey+2
              rsi_x_c(iy,2)  = ( 1.d0*iy - 0.5d0 ) 
          enddo
c
c         Cell interfaces
c
          do iy=0,irsi_frsizey+2
              rsi_x_i(iy,2)  = ( 1.d0*iy - 1.0d0 ) 
          enddo
c
c         Distance between cell centers
c
          do iy=0,irsi_frsizey+2
              rsi_dx_c(iy,2) = rsi_x_c(iy,2) - rsi_x_c(iy-1,2)
          enddo
          rsi_dx_c(irsi_frsizey+2,2) = 0.d0
c
c         Distance between cell interfaces
c
          do iy=0,irsi_frsizey+1
              rsi_dx_i(iy,2) = rsi_x_i(iy+1,2) - rsi_x_i(iy,2)
          enddo
c
c         Now the dx_cir (distance between cell center and its right cell interface)
c          
          do iy=0,irsi_frsizey+1
              rsi_dx_cir(iy,2) = rsi_x_i(iy+1,2) - rsi_x_c(iy,2)
          enddo
          rsi_dx_cir(irsi_frsizey+2,2) = 0.d0
c
c         Now the dx_cil (distance between cell center and its left cell interface)
c          
          do iy=0,irsi_frsizey+1
              rsi_dx_cil(iy,2) = rsi_x_c(iy,2) - rsi_x_i(iy,2)
          enddo
          rsi_dx_cil(irsi_frsizey+2,2) = 0.d0
      elseif(spacegrid_type_y.eq.'L') then
c
c         Regularly spaced Theta, with half spacing at North and South pole,
c         and without Theta=Pi/2 (EVEN irsi_frsizey required!)
c
          dth = pi / irsi_frsizey
c     
c         Cell centers
c          
          do iy=-1,irsi_frsizey+2
              rsi_x_c(iy,2)  = ( 1.d0*iy - 0.5d0 ) * dth
          enddo
c         
c         Cell interfaces
c         
          do iy=0,irsi_frsizey+2
              rsi_x_i(iy,2)  = ( 1.d0*iy - 1.0d0 ) * dth
          enddo
c         
c         Distance between cell centers
c         
          do iy=0,irsi_frsizey+2
              rsi_dx_c(iy,2) = rsi_x_c(iy,2) - rsi_x_c(iy-1,2)
          enddo
          rsi_dx_c(irsi_frsizey+2,2) = 0.d0
c         
c         Distance between cell interfaces
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_i(iy,2) = rsi_x_i(iy+1,2) - rsi_x_i(iy,2)
          enddo
c         
c         Now the dx_cir (distance between cell center and its right cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cir(iy,2) = rsi_x_i(iy+1,2) - rsi_x_c(iy,2)
          enddo
          rsi_dx_cir(irsi_frsizey+2,2) = 0.d0
c         
c         Now the dx_cil (distance between cell center and its left cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cil(iy,2) = rsi_x_c(iy,2) - rsi_x_i(iy,2)
          enddo
          rsi_dx_cil(irsi_frsizey+2,2) = 0.d0
      elseif(spacegrid_type_y.eq.'D') then
c
c         Regularly spaced Theta, with half spacing at North and South pole,
c         and with spacing of spacegrid_equator_eps between Theta gridline
c         and the equator at pi/2. The equator itself is not represented by
c         a coordinate (EVEN irsi_frsizey required!)
c
          if((irsi_frsizey+1)/2.ne.irsi_frsizey/2) then
              write(*,*) 'Must have even number of Theta gridpoints'
              stop 13
          endif
          if(abs(spacegrid_equator_eps).lt.1.e-7) then
              write(*,*) 'spacegrid_equator_eps not defined!'
              stop 13
          endif
          dth = ( pi - 2.d0*spacegrid_equator_eps ) / 
     %              (irsi_frsizey-1.d0)
c     
c         Cell centers
c          
          thcur = -1.5d0 * dth
          do iy=-1,irsi_frsizey/2
              rsi_x_c(iy,2)  = thcur
              thcur = thcur + dth
          enddo
          thcur = thcur - dth + 2*spacegrid_equator_eps
          do iy=irsi_frsizey/2+1,irsi_frsizey+2
              rsi_x_c(iy,2)  = thcur
              thcur = thcur + dth
          enddo
c         
c         Cell interfaces
c         
          thcur = - dth
          do iy=0,irsi_frsizey/2
              rsi_x_i(iy,2)  = thcur
              thcur = thcur + dth
          enddo
          rsi_x_i(irsi_frsizey/2+1,2)  = 0.5d0*pi
          thcur = pi - (irsi_frsizey/2-1)*dth
          do iy=irsi_frsizey/2+2,irsi_frsizey+2
              rsi_x_i(iy,2)  = thcur
              thcur = thcur + dth
          enddo
c         
c         Distance between cell centers
c         
          do iy=0,irsi_frsizey+2
              rsi_dx_c(iy,2) = rsi_x_c(iy,2) - rsi_x_c(iy-1,2)
          enddo
          rsi_dx_c(irsi_frsizey+2,2) = 0.d0
c         
c         Distance between cell interfaces
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_i(iy,2) = rsi_x_i(iy+1,2) - rsi_x_i(iy,2)
          enddo
c         
c         Now the dx_cir (distance between cell center and its right cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cir(iy,2) = rsi_x_i(iy+1,2) - rsi_x_c(iy,2)
          enddo
          rsi_dx_cir(irsi_frsizey+2,2) = 0.d0
c         
c         Now the dx_cil (distance between cell center and its left cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cil(iy,2) = rsi_x_c(iy,2) - rsi_x_i(iy,2)
          enddo
          rsi_dx_cil(irsi_frsizey+2,2) = 0.d0
c
      elseif(spacegrid_type_y.eq.'E') then
c
c         Read the Theta grid from file
c
#         ifdef RADGRID_ONEDIM
          write(*,*) 'In 1-D, no theta.inp read!'
          stop 13
#         endif
          open(unit=1,file='theta.inp',status='old',err=701)
          read(1,*) iy
          if(iy.ne.irsi_frsizey) then
              write(*,*) 'Theta size in theta.inp and in schar.inp do ',
     %          'not match'
              stop 13
          endif
          do iy=1,irsi_frsizey
              read(1,*) dummy
              rsi_x_c(iy,2) = dummy
          enddo
          close(1)
          spacegrid_theta_read = 1
c
          rsi_x_c(0,2)  = - rsi_x_c(1,2)
          rsi_x_c(-1,2) = - rsi_x_c(2,2)
          rsi_x_c(irsi_frsizey+1,2) = 2 * 3.1415926e0 
     %                  - rsi_x_c(irsi_frsizey,2)
          rsi_x_c(irsi_frsizey+2,2) = 2 * 3.1415926e0 
     %                  - rsi_x_c(irsi_frsizey-1,2)
c         
c         Cell interfaces
c         
          do iy=0,irsi_frsizey+2
              rsi_x_i(iy,2)  = 0.5d0 * ( rsi_x_c(iy-1,2) + 
     %                                   rsi_x_c(iy,2) )
          enddo
c         
c         Distance between cell centers
c         
          do iy=0,irsi_frsizey+2
              rsi_dx_c(iy,2) = rsi_x_c(iy,2) - rsi_x_c(iy-1,2)
          enddo
          rsi_dx_c(irsi_frsizey+2,2) = 0.d0
c         
c         Distance between cell interfaces
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_i(iy,2) = rsi_x_i(iy+1,2) - rsi_x_i(iy,2)
          enddo
c         
c         Now the dx_cir (distance between cell center and its right cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cir(iy,2) = rsi_x_i(iy+1,2) - rsi_x_c(iy,2)
          enddo
          rsi_dx_cir(irsi_frsizey+2,2) = 0.d0
c         
c         Now the dx_cil (distance between cell center and its left cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cil(iy,2) = rsi_x_c(iy,2) - rsi_x_i(iy,2)
          enddo
          rsi_dx_cil(irsi_frsizey+2,2) = 0.d0
c
      elseif(spacegrid_type_y.eq.'#') then
c
c         Theta grid is already specified, just recompute all
c         other data. Assume irsi_frsizey is set also.
c
          rsi_x_c(0,2)  = - rsi_x_c(1,2)
          rsi_x_c(-1,2) = - rsi_x_c(2,2)
          rsi_x_c(irsi_frsizey+1,2) = 2 * 3.1415926e0 
     %                  - rsi_x_c(irsi_frsizey,2)
          rsi_x_c(irsi_frsizey+2,2) = 2 * 3.1415926e0 
     %                  - rsi_x_c(irsi_frsizey-1,2)
c         
c         Cell interfaces
c         
          do iy=0,irsi_frsizey+2
              rsi_x_i(iy,2)  = 0.5d0 * ( rsi_x_c(iy-1,2) + 
     %                                   rsi_x_c(iy,2) )
          enddo
c         
c         Distance between cell centers
c         
          do iy=0,irsi_frsizey+2
              rsi_dx_c(iy,2) = rsi_x_c(iy,2) - rsi_x_c(iy-1,2)
          enddo
          rsi_dx_c(irsi_frsizey+2,2) = 0.d0
c         
c         Distance between cell interfaces
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_i(iy,2) = rsi_x_i(iy+1,2) - rsi_x_i(iy,2)
          enddo
c         
c         Now the dx_cir (distance between cell center and its right cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cir(iy,2) = rsi_x_i(iy+1,2) - rsi_x_c(iy,2)
          enddo
          rsi_dx_cir(irsi_frsizey+2,2) = 0.d0
c         
c         Now the dx_cil (distance between cell center and its left cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cil(iy,2) = rsi_x_c(iy,2) - rsi_x_i(iy,2)
          enddo
          rsi_dx_cil(irsi_frsizey+2,2) = 0.d0
c
      elseif(spacegrid_type_y.eq.'-') then
c
c         Read Theta grid from file, and recompute all other data. 
c         Beware of the theta mirroring option.
c
          open(unit=3,file='theta.inp',status='old')
          read(3,*) irsi_frsizey,imirt
          itmax = irsi_frsizey
#         ifdef RADGRID_ONEDIM
          imirt = 0
          if(irsi_frsizey.ne.1) then
              write(*,*) 'ERROR: RADICAL compiled for 1-D, yet the '
              write(*,*) '       theta.inp file has more than 1 ',
     %                           'grid point.'
              stop 13
          endif
          thmaxerror = 3.14159265359d0
#         else /* ifdef RADGRID_ONEDIM */
          if(imirt.ne.0) then
              irsi_frsizey = 2*irsi_frsizey
          endif
          thmaxerror = 3.14159265359d0
#         ifdef MIRROR_THETA
          thmaxerror = 1.5707963268d0
          if(imirt.eq.0) then
              write(*,*) 'ERROR: RADICAL is compiled with MIRROR_THETA'
              write(*,*) '       So I expect that the theta.inp file'
              write(*,*) '       also has Theta-mirroring, and that it'
              write(*,*) '       is specified by 1 as second number in'
              write(*,*) '       the first line. '
              stop 13
          endif
#         else
          if(imirt.eq.1) then
              write(*,*) 'ERROR: RADICAL is compiled without ',
     %                           'MIRROR_THETA'
              write(*,*) '       So I expect that the theta.inp file'
              write(*,*) '       to have no Theta-mirroring. Yet it '
              write(*,*) '       does seem to have Theta-mirroring!'
              stop 13
          endif
#         endif
          if(irsi_frsizey.gt.FRSIZE_Y) then
              write(*,*) 'ERROR: theta.inp has more Theta gridpoints'
              write(*,*) '       than FRSIZE_Y. '
              write(*,*) '       Recompile with larger FRSIZE_Y'
              stop 13
          endif
#         endif /* ifdef else RADGRID_ONEDIM */
          do iy=1,itmax
              read(3,*) rsi_x_c(iy,2)
              if(rsi_x_c(iy,2).gt.thmaxerror) then
                  write(*,*) 'ERROR: theta.inp has too large values'
                  write(*,*) '       of theta. The limit is ',thmaxerror
                  stop 13
              endif
              if(rsi_x_c(iy,2).lt.0.d0) then
                  write(*,*) 'ERROR: theta.inp has negative values.'
                  stop 13                  
              endif
          enddo
          if(imirt.ne.0) then
              do iy=1,itmax
                  rsi_x_c(2*itmax+1-iy,2) = 
cccc ERROR (09-02-00)    %                      1.5707963268d0 - rsi_x_c(iy,2)
     %                      3.14159265359d0 - rsi_x_c(iy,2)
              enddo      
          endif
          close(3)
          spacegrid_theta_read = 1
c
c         Do boundaries right
c
          rsi_x_c(0,2)  = - rsi_x_c(1,2)
          rsi_x_c(-1,2) = - rsi_x_c(2,2)
          rsi_x_c(irsi_frsizey+1,2) = 2 * 3.1415926e0 
     %                  - rsi_x_c(irsi_frsizey,2)
          rsi_x_c(irsi_frsizey+2,2) = 2 * 3.1415926e0 
     %                  - rsi_x_c(irsi_frsizey-1,2)
c         
c         Cell interfaces
c         
          do iy=0,irsi_frsizey+2
              rsi_x_i(iy,2)  = 0.5d0 * ( rsi_x_c(iy-1,2) + 
     %                                   rsi_x_c(iy,2) )
          enddo
c         
c         Distance between cell centers
c         
          do iy=0,irsi_frsizey+2
              rsi_dx_c(iy,2) = rsi_x_c(iy,2) - rsi_x_c(iy-1,2)
          enddo
          rsi_dx_c(irsi_frsizey+2,2) = 0.d0
c         
c         Distance between cell interfaces
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_i(iy,2) = rsi_x_i(iy+1,2) - rsi_x_i(iy,2)
          enddo
c         
c         Now the dx_cir (distance between cell center and its right cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cir(iy,2) = rsi_x_i(iy+1,2) - rsi_x_c(iy,2)
          enddo
          rsi_dx_cir(irsi_frsizey+2,2) = 0.d0
c         
c         Now the dx_cil (distance between cell center and its left cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cil(iy,2) = rsi_x_c(iy,2) - rsi_x_i(iy,2)
          enddo
          rsi_dx_cil(irsi_frsizey+2,2) = 0.d0
c
      elseif(spacegrid_type_y.eq.'r') then
c
c         Regularly spaced Theta, with a location of smooth grid refinement
c         The refinement is introduced using a new method: that of providing
c         a continuous function i_theta(theta) on the domain of theta. The
c         root of this-function-minus-an-integer is a gridpoint. This is a very
c         easy way to introduce all kinds of clever smooth griddings. And if
c         you want to add any kind of additional grid refinement, simply
c         add a term to the equation provided in the function grid_linear_refine()
c         et voila. The root finding is done using zbrent, a very robust 
c         Numerical Recipes algorithm for finding the root of a bracketed
c         function. 
c
          if((irsi_frsizey+1)/2.ne.irsi_frsizey/2) then
              write(*,*) 'Must have even number of Theta gridpoints'
              stop 13
          endif
          if(irsi_frsizey/2.lt.4) then
              write(*,*) 'Need at least 8 gridpoints in Theta'
              write(*,*) ' (from pole to pole, i.e. 4 from pole to eq).'
              stop 9
          endif
          if(abs(spacegrid_refine_a).lt.1.e-7) then
              write(*,*) 'spacegrid_refine_a not defined!'
              stop 13
          endif
          if(irsi_frsizey.gt.FRSIZE_Y) then
              write(*,*) 'ERROR: Too many Theta points attempted'
              stop 18
          endif
          grsol_x0   = 0.0000000000d0
          grsol_x1   = 1.5707963268d0
          grsol_xr   = spacegrid_refine_thetar
          grsol_a    = spacegrid_refine_a
          grsol_N    = irsi_frsizey/2 - spacegrid_refine_n
          if(grsol_N.le.4) then
              write(*,*) 'WARNING: Theta grid refinement eats too ',
     %                   'many gridpoints away'
              write(*,*) '         I therefore reduce the refinement'
              grsol_N = 4
              spacegrid_refine_n = irsi_frsizey/2 - 4
          endif
          grsol_M    = spacegrid_refine_n
c
c         Create the actual grid cell center points and cell interfaces
c         The functions called here will compute the gridpoints by
c         solving the equation i_theta(theta)=1,2,3,.... etc. 
c
          call solve_grid_c(grid_linear_refine,grsol_x0,grsol_x1,
     %          grsol_N+grsol_M,1,0,rsi_x_c(-1,2),-1,FRSIZE_MAX+2)
          call solve_grid_i(grid_linear_refine,grsol_x0,grsol_x1,
     %          grsol_N+grsol_M,1,0,rsi_x_i(-1,2),-1,FRSIZE_MAX+2)
c
c         Now make the copies to the other quadrant
c
          do iy=0,grsol_N+grsol_M
              rsi_x_c(2*(grsol_N+grsol_M)+1-iy,2) = 
     %           3.14159265359d0 - rsi_x_c(iy,2)
          enddo
          do iy=0,grsol_N+grsol_M
              rsi_x_i(2*(grsol_N+grsol_M)+2-iy,2) = 
     %            3.14159265359d0 - rsi_x_i(iy,2)
          enddo
cc##################################
cc     Print out:
cc
c      write(*,*) spacegrid_refine_thetar,spacegrid_refine_a
c      write(*,*) grsol_N,grsol_M
c      do iy=0,2*(grsol_N+grsol_M)+1
c          write(*,*) iy,rsi_x_c(iy,2)
c      enddo
c      write(*,*) 
c      do iy=0,2*(grsol_N+grsol_M)+2
c          write(*,101) rsi_x_i(iy,2)
c      enddo
c  101 format(1x,E13.6)
c          dth = 1d33
cc##################################
c         
c         Distance between cell centers
c         
          do iy=0,irsi_frsizey+2
              rsi_dx_c(iy,2) = rsi_x_c(iy,2) - rsi_x_c(iy-1,2)
          enddo
          rsi_dx_c(irsi_frsizey+2,2) = 0.d0
c         
c         Distance between cell interfaces
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_i(iy,2) = rsi_x_i(iy+1,2) - rsi_x_i(iy,2)
          enddo
c         
c         Now the dx_cir (distance between cell center and its right cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cir(iy,2) = rsi_x_i(iy+1,2) - rsi_x_c(iy,2)
          enddo
          rsi_dx_cir(irsi_frsizey+2,2) = 0.d0
c         
c         Now the dx_cil (distance between cell center and its left cell interface)
c         
          do iy=0,irsi_frsizey+1
              rsi_dx_cil(iy,2) = rsi_x_c(iy,2) - rsi_x_i(iy,2)
          enddo
          rsi_dx_cil(irsi_frsizey+2,2) = 0.d0
      else
          write(*,*) 'Other theta grid types than L or l ',
     %     'not yet implemented'
          write(*,*) 'ABORTING'
          stop 13
      endif      
c
      goto 710
  701 continue
      write(*,*) 'Could not open file theta.inp'
      stop 13
  710 continue
      return
      end
      


c     --------------------------------------------------------------
c                      WRITE GRID TO FILE
c     --------------------------------------------------------------
      subroutine write_grid()
      implicit none
c
#include "common_grid.h"
c
      integer is,iangset,imirt,itmax
c
      iangset = 1
c
      call set_mirror_limits_theta(itmax,imirt)
c
      open(unit=1,file='radius.dat',status='unknown')
      write(1,*) irsi_frsizex
      do is=1,irsi_frsizex
          write(1,*) rsi_x_c(is,1)
      enddo
      write(1,*) ' '
      do is=1,irsi_frsizex+1
          write(1,*) rsi_x_i(is,1)
      enddo
      close(1)
c     
      open(unit=1,file='theta.dat',status='unknown')
      write(1,*) itmax,imirt
c
c        Warning: the imirt here has been added later (01-08-00)
c
c        -------------------------------------------------------
c        WARNING: Completely changed this: now only upper 
c                 quadrant dumped in case of imirt.eq.1
c        -------------------------------------------------------
c
      do is=1,itmax
          write(1,*) rsi_x_c(is,2)
      enddo
      close(1)
c     
#ifndef NO_SHORT_CHARS
c
      open(unit=1,file='mu.dat',status='unknown')
      write(1,*) anggrid_frsizemu
      do is=-nrmu(iangset),nrmu(iangset)
          if((is.ne.0).or.(anggrid_mu_zero.eq.1)) then
              write(1,*) rmu(is,iangset)
          endif
      enddo
      write(1,*) ' '
      do is=-nrmu(iangset)-1,-1
          write(1,*) rmu_i(is,iangset)
      enddo
      if(anggrid_mu_zero.eq.1) then
          do is=1,nrmu(iangset)+1
              write(1,*) rmu_i(is,iangset)
          enddo
      else
          do is=2,nrmu(iangset)+1
              write(1,*) rmu_i(is,iangset)
          enddo
      endif
      close(1)
c     
      open(unit=1,file='phi.dat',status='unknown')
      write(1,*) anggrid_frsizephi
      do is=1,nrphi(iangset)
          write(1,*) rphi(is,iangset)
      enddo
      write(1,*) ' '
      do is=1,nrphi(iangset)
          write(1,*) rphi_i(is,iangset)
      enddo
      close(1)
c
#endif
c     
      if(freq_nr.gt.1) then
         open(unit=1,file='frequency.dat',status='unknown')
         write(1,*) freq_nr
         do is=1,freq_nr
            write(1,*) freq_nu(is)
         enddo
         close(1)
      endif
c
      end



c     --------------------------------------------------------------
c                       CREATE THE FREQUENCY GRID
c     --------------------------------------------------------------
      subroutine create_freqgrid(nu_start,nu_end,nr,type,iline)
      implicit none
      doubleprecision nu_start,nu_end
      integer nr,type,iline
      logical fex
c     
#include "common_grid.h"
c
      integer inu
      doubleprecision frnu,dfrnu
c
      if(nr.le.1) then
          cont_freq_nu(1) = nu_start
          icont_freq_nr   = 1
          return
      endif
c
      if(type.ne.-1) then
          write(*,*) 'ERROR: in create_freqgrid...'
          stop 87090
      endif
c
c         Read frequencies from file
c
      inquire(file='frequency.inp',exist=fex)
      if(fex) then 
          open(unit=1,file='frequency.inp',status='old',err=701)
          read(1,*) nr
          if(nr.gt.FRSIZE_FREQ) then
              write(*,*) 'ERROR: frequency.inp has more gridpoints'
              write(*,*) '       than FRSIZE_FREQ. '
              write(*,*) '       Recompile with larger FRSIZE_FREQ'
              stop 13
          endif          
          if(nr.eq.1) then
              write(*,*) 'Going into frequency-integrated (grey ',
     %             'opacity) mode'
          else
              do inu=1,nr
                  read(1,*) cont_freq_nu(inu)
              enddo
          endif          
          icont_freq_nr = nr
          close(1)
          freq_done = 1
      endif
ccc
ccc     NEW 28-03-06: 
ccc   ***** REMOVED BECAUSE WE NOW DO NOT DO PURE DUST ANYMORE *****
cc      if(iline.eq.0) then
cc          freq_nr = icont_freq_nr
cc          do inu=1,nr
cc              freq_nu(inu) = cont_freq_nu(inu) 
cc          enddo
cc      endif
c
      goto 710
  701 continue
      write(*,*) 'Could not open file frequency.inp'
      stop 13
  710 continue
      return
      end



c     --------------------------------------------------------------
c           ROUTINE FOR GENERATING THE CELL-CENTER GRID POINTS
c
c     The routine grfunc() is handed over to zbrent to solve for the
c     values of x \in [x0,x1] belonging to the gridpoint numbers
c      1,2,3,...,Ntot . The Nex is the amount of ghost points that
c     are to be computed as well. Note that these ghost points are
c     simply mirror copies of the real points. 
c     --------------------------------------------------------------
      subroutine solve_grid_c(grfunc,x0,x1,Ntot,Nex0,Nex1,
     %                        grid,Nmin,Nmax)
      implicit none
      integer Ntot,Nex0,Nex1,Nmin,Nmax
      doubleprecision x0,x1
      doubleprecision grid(Nmin:Nmax)
      doubleprecision zbrent
      doubleprecision grid_linear_refine
      external grfunc
      doubleprecision grfunc
c
      common/grsolve/grsol_ix
      doubleprecision grsol_ix
c
      doubleprecision tol,eps
      integer i
c
      tol = 1.d-5 * ( abs(x0) + abs(x1) )
c
c     Loop over real gridpoints
c
      do i=1,Ntot
          grsol_ix = 1.d0 * i
          grid(i)  = zbrent(grfunc,x0-eps,x1+eps,tol)
      enddo
c
c     Now mirror copy
c
      do i=0,1-Nex0
          grid(i) = 2.d0*x0 - grid(1-i)
      enddo
      do i=Ntot+1,Ntot+Nex1
          grid(i) = 2.d0*x1 - grid(2*Ntot-i+1)
      enddo
c
      end


c     --------------------------------------------------------------
c           ROUTINE FOR GENERATING THE CELL-INTERFACE GRID POINTS
c
c     The routine grfunc() is handed over to zbrent to solve for the
c     values of x \in [x0,x1] belonging to the gridpoint numbers
c      0.5,1.5,2.5,...,Ntot+0.5 . The Nex0,Nex1 are the amount of ghost 
c     points that are to be computed as well on each side. Note that 
c     these ghost points are simply mirror copies of the real points. 
c     --------------------------------------------------------------
      subroutine solve_grid_i(grfunc,x0,x1,Ntot,Nex0,Nex1,
     %                        grid,Nmin,Nmax)
      implicit none
      integer Ntot,Nex0,Nex1,Nmin,Nmax
      doubleprecision x0,x1
      doubleprecision grid(Nmin:Nmax)
      doubleprecision zbrent
      external grfunc
      doubleprecision grfunc
c
      common/grsolve/grsol_ix
      doubleprecision grsol_ix
c
      doubleprecision tol,eps
      integer i
c
      tol = 1.d-5 * ( abs(x0) + abs(x1) )
      eps = 1.d-5 * ( abs(x0) + abs(x1) )
c
c     Loop over real gridpoints
c
      do i=1,Ntot+1
          grsol_ix = 1.d0 * i - 0.5d0
          grid(i)  = zbrent(grfunc,x0-eps,x1+eps,tol)
      enddo
c
c     Now mirror copy
c
      do i=0,1-Nex0
          grid(i) = 2.d0*x0 - grid(2-i)
      enddo
      do i=Ntot+2,Ntot+Nex1+1
          grid(i) = 2.d0*x1 - grid(2*Ntot-i+2)
      enddo
c
      end


c     --------------------------------------------------------------
c               ONE OF THE POSSIBLE DISCRETIZATION FUNCTIONS
c
c     This function provides the index (in float form) ix for a 
c     given x. The solve_grid() routine then will solve for 
c     the actual locations of the gridpoints ix=1,2,3...,N+M,
c     by invoking Num Recip routine zbrent(), which finds the
c     root of the equations ix(x) = 1,2,3,....M+N .
c
c     This particular function is a linear gridding, including one
c     smooth zone of refinement. The function is:
c
c                  1     x  - x0     tanh[(x -xr)/a] - tanh[(x0-xr)/a]
c        ix(x)  =  - + N ------- + M -------------------------------
c                  2     x1 - x0     tanh[(x1-xr)/a] - tanh[(x0-xr)/a]
c
c     N    = The number of linear gridpoints
c     M    = The number of gridpoints in the refinement zone
c     x0   = Lower boundary of x
c     x1   = Upper boundary of x
c     xr   = Center of refinement region
c     a    = coefficient of refinement
c
c     --------------------------------------------------------------
      function grid_linear_refine(x)
      implicit none
      doubleprecision grid_linear_refine,x
c
      common/gridlinref/grsol_x0,grsol_x1,grsol_xr,grsol_a
      doubleprecision grsol_x0,grsol_x1,grsol_xr,grsol_a
      common/igridlinref/grsol_N,grsol_M
      integer grsol_N,grsol_M
c
      common/grsolve/grsol_ix
      doubleprecision grsol_ix
c
      grid_linear_refine = 0.5d0 
     %  + grsol_N * ( x - grsol_x0 ) / ( grsol_x1 - grsol_x0 )
     %  + grsol_M * ( tanh((x-grsol_xr)/grsol_a) -
     %                tanh((grsol_x0-grsol_xr)/grsol_a) ) /
     %              ( tanh((grsol_x1-grsol_xr)/grsol_a) -
     %                tanh((grsol_x0-grsol_xr)/grsol_a) )
     %  - grsol_ix
c
      return
      end


c     --------------------------------------------------------------
c              HELP ROUTINE FOR AT THE START OF A MAIN LOOP
c
c     This routine will automatically set the correct limits of the
c     big spatial arrays, taking into account the mirror symmetries.
c     In principle this is stupid way of doing things. These limits
c     should be set in global variables. But for historic reasons
c     I re-computed them at the start of each loop. To remove (some
c     of) these disturbing 10 or 12 lines, they can be replaced with
c     a call to this routine.
c     --------------------------------------------------------------
      subroutine set_mirror_limits_theta(itmax,imirt)
      implicit none
      integer itmax,imirt
#     include "common_grid.h"
      itmax   = irsi_frsizey
      imirt   = 0
#     ifdef MIRROR_THETA
      itmax = (irsi_frsizey+1)/2
      imirt = 1
#     endif
#     ifdef RADGRID_ONEDIM
      itmax = 1
#     endif
      return
      end


c     --------------------------------------------------------------
c               FIND LARGEST GRID SPACING OF RADIAL GRID
c     --------------------------------------------------------------
      subroutine find_largest_rgrid_drr(drr)
      implicit none
c
#     include "common_grid.h"
c
      doubleprecision ratio,drr
      integer ir
c
      ratio = rsi_x_c(2,1)/rsi_x_c(1,1)
      do ir=2,irsi_frsizex-1
          if(rsi_x_c(ir+1,1)/rsi_x_c(ir,1).gt.ratio) then
             ratio = rsi_x_c(ir+1,1)/rsi_x_c(ir,1)
          endif
      enddo
c
      drr=ratio-1.d0
      return
      end
