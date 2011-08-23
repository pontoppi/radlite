#include "main.h"

c     --------------------------------------------------------------
c                 COPY THE PRECOMPUTED PHYSICAL VARIABLES
c
c     This routine is actually a left-over from an inefficient way
c     of doing things. Perhaps it can be eliminated once.
c
c     REMARK: It now automatically checks if it remains within its
c             bounds. This is especially relevant for mirror
c             symmetries.
c     --------------------------------------------------------------
      subroutine get_phys_var(ir,it)
      implicit none
c
      integer ir,it
c
#include "common_grid.h"
#include "common_setup.h"
c
#     ifdef SAFETY_CHECKS_ACTIVE
#     ifndef RADGRID_ONEDIM
#     ifdef MIRROR_THETA
      if(it.gt.(irsi_frsizey+1)/2) then
          write(*,*) 'ERROR: medium.F/get_phys_var()'
          write(*,*) '       it in mirror quadrant!'
          stop 13
      endif
#     endif
#     endif
#     endif
c
      medium_rho = medium_arr_rho(it,ir)
      medium_el_temp = medium_arr_el_temp(it,ir)
      medium_bmag = medium_arr_bmag(it,ir)
c
      end


c     --------------------------------------------------------------
c                 WRITE THE PHYSICAL VARIABLES TO FILE
c     --------------------------------------------------------------
      subroutine write_phys_var()
      implicit none
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
#include "common_dust.h"
c
      integer ir,it,imirt,itmax,ispec
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
c     First the backwards consistency file: the archaic 
c     physvar.dat file. 
c
      open(unit=1,file='physvar.dat',status='unknown')
      write(1,*) irsi_frsizex,itmax,imirt
      write(1,*) ' '
      do ir=1,irsi_frsizex
          do it=1,itmax
              call get_phys_var(ir,it)
              write(1,*) medium_rho,medium_el_temp
          enddo
          write(1,*) ' '
      enddo
#     ifdef LINE_VELOCITIES
      write(1,*) ' '
      do ir=1,irsi_frsizex
          do it=1,itmax
              write(1,*) medium_arr_velocity(1,it,ir),
     %                   medium_arr_velocity(2,it,ir),
     %                   medium_arr_velocity(3,it,ir)
          enddo
          write(1,*) ' '
      enddo
#     endif
      close(1)
c
c     If the main density is useful, then write it in a 
c     separate file
c
      if(iradproc_dust.eq.0) then
          open(unit=1,file='density.dat',status='unknown')
          write(1,*) irsi_frsizex,itmax,imirt
          write(1,*) ' '
          do ir=1,irsi_frsizex
              do it=1,itmax
                  call get_phys_var(ir,it)
                  write(1,*) medium_rho
              enddo
          enddo
          write(1,*) ' '
          close(1)
      endif
c
c     If the main temperature is useful, then write it in a 
c     separate file
c
      if((iradproc_dust.eq.0).and.(iradproc_solve_temp.eq.0)) then
          open(unit=1,file='temperature.dat',status='unknown')
          write(1,*) irsi_frsizex,itmax,imirt
          write(1,*) ' '
          do ir=1,irsi_frsizex
              do it=1,itmax
                  call get_phys_var(ir,it)
                  write(1,*) medium_el_temp
              enddo
          enddo
          write(1,*) ' '
          close(1)
      endif
c
c     Now, if there are any dust species, then
c
#ifdef INCLUDE_DUST
      if((iradproc_dust.ne.0).and.
     %   (dust_setup_nrspecies.ge.1)) then
          open(unit=1,file='dustdens.dat',status='unknown')
          write(1,*) dust_setup_nrspecies,irsi_frsizex,itmax,imirt
          write(1,*) ' '
          do ispec=1,dust_setup_nrspecies
              do ir=1,irsi_frsizex
                  do it=1,itmax
                      write(1,*) dust_rho(ispec,it,ir)
                  enddo
              enddo
              write(1,*) ' '
          enddo
          write(1,*) ' '
          close(1)
      endif
#endif
      
      end


c     --------------------------------------------------------------
c          READ EXTERNALLY PROVIDED DUST TEMPERATURES FROM FILE
c
c     For certain setup styles not only the dust density, but also
c     the dust temperatures are being read. No iteration is then
c     needed, and a direct spectrum and images can be calculated 
c     by ray tracing. Here is the routine for importing the dust
c     temperatures.
c     --------------------------------------------------------------
      subroutine read_dust_temperatures(inpfile)
      implicit none
      logical inpfile
#ifdef INCLUDE_DUST
c
      integer ir,it
      logical fex
      character*80 string
      integer nr,nt,imirt,i,idum
      integer ispec,flnr
      doubleprecision r(FRSIZE_X),theta(FRSIZE_Y)
      doubleprecision dum,xi,dxi,dxirel,dumar(10),dummy
      character*80 comstring,filename,base,ext
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
#include "common_dust.h"
c
c     Find out which filename
c
      if(inpfile) then
          filename = 'dusttemp.inp'
          write(*,*) 'Reading dust temperature data from dusttemp.inp'
          inquire(file=filename,EXIST=fex)
          if(.not.fex) then 
              write(*,*) 'ERROR: dusttemp.inp not found'
              stop 1
          endif
      else
          inquire(file='dusttemp.info',EXIST=fex)
          if(fex) then
              open(unit=1,file='dusttemp.info',status='old')
              flnr=0
              read(1,*) flnr
              close(1)
              if(flnr.eq.-2) then
                  filename='dusttemp_final.dat'
              else
                  base='dusttemp_'
                  ext ='.dat'
                  call make_indexed_filename(base,flnr,ext,filename)
              endif
          else
              write(*,*) 'Could not find dusttemp.info'
              write(*,*) 'Dont know which dusttemp file is latest'
              stop 13
          endif
          write(*,*) 'Reading dust temperature data from ',filename
      endif
c
c     Read the dust temperatures. The dust temp arrays have big Theta-size,
c     as opposed to the usual arrays. Sorry for this inconsistency. In
c     any case: By passing 1 as second argument, we make this clear to the
c     scalar field reader.
c
cc      call read_vecvector_field(filename,1,DUST_SPECIES_MAX,
cc     %                       DUST_SIZE_MAX,
cc     %                       0,FRSIZE_X,0,FRSIZE_Y,dust_temp,
cc     %                       dust_setup_nrspecies,dust_setup_nrsizes,
cc     %                       nr,nt,imirt)
cc
      call read_vecvector_field(filename,0,DUST_SPECIES_MAX,
     %                       DUST_SIZE_MAX,
     %                       0,FRSIZE_X,0,FRSIZE_Y,dust_temp,
     %                       dust_setup_nrspecies,dust_setup_nrsizes,
     %                       nr,nt,imirt)
c
c     Return some values
c
      setup_dusttemp_read  = 1
c
c     Check if the dimensions are consistent with the dimensions
c     of the radlite.inp file
c
      if((nr.ne.irsi_frsizex).or.(nt.ne.irsi_frsizey)) then
          write(*,*) 'ERROR: The dimensions of dusttemp.inp'
          write(*,*) '       are inconsistent with the spatial '
          write(*,*) '       dimensions specified in radlite.inp'
          write(*,*) '       or radius.inp / theta.inp. '
          write(*,*) '       Modify radlite.inp to make it consistent.'
          stop 13
      endif
c
#else 
      stop 66
#endif
      end



c     --------------------------------------------------------------
c             READ SPATIAL VECTOR-VECTOR ARRAY FROM FILE
c
c     This routine is used by read_dust_temperature(). 
c     It is a general-purpose routine for reading a vector field
c     containing vectors of varying length as its elements.
c     The grid in R and Theta is found in the input
c     files radius.inp and theta.inp. 
c
c     A vector-vector may look like this (left a vector-vector,
c     right a normal vector):
c
c     Vector-vector (4 elements):     Normal vector (4 elements):
c
c       -----                              -
c       ---                                -
c       -------                            -
c       --                                 -
c
c     This is typical for dust temperatures, since each species
c     can have a different number of dust sizes.
c     
c     Note: the vectors must have starting index 1. 
c
c     ARGUMENTS:
c       filename             The file to read from
c       nvmax                The vector maximum size
c       nvvmax               The element maximum vector size
c       nrmin                The starting index of radial dimension
c       ntmin                The starting index of theta dimension
c       nrmax                The radial maximum size of the array
c       ntmax                The theta maximum size of the array
c     
c     RESULT:
c       vector               The spatial vector-vector field array to read
c       nv                   The actual number of vector components
c       nvv                  An array of element vector lengths
c       nr                   The actual number of radial points
c       nt                   The actual number of theta points
c       imirt                Whether or not mirror symmetry in Theta
c                            was present (nt is then twice the actual
c                            theta points read).
c       isep                 0=vector is inner index, 1=vector is outer index
c
c     --------------------------------------------------------------
      subroutine read_vecvector_field(filename,icopymirr,nvmax,
     %                             nvvmax,nrmin,nrmax,ntmin,ntmax,
     %                             vecvector,nv,nvv,nr,nt,imirt)
      implicit none
      integer nrmin,ntmin,nrmax,ntmax,nvmax,nvvmax
      integer nv,nr,nt,imirt,icopymirr,nvv(nvvmax)
      doubleprecision vecvector(nvvmax,nvmax,ntmin:ntmax,nrmin:nrmax)
      integer itmax,it,ir,iv,ivv
      character*80 filename
c
#include "common_grid.h"
c
c     First check if the radius.inp and theta.inp have been read.
c
#     ifdef RADGRID_ONEDIM
      spacegrid_theta_read = 1
#     endif      
      if((spacegrid_radius_read.eq.0).or.
     %     (spacegrid_theta_read.eq.0)) then
          write(*,*) 'Sorry, before reading a vector-vector field from '
          write(*,*) '       a file, the radius.inp and theta.inp'
          write(*,*) '       files must have been read. This is'
          write(*,*) '       enabled by specifying - as the grid '
          write(*,*) '       type on the radlite.inp file'
         stop 13
      endif
c
c     Open file
c
      open(unit=1,file=filename,status='old',err=10)
c
c     First read the dimensions of the array
c
      read(1,*) nv,nr,nt,imirt
      itmax = nt
#     ifndef RADGRID_ONEDIM
      if(imirt.ne.0) then
          nt = 2*nt
      endif
#     else
      imirt = 0
      if(nt.ne.1) then
          write(*,*) 'ERROR while reading dustdensity.inp'
          write(*,*) '      RADICAL is compiled for 1-D, but'
          write(*,*) '      ',filename,' seems to be 2-D'
          stop 13
      endif
#     endif
c
c     Check if the dimensions do not exceed maximum ones
c
      if(nv.gt.nvmax) then
          write(*,*) 'Too many vector components!'
          write(*,*) 'File has nv=',nv,' while nvmax=',nvmax
          stop 13
      endif
      if(nr.gt.FRSIZE_X) then
          write(*,*) 'R-size too large!'
          write(*,*) 'File has nr=',nr,' while FRSIZE_X=',FRSIZE_X
          write(*,*) 'Recompile with larger array size FRSIZE_X'
          stop 13
      endif
      if(nt.gt.FRSIZE_Y) then
          write(*,*) 'Theta-size too large!'
          write(*,*) 'File has nt=',nt,' while FRSIZE_Y=',FRSIZE_Y
          if(imirt.ne.0) then
              write(*,*) 'NOTE: With mirror symmetry nt = twice that'
              write(*,*) '      which is specified in the file.'
          endif
          write(*,*) 'Recompile with larger array size FRSIZE_Y'
          stop 13
      endif
#     ifdef MIRROR_THETA
      if(imirt.eq.0) then
          write(*,*) 'ERROR: RADICAL is compiled with MIRROR_THETA'
          write(*,*) '       So I expect that the input files'
          write(*,*) '       also has Theta-mirroring, and that it'
          write(*,*) '       is specified by 1 as third number in'
          write(*,*) '       the first line. '
          stop 13
      endif
#     endif
c
c     Now read the vector field 
c
      do iv=1,nv
          read(1,*,err=30,end=20) nvv(iv)
          if(nvv(iv).gt.nvvmax) then
              write(*,*) 'Too many vecvector components at iv=',iv
              write(*,*) 'nvv=',nvv(iv),' while nvvmax=',nvvmax
              stop 13
          endif
          do ivv=1,nvv(iv)
              do ir=1,nr
                  do it=1,itmax
                      read(1,*,err=30,end=20) vecvector(ivv,iv,it,ir)
                  enddo
              enddo
#             ifndef RADGRID_ONEDIM
              if((imirt.ne.0).and.(icopymirr.eq.1)) then
c                 BUGFIX (16-09-00): Forgot do ir=1,nr
                  do ir=1,nr
                      do it=1,itmax
                          vecvector(ivv,iv,2*itmax+1-it,ir) = 
     %                         vecvector(ivv,iv,it,ir)
                      enddo
                  enddo
              endif
#             endif
          enddo
      enddo
c
c     Close file
c
      close(1)
c
      goto 100
   10 continue
      write(*,*) 'Could not open file ',filename
      stop 13      
   20 continue
      write(*,*) 'Prematurely reached end of file ',filename
      stop 13
   30 continue
      write(*,*) 'Could not understand contents of ',filename
      stop 13
  100 continue
c
      end


c     --------------------------------------------------------------
c              COMPUTE ONLY OUTER FREQ-DEPENDENT LUMINOSITY
c
c     This routine is much faster than compute_lumfunc_freq(), and
c     is therefore more appropriate for use after each iteration.
c     But it only computes the flux at the very outer edge. 
c     --------------------------------------------------------------
      subroutine outer_luminosity_freq(inu)
      implicit none
      integer inu
c
#include "common_source.h"
#include "common_radfield.h"
#include "common_grid.h"
c
      doubleprecision pi,r,cost1,cost2,surf,dum,dumm
      integer iangset,itmax,imumax,iphimax
      integer ir,it,imu,iphi
      parameter(pi=PICONST)
c
c     At the moment allow only one angular set
c
      iangset   = 1
c
c     Get some important integer values...
c
      itmax   = irsi_frsizey
      imumax  = nrmu(iangset)
      iphimax = nrphi(iangset)
c
c     Now the mirror things...
c
#ifdef MIRROR_THETA
      itmax   = (irsi_frsizey+1)/2
#endif
#ifdef MIRROR_PHI
      iphimax = nrphi(iangset) / 2
#endif
#ifdef RADGRID_ONEDIM
      iphimax = 1
      itmax   = 1
#endif
c
      ir=irsi_frsizex
      tlum_lumnu(inu,ir) = 0.d0
      r = rsi_x_c(ir,1)
      do it=1,itmax
          cost1 = cos(rsi_x_i(it,2))
          cost2 = cos(rsi_x_i(it+1,2))
          surf  = 2.d0*pi*r*r*abs(cost1-cost2)
          dum   = 0.d0
          do imu=1,imumax
              dumm = 0.d0
              do iphi=1,iphimax
                  dumm = dumm + intens(iphi,imu,it,ir)
                  dumm = dumm - intens(iphi,-imu,it,ir)
              enddo
              dum = dum + dumm * rmu(imu,iangset) *
     %             ( rmu_i(imu+1,iangset) -
     %               rmu_i(imu,iangset) )
          enddo
          dum = 2.d0 * pi * dum / ( 1.d0*iphimax )
          tlum_lumnu(inu,ir) = tlum_lumnu(inu,ir) + dum * surf
      enddo
#     ifdef MIRROR_THETA
      tlum_lumnu(inu,ir) = 2.d0 * tlum_lumnu(inu,ir)
#     endif
c    
      end


c     --------------------------------------------------------------
c                READ SPATIAL VECTOR ARRAY FROM FILE
c
c     This routine is used by read_dust_density() etc.
c     It is a general-purpose routine for reading vector field
c     input files. The grid in R and Theta is found in the input
c     files radius.inp and theta.inp. 
c
c     ARGUMENTS:
c       filename             The file to read from
c       nvmin                The starting index of vector dimension
c       nrmin                The starting index of radial dimension
c       ntmin                The starting index of theta dimension
c       nvmax                The vector maximum size
c       nrmax                The radial maximum size of the array
c       ntmax                The theta maximum size of the array
c     
c     RESULT:
c       vector               The spatial vector field array to read
c       nv                   The actual number of vector components
c       nr                   The actual number of radial points
c       nt                   The actual number of theta points
c       imirt                Whether or not mirror symmetry in Theta
c                            was present (nt is then twice the actual
c                            theta points read).
c       isep                 0=vector is inner index, 1=vector is outer index
c
c     --------------------------------------------------------------
      subroutine read_vector_field(filename,icopymirr,nvmin,nvmax,
     %                             nrmin,nrmax,ntmin,ntmax,
     %                             vector,nv,nr,nt,imirt,isep)
      implicit none
      integer nvmin,nrmin,ntmin,nvmax,nrmax,ntmax
      integer nv,nr,nt,imirt,icopymirr,isep
      doubleprecision vector(nvmin:nvmax,ntmin:ntmax,nrmin:nrmax)
      integer itmax,it,ir,iv
      character*80 filename
c
#include "common_grid.h"
c
c     First check if the radius.inp and theta.inp have been read.
c
#     ifdef RADGRID_ONEDIM
      spacegrid_theta_read = 1
#     endif      
      if((spacegrid_radius_read.eq.0).or.
     %     (spacegrid_theta_read.eq.0)) then
          write(*,*) 'Sorry, before reading a vector field from '
          write(*,*) '       a file, the radius.inp and theta.inp'
          write(*,*) '       files must have been read. This is'
          write(*,*) '       enabled by specifying - as the grid '
          write(*,*) '       type on the radlite.inp file'
         stop 13
      endif
c
c     Open file
c
      open(unit=1,file=filename,status='old',err=10)
c
c     First read the dimensions of the array
c
      read(1,*) nv,nr,nt,imirt
      itmax = nt
#     ifndef RADGRID_ONEDIM
      if(imirt.ne.0) then
          nt = 2*nt
      endif
#     else
      imirt = 0
      if(nt.ne.1) then
          write(*,*) 'ERROR while reading dustdensity.inp'
          write(*,*) '      RADICAL is compiled for 1-D, but'
          write(*,*) '      ',filename,' seems to be 2-D'
          stop 13
      endif
#     endif
c
c     Check if the dimensions do not exceed maximum ones
c
      if(nv.gt.nvmax) then
          write(*,*) 'Too many vector components!'
          write(*,*) 'File has nv=',nv,' while nvmax=',nvmax
          stop 13
      endif
      if(nr.gt.FRSIZE_X) then
          write(*,*) 'R-size too large!'
          write(*,*) 'File has nr=',nr,' while FRSIZE_X=',FRSIZE_X
          write(*,*) 'Recompile with larger array size FRSIZE_X'
          stop 13
      endif
      if(nt.gt.FRSIZE_Y) then
          write(*,*) 'Theta-size too large!'
          write(*,*) 'File has nt=',nt,' while FRSIZE_Y=',FRSIZE_Y
          if(imirt.ne.0) then
              write(*,*) 'NOTE: With mirror symmetry nt = twice that'
              write(*,*) '      which is specified in the file.'
          endif
          write(*,*) 'Recompile with larger array size FRSIZE_Y'
          stop 13
      endif
#     ifdef MIRROR_THETA
      if(imirt.eq.0) then
          write(*,*) 'ERROR: RADICAL is compiled with MIRROR_THETA'
          write(*,*) '       So I expect that the input files'
          write(*,*) '       also has Theta-mirroring, and that it'
          write(*,*) '       is specified by 1 as third number in'
          write(*,*) '       the first line. '
          stop 13
      endif
#     endif
c
c     Now read the vector field 
c
      if(isep.eq.1) then
c
c         Each component is separate
c
          do iv=1,nv
              do ir=1,nr
                  do it=1,itmax
                      read(1,*,err=30,end=20) vector(iv,it,ir)
                  enddo
              enddo
#             ifndef RADGRID_ONEDIM
              if((imirt.ne.0).and.(icopymirr.eq.1)) then
c                 BUGFIX (16-09-00): Forgot do ir=1,nr
                  do ir=1,nr
                      do it=1,itmax
                          vector(iv,2*itmax+1-it,ir) = vector(iv,it,ir)
                      enddo
                  enddo
              endif
#             endif
          enddo
      else
c
c         Components are collected together on a line
c
          do ir=1,nr
              do it=1,itmax
                  read(1,*,err=30,end=20) (vector(iv,it,ir),iv=1,nv)
              enddo
          enddo
#         ifndef RADGRID_ONEDIM
          if((imirt.ne.0).and.(icopymirr.eq.1)) then
c             BUGFIX (16-09-00): Forgot do ir=1,nr
              do ir=1,nr
                  do it=1,itmax
                      do iv=1,nv
                          vector(iv,2*itmax+1-it,ir) = vector(iv,it,ir)
                      enddo     
                  enddo
              enddo
          endif
#         endif
      endif
c
c     Close file
c
      close(1)
c
      goto 100
   10 continue
      write(*,*) 'Could not open file ',filename
      stop 13      
   20 continue
      write(*,*) 'Prematurely reached end of file ',filename
      stop 13
   30 continue
      write(*,*) 'Could not understand contents of ',filename
      stop 13
  100 continue
c
      end




c     --------------------------------------------------------------
c            READ EXTERNALLY PROVIDED DUST DENSITY FROM FILE
c
c     For certain setup styles (e.g. -500) the density of each of
c     the dust species is given in a file named "dustdensity.inp". 
c     --------------------------------------------------------------
      subroutine read_dust_density()
#ifdef INCLUDE_DUST
      implicit none
c
      integer ir,it,nrspecies
      logical fex
      character*80 string
      integer nr,nt,imirt,i,idum
      integer ispec
      doubleprecision r(FRSIZE_X),theta(FRSIZE_Y)
      doubleprecision dum,xi,dxi,dxirel,dumar(10),dummy
      character*80 filename
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
#include "common_dust.h"
c
c     Check the presence of the dustdensity.inp file
c
      filename = 'dustdensity.inp'
      write(*,*) 'Reading dust density data from dustdensity.inp'
      inquire(file=filename,EXIST=fex)
      if(.not.fex) then 
          filename = 'dustdens.inp'
          inquire(file=filename,EXIST=fex)
          if(.not.fex) then 
              write(*,*) 'ERROR: dustdensity.inp nor '
              write(*,*) '       dustdens.inp found'
              stop 1
          endif
      endif
c
c     Read the dust densities. The dust density arrays have big Theta-size,
c     as opposed to the usual density. Sorry for this inconsistency. In
c     any case: By passing 1 as second argument, we make this clear to the
c     scalar field reader.
c
cc      call read_vector_field(filename,1,1,DUST_SPECIES_MAX,
cc     %                       0,FRSIZE_X,0,FRSIZE_Y,dust_rho,
cc     %                       nrspecies,nr,nt,imirt,1)
cc
      call read_vector_field(filename,0,1,DUST_SPECIES_MAX,
     %                       0,FRSIZE_X,0,FRSIZE_Y,dust_rho,
     %                       nrspecies,nr,nt,imirt,1)
c
c     Return some values
c
      dust_setup_nrspecies = nrspecies
      setup_dustdens_read  = 1
c
c     Check if the dimensions are consistent with the dimensions
c     of the radlite.inp file
c
      if((nr.ne.irsi_frsizex).or.(nt.ne.irsi_frsizey)) then
          write(*,*) 'ERROR: The dimensions of dustdensity.inp'
          write(*,*) '       are inconsistent with the spatial '
          write(*,*) '       dimensions specified in radlite.inp'
          write(*,*) '       or radius.inp / theta.inp. '
          write(*,*) '       Modify radlite.inp to make it consistent.'
          write(*,*) '       nr=',nr,' nt=',nt
          stop 13
      endif
c
#endif
      end


c     --------------------------------------------------------------
c                          READ VELOCITY FIELD
c
c     This is necessary for the line transfer
c     --------------------------------------------------------------
      subroutine read_velocity()
      implicit none
c
      integer ir,it
      logical fex
      integer nr,nt
      doubleprecision vr,vt,vp
      character*80 filename
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
#include "common_lines.h"
c
c     Check the presence of the velocity.inp file
c
      filename = 'velocity.inp'
      inquire(file=filename,EXIST=fex)
      if(fex) then 
          write(*,*) 'Reading velocity data from velocity.inp'
c
          open(unit=1,file=filename)
          read(1,*) nr,nt
          if((nr.ne.irsi_frsizex).or.(nt.ne.irsi_frsizey/2)) stop 325
          do ir=1,nr
              do it=1,nt
                  read(1,*) vr,vt,vp
                  medium_arr_velocity(1,it,ir) = vr
                  medium_arr_velocity(2,it,ir) = vt
                  medium_arr_velocity(3,it,ir) = vp
              enddo
          enddo
          close(1)
      endif
c
      end



c     --------------------------------------------------------------
c                        READ TURBULENT VELOCITIES
c
c     This is necessary for the line transfer
c
c     NOTE: At the moment only a Gaussian TOTAL line width is used.
c           This includes also the thermal one, and is read here.
c           So the user will have to compute the thermal width
c           him/herself. 31-10-05.
c     --------------------------------------------------------------
      subroutine read_turbulence()
      implicit none
c
      integer ir,it,iformat,iline
      logical fex
      integer nr,nt
      doubleprecision vr,vt,vp
      character*80 filename
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
#include "common_lines.h"
c
c     Check the presence of the velocity.inp file
c
      filename = 'turbulence.inp'
      inquire(file=filename,EXIST=fex)
      if(fex) then 
          write(*,*) 'Reading turbulence data from turbulence.inp'
c
          open(unit=1,file=filename)
          read(1,*) iformat
          read(1,*) nr,nt
          if((nr.ne.irsi_frsizex).or.(nt.ne.irsi_frsizey/2)) stop 325
          do ir=1,nr
              do it=1,nt
c
c                 Read the local line width
c
                  read(1,*) locprof_linewidth(1,it,ir)
c
c                 Copy this to all (originally the linewidth could be
c                 line-dependent, but that is silly and unphysical anyway)
c
                  do iline=2,SZ_NLINES 
                      locprof_linewidth(iline,it,ir) = 
     %                         locprof_linewidth(1,it,ir)
                  enddo
              enddo
          enddo
          close(1)
      else
          write(*,*) 'No file turbulence.inp, so assuming the global'
          write(*,*) 'turbulenct width from the line input files.'
      endif
c
      end



c     --------------------------------------------------------------
c                        READ MOLECULAR ABUNDANCES
c
c     This is necessary for the line transfer. Note: Only a single
c     molecule is possible for now. The abundance is that of the 
c     molecule itself and that of its most important collision 
c     partner. The abundances are `number per as molecule', where
c     umass_av (typically 2.3) is the mean molecular weight to
c     convert densitry into number density.
c     --------------------------------------------------------------
      subroutine read_abundance()
      implicit none
c
      integer ir,it,ithetamax
      logical fex
      integer nr,nt
c      doubleprecision abund_molec,abund_colpartner
      character*80 filename
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
#include "common_lines.h"
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
c     Check the presence of the abundance.inp file
c
      filename = 'abundance.inp'
      inquire(file=filename,EXIST=fex)
      if(fex) then 
#         ifndef LINE_LOCAL_ABUNDANCE
c
          write(*,*) 'ERROR: Found file abundance.inp for local'
          write(*,*) '       molecular abundance. But code is not'
          write(*,*) '       compiled with LINE_LOCAL_ABUNDANCE.'
          stop
c
#         else /* ifndef LINE_LOCAL_ABUNDANCE */
c
          write(*,*) 'Reading abundance data from abundance.inp'
          open(unit=1,file=filename)
          read(1,*) nr,nt
          if((nr.ne.irsi_frsizex).or.(nt.ne.irsi_frsizey/2)) stop 325
          do ir=1,irsi_frsizex
              do it=1,ithetamax
                  read(1,*) abund_molec,abund_colpartner
                  locabun_abund_mol(it,ir)  = abund_molec
                  locabun_abund_coll(it,ir) = abund_colpartner
              enddo
          enddo
          locabun_done = 1
c
#         endif /* ifndef LINE_LOCAL_ABUNDANCE */
      endif
      abund_molec      = 0.d0
      abund_colpartner = 0.d0
c
      end




c     --------------------------------------------------------------
c                THE BLACKBODY PLANCK FUNCTION B_nu(T)
c
c     This function computes the Blackbody function 
c
c                    2 h nu^3 / c^2
c        B_nu(T)  = -------------------    [ erg / cm^2 s ster Hz ]
c                   exp(-h nu / kT) - 1
c
c     ARGUMENTS:
c        nu    [Hz]            = Frequency
c        temp  [K]             = Electron temperature
c     --------------------------------------------------------------
      function bplanck(temp,nu)
      implicit none
      doubleprecision temp
      doubleprecision nu
      doubleprecision bplanck
c
      if(temp.eq.0.d0) then 
          bplanck = 0.d0
          return
      endif
c
      bplanck = 1.47455d-47 * nu * nu * nu / 
     %      (exp(4.7989d-11 * nu / temp)-1.d0) + 1.d-290
c
      return
      end


c     =============================================================
c     This is a very scaled down version of the original setup
c     file. It is no longer intended to be a user-definable file and
c     will hard code just one type of setup (the old 551 setup). 
c     I retain Kees comments for now.
c
c     Klaus, 23/03/07
c
c     =============================================================



c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c                 ( MODULE: SET UP THE UNDERLYING MEDIUM )
c
c                            Leiden, May 1999
c                             C.P. Dullemond   
c
c     This module is responsible for the setup of the problem. 
c     The parameters for the analytical functions for the problem 
c     at hand will be read in and the functions will be put on 
c     the grid as well. The routines may also decide to modify
c     some of the default settings, and to in/activate certain
c     radiative processes. 
c
c     =============================================================
c                LIST OF FUNDAMENTAL CONSTANTS OF RADICAL
c
c     Parameters:
c        emisdum_kappa        The kappa for the dummy emission 
c        scati_kappa          The kappa for the dummy scattering
c        .....
c
c     Switches:
c        iradproc_dum_emis    Activate dummy emission
c        iradproc_dum_scat    Activate dummy scattering
c        iradproc_dust        Activate dust transfer
c        iradproc_lines       Activate line transfer
c        .....
c 
c     =============================================================

c     -551      Line envelope: read gas density and temperature for the lines
c     =============================================================



c     --------------------------------------------------------------
c                     READ USER-SETUP PARAMETERS
c
c     This routine reads the data from the end of the input file. 
c     It presumes that the file counter is at the appropriate 
c     position, which should automatically be the case since this
c     routine is called after the reading of the general parameters
c     by the read_radical_input() routine. The user can add her/his
c     own problem-setup by choosing an unused number (say e.g. 777)
c     for identification, and adding an "elseif(wp_type.eq.777) then"
c     line and some read statements. A new common block for private
c     use can be defined in the file "common_setup.h", so that the
c     user can specify her/his own setup parameters (e.g. 
c     common/mydisk/md_thickness,md_temp). The actual setup formulae
c     should be programmed in the routine calc_phys_var(), below. 
c
c     After the parameters have been read, the read_compute_medium()
c     routine will automatically successively call calc_phys_var()
c     to put the relevant variables on the spatial grid.
c
c     IMPORTANT: 
c     Each setup must specify which of the radiative processes
c     or radiation modes are activated. This can be done by setting
c     to 1 one of the following global variables. Note that some 
c     may require hard-switches (compiler switches) to be enabled as 
c     well (see configure.h file). This is noted in the description
c     below. Note that not setting either of the modes below produces
c     an error, since RADICAL understands that it has no reasonable
c     duty to do: no radiative processes are present. 
c
c      MAJOR MODES (EXCLUSIVE):
c
c       iradproc_compt            Compton/Thompson mode
c                                 =1 Thompson scattering
c                                 =2 Isotr Compton in diffusion limit
c                                 =3 Isotr Compton in CSK approach
c                                 (define INCLUDE_COMPTON in configure.h)
c       iradproc_line             Line (molec/atom) mode
c                                 (define INCLUDE_LINES in configure.h)
c       iradproc_dust             Dust continuum mode
c                                 =1  Full dust mode, incl dust temp solving
c                                 =-1 No dust temp solving, only emis/abs/scat
c                                 (define INCLUDE_DUST in configure.h)
c                                 NOTE: You must specify the densities 
c                                       for each of the dust species 
c                                       separately in the "dust_rho()" array.
c                                       The "density" value has no meaning
c                                       for dust continuum transfer.
c                                       
c
c      MINOR MODES (ADDITIVE): 
c                                 
c       iradproc_brems            Bremsstrahlung emis/abs added 
c       iradproc_synch            Synchrotron emis/abs added
c       iradproc_dum_emis         Simple emission/absorption where
c                                 the source terms are specified
c                                 in advance.
c       iradproc_dum_scat         Simple isotropic scattering added
c       iradproc_line_dust        (not active yet)
c       iradproc_pseudo_disk      Disk at equatorial boundary condition
c                                 (define EQUATOR_DISK)
c
c     --------------------------------------------------------------
      subroutine read_compute_medium_lean(ifile)
      implicit none
      integer ifile
c
#include "common_grid.h"
#include "common_source.h"
#include "common_setup.h"
#include "common_dust.h"
#include "common_lines.h"
#include "common_boundary.h"
c
      integer ir,it,imirt,itmax,ispec,isize,idoradtrans
      character*80 string
      logical fex
c
c     ======================================================
c           FIRST SOME THINGS THAT HAVE TO BE HANDLED
c            >>>>>> DO NOT CHANGE THESE LINES <<<<<<<
c     ======================================================
c
      call set_mirror_limits_theta(itmax,imirt)
#     ifndef INCLUDE_LINES
      if(umass_av.eq.0.d0) then 
          umass_av=-1.d66
      endif
#     endif
      ipsrc_ir     = -10
      ipsrc_itheta = -10
      ipsrc_imu    = -1000
      ipsrc_iphi   = -10
      dust_setup_nrspecies = -1
#     ifdef INCLUDE_LINES
      call read_linedata(1)
#     endif
#     ifdef INCLUDE_DUST
#     ifdef INCLUDE_DUST_ARRAYS
      do ir=1,irsi_frsizex
          do it=1,itmax
              do ispec=1,DUST_SPECIES_MAX
                  dust_rho(ispec,it,ir) = 0.d0
                  do isize=1,DUST_SIZE_MAX
                      dust_temp(isize,ispec,it,ir) = 0.d0
                  enddo
              enddo
          enddo
      enddo
#     endif
#     endif
c
      read(ifile,*,err=1001,end=1002) wp_type
c
c
c
c     ======================================================
c                USER EDITABLE PART COMES NOW:
c                PARSING THE SETUP PARAMETERS
c                 THESE LINES MAY BE CHANGED 
c     ======================================================
       if(wp_type.eq.-551) then
c         --------------------------------------------------
c         Line envelope, with externally specified density
c         and temperature.
c         NEW: 31-10-05
c         --------------------------------------------------
          idoradtrans = -1
          iradbnd_in_itype = 2
          read(ifile,*,err=1001,end=1002) iradproc_dust
          read(ifile,*,err=1001,end=1002) iradproc_line
          read(ifile,*,err=1001,end=1002) iradproc_line_dust
          read(ifile,*,err=1001,end=1002) iradproc_line_starpump
          if((iradproc_dust.ne.0).and.(iradproc_line.ne.0)) then
              write(*,*) 'ERROR: Cannot do dust and lines together...'
              stop
          endif
          if(iradproc_line.eq.0) iradproc_line_dust=0
c
c         If required, read the dust
c
          if((iradproc_dust.ne.0).or.(iradproc_line_dust.ne.0)) then
              call read_dust_density()
              call compute_sum_dust_density()
              call read_scatsource()
              if(idoradtrans.eq.0) then
c
c                 Read from dusttemp.inp file
c
                  call read_dust_temperatures(.true.)
                  if(iradproc_line_dust.eq.0) then
                      iscat_done_nonlte       = 1 ! Don't read nonlte.dat
                  else
                      iscat_done_nonlte       = 0 ! Read the line data!!
                  endif
                  iscat_update_dust_src   = 1 ! But update src from dust temp
                  iscat_no_radtrans       = 1 ! Don't allow rad trans
              elseif(idoradtrans.eq.-1) then
c
c                 Read from dusttemp_***.dat file (find out which 
c                 file from dusttemp.info)
c
                  call read_dust_temperatures(.false.)
                  if(iradproc_line_dust.eq.0) then
                      iscat_done_nonlte       = 1 ! Don't read nonlte.dat
                  else
                      iscat_done_nonlte       = 0 ! Read the line data!!
                  endif
                  iscat_update_dust_src   = 1 ! But update src from dust temp
                  if(iradproc_line.ne.0) then
                      iscat_no_radtrans   = 0 ! Allow rad trans for lines
                  else
                      iscat_no_radtrans   = 1 ! Don't allow rad trans for dust
                                              ! Use RADMC for that.
                  endif
              endif
          endif
c
c         If lines active, read mass of molecules/atoms and the 
c         gas density, temperature, as well as the abundances of 
c         the molecules
c
          if(iradproc_line.ne.0) then
              call read_line_umass()
              call read_density()
c              call read_temperature()
              call read_abundance()
              call read_velocity()
              call read_turbulence()
          endif
      else
          write(*,*) 'Dont know type ',wp_type
          stop 1
      endif
c
c     ======================================================
c            END OF PARSING THE SETUP PARAMETERS
c     ======================================================
c
c     Next a loop over the ir and it, to fill all the physical
c     variables into the big arrays
c
c     NOTE: The dust densities have been done already at the
c           end of the calc_phys_var(). Okay, it's messy; it
c           will be done more systematically in the future...
c
      do ir=1,irsi_frsizex
          do it=1,itmax
              call calc_phys_var(ir,it)
              if(setup_dens_read.eq.0) then
                  medium_arr_rho(it,ir) = medium_rho
              endif
              if(setup_temp_read.eq.0) then
                  medium_arr_el_temp(it,ir) = medium_el_temp
              endif
              if(setup_bmag_read.eq.0) then
                  medium_arr_bmag(it,ir) = medium_bmag
              endif
              medium_arr_velocity(1,it,ir) = medium_velocity(1)
              medium_arr_velocity(2,it,ir) = medium_velocity(2)
              medium_arr_velocity(3,it,ir) = medium_velocity(3)
          enddo
      enddo
c
c     Check if dust density is set properly.
c
#     ifdef INCLUDE_DUST
      if(iradproc_dust.ne.0) then
          call dust_check_grid()
      endif
#     endif
c
c     Error handling
c
      goto 1010
 1001 continue
c
c     Error: could not read line in input file
c     
      read(2,*,err=1003,end=1003) string
      write(*,*) 'ERROR in reading input file. Last line:'
      write(*,*) string
      write(*,*) 'Next few lines:'
      read(2,*,err=1004,end=1004) string
      write(*,*) string
      read(2,*,err=1004,end=1004) string
      write(*,*) string
      read(2,*,err=1004,end=1004) string
      write(*,*) string      
      read(2,*,err=1004,end=1004) string
      write(*,*) string      
      read(2,*,err=1004,end=1004) string
      write(*,*) string      
      read(2,*,err=1004,end=1004) string
      write(*,*) string      
      read(2,*,err=1004,end=1004) string
      write(*,*) string      
      stop 13
 1002 continue
      write(*,*) 'ERROR in reading input file while reading setup. ',
     %              'EOF detected'
      stop 13
 1003 continue
      write(*,*) 'FATAL ERROR in reading input file'
      stop 13
 1004 continue
      stop 13
 1010 continue      
c
      end

c     --------------------------------------------------------------
c                     RETURN THE PHYSICAL VARIABLES 
c
c     This is the function that returns the density and temperature
c     (and possibly some other local physical variables) at the 
c     position R=R(ir), Theta=Theta(it). 
c
c     It is here that the actual setup is built. Various setups are
c     possible, and one can select them by specifying the variable
c     wp_type from the common_setup.h file. These are selected in 
c     an if()-elseif()-elseif()-endif switch. You can add your own
c     setup by adding another elseif() block, and assigning a new
c     (and yet unused) number to identify your own setup. In this
c     block you must specify the following 2 variables:
c
c       - density                 = The density         [g/cm^3]
c       - temperature             = The temperature     [K]
c
c     Depending on the case at hand the other variables that must
c     be specified are:
c
c       - bgauss                  = The magnetic field  [Gauss]
c
c     If you have enabled iradproc_dust (and of course enabled
c     the define INCLUDE_DUST in the configure.h file), then
c     you must also specify the density of the dust species:
c
c       - dust_density(ispec)     = The dust density of 
c                                   species ispec       [g/cm^3]
c
c     If you want to let the number of species be determined
c     by the dust opacity input files, then simply compute density
c     and copy this into all dust_density() slots (with a do-loop
c     from 1 to DUST_SPECIES_MAX). In that case, the opacity files 
c     for the species should then be such that they are already 
c     normalized with their abundance to the gas, since their
c     density is taken equal to the gas density, while this is
c     obviously not so. But if you simply scale (in the opacity
c     files) the opacity by the abundance of the species, the
c     effect (from the radiative transfer point of view) is the
c     same as if the true opacity was taken, with the true abundance.
c     So that is why this is implemented like this here. But if
c     you wish to control the opacities directly, then simply
c     set the density of each of the species yourself. 
c
c     Then, you must switch on the radiative processes you wish.
c     They are switched by the flags in the common_source.h file
c     by the flags in the /radprocesses/ common block. These
c     flags all begin with the "iradproc_" precursor. 
c
c     Of course, your setup may wish to receive certain input 
c     parameters from the schar.inp file. You can do so by making
c     your own common block in common_setup.h, containing the
c     parameters you need. Then you must add a block in the
c     read_compute_medium() routine (above) which reads in your
c     set of parameters from the schar.inp file. 
c
c     --------------------------------------------------------------
      subroutine calc_phys_var(ir,it)
      implicit none
      integer ir,it
      integer onedim
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
#include "common_lines.h"
#include "common_boundary.h"
#include "common_dust.h"
c
c     The usual definitions
c
      doubleprecision r,theta,sint2,cost2,aa,dtheta
      doubleprecision rhobelt,rhomed,rhoblob,rhowall
      doubleprecision density,temperature,bgauss,iontemp
      doubleprecision ficke,wallexp,beltexp,blobexp,omkr
      doubleprecision hcbmu,hcbmu0,hcbmdot,t0,sonfact
      doubleprecision hscale,z,rcyl,slope,sigma
      doubleprecision dust_density(1:DUST_SPECIES_MAX)
      doubleprecision dust_temperature(1:DUST_SIZE_MAX,
     %                                 1:DUST_SPECIES_MAX)
      integer ispec,isize,ishell
c
#ifndef INCLUDE_LINES
      if(umass_av.eq.0.d0) then 
          umass_av=-1.d66
      endif
#endif
c
      r     = rsi_x_c(ir,1)
      theta = rsi_x_c(it,2)
c
c     Check which dimensions: 1-D or 2-D?
c
#     ifdef RADGRID_ONEDIM
      onedim = 1
      theta  = 1.5707963268d0    ! 1-D = Theta = pi/2 = equator
#     else
      onedim = 0
#     endif
c
c     Reset physical variables
c
      density     = 0.d0
      temperature = 0.d0
      bgauss      = 0.d0
      iontemp     = 0.d0
#     ifdef INCLUDE_DUST
#     ifdef INCLUDE_DUST_ARRAYS
      do ispec=1,DUST_SPECIES_MAX
          dust_density(ispec)=0.d0
          do isize=1,DUST_SIZE_MAX
              dust_temperature(isize,ispec)=0.d0
          enddo
      enddo
#     endif
#     endif
c
c     Some geometric variables
c
      sint2 = sin(theta)**2
      cost2 = 1.d0 - sint2
c
c     This is the point where there used to be a 'section' for each setup.
c     I'm leaving only the 551 setup.
c
c     --------------------------------------------------
c     Dusty envelope, with externally specified density
c     The density of all the species is specified in
c     the file dustdensity.inp
c     --------------------------------------------------
c     
      if(wp_type.eq.-551) then
         density            = medium_arr_rho(it,ir)
         temperature        = medium_arr_el_temp(it,ir)
         medium_velocity(1) = medium_arr_velocity(1,it,ir)
         medium_velocity(2) = medium_arr_velocity(2,it,ir)
         medium_velocity(3) = medium_arr_velocity(3,it,ir)
      else
c     --------------------------------------------------
c     >> NOT RECOGNIZED THE SETUP NUMBER <<
c     --------------------------------------------------
         write(*,*) 'Sorry, I do not know setup type',wp_type
         stop 1         
      endif
c
c     Copy the local variables "density" and "temperature" in 
c     the global common blocks
c
      medium_rho     = density
      medium_el_temp = temperature
      medium_bmag    = bgauss
#     ifdef INCLUDE_DUST
      if(setup_dustdens_read.eq.0) then
          do ispec=1,DUST_SPECIES_MAX
              dust_rho(ispec,it,ir) = dust_density(ispec)
          enddo
      endif
      if(setup_dusttemp_read.eq.0) then
          do ispec=1,DUST_SPECIES_MAX
              do isize=1,DUST_SIZE_MAX
                  dust_temp(isize,ispec,it,ir) = 
     %                 dust_temperature(isize,ispec)
              enddo
          enddo
      endif
#     endif
c

      end

c     --------------------------------------------------------------
c               READ EXTERNALLY PROVIDED DENSITY FROM FILE
c
c     For certain setup styles the density is given in a file named 
c     "density.inp". 
c     --------------------------------------------------------------
      subroutine read_density()
      implicit none
c
      integer ir,it
      logical fex
      character*80 string
      integer nr,nt,imirt,i,idum
      integer ispec
      doubleprecision r(FRSIZE_X),theta(FRSIZE_Y)
      doubleprecision dum,xi,dxi,dxirel,dumar(10),dummy
      character*80 filename
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
c
c     Check the presence of the density.inp file
c
      filename = 'density.inp'
      write(*,*) 'Reading density data from density.inp'
      inquire(file=filename,EXIST=fex)
      if(.not.fex) then 
          write(*,*) 'ERROR: density.inp not found'
          stop 1
      endif
c
c     Read the density. The density array has small Theta-size, i.e.
c     in case of MIRROR_THETA, the pi/2 to pi part is not stored. 
c     By passing 0 as second argument, we make this clear to the
c     scalar field reader.
c
      call read_scalar_field(filename,0,0,FRSIZE_X,0,FRSIZE_Y_SMALL,
     %                       medium_arr_rho,nr,nt,imirt)
      setup_dens_read  = 1
c
c     Check if the dimensions are consistent with the dimensions
c     of the radlite.inp file
c
      if((nr.ne.irsi_frsizex).or.(nt.ne.irsi_frsizey)) then
          write(*,*) 'ERROR: The dimensions of the density.inp array'
          write(*,*) '       are inconsistent with the spatial '
          write(*,*) '       dimensions specified in radlite.inp'
          write(*,*) '       or radius.inp / theta.inp. '
          write(*,*) '       Modify radlite.inp to make it consistent.'
          stop 13
      endif
c
      end


c     --------------------------------------------------------------
c               READ EXTERNALLY PROVIDED TEMPERATURE FROM FILE
c
c     For certain setup styles the density is given in a file named 
c     "temperature.inp". 
c     --------------------------------------------------------------
      subroutine read_temperature()
      implicit none
c
      integer ir,it
      logical fex
      character*80 string
      integer nr,nt,imirt,i,idum
      integer ispec
      doubleprecision r(FRSIZE_X),theta(FRSIZE_Y)
      doubleprecision dum,xi,dxi,dxirel,dumar(10),dummy
      character*80 filename
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
c
c     Check the presence of the temperature.inp file
c
      filename = 'temperature.inp'
      write(*,*) 'Reading temperature data from temperature.inp'
      inquire(file=filename,EXIST=fex)
      if(.not.fex) then 
          write(*,*) 'ERROR: temperature.inp not found'
          stop 1
      endif
c
c     Read the temperature. The temperature array has small Theta-size, i.e.
c     in case of MIRROR_THETA, the pi/2 to pi part is not stored. 
c     By passing 0 as second argument, we make this clear to the
c     scalar field reader.
c
      call read_scalar_field(filename,0,0,FRSIZE_X,0,FRSIZE_Y_SMALL,
     %                       medium_arr_el_temp,nr,nt,imirt)
      setup_temp_read = 1
c
c     Check if the dimensions are consistent with the dimensions
c     of the radlite.inp file
c
      if((nr.ne.irsi_frsizex).or.(nt.ne.irsi_frsizey)) then
          write(*,*) 'ERROR: The dimensions of temperature.inp array'
          write(*,*) '       are inconsistent with the spatial '
          write(*,*) '       dimensions specified in radlite.inp'
          write(*,*) '       or radius.inp / theta.inp. '
          write(*,*) '       Modify radlite.inp to make it consistent.'
          stop 13
      endif
c
      end


c     --------------------------------------------------------------
c                READ SPATIAL SCALAR ARRAY FROM FILE
c
c     This routine is used by read_density() etc. 
c     It is a general-purpose routine for reading scalar field
c     input files. The grid in R and Theta is found in the input
c     files radius.inp and theta.inp. 
c
c     ARGUMENTS:
c       filename             The file to read from
c       icopymirr            =1 Mirror part is copied, =0 No copy of mirror
c       nrmin                The starting index of radial dimension
c       ntmin                The starting index of theta dimension
c       nrmax                The radial maximum size of the array
c       ntmax                The theta maximum size of the array
c     
c     RESULT:
c       scalar               The spatial scalar field array to read
c       nr                   The actual number of radial points
c       nt                   The actual number of theta points
c       imirt                Whether or not mirror symmetry in Theta
c                            was present (nt is then twice the actual
c                            theta points read).
c
c     --------------------------------------------------------------
      subroutine read_scalar_field(filename,icopymirr,nrmin,nrmax,
     %                             ntmin,ntmax,scalar,
     %                             nr,nt,imirt)
      implicit none
      integer nrmin,nrmax,ntmin,ntmax,nr,nt,imirt,icopymirr
      doubleprecision scalar(ntmin:ntmax,nrmin:nrmax)
      integer itmax,it,ir
      character*80 filename
c
#include "common_grid.h"
c
c     First check if the radius.inp and theta.inp have been read.
c
#     ifdef RADGRID_ONEDIM
      spacegrid_theta_read = 1
#     endif      
      if((spacegrid_radius_read.eq.0).or.
     %     (spacegrid_theta_read.eq.0)) then
          write(*,*) 'Sorry, before reading a scalar field from '
          write(*,*) '       a file, the radius.inp and theta.inp'
          write(*,*) '       files must have been read. This is'
          write(*,*) '       enabled by specifying - as the grid '
          write(*,*) '       type on the radlite.inp file'
         stop 13
      endif
c
c     Open file
c
      open(unit=1,file=filename,status='old',err=10)
c
c     First read the dimensions of the array
c
      read(1,*) nr,nt,imirt
      itmax = nt
#     ifndef RADGRID_ONEDIM
      if(imirt.ne.0) then
          nt = 2*nt
      endif
#     else
      imirt = 0
      if(nt.ne.1) then
          write(*,*) 'ERROR while reading dustdensity.inp'
          write(*,*) '      RADICAL is compiled for 1-D, but'
          write(*,*) '      ',filename,' seems to be 2-D'
          stop 13
      endif
#     endif
c
c     Check if the dimensions do not exceed maximum ones
c
      if(nr.gt.FRSIZE_X) then
          write(*,*) 'R-size too large!'
          write(*,*) 'File has nr=',nr,' while FRSIZE_X=',FRSIZE_X
          write(*,*) 'Recompile with larger array size FRSIZE_X'
          stop 13
      endif
      if(nt.gt.FRSIZE_Y) then
          write(*,*) 'Theta-size too large!'
          write(*,*) 'File has nt=',nt,' while FRSIZE_Y=',FRSIZE_Y
          if(imirt.ne.0) then
              write(*,*) 'NOTE: With mirror symmetry nt = twice that'
              write(*,*) '      which is specified in the file.'
          endif
          write(*,*) 'Recompile with larger array size FRSIZE_Y'
          stop 13
      endif
#     ifdef MIRROR_THETA
      if(imirt.eq.0) then
          write(*,*) 'ERROR: RADICAL is compiled with MIRROR_THETA'
          write(*,*) '       So I expect that the input files'
          write(*,*) '       also has Theta-mirroring, and that it'
          write(*,*) '       is specified by 1 as third number in'
          write(*,*) '       the first line. '
          stop 13
      endif
#     endif
c
c     Now read the scalar field
c
      do ir=1,nr
          do it=1,itmax
              read(1,*,err=30,end=20) scalar(it,ir)
          enddo
          if((imirt.ne.0).and.(icopymirr.eq.1)) then
              do it=1,itmax
                  scalar(2*itmax+1-it,ir) = scalar(it,ir)
              enddo      
          endif
      enddo
c
c     Close file
c
      close(1)
c
      goto 100
   10 continue
      write(*,*) 'Could not open file ',filename
      stop 13      
   20 continue
      write(*,*) 'Prematurely reached end of file ',filename
      stop 13
   30 continue
      write(*,*) 'Could not understand contents of ',filename
      stop 13
  100 continue
c
      end



c     --------------------------------------------------------------
c          COMPUTE THE TOTAL DENSITY BY SUMMING DUST SPECIES
c
c     This function is there only for consistency with the setup.pro
c     IDL routine, which reads in BOTH the dust density and the 
c     total density file. This function has no use for the computation
c     of the dust continuum itself, which is only dependent on the
c     dust densities.
c
c     NOTE: If this routine be consistent, then make sure that the
c           dust densities represent the TRUE densities of the various
c           dust species. Not the total dust density, with the opacities
c           containing the correction for the relative abundance. It
c           is of course no problem if you do it the latter way; it's
c           just that the computed density here is then not the correct
c           one.
c
c     NOTE: Make sure that dust_setup_nrspecies is set to the nr of
c           species expected
c     --------------------------------------------------------------
      subroutine compute_sum_dust_density()
#ifdef INCLUDE_DUST
      implicit none
c
      integer ir,it,nrspecies
      integer nr,nt,itmax,i,idum
      integer ispec,nspec
c
#include "common_grid.h"
#include "common_setup.h"
#include "common_source.h"
#include "common_dust.h"
c
c     Some computations
c
      nspec = dust_setup_nrspecies
      if(nspec.eq.0) then 
          nspec = DUST_SPECIES_MAX
      endif
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
      do ir=1,irsi_frsizex
          do it=1,itmax
              medium_arr_rho(it,ir) = 0.d0
              do ispec=1,nspec
                  medium_arr_rho(it,ir) = medium_arr_rho(it,ir) +
     %                        dust_rho(ispec,it,ir)
              enddo
          enddo
      enddo
c
#endif
c     
      end

