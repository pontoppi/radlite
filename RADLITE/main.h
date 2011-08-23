c===================================================================
c                        RADICAL 2-D 
c
c                          (C) 1999
c
c               C.P. Dullemond  and  R. Turolla
c
c
c This is the configuration file for RADICAL. This file is made such 
c that by reading carefully through the file from top to bottom, and 
c following the directives, one should be able to configurate the 
c code to the desired form. 
c===================================================================
c
#include "configure.h"
c
c===================================================================
c                        SYSTEM SPECIFIC STUFF
c===================================================================
c
#define OPENACC access
c
c-------------------------------------------------------------------
c Specific to the XL Fortran compiler for AIX
c-------------------------------------------------------------------
c
#ifdef COMPILER_XLAIX
#undef OPENACC
#define OPENACC position
#define sign sgnf
#endif
c
c===================================================================
c                         VERSION PARAMETER
c===================================================================
c
c Date: 14-Oct-1999 
c#define VERSION_1_1
c
c===================================================================
c                           GENERAL CHECKS
c===================================================================
c
c  THIS IS A LINE-ONLY CODE (KEES & KLAUS 22.03.07)
c
#ifndef INCLUDE_LINES
# ERROR
#endif
#ifndef INCLUDE_DUST
# ERROR
#endif
#ifndef LINE_DUSTCONT
# ERROR
#endif
c
c-------------------------------------------------------------------
c First check if the most important defines are indeed defined
c-------------------------------------------------------------------
c
#ifndef FRSIZE_X
#  error : Need to define FRSIZE_X
#endif
#ifndef FRSIZE_Y
#  error : Need to define FRSIZE_Y
#endif
c#ifndef COORD_SPHERICAL
c#ifndef COORD_CARTESIAN
c#ifndef COORD_CYLINDRICAL
c#   ERROR: Must specify either COORD_SPHERICAL, COORD_CARTESIAN or COORD_CYLINDRICAL
c#endif
c#endif
c#endif
#ifndef COORD_SPHERICAL
#   ERROR: Must specify COORD_SPHERICAL for now
#endif
c
c-------------------------------------------------------------------
c If 1-D, then the y-size of the array should be 1
c-------------------------------------------------------------------
c
#ifdef RADGRID_ONEDIM
#undef FRSIZE_Y
#define FRSIZE_Y 1
#endif
c
c-------------------------------------------------------------------
c Use new 1-D SC option: add point at s=0
c-------------------------------------------------------------------
c
c#define SCONEDIM_ADDPOINT
c
c===================================================================
c                     CONSISTENCY CHECKS
c===================================================================
c
c-------------------------------------------------------------------
c Cartesian coordinates
c-------------------------------------------------------------------
c
#ifdef COORD_CARTESIAN
#ifndef RADGRID_ONEDIM
#     ERROR: CARTESIAN ONLY ALLOWED IN 1-D SO FAR
#endif
#ifdef TELESCOPE_CAMERA
#     ERROR: IF CARTESIAN COORDINATES, THEN NO TELESCOPE_CAMERA!!
#endif
#endif
c
c-------------------------------------------------------------------
c Cylindrical coordinates -----> not active
c-------------------------------------------------------------------
c
#ifdef COORD_CYLINDRICAL
#     ERROR: CYLINDRICAL COORDINATES NOT YET ALLOWED
#endif
c
c-------------------------------------------------------------------
c Check the dimensions 
c-------------------------------------------------------------------
c
#ifdef RADGRID_ONEDIM
#ifdef RADGRID_TWODIM
#       ERROR: Cannot specify both RADGRID_ONEDIM and RADGRID_TWODIM
#endif
#endif
#ifndef RADGRID_ONEDIM
#ifndef RADGRID_TWODIM
#       ERROR: Must specify either RADGRID_ONEDIM or RADGRID_TWODIM
#endif
#if (FRSIZE_Y+2)/2==(FRSIZE_Y+1)/2
#       ERROR:  FRSIZE_Y SHOULD BE EVEN
#endif
#endif
#if (FRSIZE_PHI+2)/2==(FRSIZE_PHI+1)/2
#       ERROR: FRSIZE_PHI SHOULD BE EVEN
#else
#if (FRSIZE_PHI/2+2)/2==(FRSIZE_PHI/2+1)/2
#       ERROR: FRSIZE_PHI NOT DIVISIBLE BY 4
#endif
#endif
c
c-------------------------------------------------------------------
c Check interpolation schemes
c-------------------------------------------------------------------
c
#ifdef INTERPOL_ORDER_1
#ifdef INTERPOL_THETA_3
#error : FIRST ORDER, SO NO THETA_3
#endif
#ifdef INTERPOL_PHI_3
#error : FIRST ORDER, SO NO PHI_3
#endif
#endif
c
c-------------------------------------------------------------------
c If onedimensional --> Only ESC_COMPLETE
c-------------------------------------------------------------------
c
#ifdef RADGRID_ONEDIM
#undef ESC_MINIMAL_EXTENSION 
#undef ESC_COMPLETE
#define ESC_COMPLETE
#endif
c
c===================================================================
c                end of consistency checks
c===================================================================
c
c
c
c
c
c===================================================================
c                 MEMORY SAVING DEFINES
c===================================================================
c
c-------------------------------------------------------------------
c If we want RADICAL only to do the raytracing stuff (i.e. no
c ESC/non-LTE/scattering stuff to be computed), then defining
c this one will remove all that stuff and save ENORMOUS amounts
c of memory.
c-------------------------------------------------------------------
c
#ifdef ONLY_RAY_TRACING
#define NO_SHORT_CHARS
#define NO_INTENSITY_STORAGE
#undef LC_INTEGRATION
#endif
c
c-------------------------------------------------------------------
c If the Long characteristics option is active, then make sure that
c the short characteristics are off.
c-------------------------------------------------------------------
#ifdef LC_INTEGRATION
#ifdef ESC_MINIMAL_EXTENSION 
#error LC and MESC are mutually exclusive
#endif
#ifdef ESC_COMPLETE
#error LC and ESC are mutually exclusive
#endif
#define NO_SHORT_CHARS
#endif
c
c-------------------------------------------------------------------
c Memory saving compact shortchar storage stuff
c-------------------------------------------------------------------
c
#ifndef NO_SHORT_CHARS
c
c Activate?
c
#define SHORCHAR_COMPACT_ARRAYS
c
c If active...
c
#ifdef SHORCHAR_COMPACT_ARRAYS
#define FRSIZE_CHR1 1
#ifndef LONGESC_FRACTION
#define LONGESC_FRACTION 0.01
#endif
#else
#define FRSIZE_CHR1 FRSIZE_CHAR
#endif
c
#define FRSIZE_CHAR_BIG (10*FRSIZE_CHAR)
#endif
c
c-------------------------------------------------------------------
c Then we must take care of the mirroring and stuff.
c-------------------------------------------------------------------
c
#ifndef RADGRID_ONEDIM
c
c-------------------------------------------------------------------
c 2-D case
c-------------------------------------------------------------------
c
#ifndef MIRROR_THETA
#define FRSIZE_Y_SMALL FRSIZE_Y
#else
#define FRSIZE_Y_SMALL (FRSIZE_Y/2)
#endif
#ifndef MIRROR_PHI
#define FRSIZE_PHI_SMALL FRSIZE_PHI
#else
#define FRSIZE_PHI_SMALL (FRSIZE_PHI/2)
#endif
#define FRSIZE_INF_PHI_SMALL FRSIZE_INF_PHI
#else
c
c-------------------------------------------------------------------
c 1-D case
c-------------------------------------------------------------------
c
#undef MIRROR_THETA
#undef MIRROR_PHI
#define FRSIZE_Y_SMALL 1     
#define FRSIZE_PHI_SMALL 1   
#define FRSIZE_INF_PHI_SMALL 1  
c
#endif
c
c-------------------------------------------------------------------
c For small memory: all frequencies use same memory
c-------------------------------------------------------------------
c 
#ifdef SMALL_MEMORY
# define FRSIZE_INF_FREQ 1
#else
# define FRSIZE_INF_FREQ FRSIZE_FREQ
#endif
c
c===================================================================
c               end of memory saving defines
c===================================================================
c
c
c
c
c===================================================================
c                   LINE TRANSFER DEFINES 
c===================================================================
c
#ifdef INCLUDE_LINES
c
c-------------------------------------------------------------------
c From now on, always put on the LINE_VELOCITIES
c-------------------------------------------------------------------
c
#define LINE_VELOCITIES
c
c-------------------------------------------------------------------
c If line transfer, then sometimes the dust arrays are useful
c anyway (for instance if dust opacity is included in line trans)
c-------------------------------------------------------------------
c
#define INCLUDE_DUST_ARRAYS
c
c-------------------------------------------------------------------
c Check if the size of the line profiles are consistent
c-------------------------------------------------------------------
c
#ifndef SZ_LINEPROFILE_SMALL
#define SZ_LINEPROFILE_SMALL SZ_LINEPROFILE
#endif
#if SZ_LINEPROFILE<SZ_LINEPROFILE_SMALL 
#error SZ_LINEPROFILE_SMALL must be smaller than SZ_LINEPROFILE
#endif
c
c-------------------------------------------------------------------
c Switch on the line-variant of Ng acceleration, which acts on 
c the Jbar, rather than the J or the S itself.
c-------------------------------------------------------------------
c
#define LINE_NG_ACCEL
c
c-------------------------------------------------------------------
c Line transfer can in principle be done in two ways. One way
c is to make use of the usual scat_src arrays etc. This excludes
c systematic velocities. The better way (which should be the
c default!!) is to define the macro LINE_VELOCITIES, which enables
c the real line transfer: including R,Theta-dependent systematic
c velocities and much more fancy stuff. This mode also EXCLUDES
c per default the scat_src, intmom_0 etc arrays, in order to save
c memory. This is done here: by setting the macro JSRC_NO_ARRAYS.
c You can undo this line, which has no serious consequences. But
c you may run into memory problems in 2-D line transfer calculations.
c
c MODIFIED 21-10-05: If LINE_DUSTCONT is defined, then I do still
c    include the arrays, because then I need to include the dust
c    data to compute the continuum...
c-------------------------------------------------------------------
c
#ifdef LINE_VELOCITIES
# ifndef LINE_DUSTCONT
# define JSRC_NO_ARRAYS
# endif
# define LONGCHAR_EXTRA_POINTS
# define INTENSITY_ARRAY
# define PHYSVAR_VELOCITIES
# define LINE_LOCAL_ABUNDANCE
# define LINE_LOCAL_WIDTH 
cc# define LINE_PROFILE_ARRAY
# ifdef INTARRAY_FREQ
#  error : The intensity array is already used by something else
# endif
# define INTARRAY_FREQ SZ_LINEPROFILE_SMALL
# ifdef NO_INTENSITY_STORAGE
#  error : for line transfer intensity storage is crucial
# endif
#endif
c
#endif 
c===================================================================
c               end of line transfer defines
c===================================================================
c
c-------------------------------------------------------------------
c If line transfer is not activated, be sure that everything
c involving line transfer is indeed de-activated
c-------------------------------------------------------------------
c
#ifndef INCLUDE_LINES
# undef LINE_NG_ACCEL
# undef LINE_LOCAL_ABUNDANCE
# undef LINE_LOCAL_WIDTH
# undef LINE_VELOCITIES
#endif
c
c-------------------------------------------------------------------
c Dust and lines are incompatible
c-------------------------------------------------------------------
c
c#ifdef INCLUDE_LINES
c#ifdef INCLUDE_DUST
c#error : DUST and LINES are incompatible      
c#endif
c#endif
c
c-------------------------------------------------------------------
c Dust defines
c-------------------------------------------------------------------
c
#define DUST_RADICAL
#ifndef DUST_TRANGE_MAX
#define DUST_TRANGE_MAX 1
#endif
#define DUST_RHOMIN 1d-99
c
c-------------------------------------------------------------------
c Defines for the CSK stuff
c-------------------------------------------------------------------
c
#ifndef NGAM1
#   define NGAM1 ISOCSK_NGAM1
#endif
#ifdef CSK_TAU_INTERPOL
#define NRTAUMAX CSK_TAU_INTERPOL
#else
#define NRTAUMAX ((FRSIZE_X+1)*(FRSIZE_Y_SMALL+1))
#endif
c
c
c-------------------------------------------------------------------
c Define for grid setup stuff
c-------------------------------------------------------------------
c
#define MUTHRES 1.d-6
#define FRSIZE_MU_HALF (FRSIZE_MU/2)
c
c-------------------------------------------------------------------
c Frequency grid define (for dust tridiagonal ALI, future).
c-------------------------------------------------------------------
c
#ifndef DUST_TRIDAG_ALI
#define FREQ_START 1
#else
#define FREQ_START 0
#endif
#ifndef DUST_TEMP_MIN 
#define DUST_TEMP_MIN 0.02d0
#endif
#ifndef DUST_TEMP_MAX 
#define DUST_TEMP_MAX 3000.d0
#endif
#ifndef DUST_TEMP_ACCUR
#define DUST_TEMP_ACCUR 0.01d0
#endif
c
c-------------------------------------------------------------------
c Version number
c-------------------------------------------------------------------
c
#define VERSION "V. 0.3"
c
c-------------------------------------------------------------------
c Find the maximum dimension
c-------------------------------------------------------------------
c
#if FRSIZE_X > FRSIZE_Y           /* Be sure to have the correct one... */
#   define FRSIZE_MAX FRSIZE_X    /* This is the biggest mesh size. */
#else
#   define FRSIZE_MAX FRSIZE_Y    /* This is the biggest mesh size. */
#endif
c
c-------------------------------------------------------------------
c Some constants
c-------------------------------------------------------------------
c
#define LOWER_DOUBLE (1e-33)
#define E_LOW       (1e-12)
#define PICONST 3.1415926535897932385d0
#define TEMPCMB 2.728d0
#define SYNCH_LOWER_T_CUTOFF (5.d8)
c
c-------------------------------------------------------------------
c Stuff for the interpolations
c-------------------------------------------------------------------
c
#ifdef MIRROR_THETA
#define IDX0 0
#define IDX1 1
#define IDX2 2
#define IDX3 3
#else
#define IDX0 0
#define IDX1 0
#define IDX2 0
#define IDX3 0
#endif
c
c-------------------------------------------------------------------
c TELESCOPE/LONGCHAR settings
c-------------------------------------------------------------------
c
#ifndef TELESC_EPS
#define TELESC_EPS 1.d-6
#endif
