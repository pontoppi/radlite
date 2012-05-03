#==========================================================================
#                         RADLite Makefile
#
#                              (C) 2009
#
#               K.M. Pontoppidan & C.P. Dullemond
#
#==========================================================================
#
#
#--------------------------------------------------------------------------
#                           GENERAL STUFF:
#--------------------------------------------------------------------------

# Name of executable and main source
EXENAME  = RADlite

# Optimize
OPTIM = -O3

# Which compiler to use
#FF = ifort -132
FF = gfortran -m64 -ffixed-line-length-132 -fpic 


# Glueing all the things together
FFLAGS   = $(OPTIM) 
#-Wall
LFLAGS   = $(OPTIM) 
#-Wall

# The lists of the files 
COMMONS = common_main.h common_boundary.h common_grid.h common_setup.h common_radfield.h common_shortchar.h common_source.h common_telescope.h common_lines.h common_dust.h common_check.h common_isave.h
SRCMAIN = main.F main.h configure.h
OBJMAIN = main.o
SRCGRID = grid.F main.h configure.h
OBJGRID = grid.o
SRCSHORTCHAR = shortchar.F main.h configure.h
OBJSHORTCHAR = shortchar.o
SRCTRANSFER = transfer.F main.h configure.h
OBJTRANSFER = transfer.o
SRCINTERPOL = interpol.F main.h configure.h
OBJINTERPOL = interpol.o
SRCTELESCOPE = telescope.F main.h configure.h
OBJTELESCOPE = telescope.o
SRCSOURCE = source.F main.h configure.h
OBJSOURCE = source.o
SRCSTAR = star.F main.h configure.h
OBJSTAR = star.o
SRCSETUP = setup.F main.h configure.h
OBJSETUP = setup.o
SRCLINE = line.F main.h configure.h
OBJLINE = line.o
SRCNRECIP = nrecip.F main.h configure.h
OBJNRECIP = nrecip.o
SRCDUST = dust.F main.h configure.h
OBJDUST = dust.o
SRCTOOLS = tools.F
OBJTOOLS = tools.o
OBJ      = $(OBJMAIN) $(OBJGRID) $(OBJSHORTCHAR) $(OBJTRANSFER) $(OBJSOURCE) $(OBJEMISSRC) $(OBJINTERPOL) $(OBJTELESCOPE) $(OBJSCATSRC) $(OBJITERATION) $(OBJSTAR) $(OBJSETUP) $(OBJVETMETHOD) $(OBJMONTECARLO) $(OBJLINE) $(OBJDUST) $(OBJCOMPTON) $(OBJNRECIP) $(OBJSOLVETE)  $(OBJAUTOINSPECT) $(OBJTOOLS)
#################################################
#                   RULES                       #
#################################################

#all:	$(EXENAME) Makefile
#	./install.perl

RADlite: $(OBJ) Makefile
	$(FF) $(LFLAGS) $(OBJ) $(LIBS) -o $@ 

main.o: $(SRCMAIN) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  main.F -o $@

grid.o: $(SRCGRID) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  grid.F -o $@

shortchar.o: $(SRCSHORTCHAR) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  shortchar.F -o $@

transfer.o:  $(SRCTRANSFER) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  transfer.F -o $@

interpol.o:  $(SRCINTERPOL) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  interpol.F -o $@

telescope.o:  $(SRCTELESCOPE) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  telescope.F -o $@

source.o: $(SRCSOURCE) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  source.F -o $@

thermsrc.o: $(SRCEMISSRC) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  thermsrc.F -o $@

scatnonlte.o: $(SRCSCATSRC) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  scatnonlte.F -o $@

star.o: $(SRCSTAR) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  star.F -o $@

setup.o: $(SRCSETUP) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  setup.F -o $@

line.o:  $(SRCLINE) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  line.F -o $@

dust.o:  $(SRCDUST) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  dust.F -o $@

nrecip.o: $(SRCNRECIP) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  nrecip.F -o $@

tools.o: $(SRCTOOLS) Makefile $(COMMONS)
	$(FF) -c $(FFLAGS)  tools.F -o $@
