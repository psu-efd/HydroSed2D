.SUFFIXES: .o .f90  

COMMAND=HydroSed2D

FC=ifort
F90=ifort

HydroSed2D_DIR = $(HOME)/research/HydroSed2D_public

BOPT = O2

SDIR = $(HydroSed2D_DIR)/src
ODIR = $(HydroSed2D_DIR)/lib
XDIR = $(HydroSed2D_DIR)/bin
MODULES_DIR=$(HydroSed2D_DIR)/lib
#FFLAGS= -${BOPT} -c -openmp -extend-source 132 -check bounds
FFLAGS= -${BOPT} -c -openmp -extend-source 132 

LIBPATH=-L/usr/local/lib
SYSLIBS=-lm
 

# Fortran source code :
SOURCES= \
$(SDIR)/module.f90\
$(SDIR)/auxfuncs.f90\
$(SDIR)/auxsubs.f90\
$(SDIR)/buildMeshData.f90\
$(SDIR)/initialization.f90\
$(SDIR)/readGambitMesh.f90\
$(SDIR)/readGMSHMesh.f90\
$(SDIR)/restart_output.f90\
$(SDIR)/results_output.f90\
$(SDIR)/sediment.f90\
$(SDIR)/swe.f90\
$(SDIR)/swefvm.f90

# Object code (like F90SRCS but with .o suffix):
OBJS= \
$(ODIR)/module.o\
$(ODIR)/auxfuncs.o\
$(ODIR)/auxsubs.o\
$(ODIR)/buildMeshData.o\
$(ODIR)/initialization.o\
$(ODIR)/readGambitMesh.o\
$(ODIR)/readGMSHMesh.o\
$(ODIR)/restart_output.o\
$(ODIR)/results_output.o\
$(ODIR)/sediment.o\
$(ODIR)/swe.o\
$(ODIR)/swefvm.o

$(COMMAND):	$(OBJS)
#	$(FC) -$(BOPT) -openmp -check bounds -o $(XDIR)/$(COMMAND) $(OBJS)
	$(FC) -$(BOPT) -openmp -check -o $(XDIR)/$(COMMAND) $(OBJS)
#$(LIBPATH) $(SYSLIBS)

clean:
	-rm $(ODIR)/*
	-rm $(XDIR)/* 

.f90.o: 
	$(F90) -c -$(BOPT) $(FFLAGS) $(SDIR)/$<

$(ODIR)/auxfuncs.o : $(SDIR)/auxfuncs.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/auxsubs.o : $(SDIR)/auxsubs.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/module.o : $(SDIR)/module.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/buildMeshData.o : $(SDIR)/buildMeshData.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/initialization.o : $(SDIR)/initialization.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/readGambitMesh.o : $(SDIR)/readGambitMesh.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/readGMSHMesh.o : $(SDIR)/readGMSHMesh.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/restart_output.o : $(SDIR)/restart_output.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/results_output.o : $(SDIR)/results_output.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/sediment.o : $(SDIR)/sediment.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/swe.o : $(SDIR)/swe.f90
	$(F90) $(FFLAGS) $< -o $@
$(ODIR)/swefvm.o : $(SDIR)/swefvm.f90
	$(F90) $(FFLAGS) $< -o $@
 
