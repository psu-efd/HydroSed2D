.SUFFIXES: .o .f90  

COMMAND=DEMInterpolation

FC=gfortran
F90=gfortran

HydroSed2D_DIR = $(HOME)/research/HydroSed2D_public

SDIR = $(HydroSed2D_DIR)/tools/DEMInterpolation
ODIR = $(SDIR)
XDIR = $(SDIR)
MODULES_DIR=$(SDIR)
FFLAGS= -M $(MODULES_DIR) -I $(MODULES_DIR) -ffree-line-length-none -fopenmp -fbounds-check
#FFLAGS= -I $(MODULES_DIR) -ffree-line-length-none -openmp
#FFLAGS= -M $(MODULES_DIR) -I $(MODULES_DIR) -ffree-line-length-none

BOPT = O2

LIBPATH=-L/usr/local/lib
SYSLIBS=-lm
 

# Fortran source code :
SOURCES= \
$(SDIR)/module.f90\
$(SDIR)/readGambitMesh.f90\
$(SDIR)/readGMSHMesh.f90\
$(SDIR)/DEMInterpolation.f90

# Object code (like F90SRCS but with .o suffix):
OBJS= \
$(ODIR)/module.o\
$(ODIR)/readGambitMesh.o\
$(ODIR)/readGMSHMesh.o\
$(ODIR)/DEMInterpolation.o

$(COMMAND):	$(OBJS)
	$(FC) -${BOPT} -fopenmp -o $(XDIR)/$(COMMAND) $(OBJS)
#	$(FC) $(FFLAGS) -o $(XDIR)/$(COMMAND) $(OBJS)
#$(LIBPATH) $(SYSLIBS)

clean:
	-rm $(ODIR)/*.o 
	-rm $(ODIR)/*.mod
	-rm $(XDIR)/$(COMMAND)

.f90.o: 
	$(F90) -c -$(BOPT) $(FFLAGS) $(SDIR)/$<

$(ODIR)/module.o : $(SDIR)/module.f90
	$(F90) -c -$(BOPT) $(FFLAGS) $<
$(ODIR)/readGambitMesh.o : $(SDIR)/readGambitMesh.f90
	$(F90) -c -$(BOPT) $(FFLAGS) $<
$(ODIR)/readGMSHMesh.o : $(SDIR)/readGMSHMesh.f90
	$(F90) -c -$(BOPT) $(FFLAGS) $<
$(ODIR)/DEMInterpolation.o : $(SDIR)/DEMInterpolation.f90
	$(F90) -c -$(BOPT) $(FFLAGS) $<
 
