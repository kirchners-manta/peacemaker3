SOURCES = src/kinds.f90 \
		  src/lengths.f90 \
		  src/constants.f90 \
		  src/iso_varying_string.f90 \
		  src/auxiliary.f90 \
		  src/error.f90 \
		  src/atomic_data.f90 \
		  src/config.f90 \
		  src/input.f90 \
		  src/cluster.f90 \
		  src/info.f90 \
		  src/polynomial.f90 \
		  src/shared_data.f90 \
		  src/partition_functions.f90 \
		  src/thermo.f90 \
		  src/qce.f90 \
		  src/main.f90

all: release

FC = gfortran
FFLAGS = -cpp
OPENMPFLAGS = -fopenmp
DEBUGFLAGS = -Og -std=f2008 -pedantic -Wall -Wextra -fcheck=all -ggdb -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal
PROFILEFLAGS = -pg
RELEASEFLAGS = -O3 -flto

#FC = ifort
#FFLAGS = -fpp
#OPENMPFLAGS = -qopenmp
#DEBUGFLAGS = -O0 -g -traceback -std08 -warn all -check all -debug all -ftrapuv
#PROFILEFLAGS = -p
#RELEASEFLAGS = -O3 -xHost -ipo

debug:
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(DEBUGFLAGS) ${SOURCES} -o peacemaker
debug_serial:
	$(FC) $(FFLAGS) $(DEBUGFLAGS) ${SOURCES} -o peacemaker
release:
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(RELEASEFLAGS) ${SOURCES} -o peacemaker
release_serial:
	$(FC) $(FFLAGS) $(RELEASEFLAGS) ${SOURCES} -o peacemaker
profile:
	$(FC) $(FFLAGS) $(PROFILEFLAGS) ${SOURCES} -o peacemaker
clean:
	rm -rvf *.o *.mod
