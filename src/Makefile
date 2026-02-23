# Fortran compiler
FC = mpifort
# Compiler flags
# FFLAGS := -O3 -cpp -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -mcmodel=medium -fcheck=all -fbacktrace -fno-range-check -g
FFLAGS := -cpp -ffree-line-length-none -fdefault-real-8 -mcmodel=medium -fcheck=all -fbacktrace -g
# FFLAGS += -DDOUBLE_PREC
BIG :=
DBG :=
OPT :=
PROF :=
OMP :=

# 2decomp&fft
include libs/2decomp_fft/src/Makefile.inc
INCLUDE = -I libs/2decomp_fft/include

# FFTW (local)
FFTW_DIR = /usr/local
LIBFFTW = -L$(FFTW_DIR)/lib -lfftw3 -lm
INCFFTW = -I$(FFTW_DIR)/include

# Libraries
LIB = -Llibs/2decomp_fft/lib -l2decomp_fft -Llibs/fft/lib -lfft $(LIBFFTW)

TARGET = iparts

SRC = gridsize.f90 param.f90 common.f90 init.f90 initmpi.f90 initsolver.f90 bound.f90 chkdiv.f90 chkdt.f90 loadd.f90 mom.f90 scal.f90 test_sources.f90 rk3.f90 tests.f90 fillps.f90 zredistribute.f90 solver.f90 correc.f90 initparticles.f90 coordsfp.f90 kernel.f90 intgr_over_sphere.f90 forcing.f90 interp_spread.f90 collisions.f90 intgr_nwtn_eulr.f90 intgr_nwtn_eulr_subc.f90 phase_indicator.f90 post.f90 output.f90 fftw.f90 autocorr.f90 main.f90

OBJ = $(SRC:.f90=.o)

# Build rules
all: libraries $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(BIG) $(DBG) $(OPT) $(PROF) $(OMP) -o $@ $(OBJ) $(LIB)

%.o: %.f90
	$(FC) $(INCLUDE) $(INCFFTW) $(FFLAGS) $(BIG) $(DBG) $(OPT) $(PROF) $(OMP) -c $<

# Utilities
.PHONY: clean veryclean libraries

clean:
	rm -f *.o *.mod $(TARGET)

veryclean:
	rm -f *.o *.mod $(TARGET) .depend
	cd libs/fft && rm -f fft.o lib/libfft.a
	cd libs/2decomp_fft/src && make clean

libraries:
	cd libs/fft && \
	$(FC) $(FFLAGS) -c fft.f && \
	ar rcs libfft.a fft.o && mkdir -p lib && mv libfft.a lib && cd ../../ && \
	cd libs/2decomp_fft/src && make clean && make FC=$(FC) && cd ../../../
# 	cd libs/2decomp_fft/src && make clean && make FC=$(FC) OPTIONS="-DDOUBLE_PREC" && cd ../../../

.depend dep:
	./.makedepo $(SRC) > .depend

# -include .depend
