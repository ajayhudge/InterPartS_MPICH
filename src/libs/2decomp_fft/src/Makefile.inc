#-------------------- Preprocessor options ----------------
OPTIONS = -DDOUBLE_PREC

# Choose FFT engine
FFT = fftw3

#-------------------- FFTW paths ---------------------------
# Adjust if FFTW is installed elsewhere
FFTW_PATH = /usr/local
INC = -I$(FFTW_PATH)/include
LIB_FFTW3 = -L$(FFTW_PATH)/lib -lfftw3 -lfftw3f

#-------------------- MPI Compilers ------------------------
F90 = mpif90
CC  = mpicc

#-------------------- Compiler flags -----------------------
CPPFLAGS = -cpp
OPTIM    = -O3

F90FLAGS = $(OPTIM) $(CPPFLAGS) $(INC)
CFLAGS   = $(OPTIM)
LDFLAGS  = $(OPTIM)

#-------------------- External libraries ------------------
ifeq ($(FFT),generic)
  LIBFFT =
else ifeq ($(FFT),fftw3)
  LIBFFT = $(LIB_FFTW3)
else ifeq ($(FFT),fftw3_f03)
  LIBFFT = $(LIB_FFTW3)
endif
