PROGRAM=distrlfp

FC   = gfortran 
###FC   = ifort #gfortran 
#####version := $(shell $(FC) --version )
#####ifeq ($(findstring ifort,$(version)),ifort)
#####  ftype := ifort
#####endif
#####ifeq ($(findstring GNU,$(version)),GNU)
#####  ftype := gfortran
#####endif
#####
#####ifndef ftype
#####nofortran:
#####	@echo no fortran compiler ifort or gfortran found
#####	@exit 1
#####endif

#ifeq ($(ftype),gfortran)
  FFLAGS := -g -O2 -ffixed-line-length-none -fdefault-real-8 -fopenmp
#endif

####ifeq ($(ftype),ifort)
####  FFLAGS := -g -r8 -O3 -132 -traceback -fopenmp#-CB-
####endif

all: $(PROGRAM)

OBJS= lfp.o loadd.o param.o

SRCS = $(OBJS:.o=.f)

$(PROGRAM): $(OBJS)
	$(FC) -o $@ $(OBJS) $(FFLAGS)

clean:
	rm -f core $(OBJS) *.mod $(PROGRAM)

.depend dep:
	./makedepo $(SRCS) > .depend
	
%.o: %.f
	$(FC) -c -o $@ $(FFLAGS) $<

include .depend
