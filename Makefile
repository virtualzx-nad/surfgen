# makefile for surfgen, the coupled potential energy surface generator
# 2013, Xiaolei Zhu, Yarkony Group, Johns Hopkins University
# currently tested on Mac OSX and linux ONLY
# currently tested for gfortran and ifort ONLY
# need LAPACK and BLAS.  tested for OSX vecLib and Intel MKL ONLY

#note: it works on Mac, but you might as well want to use the XCode build
#      rules that comes with the repository.

# Objects needed for standalone fitting program
OBJSf   = hddata.o diis.o rdleclse.o combinatorial.o progdata.o libutil.o \
            libsym.o libinternal.o localcoord.o makesurf.o io.o surfgen.o 

# Objects needed for runtime interface library
OBJSLf  = hddata.o combinatorial.o progdata.o libutil.o libsym.o libinternal.o\
            io.o potlib.o

# Set surfgen vesion
SGENVER := 2.2.2

# Get the OS name and version
UNAME := $(shell uname -a)
OS    := $(word 1,$(UNAME))
OSV   := $(word 3,$(UNAME))
ARC   := $(shell uname -m)
$(info Operating system: $(OS))
$(info OS Version: $(OSV))

# Set up directories
bindir:= bin
libdir:= lib
JDIR  := $(shell pwd)
BDIR  := $(JDIR)/$(bindir)
LDIR  := $(JDIR)/$(libdir)
SDIR  := $(JDIR)/source
CDS   := cd $(SDIR);

OBJS  := $(addprefix $(SDIR)/,$(OBJSf))
OBJSL := $(addprefix $(SDIR)/,$(OBJSLf))

# Set compiler
CPLIST := g95 pgf90 /usr/bin/gfortran gfortran 
ifdef FC  #use predefined fortran compiler
  COMPILER := $(FC)
  $(info Using Fortran compiler : $(FC))
else       #find default compilers
  $(error Compiler NOT defined! Please set it with variable FC)
endif

# set up product name
EXEC  := $(BDIR)/surfgen-$(SGENVER)-$(OS)-$(ARC)
LIBF  := $(LDIR)/libsurfgen-$(SGENVER)-$(OS)-$(ARC).a

# set debugging flags
ifeq ($(DEBUGGING_SYMBOLS),YES)
  ifndef DEBUGFLAG
   ifneq ($(findstring gfortran,$(COMPILER)),)
    DEBUGFLAG = -g -fbounds-check -fbacktrace -Wall -Wextra
   else
    ifneq ($(findstring ifort,$(COMPILER)),)
        DEBUGFLAG = -g -check uninit -check bound -check pointers -traceback -debug
    else
        DEBUGFLAG = -g
    endif
   endif
  endif
else
  DEBUGFLAG := 
endif

# set BLAS and LAPACK libraries.  
# Please set variable $BLAS_LIB or $LIBS to enable the program to find the link line.
# On Mac the default is to use the vecLib framework
ifndef LIBS
 ifndef BLAS_LIB
  ifeq ($(OS),Darwin)  
     #on mac, use frameworks
     LIBS := -framework vecLib
     $(info Using vecLib framework for Mac OS X)
  else
     $(info BLAS_LIB not set.  Trying to determine LAPACK link options...)
     #BLAS_LIB is not set.  check LD_LIBRARY_PATH for mkl
     ifneq ($(findstring mkl,$(LD_LIBRARY_PATH)),)
        $(info Found mkl in LD_LIBRARY_PATH. Using dynamic link to MKL.)
        LIBS := -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
     endif #ifneq mkl,LD_...
  endif #OS==Darwin
 else
    LIBS := $(BLAS_LIB)
 endif #BLAS_LIB
endif #ifndef $LIBS
ifndef LIBS
  $(info Warning:  Lapack library link line options not determined. \
     Use variable LIBS or BLAS_LIB to define these libraries.)
endif

    RM := rm -rf

ifndef AR
    AR := ar
endif

# set default compiler flags
ifneq ($(findstring ifort,$(COMPILER)),)
  ifdef NO_I8
    CPOPT    := -auto -c -assume byterecl -parallel -O3 -lpthread -openmp
  else
    CPOPT    := -auto -c -assume byterecl -parallel -O3 -lpthread -openmp -xHost -no-opt-matmul -i8
  endif
  LKOPT    := -auto -lpthread -parallel
else
 ifneq ($(findstring gfortran,$(COMPILER)),)
  ifneq (,$(filter $(NO_I8),YES yes))
    CPOPT    := -fopenmp -O3 
  else
    CPOPT    := -fopenmp -O3 -m64
  endif
  LKOPT    :=
 else
  CPOPT    := 
  LKOPT    :=
 endif
endif

# build everything
all  :  surfgen libs
	@echo 'Finished building surfgen.'
	@echo ''

#target for standalone fitting and analysis program
surfgen  :  $(OBJS) | $(BDIR)
	@echo '-----------------------------------------'
	@echo '   SURFGEN FITTING PROGRAM '
	@echo 'Program Version:     $(SGENVER)'
	@echo 'Executable Name:     $(EXEC)'
	@echo 'Executable Saved to: $(BDIR)'
	@echo 'BLAS/LAPACK LIB:     $(LIBS)'
	@echo 'Debug Flag:          $(DEBUGFLAG)'
	@echo 'Compiler options:    $(CPOPT)'
	@echo 'Linking options:     $(LKOPT)'
	@echo '-----------------------------------------'
	@echo 'Building target:     $@'
	@echo 'Invoking: Linker'
	$(CDS) $(COMPILER) -o $(EXEC) $(OBJS) $(LIBS) $(LKOPT) $(LDFLAGS)
	@echo 'Finished building target: $@'
	@echo '-----------------------------------------'
	@echo ''
	@echo ''

#target for runtime interface library
libs  :  $(OBJSL) | $(LDIR)
	@echo '-----------------------------------------'
	@echo '  SURFGEN EVALUATION LIBRARY'
	@echo 'Program Version:     $(SGENVER)'
	@echo 'Archiver (AR):       $(AR)'
	@echo 'Library File Name:   $(LIBF)'
	@echo '-----------------------------------------'
	@echo 'Building target: $@'
	@echo 'Archiving the object into library '
	$(CDS) $(AR) -r -v  $(LIBF) $(OBJSL)
	@echo '-----------------------------------------'
	@echo ''
	@echo ''

#
clean:
	-$(RM) $(OBJS) $(SDIR)/*.mod $(SDIR)/potlib.o 
	@echo 'Finished cleaning'

#
$(SDIR)/%.o : $(SDIR)/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Compiler'
	$(CDS) $(COMPILER) -c -o $@ $< $(CPOPT) $(DEBUGFLAG) $(FFLAGS)
	@echo 'Finished building: $<'
	@echo ' '

$(BDIR) $(LDIR) :
	@echo 'Creating directory $@'
	@mkdir -p $@
