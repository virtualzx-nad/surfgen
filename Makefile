# makefile for surfgen, the coupled potential energy surface generator
# 2013, Xiaolei Zhu, Yarkony Group, Johns Hopkins University
# currently tested on Mac OSX and linux ONLY
# currently tested for gfortran and ifort ONLY
# need LAPACK and BLAS.  tested for OSX vecLib and Intel MKL ONLY

#note: it works on Mac, but you might as well want to use the XCode build
#      rules that comes with the repository.

# Objects needed for standalone fitting program
OBJSf   = hddata.o diis.o rdleclse.o combinatorial.o progdata.o libutil.o \
            libsym.o libinternal.o localcoord.o makesurf.o io.o linsteps.o surfgen.o 

# Objects needed for runtime interface library
OBJSLf  = hddata.o combinatorial.o progdata.o libutil.o libsym.o libinternal.o\
            io.o potlib.o

# Objects needed for test programs 
OBJTf   =  hddata.o diis.o rdleclse.o combinatorial.o progdata.o libutil.o \
           libsym.o libinternal.o localcoord.o makesurf.o linsteps.o io.o testsurfgen.o

# Modules stored in the surfgen libraries
#MODLf   =  hddata.mod potdata.mod potdata.mod cnpi.mod  

# Objects for version labeling
OBJVf   =  getver.o

# manual pages
PDFfl   =  surfgen.pdf surfgen.in.pdf points.in.pdf coord.in.pdf 


# Set surfgen vesion
SGENVER := 2.5.15

# Get the OS name and version
UNAME := $(shell uname -a)
OS    := $(word 1,$(UNAME))
OSV   := $(word 3,$(UNAME))
ARC   := $(shell uname -m)
PS2PDF:= ps2pdf
$(info Operating system: $(OS))
$(info OS Version: $(OSV))

# Set up directories
bindir:= bin
libdir:= lib
tstdir:= test
srcdir:= source
docdir:= pdf
JDIR  := $(shell pwd)
BDIR  := $(JDIR)/$(bindir)
LDIR  := $(JDIR)/$(libdir)
TDIR  := $(JDIR)/$(tstdir)
SDIR  := $(JDIR)/$(srcdir)
DDIR  := $(JDIR)/$(docdir)
MANDIR:= $(JDIR)/man
CDS   := cd $(SDIR);

OBJS  := $(addprefix $(SDIR)/,$(OBJSf))
OBJSL := $(addprefix $(SDIR)/,$(OBJSLf))
OBJT  := $(addprefix $(SDIR)/,$(OBJTf))
OBJV  := $(addprefix $(SDIR)/,$(OBJVf))
PDFPG := $(addprefix $(DDIR)/,$(PDFfl))
#MODL  := $(addprefix $(SDIR)/,$(MODLf))

PDFsrc := $(PDFfl:.pdf=.1)
PDFMN  := $(addprefix $(MANDIR)/man1/,$(PDFsrc))

# Set compiler
CPLIST := g95 pgf90 /usr/bin/gfortran gfortran 
ifdef FC  #use predefined fortran compiler
  COMPILER := $(FC)
  $(info Using Fortran compiler : $(FC))
else       #find default compilers
  ifdef NERSC_HOST
    FC = ifort
    COMPILER = ifort
  else
    $(error Compiler NOT defined! Please set it with variable FC)
  endif
endif

# set up product name
EXEC  := $(BDIR)/surfgen-$(SGENVER)-$(OS)-$(ARC)
LIBF  := $(LDIR)/libsurfgen-$(SGENVER)-$(OS)-$(ARC).a
TSTX  := $(TDIR)/test

# set debugging flags
ifdef DEBUGGING_SYMBOLS
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
  ifdef NO_I8
    CPOPT = 
  else
    CPOPT =  -i8
  endif
  LKOPT = 
else
  DEBUGFLAG := 
  # set default compiler flags
  ifneq ($(findstring ifort,$(COMPILER)),)
    ifdef NO_I8
      CPOPT    := -auto -assume byterecl -parallel -O3 -lpthread -openmp
    else
      CPOPT    := -auto  -assume byterecl -parallel -O3 -lpthread -openmp -no-opt-matmul -i8
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

endif

# set BLAS and LAPACK libraries.  
# Please set variable $BLAS_LIB or $LIBS to enable the program to find the link line.
# On Mac the default is to use the vecLib framework
ifndef LIBS
 ifndef BLAS_LIB
  ifdef NERSC_HOST  
    # On NERSC.   Define proper MKL library pathes for hopper and carver
    ifeq ($(NERSC_HOST),hopper)
      MKLROOT=/opt/intel/composer_xe_2013.1.117/mkl
    else
    # i am assuming carver here.   
      MKLROOT=/usr/common/usg/intel/13.0.028/composer_xe_2013.1.117/mkl
    endif
    BLAS_LIB:=-Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
          $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
    LIBS:=$(BLAS_LIB)
    FFLAGS:=-openmp -i8 -I$(MKLROOT)/include
    LDFLAGS:=-openmp
  else
    ifeq ($(OS),Darwin)  
       #on mac, use frameworks
       LIBS := -framework vecLib
       $(info Using vecLib framework for Mac OS X)
       PS2PDF := /sw/bin/ps2pdf
    else
       $(info BLAS_LIB not set.  Trying to determine LAPACK link options...)
       #BLAS_LIB is not set.  check LD_LIBRARY_PATH for mkl
       ifneq ($(findstring mkl,$(LD_LIBRARY_PATH)),)
          $(info Found mkl in LD_LIBRARY_PATH. Using dynamic link to MKL.)
        LIBS := -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
       endif #ifneq mkl,LD_...
    endif #OS==Darwin
  endif #is on hopper.nersc.gov
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


# build everything
all  :  surfgen libs man tests
	@echo 'Finished building surfgen.'
	@echo ''

#target for standalone fitting and analysis program
surfgen  : $(OBJV)  $(OBJS) | $(BDIR)
	@echo ''
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
	$(CDS) $(COMPILER) -o $(EXEC) $(OBJS) $(OBJV) $(LIBS) $(LKOPT) $(LDFLAGS)
	@echo 'Finished building target: $@'
	@echo '-----------------------------------------'
	@echo 'Creating symbolic link to the new binary'
	ln -sf $(EXEC) $(BDIR)/surfgen
	@echo ''

#target for runtime interface library
libs  :  $(OBJV) $(OBJSL) | $(LDIR)
	@echo ''
	@echo '-----------------------------------------'
	@echo '  SURFGEN EVALUATION LIBRARY'
	@echo 'Program Version:     $(SGENVER)'
	@echo 'Archiver (AR):       $(AR)'
	@echo 'Library File Name:   $(LIBF)'
	@echo '-----------------------------------------'
	@echo 'Building target: $@'
	@echo 'Archiving the object into library '
	$(CDS) $(AR) -r -v  $(LIBF) $(OBJSL)  $(OBJV) 
	@echo '-----------------------------------------'
	@echo 'Creating symbolic link to new library'
	ln -sf $(LIBF) $(LDIR)/libsurfgen.a
	@echo ''

#
clean:
	-$(RM) $(OBJS) $(SDIR)/*.mod $(OBJSL) $(OBJT) $(OBJV) 
	@echo 'Finished cleaning'

.PHONY : clean

.FORCE :

# Compile and run test code 
tests : $(OBJV) $(OBJT) | $(TDIR) 
	@echo '-----------------------------------------'
	@echo '  SURFGEN TESTING PROGRAM'
	@echo 'Program Version:     $(SGENVER)'
	@echo 'Test Program:        $(TSTX)'
	@echo 'Execution Directory: $(TDIR)'
	@echo 'BLAS/LAPACK LIB:     $(LIBS)'
	@echo 'Debug Flag:          $(DEBUGFLAG)'
	@echo 'Linking Options:     $(LKOPT)'
	@echo '-----------------------------------------'
	@echo 'Building Target: $@'
	@echo 'Invoking Linkier'
	cd $(SDIR); $(COMPILER) -o $(TSTX) $(OBJT) $(OBJV) $(LIBS) $(LKOPT) $(LDFLAGS)
	@echo '-----------------------------------------'
	@echo 'Performing Tests.  Please see test.log for details'
	@echo ''
	@cd $(TDIR); $(TSTX) > $(TDIR)/test.log
	@echo ''
	@echo 'All tests finished. '
	@echo '-----------------------------------------'

# make sure version subroutines are always updated
$(OBJV) : .FORCE

.PHONY : man install

# compile man pages into pdf files, and copy man pages to system man page directory
install : $(PDFMN)
	@echo 'Copying man pages to /usr/share/man'
	-cp $(PDFMN) /usr/share/man/man1

$(DDIR)/%.pdf : $(MANDIR)/man1/%.1 | $(DDIR)
	@echo 'Constructing pdf manual from man page of $(notdir $(basename $<))'
	@man -t -M $(MANDIR) $(notdir $(basename $@)) > $(notdir $(basename $@)).tmp.ps
	@$(PS2PDF) $(notdir $(basename $@)).tmp.ps $@
	@rm $(notdir $(basename $@)).tmp.ps

man : $(PDFPG) 

# Compile source files
$(SDIR)/%.o : $(SDIR)/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Compiler'
	$(CDS) $(COMPILER) -c -o $@ $< $(CPOPT) $(DEBUGFLAG) $(FFLAGS)
	@echo 'Finished building: $<'
	@echo ' '

# Compile source file that needs preprocessing
$(SDIR)/%.o : $(SDIR)/%.F90
	@echo 'Building file: $<'
	@echo 'Invoking: Compiler'
	$(CDS) $(COMPILER) -c -o $@ $< $(CPOPT) $(DEBUGFLAG) $(FFLAGS) -DSGENVER=\"$(SGENVER)\"
	@echo 'Finished building: $<'
	@echo ' '
# Compile version string subroutine

$(BDIR) $(LDIR) $(TDIR) $(DDIR):
	@echo 'Creating directory $@'
	@mkdir -p $@
