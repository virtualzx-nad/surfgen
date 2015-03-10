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
SGENVER :=2.7.2

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
  $(info Using supplied Fortran compiler : $(FC))
else       #find default compilers
  ifdef NERSC_HOST
    $(info Running on NERSC)
    ifeq ($(NERSC_HOST),hopper)
      $(info System : hopper) 
      FC=ftn
      COMPILER=ftn
    else
      FC = ifort
      COMPILER = ifort
    endif
  else
    $(error Compiler NOT defined! Please set it with variable FC)
  endif
endif

# set up product name
EXEC  := $(BDIR)/surfgen-$(SGENVER)-$(OS)-$(ARC)
LIBF  := $(LDIR)/libsurfgen-$(SGENVER)-$(OS)-$(ARC).a
DYLIBF:= $(LDIR)/surfgen.dylib
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
    ifdef DYNAMIC
      CPOPT := $(CPOPT) -dynamic
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
      $(info Setting LAPACK flags with Hopper)
      MKLROOT=
      BLAS_LIB=
      LIBS=
      FC=ftn
      CPOPT=
      FFLAGS=
      LDFLAGS=
    else
      ifeq ($(NERSC_HOST),edison)
        BLAS_LIB:=-Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
            $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
        LIBS:=$(BLAS_LIB)
        FFLAGS:=-openmp -i8 -I$(MKLROOT)/include
        LDFLAGS:=-openmp
      else
    # i am assuming carver here.   
        MKLROOT=/usr/common/usg/intel/13.0.028/composer_xe_2013.1.117/mkl
        BLAS_LIB:=-Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
            $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
        LIBS:=$(BLAS_LIB)
        FFLAGS:=-openmp -i8 -I$(MKLROOT)/include
        LDFLAGS:=-openmp
      endif
    endif
  else
    ifeq ($(OS),Darwin)  
       #on mac, use frameworks
       LIBS := -framework vecLib
       $(info Using vecLib framework for Mac OS X)
       PS2PDF := /sw/bin/ps2pdf
    else
       $(info BLAS_LIB not set.  Trying to determine LAPACK link options...)
       #BLAS_LIB is not set.  check MKL_HOME and LD_LIBRARY_PATH for mkl
       ifdef MKL_HOME
         $(info Using path specified by MKL_HOME variable)
         MKLROOT=$(MKL_HOME)
         BLAS_LIB:=-Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
                 $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
         LIBS:=$(BLAS_LIB)
       else
         ifneq ($(findstring mkl,$(LD_LIBRARY_PATH)),)
            $(info Found mkl in LD_LIBRARY_PATH. Using dynamic link to MKL.)
            LIBS := -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
          endif #ifneq mkl,LD_...
       endif #ifdef MKL_HOME
    endif #OS==Darwin
  endif #is on hopper.nersc.gov
 else
    LIBS := $(BLAS_LIB)
 endif #BLAS_LIB
endif #ifndef $LIBS

    RM := rm -rf

ifndef AR
    AR := ar
endif

LTOOL=libtool -dynamic -install_name libsurfgen -current_version $(SGENVER)  -undefined suppress -flat_namespace

SURFGENLIB= $(LDIR)/libsurfgen.a $(LIBS) $(LKOPT) $(LDFLAGS)
# build everything
all  :  surfgen libs tests
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
	@echo '-----------------------------------------'
	@echo 'Compile and linking options:'
	@echo 'To compile a program that use surfgen subroutines,'
	@echo '* Compile with these options:' 
	@echo ''
	@echo '  $(CPOPT)'
	@echo '' 
	@echo '* Link link with these options :'
	@echo ''
	@echo '  $(SURFGENLIB)'
	@echo ''
	@echo '* You may also source the setsgenvar script in bin directory for your current shell.' 
	@echo '  It will set up compilation options as $$SGENFLAG and link options as $$SGENLIB.'
	@echo '  This shell script will also raise the stack limit.'
	@echo '' #sh/bash source file for variable settings 
	@echo '#!/bin/bash' > $(BDIR)/setsgenvars.sh 
	@echo 'ulimit -s unlimited' >> $(BDIR)/setsgenvars.sh 
	@echo '#Set surfgen compilation flags and link line options' >> $(BDIR)/setsgenvars.sh 
	@echo "export SGENFLAG='$(CPOPT)'" >> $(BDIR)/setsgenvars.sh
	@echo "export SGENLIB='$(SURFGENLIB)'" >> $(BDIR)/setsgenvars.sh
	@echo "export SGENFC='$(COMPILER)'" >> $(BDIR)/setsgenvars.sh
	@echo "export SGENDIR='$(BDIR)'" >> $(BDIR)/setsgenvars.sh
	@echo 'if [[ '\"':$$PATH:'\"' =~ '\"':$$SGENDIR:'\"' ]]; then' >> $(BDIR)/setsgenvars.sh
	@echo "  echo surfgen found in PATH variable">> $(BDIR)/setsgenvars.sh
	@echo "else">> $(BDIR)/setsgenvars.sh
	@echo '  export PATH=$$PATH:$$SGENDIR'>> $(BDIR)/setsgenvars.sh
	@echo "fi" >> $(BDIR)/setsgenvars.sh
	@echo '' #csh/tcsh source file for variable settings 
	@echo '#!/bin/tcsh' > $(BDIR)/setsgenvars.csh 
	@echo 'set stacksize unlimited' >> $(BDIR)/setsgenvars.csh 
	@echo '#Set surfgen compilation flags and link line options' >> $(BDIR)/setsgenvars.csh 
	@echo "setenv SGENFLAG '$(CPOPT)'" >> $(BDIR)/setsgenvars.csh
	@echo "setenv SGENLIB '$(SURFGENLIB)'" >> $(BDIR)/setsgenvars.csh
	@echo "setenv SGENFC '$(COMPILER)'" >> $(BDIR)/setsgenvars.csh
	@echo "setenv SGENDIR '$(BDIR)'" >> $(BDIR)/setsgenvars.csh
	@echo 'set found=`echo '\"'$$PATH'\"' | tr '\"':'\"' '\"'\n'\"' | grep -x '\"'$$SGENDIR'\"'`'>> $(BDIR)/setsgenvars.csh 
	@echo 'if ( $${?found} == 0 ) then'>> $(BDIR)/setsgenvars.csh 
	@echo '   setenv PATH '\"'$${PATH}:$${SGENDIR}'\" >> $(BDIR)/setsgenvars.csh
	@echo 'else ' >> $(BDIR)/setsgenvars.csh
	@echo '   echo surfgen found in PATH ' >> $(BDIR)/setsgenvars.csh
	@echo 'endif' >> $(BDIR)/setsgenvars.csh

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
	$(CDS) $(COMPILER) -c -o $@  $(CPOPT) $(DEBUGFLAG) $(FFLAGS) $<
	@echo 'Finished building: $<'
	@echo ' '

# Compile source file that needs preprocessing
$(SDIR)/%.o : $(SDIR)/%.F90
	@echo 'Building file: $<'
	@echo 'Invoking: Compiler'
	$(CDS) $(COMPILER) -c -o $@ $(CPOPT) $(DEBUGFLAG) $(FFLAGS) -DSGENVER=\"$(SGENVER)\" $<
	@echo 'Finished building: $<'
	@echo ' '
# Compile version string subroutine

$(BDIR) $(LDIR) $(TDIR) $(DDIR):
	@echo 'Creating directory $@'
	@mkdir -p $@

# Create Mac OS X dylib library
dylib : $(OBJV) $(OBJSL) | $(LDIR)
	@echo ''
	@echo '-----------------------------------------'
	@echo '  SURFGEN dylib LIBRARY for Max OS X'
	@echo 'Program Version:     $(SGENVER)'
	@echo 'Library Tool:        $(LTOOL)'
	@echo 'Library File Name:   $(DYLIBF)'
	@echo '-----------------------------------------'
	@echo 'Building target: $@'
	@echo 'Creating dynamic linked library'
	$(CDS) $(LTOOL) -o $(DYLIBF) $(OBJSL)  $(OBJV) 
	@echo ''

