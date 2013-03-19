
# Objects needed for standalone fitting program
OBJS     := hddata.o diis.o rdleclse.o combinatorial.o progdata.o libutil.o \
            libsym.o libinternal.o localcoord.o makesurf.o minmex.o io.o surfgen.o 

# Objects needed for runtime interface library
OBJSL    := hddata.o combinatorial.o progdata.o libutil.o libsym.o libinternal.o\
            localcoord.o io.o potlib.o minmex.o

HNAME := $(shell hostname)
JDIR  := $(shell pwd)

# set debugging flags
ifeq ($(DEBUGGING_SYMBOLS),YES)
  DEBUGFLAG := -g -fbounds-check -fbacktrace -Wall -Wextra
else
  DEBUGFLAG := 
endif

# setting up the parameters to work from both local (mac) and on the cluster
# on the cluster, ifort needs to be in the seek path and MKL libraries 
# needs to be properly defiled with variable $BLAS_LIB 
ifneq ($(findstring Mac,$(HNAME)),)
  Msg:=Building on local Mac
  COMPILER := gfortran
  AR := ar -rv
  # location of LAPACK and BLAS library for Mac
  LIBS     := -framework vecLib
  #compile options.  
  CPOPT    := -c -fopenmp
  #link options
  LKOPT    :=     
  #executable naming
  EXEC     :=surfgen
else
  Msg:=Building on Linux
  COMPILER := ifort
  AR := xiar -rv
  # location of LAPACK and BLAS library for 64 bit
  LIBS     := $(BLAS_LIB) 
  #compile options
  CPOPT    := -auto -c -assume byterecl -parallel -O3 -lpthread -openmp -xHost -no-opt-matmul -i8 
  #link options
  LKOPT    := -auto -lpthread -parallel
  #executable naming
  EXEC     :=surfgen  
endif

RM := rm -rf

#target for standalone fitting and analysis program
all  :  $(OBJS)    
	@echo 'MESSAGE : $(Msg), Debug flag=$(DEBUGFLAG)'
	@echo 'BLAS/LAPACK LIB:  $(LIBS)'
	@echo 'Compiler options: $(CPOPT)'
	@echo 'Linking options:  $(LKOPT)'
	@echo 'Building target:  $@'
	@echo 'Invoking: Linker'
	$(COMPILER) -o $(prefix)$(EXEC) $(OBJS) $(LIBS) $(LKOPT)
	@echo 'Finished building target: $@'
	@echo ' '

#target for runtime interface library
lib  :  $(OBJSL)
	@echo 'Building target: $@'
	@echo 'Archiving the object into library libsurfgen.a'
	$(AR) $(prefix)libsurfgen.a $(OBJSL)
	@echo ' '

#
clean:
	-$(RM) $(OBJS) *.mod
	-@echo 'Finished cleaning'

#
%.o : ./%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Compiler'
	$(COMPILER) -o $@ $< $(CPOPT) $(DEBUGFLAG)
	@echo 'Finished building: $<'
	@echo ' '



