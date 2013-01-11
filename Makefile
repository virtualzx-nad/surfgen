
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
  DEBUGFLAG := -g -fbounds-check -fbacktrace
else
  DEBUGFLAG :=
endif

# setting up the parameters to work from both local (mac) and on the cluster
# on the cluster, for both running from node056 (64 bit node with all libraries)
# from other nodes, log into node056 and compile
ifneq ($(findstring Mac,$(HNAME)),)
  Msg:=Building on local Mac
  PREFIX := sh -c
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
  ifeq ($(HNAME),node056)
    PREFIX := sh -c 
    COMPILER := ifort
    AR := xiar -rv
  else
    PREFIX := ssh node056 
    COMPILER := cd $(JDIR) ; ifort
    AR := cd $(JDIR) ; xiar -rv
  endif
  # location of LAPACK and BLAS library for 64 bit
  LIBS     := -Wl,--start-group  $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a \
            $(MKLROOT)/lib/em64t/libmkl_intel_thread.a  \
            $(MKLROOT)/lib/em64t/libmkl_core.a -Wl,--end-group -lpthread -lm 
  #compile options
  CPOPT    := -auto -c -assume byterecl -parallel -O3 -lpthread -openmp -axS
  #link options
  LKOPT    := -auto -lpthread -parallel
  #executable naming
  EXEC     :=surfgen.x  
endif

RM := rm -rf

#target for standalone fitting and analysis program
all  :  $(OBJS)    
	@echo 'MESSAGE : $(Msg), Debug flag=$(DEBUGFLAG)'
	@echo 'Building target: $@'
	@echo 'Invoking: Linker'
	$(PREFIX) "$(COMPILER) -o $(EXEC) $(OBJS) $(LIBS) $(LKOPT)"
	@echo 'Finished building target: $@'
	@echo ' '

#target for runtime interface library
lib  :  $(OBJSL)
	@echo 'Building target: $@'
	@echo 'Archiving the object into library libsurfgen.a'
	$(PREFIX) "$(AR) libsurfgen.a $(OBJSL)"
	@echo ' '

#
clean:
	-$(RM) $(OBJS) *.mod
	-@echo 'Finished cleaning'

#
%.o : ./%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Compiler'
	$(PREFIX) "$(COMPILER) -o $@ $< $(CPOPT) $(DEBUGFLAG)"
	@echo 'Finished building: $<'
	@echo ' '



