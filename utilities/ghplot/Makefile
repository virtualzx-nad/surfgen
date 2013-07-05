OBJS  :=  ghplot.o 

COMPILER := ifort -i8

SURFGEN_LIB :=  /home2/xzhu/surfgen/lib/libsurfgen.a

all  :  $(OBJS)
	@echo '---------------------------------'
	@echo '   surfgen g-h plot utility'
	@echo '---------------------------------'
	@echo 'Building target: $@'
	@echo 'Invoking: Linker'
	$(COMPILER) -o ghplot.x $(OBJS) $(SURFGEN_LIB) $(LAPACK_LIB) -lpthread -auto
	@echo 'Finished building target: $@'
	@echo ' '
%.o : ./%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Compiler'
	$(COMPILER) -o"$@" "$<" -c -auto
	@echo 'Finished building: $<'
	@echo ' '

