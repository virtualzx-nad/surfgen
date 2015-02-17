# Makefile for pointDist.x
OBJS	:= pointDist.o
CC	:= icc
COPT	:= 
# Executable
EXE	:= pointDist.x

%.o: %.c
	sh -c "$(CC) -c -o $(COPT) $@ $< "

$(EXE):	$(OBJS)
	sh -c "$(CC) -o $@ $< "

install: 
	@echo "Installing..."
	sh -c "mv pointDist.x ../../../bin/ "

clean:
	@echo "Cleaning..."
	rm -f ./*.o