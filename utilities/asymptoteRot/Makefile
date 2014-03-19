# Can use gcc or icc, but we chose gcc
CC		:= gcc
CFLAGS 	:= -lm -O2
DEBUG 	:= -DDEBUG

EXE		:= asymptoteRot.x

%.o: %.c 
	sh -c "$(CC) -c -o $@ $< $(CFLAGS)"

$(EXE): asymptoteRot.o
	sh -c "$(CC) -o $@ $< $(CFLAGS)"

install:
	@echo "Installing..."
	sh -c "mv asymptoteRot.x ../../bin/ "

clean:
	@echo "Cleaning..."
	rm -f ./*.o

