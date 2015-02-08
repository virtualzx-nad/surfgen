 
 LIB=memmac.o cgetmem.o cdatime.o


 memlib.a: $(LIB)
	-ar rvu memlib.a $(LIB)
	ranlib memlib.a
		

 
 memmac.o: memmac.f
	xlf -c -O memmac.f

 cgetmem.o: cgetmem.c
	cc -c cgetmem.c

 cdatime.o: cdatime.c
	cc -c cdatime.c
