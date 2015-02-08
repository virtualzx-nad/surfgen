#  fortran1 = rev3 
 LIB=iolib1.o iolib2.o pack77.o


iolib.a: $(LIB)
	-ar rvu iolib.a $(LIB)
	ranlib iolib.a
		

 
 iolib1.o: iolib1.f
	fortran -c -RD -Ogv -AS iolib1.f

 iolib2.o: iolib2.f
	fortran -c -RD iolib2.f


 pack77.o: pack77.f
	fortran -c -RD -Ogv -AS pack77.f
