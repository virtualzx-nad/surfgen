 polyhes.vnew.e: polyhes.vnew.o 
	xlf -g -o polyhes.vnew.e polyhes.vnew.o memlib.a  iolib.a 

 polyhes.vnew.o: polyhes.vnew.f
	xlf -c -g  -w polyhes.vnew.f


