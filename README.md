surfgen
=======

surfgen use *ab initio* data to generate analytical coupled [potential energy surfaces(PES)]
(http://en.wikipedia.org/wiki/Potential_energy_surface) using quasi-diabatic Hamiltonians(**H**d).  
The primary purpose of this project is to provide a way to accelerate quantum mechanical 
simulations of non-adiabatic chemical processes without compromising the accuracy of 
expensive quantum chemical electronic structure methods required for non-adiabatic processes.

This project is coded with [Fortran 90](http://en.wikipedia.org/wiki/Fortran_90#Fortran_90) 
and contains both fitting programs used to generate coupled PESs and evaluation libraries 
that can be conveniently used to utilize the fit surface in any simulation program.

Features
--------

**Fitting Ab Initio Data**

The program is capable of fitting energy, energy gradient and derivative coupling (AKA nonadiabatic coupling or
vibronic coupling) data obtained from ab initio calculations.  A weighed least squares procedure is used to generate
the fit to simultaneously optimize the reproduction of adiabatic energies and energy gradients and the diabaticity of
the quasi-diabaticrepresentation, defined by the residual coupling between diabatic states.  Lagrange multipliers are
used to enable the exact reproduction of arbitrary selected set of data, such as energy and gradients at critical
points on the potential or energy difference and derivative couplings at points of conical intersections.

**Nonadiabatic couplings and Seams of Conical Intersections**

The **H**d approach is capable of extremely accurate description of nonadiabatic interactions.  It has been used to 
successfully describe large portions of the seam of conical intersections, spanning completely different 
geometrically structures.  This is enabled by the application of intersection adapted coordinates.

**Flexibility**
              
The program allows the user to define the blocks of the matrix **H**d with a large set of customizable basis functions,
enabling the flexibility to describe complex features on the surface to a high level of accuracy.

**Global Symmetry Treatment**

Projection operator method is used to allow the program to use an arbitrary subgroup of the Complete Nuclear 
Permutation Inversion group to construct symmetry adapted basis for the fitting procedure.  With the help of 
such feature, the program can correctly treat the symmetry in problems that involve large amplitude motions, 
as well as vibrational problems.

**Efficiency**

With the fully analytical form of Hd, the evaluation time for a single point is usually within 50ms. Future 
update is planned to allow vectorized evaluation of large number of data points.  The fitting program achieve 
high efficiency through the extensive use of optimized and threaded LAPACK and BLAS libraries.   Other 
functionalities such as automatic local internal coordinate construction, automatic null space removal and 
GDIIS extrapolations provide tools to enhance performance.

Documentation
-------------

Chapter 1 manual pages are provided for surfgen and several input files are provided in `man/man1`.
To view these manuals, add the man directory to `$MANPATH` or copy the manual files to `/usr/shared/man/man1`.

    $ cp man/man1/* /usr/shared/man/man1

Then simply use `man` to view the pages.  

    $ man surfgen

You can also use `-M` flag to explicitly specify the manual directory or use `-t` flag to generate PostScript
version of the manual page.  For example, use the following command to view PDF versions on Mac OS X

    $ man -M man -t surfgen | open -f -a /Applications/Preview

PDF versions of the documents can also be found in `pdf` directory.  However, they may not be as up-to-date
as the man files.

Installation
------------

**Please be advised that we have only tested installation on Mac OS X Mountain Lion, RHEL 4.x
and CentOS 6.x, with Intel Fortran compiler or gfortran.**  Other set up should also work but
might require modifications of code or Makefile.   A proper Fortran 90 compiler, GNU make 
and implementations of BLAS and LAPACK are also required.   They can all be obtained free of charge
([gfortran](http://gcc.gnu.org/wiki/GFortran), [GNU make](http://www.gnu.org/software/make/), 
[ATLAS](http://math-atlas.sourceforge.net))

A number of precompiled binaries fitting programs can be found in `bin` directory.

A library is built with `ifort` with static link to LAPACK/BLAS on CentOS 6.  It can be found in
`lib` directory.

If you cannot file a working binary or library, it will only take you a few minutes to build them
from source code, located in `source` directory.  See the following sections for compilation guides.

Compiling the Fitting Program
-----------------------------

If you cloned the repository in [Xcode](https://developer.apple.com/xcode/) on Mac OS X, 
you should be able to use Xcode to run or build for any type immediately after download 
if you have `gfortran` installed in `/usr/bin`.

To compile the fitting program, you have to set up environmental variables `$FC` to point to 
your fortran compiler and `$LIBS`(or `$BLAS_LIB`) to be your LAPACK/BLAS link line flags.  
If you included MKL path in `$LD_LIBRARY_PATH` or if you have included LAPACK/BLAS link options 
in `LDFLAGS`, then you can skip setting up the `$LIBS`(or `$BLAS_LIB`) variable.

After the variables are set, simply do

    $ make 

The compiled program can be found in `bin` directory with the name 
`surfgen-{version}-{OS}-{OS version}-{compiler}`

You can also pass any unset variables through arguments of `make`

    $ make FC=ifort BLAS_LIB="-Wl,--start-group  $MKLROOT/lib/intel64/libmkl_intel_ilp64.a $MKLROOT/lib/intel64/libmkl_intel_thread.a $MKLROOT/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm"

If for any reason you suspect a bug, you can include debug symbols or check flags by setting `$DEBUGFLAG` variable 
or set `DEBUGGING_SYMBOLS` to `YES`.  Note that if you click the Run button in Xcode or compile for Debugging, this 
option is automatically enabled.

Extra compilation flags can be added through `$CPOPT` or `$FFLAGS`.  Extra linking flags can be added through
`$LKOPT` or `$LDFLAGS`.

Compiling the Evaluation Library
--------------------------------

To compile the library, you also need to supply environmental variable `$FC`, 
either in the shell or through arguments.

Example:

    $ make lib FC=gfortran

Extra compilation flags can be added through `$CPOPT` or `$FFLAGS`.  Extra linking flags can be added through
`$LKOPT` or `$LDFLAGS`.  You can also change the default archiver by setting `$AR`.  

For example, if you want to enable ipo:

    $ make lib FC=ifort FFLAGS="-parallel -O3 -xHost -i8 -ip -ipo" AR="xiar -rv"

Library can be found in `lib` directory.

Authors
-------

surfgen is created by [Xiaolei Zhu](http://www.linkedin.com/in/virtualzx) and 
[other awesome members](http://www.jhu.edu/~chem/yarkony/group.html)
from [Yarkony Group](http://www.jhu.edu/~chem/yarkony/), Department of Chemistry, Johns Hopkins University

This program is based on, and would not have been possible without the quadratic Hd fitting program 
for vibrational spectroscopy simulation, developed by 
[Michael Schuurman](ca.linkedin.com/pub/michael-schuurman/7/996/32), then also a member of Yarkony group.

References
----------
  
Xiaolei Zhu and David Yarkony, *Toward eliminating the electronic structure bottleneck in* 
*nonadiabatic dynamics on the fly: An algorithm to fit nonlocal, quasidiabatic, coupled* 
*electronic state Hamiltonians based on ab initio electronic structure data* ,
[*J. Chem. Phys.*, **132**, 104101 (2010)](http://dx.doi.org/10.1063/1.3324982)

Copyright
---------

Copyright (c) 2013, Yarkony Group, Johns Hopkins University

Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is 
hereby granted, provided that the above copyright notice and this permission notice appear in all copies
and that any documentation, advertising materials, and other materials related to such distribution and 
use are provided with proper reference to 
[surfgen publications](https://github.com/virtualzx-nad/surfgen/blob/master/README.md#references).

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE 
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE 
FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT 
OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.