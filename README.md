surfgen
=======

surfgen use *ab initio* data to generate analytical coupled [potential energy surfaces(PES)]
(http://en.wikipedia.org/wiki/Potential_energy_surface) the quasi-diabatic Hamiltonian(**H**<sup>d</sup>) approach.  
The primary purpose of this project is to provide a way to accelerate dynamics
simulations of non-adiabatic processes without compromising the accuracy of the representation of 
electronic states, which requires extremely expensive correlated multi-reference electronic structure methods 
such as MRCI.

The form of expansion and symmetry treatment is general and versatile, making it essentially capable of
fitting almost any system as long as the underlying *ab initio* data can be supplied.  The method has been used
to generate coupled PESs for systems as large as phenol ([*J. Chem. Phys.* **140**, 024112 (2014)](http://dx.doi.org/10.1063/1.4857335)), where the coupled PESs accurately describe the full 33 
dimensionalities in their entire dynamically relevant region, to a very high energy range of 50,000cm<sup>-1</sup>,
including the 4 lowest singlet states.  

The program uses energy, energy gradient and derivative coupling data from *ab initio* data to generate the fit. 
The use of energy gradient drastically reduce the number of data points needed to construct a converged fit, 
especially in higher dimensional cases, and at the same time significantly reduces oscilation in the fit surfaces.
The program use derivative coupling data to automatically generate the most diabatic representation for the system,
in a least-squares sense.   Gradients and couplings are not required for all data points, but are crucial
for the regions with strong couplings. We have found that obtaining high level *ab initio* data that contains gradient
and couplings data and consistent in the entire domain is the most difficult step in the procedure.

Molecular symmetry, or subgroups of Complete Permutation Inversion Symmetry (CNPI) group is fully implemented 
in a general manner through a projection operators approach.  CNPI group has an almost infinite number of 
group structures and irreducible representations, and it is in many aspect very different from point group 
symmetry. The user is advised to carefully study the concept of CNPI symmetry and the problem at hand to 
correctly choose the right symmetry for your system.  The program can in principle treat any symmetry, but the
user will have to supply the representation matrices.

Fortran 90 and 70 interfaces are supplied to facilitate the use of fit PESs in other programs.  

This project is coded with [Fortran 90](http://en.wikipedia.org/wiki/Fortran_90#Fortran_90) 
and contains both fitting programs used to generate coupled PESs and evaluation libraries 
that can be conveniently used to utilize the fit surface in any simulation program.

Features
--------

**Fitting Ab Initio Data**

The program is capable of fitting energy, energy gradient and derivative coupling (AKA nonadiabatic coupling or
vibronic coupling) data obtained from ab initio calculations.  A weighed least squares procedure is used to generate
the fit to simultaneously optimize the reproduction of adiabatic energies and energy gradients and the diabaticity of
the quasi-diabaticrepresentation, defined by the residual coupling between diabatic states.  Lagrange multipliers 
can be used to enable the exact reproduction of arbitrary selected set of data, such as energy and gradients at critical
points on the potential or energy difference and derivative couplings at points of conical intersections.

Geometries, energies, gradients and couplings should all be prepared in COLUMBUS format.  The program does not require 
all data to be present at all data points, and can operate even when large amount of data are missing.  However, 
missing data makes the fitting procedure less efficient and less reliable.  The availability of energy gradients and 
couplings will drastically reduce the number of data points need to construct a fit that describe the entire desired
region.   Gradients and couplings are also crucial in areas where states are strongly coupled.

The algorithm, through the use of gradients and couplings, is highly tolerant of discontinuities and numeric noises in 
the underlying *ab initio* data, and is found to be able to automatically smooth out these unphysical effects. 
However, if the discontinuities are too large the program can have trouble identifying and matching the electronic states.

**Nonadiabatic couplings and Seams of Conical Intersections**

The **H**<sup>d</sup> approach is capable of extremely accurate description of nonadiabatic interactions.  It has been used to 
successfully describe large portions of the seam of conical intersections, spanning completely different 
geometrically structures.  This is enabled by the application of intersection adapted representation.  The adiabatic
representation changes drastically near avoided crossings, reaching singularity at the seam of conical intersections, 
which complicates the fitting procedure.   In order to facilitate a stable fit, we use gradients and couplings to 
construct a stable intersection adapted representation (which is effectively the joint diagonalization of the matrix 
representation of vector operator ∇H) that is continuous near the seam.  *Ab initio* data and the fit Hamiltonian are both
rotated to this representation to facilitate the fit.  

Arbitrary energy difference range and arbitrary number of states can be treated, provided that all gradient and coupling
data are present. (see [*J. Chem. Phys.* **141**, 174109 (2014)](http://dx.doi.org/10.1063/1.4900631))

**Flexibility**
              
The program allows the user to define the blocks of the matrix **H**<sup>d</sup> with a large set of customizable 
basis functions, enabling the flexibility to describe complex features on the surface to a high level of accuracy.
The function basis are constructed from symmetric projections of the monomials of single-coordinate functions.  
A number of forms of single coordinate function forms are provided, including internuclear distances, bond angles 
and different types of out-of-plane bending motions.  Many scaling options are availble for each type of coordiantes.
The user can control the number, parameters and types of the coordinates used as well as their combination to generate
the monomials through `coords.in` file.

This in principle allows the user to achieve any precision.  However, in practice the memory use scales with _n_<sup>2</sup>,
the fitting time scales with _n_<sup>3</sup> and the evaluation time when using the fit H<sup>d</sup> scales with _n_,
where _n_ is the size of funciton basis used.  This can quickly render the fit impractical if the basis is too large.  
A basis size below 10,000 is prefered, but a size around 20,000 is usually still managable provided that sufficient memory is
available.  We recommend fine tuning the parameters and type of coordinates instead of simply increase the order of monomial.
Adding different types and scaling of coordinates and coordinates between different group of interacting atoms is also found 
to be much more efficient than increasing the order.  It is recommended that Morse functions for most interatomic distances 
be included, as well as sufficient out-of-plane coordinates, to ensure that the coordinates is capable of describing all motions. Hyperbolic tangent, gaussian functions and cosines of bond angles are found to be very efficient in improving 
quality of the fit when many different regions are involved, since they behave like smooth window functions.

**Global Symmetry Treatment**

Standard group-theoretic projection operator method is used by the program to facilitate an arbitrary subgroup of the 
Complete Nuclear Permutation Inversion(CNPI) group to construct symmetry adapted basis for the fitting procedure.  
With the help of such feature, the program can correctly treat the symmetry in problems that involve large amplitude motions, 
as well as vibrational problems.   

CNPI group is a global symmetry group, which provides symmetry relation between different parts of the surfaces.  This is
to be compared to point group symmetry, which relates different vibrational component of the surface at the same data point.
In fact, with few exceptions, CNPI group symmetry automatically generates the corresponding point group symmetry when the
geometry carries any specific point group.  The point group is homomorphic to a subgroup of the CNPI group induced by 
the subset of CNPI operations and/or rotations that keeps the geometry invariant.  The only exception being the C<sub>∞</sub>
axis of linear molecules, which is a pure rotation that cannot be achieved by any permutation/inversion.

Relative signs between different parts of the potential can be very difficult to access, making the assignment of CNPI 
symmetry difficult.  One should therefore take advantage of the induced subgroup at high symmetry points to determine 
the symmetry of diabatic states.  

The projection operator approach requires representation matrices.  Due to the large amount of possible combination of 
groups and irreducible representations, we currently do not supply these representation matrices.  The user needs to
work out the effect of symmetry operations and representation matrices and supply them in the `irrep.in` file.

A problem that very often come up is the change of symmetry of states in the process of a reaction.  Since only a
limited number of states are treated, states with different symmetry are often found to enter and exit the set of 
treated states in the course of a reaction, changing the symmetry.  This means whatever symmetry that we set the 
diabatic states to be, they are likely to be incorrect in some regions where states with other symmetries intrude
into the problem.  One can still converge a fit in this case, but the couplings and the topography of conical intersections
will be incorrect.  Usually this can be overcome by including 1 or 2 more states than the problem physically require, 
and allowing these states to have incorrect symmetry, so that the problematic crossings happen above the energy limit of 
the problem.  

**Efficiency**

With the fully analytical form of H<sup>d</sup>, the evaluation time for a single point is usually within 50ms. Future 
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

* Stacksize

Please set stacksize to unlimited to allow execution for larger systems.

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

Help, Problem Reports and Wiki
-------
Please email [Xiaolei Zhu](mailto:virtualzx@gmail.com) if you need help to compile or use the program. Please use the
Issues page on github to report bugs and other problems.  

We are planning on implementing an wiki site but have yet to realize that thought.  Anyone is invited to help to 
construct the wiki.

Authors
-------

surfgen is created by [Xiaolei Zhu](http://www.linkedin.com/in/virtualzx) and 
[other awesome members](http://www.jhu.edu/~chem/yarkony/group.html)
from [Yarkony Group](http://www.jhu.edu/~chem/yarkony/), Department of Chemistry, Johns Hopkins University

This program is based on, and would not have been possible without the quadratic Hd fitting program 
for vibrational spectroscopy simulation, developed by 
[Michael Schuurman](ca.linkedin.com/pub/michael-schuurman/7/996/32), then also a member of Yarkony group.

Xiaolei Zhu's Research Profile:
<a title="Follow me on ResearchGate" href="https://www.researchgate.net/profile/Xiaolei_Zhu4/?cp=shp"><img src="https://www.researchgate.net/images/public/profile_share_badge.png" alt="Follow me on ResearchGate" /></a>
<a title="My ORCID profile" href="http://orcid.org/0000-0002-7825-4815">
<img src="http://orcid.org/sites/all/themes/orcid/img/orcid-logo.png" alt="My ORCID profile" border="5"/>
</a>


References
----------

Xiaolei Zhu's Ph.D. dissertation, [The Quasi-Diabatic Hamiltonian Approach to Accurate and Efficient Nonadiabatic Dynamics with Correct Treatment of Conical Intersection Seams](https://jscholarship.library.jhu.edu/handle/1774.2/37011)(2014)

Xiaolei Zhu and David Yarkony, *Quasi-diabatic representations of adiabatic potential energy surfaces coupled by conical intersections including bond breaking: A more general construction procedure and an analysis of the diabatic representation* ,
[*J. Chem. Phys.*, **137**, 22A511 (2012)](http://dx.doi.org/10.1063/1.4734315)
  
Xiaolei Zhu and David Yarkony, *Toward eliminating the electronic structure bottleneck in* 
*nonadiabatic dynamics on the fly: An algorithm to fit nonlocal, quasidiabatic, coupled* 
*electronic state Hamiltonians based on ab initio electronic structure data* ,
[*J. Chem. Phys.*, **132**, 104101 (2010)](http://dx.doi.org/10.1063/1.3324982)

Copyright Notice
---------
Copyright 2011-2015 Yarkony Group, The Johns Hopkins University.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


<a rel="license" href="http://creativecommons.org/licenses/by/3.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by/3.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">The texts in Surfgen Program Suite</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="http://www.jhu.edu/~chem/yarkony/" property="cc:attributionName" rel="cc:attributionURL">Yarkony Group, Johns Hopkins University</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/3.0/deed.en_US">Creative Commons Attribution 3.0 Unported License</a>.
