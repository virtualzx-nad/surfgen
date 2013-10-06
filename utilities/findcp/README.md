findcp
======

`findcp` is a surfgen utility used to search for critical points on fit potential energy surfaces 

2013 Yarkony Group, Department of Chemistry, Johns Hopkins University

Synopsis
--------

findcp.x [geom_file_name [surface_index]]

Arguments
---------

* geom_file_name

Name of the file that contains the initial geometry for critical point search.   The geometries are given in cartesian coordinates,
in atomic units, with _COLUMBUS_ format.  The default is `geom`. 

* surface_index

The surface on which critical point search will be performed.

Input Files
-----------
* hd.data

File that contains the expansion coefficients of the coupled potential surfaces.  The default is `hd.data`,
and the filename can be changed in `surfgen.in`.

* surfgen.in

Standard surfgen input file that has the `GENERAL` namelist that defines the molecule and expansion.  

* coord.in

File that defines the coordinates used to construct Hd.

* connect.in [optional]

This file defines the linking conditions that generates the  feasible permutations of the CNPI group. The name of this
file can be changed in `surfgen.in`

* irrep.in

This file defiles the irreducible representation matrices of the CNPI group carried by the quasi-diabatic states used to 
construct Hd.

Method
------

The program will first display geometry information, including bond-lengths, angles and torsion angles for connected atoms,
then standard Newton-Raphson search is performed on designated surface, with Hessians calculated every iteration.  
The utility use surfgen evaluation library `libsurfgen.a` to extract Hd data in both diabatic and adiabatic representations.

Currently certain parameters are hardwared in the code:   Maximum number of iterations is 100, maximum allowed step size is
0.1, convergence tolerance for gradients and displacements are 10^-7 and 10^-5 respectively.   

Installation
-----------

This utility program uses the subroutines found in the surfgen evaluation libraries.  A static library is generated when
installing the main program through `make` or `make libs`.  When the library is properly compiled and archived, use
the `install.sh` script to compile the utility program.
