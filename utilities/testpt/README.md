testpt
======

testpt is a surfgen utility used to extract surfgen predicted data for specific geometries.

2013 Yarkony Group, Department of Chemistry, Johns Hopkins University

Synopsis
--------

testpt [geom_file_name] 

Arguments
---------

* geom_file_name

Name of the file that contains all the geometries to be evaluated.   The geometries are given in cartesian coordinates,
in atomic units, with _COLUMBUS_ format.  The default is `geom.all`.  If more than more geometires are present, simply 
append all geometry together.

* hd_file_name

Name of the file that contains the expansion coefficients of the coupled potential surfaces.  The default is `hd.data`,
and the filename can be changed in `surfgen.in`.

Input Files
-----------

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

The utility use surfgen evaluation library `libsurfgen.a` to extract Hd data in both diabatic and adiabatic representations.
Evaluations are done for geometries specified in file defiled by variable `geom_file_name`.  


