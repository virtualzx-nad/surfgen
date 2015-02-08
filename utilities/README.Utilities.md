Surfgen Utility Programs
========================

Descriptions
------------

In this directory you can find a number of small programs that use the surfgen library to do a few simple tasks.  
These are meant to serve as simple example programs that illustrates how to incorporate coupled potential energy
surfaces created by surfgen into your program.  They can also perform simple tasks for analysis purposes.

Structure
---------

Each sub-directory contains one program.  This includes the source code, a README file that describes the program,
and an installation script.  Below is a list of all utility programs and their usage

|  Program Name     |  Description                            
| ----------------- | --------------------------------------- 
|  testpoints       | Print out energies in adiabatic and diabatic representation at given geometries.
		        |  testpoints.x also prints orthogonal G and H vectors for ab_initio data. 
|  ghplot           | Output energy data on a grid along the branching coordinates from a given point  
|  findcp           | Critical point (minima or saddle points) search with frequency analysis         

Installation
------------

To install the utility programs, first compile the surfgen evaluation libraries by executing `make` or `make libs`
in surfgen main directory.   After successful compilation, go to the utility program that you need and execute the
installation script `install.sh` from the directory of that program.  The compiled executable bears the same name as the directory name.

Usage
-----

Please consult README files in individual programs.
