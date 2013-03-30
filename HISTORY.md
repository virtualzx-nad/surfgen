## 2.1.4

* phenol input made runnable
* total number of coefs before null space removal is now printed to the surfgen.out file
* OS version dropped from executable name since they seem to be compatible
* Large gradient error warnings are now more reasonable and will no long complain about errors that 
are infinitesimal
* minor fixes to output formatting

## 2.1.3 (2013-03-29)

* New coordinate 4-center dot product added to describe angular motions that are anti-symmetry with respect to 
atom permutations but symmetry with respect to inversion.
* Changed out-of-plane coordinate(type=-1,mode=0) to take a second parameter which scales the coordinate linearly.
* Fixed genEnerGroup problem where states with no input data are being classified as degenerate.
* Test jobs for phenol and hydroxymethyl added.

## 2.1.2 (2013-03-27)

* Fixed `Makefile` to properly create `bin` and `lib` directories and prevent unnecessary building
* Fixed up messy output search code in `makesurf.f90`
* Fixed minor output formatting problems that cause warnings in newer versions of `ifort`

## 2.1.1 (2013-03-25)

* Updated makefile and project settings to allow compilation of evaluation libraries.
* Removed deprecated options disscprate and disscpmidpt.  Since we use OOP angles
with correct assyptotic behavior now, coupling scaling is no longer needed.
* Binaries are libraries are for the moment not supplied with the repository.  
You will have to build it yourself.

## 2.1.0 (2013-03-22)

* SearchPaths feature: the program to search a number of directories for
input files, and files that are not present will automatically be marked
as nonexistent and excluded from fit.  This eliminates the need to make 
up placeholder data and then set then to be not used in _points.in_ .
* Input file naming patterns: You can now set a naming pattern for your input
files, just in case they are named in a non-conventional manner.
* Removed deprecated subroutines makeXCoords and removeTransRot

## 2.0.0 (2013-03-21)

* First stable release of mix-global-local **surfgen** tested for NH3 case
