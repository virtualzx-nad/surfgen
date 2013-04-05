## 2.2.2 (2013-04-05)

* The print out for errors in LSE block before every iteration were not weighted and therefore does not reflect the true quality of
fit.  Now both weighted and unweighted results are reported.   

## 2.2.1 (2013-04-03)

* Direct construction of normal equations enables the program to skip explicit construction and storage of W matrix, which becomes
memory intensive when number of equations grow close to 100,000.  
* MINMEX mode is deprecated.  Use the evaluation library with another program for the purpose.

## 2.2.0 (2013-04-02)  

This is a major update.  Part of your input will have to be redone in order to work properly.
* You can now specify a list of allowed symmetries in stead of one specific symmetry for any state group.  This allows the treatment 
of a changing symmetry.  All specified symmetry has to share the same dimensionalities, otherwise the result will be unpredictable.
* Options `groupsym` and `groupprty` are now two dimensional arrays, with the first dimension being index of state group, second 
being the index of allowed symmetry.   
* Merged matrix basis forward and backward transformation matrix for basis reconstruction to improve memory efficiency further.
* StepMethod/=0 is now depreciated since it does not really work.

## 2.1.7 (2013-03-31)

* For points that has large coupling errors, the program will now also print out the percentage error for derivative coupling
times energy difference.  This is particularly useful for diagnosis purpose in order to find out if it is gradients of Hd or
energy difference that is causing the large gradient error.   For points with small energy separations, large relative error
is sometimes unavoidable unless you force it with Lagrange multipliers using the *LD* option in `points.in`
* Fixed a bug that sometimes corrupt the permutational sign of OOP angles and 4C dotproducts, resulting in problematic symmetry properties.

## 2.1.6 (2013-03-31)

* New input option `restartdir`:  the program will save Hd coefficients to `$restartdir/hd.data.$iter` every 
iteration if this option is nonempty
* Fixed a problem where the point indices are not correctly printed in `surfgen.out` when some ab initio are absent 

## 2.1.5 (2013-03-30)

* Fixed a problem where the construction of intersection adapted coordinate at points with partially avaiable 
data will cause all ab initio data to be corrupted.
* Fixed an error in the print out where the state label for the table of equation inclusion does not line up
with the content 
* Increased file index searching limit in getFLUnit to avoid file unit ID overflow.

## 2.1.4 (2013-03-29)

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
