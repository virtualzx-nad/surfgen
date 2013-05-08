## 2.4.6 (2013-05-06)

* A new option CpOrder now allows the user to specify the maximum order of off-diagonal coupling blocks to be different from that of the diagonal blocks.
* Added g and h vector output for intersection points.

## 2.4.5 (2013-05-02)

* Added FORTRAN77 interface to `libsurfgen` which allows arrays with arbitrary fixed sizes to be passed as arguments to 
evaluation subroutine.
* Fixed a problem with matrix printing that caused the gradients and couplings output in `surfgen.out` not to align correctly.

## 2.4.4 (2013-05-01)

* Added utility program `testpt`.   This, and any future, utility programs can be found in the `utilities` directory.
* Removed some obselete features from `potlib`, such as coupling flattening, which can now be properly prevented by correct choice of coordinates.

## 2.4.3 (2013-04-30)

* Fixed the `enfDiab` not working bug.  

## 2.4.2 (2013-04-30)

**IMPORTANT** : You have to updated your `coord.in` file after updating to this version!

* _All_ coordinate types now require the two scaling coefficients.  Put them there regardless of coordinate you are defining.   
With this change, the format is the same for all coordinate definitions and you do not have to check if coefficients are expected for a specific type.
* Plain bond distance and bond angle coordinates are not shifted by coefficient 2 and then scaled by coefficient 1.   Old definitions would be equivalent
to defining with

    1.000     0.000

* Added new rij scaling type (tanh) that can be used as a smooth step function.

## 2.4.1 (2013-04-27)

* Fixed an array deallocation problem upon termination.

## 2.4.0 (2013-04-27)

* Block symmetry analysis: the program now try to detect blocks that exhibit the same symmetry properties.   
Irrep matrices for *blocks* between all possible pairs of state irreps are constructed and compared.   
Expansions(`maptab`) and basis generation (AKA null space removal, performed in subroutine `genBasis`) are only performed for blocks with unique symmetries.  
## 2.3.7 (2013-04-26)

* Fixed a bug that off-diagonal are not being properly added during the normal equation construction when partial diagonalization is performed. 

## 2.3.6.1 (2013-04-25)

* Updated documentation for the changes in internal coordinate definition for type=0 mode=4

## 2.3.6 (2013-04-17)

* Changed form of long range distance scalings (mode=4) to make it more useful.
* Fixed redundant gradient/couplings error output.

## 2.3.5 (2013-04-15)

* Cleaning up fitting output

## 2.3.4 (2013-04-14)

* Cleaned up `testsurfgen`.
* Updated project settings for XCode

## 2.3.3 (2013-04-14)

* Gradient following method implemented.  A diagonal approximation is used for the coef-coef block of Lagrangian hessian.
This prevents the construction and factorization of normal equations, which is the most memory and CPU intensive step.
* Linear search algorithm reworked.  The program now use a more robust algorithm to sample along the step direction and 
use linear gradient extrapolation to find better step length.  Each pair of adjacent points are chosen and the pair that
give the best predicted improvements are used. 
* Fixed a flag issue which made version subroutine fail to compile with `gfortran`

## 2.3.2 (2013-04-12)

* Improved linear search algorithm
* LinNegSteps is now deprecated

## 2.3.1 (2013-04-12)

* Performance improvement for Lagrangian gradient and normal equations constructions. 
* Examples updated for the newest version.

## 2.3.0 (2013-04-11)

* The energies are structured as matrices with off-diagonal elements are explicitly stored for situations where intersection 
adapted coordinates instead of Schrodinger equations are used to determine some of the states.  In these cases, the off-diagonal
will be a small but non-zero number.   The treatment is now exact even when intersection adapted coordinates are used for 
points that are not strictly degenerate (which will always be the case)
* The evaluation of Dij now use the same degeneracy groups that are used to construct intersectin adapted coordinates
* The Dij are solved with a set of linear equations which should exactly reflect the derivatives in an arbitrary groups of 
degeneracies among arbitrary number of states.  

## 2.2.10 (2013-04-10)

* Symbolic link `libsurfgen` to most recent library created in `lib` directory when compiling libraries.
* Parameter `DijScale2` is now deprecated.  Its purpose is now achieved by `DijScale`.  You can still suppress the Dij 
contributions to normal equations, but you can no longer disable Dij contributions to Lagrangians.
* pdf documentations are no longer distributed with the repository.  Instead, they can be generated from man pages with
command `$ make man`.  The man pages can be installed into `/usr/share/man` with command `$ sudo make install`. 

## 2.2.9 (2013-04-10)

* Added version numbers to `surfgen.out`, standard output, and `hd.data` files
* The eigenvector input file is now processed on a point by point basis.  For undefined points, `updateEigenVec` will be
used to generate the initial vectors.  This allows Ckl data generated from a smaller point set to be used.  It can also
be used to manually modify the adiabatic-diabatic transformation at a specific point.
* When compiling with `make surfgen`, a symbolic link `surfgen` will be created in `lib` directory, pointing to the 
newest executable.

## 2.2.8 (2013-04-10)

* Fixed the input file searching problem when printlvl<2 (Issue #8)
* Code clean up

## 2.2.7 (2013-04-10)

* Fixed the DIJ problem in normal equations

## 2.2.6 (2013-04-09)

* Fixed a bug in Lagrangian gradient evaluation where some elements of Dij are used before constructed, causing the program to
fail the Lagrangian test.  (Issue #7)

## 2.2.5 (2013-04-08)

* Test program for the gradients of Hd and the Largrangian added to the test program
* Increased the length of all filename strings from 72 to 255
* Code clean up and structural optimizations

## 2.2.4 (2013-04-08)

* Test program can now be compiled and run with `make` by doing

    $ make tests

* Another fix to error print out

## 2.2.3 (2013-04-07)

* A test program created to test the consistency and detect problems within subroutines
* Fixed a bug that ncoords was not initialized to 0, which caused 4-center dot product
 definition to behave erroneously with `gfortran` compiler.

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
