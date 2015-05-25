## 2.8.3 (2015-05-25)

* Improved findmex algorithms
* On can use hd coefficients to weight the basis selection procedure with scalebycoef option

## 2.8.3 (2015-05-20)

* Levenberg-Marquardt implementation where both diagonal and Jacobian shifts are manually set.  Diagonal shift through
LSETol and Jacobian shift through jshift parameter. 
* Set `autoshrink` to `.true.` to force the program to shrink step sizes when fitting error raises.

## 2.8.2 (2015-05-14)

* Gradients of Lagrangian after the laster iteration is skipped to avoid unnessary time cost.

## 2.8.1 (2015-05-12)

* Potlib subroutines can now use truncated expansions defined in an external basis file.  
* findmex utility has been updated to allow the user to set a number of input parameters through input file `mexopt.in`,
in namelist _mexopt_ .  One can now also set distance and angle constraints using the same input file. 

## 2.8.0 (2015-05-10)

* One can now use a term selection mechanism to select most useful terms from an extremely large expansion, in order
to limit the size of expansion while obtaining a more comprehensive basis.   Use jobtype=2 to generate an expansion,
and use `basisfl` variable to specify an expansion file when you use such expansion.  

## 2.7.7 (2015-05-08)

* A test term selection routine for jobtype=2 would determine a subset of independent terms in an expansion that spans
the same space as the full expansion.  This will be used in the future to automatically select terms from a very large
expansion.  

## 2.7.6 (2015-05-04)

* Use previously save set of diabat data to reconstruct Hd as the starting guess of the fit using `loadDiabats` option. 
* Updated manual page for `surfgen.in`
* Fixed a problem where redundant blocks were not correctly linked, causing memory consumption to be higher than necessary.

## 2.7.5 (2015-05-01)

* Added a procedure to record the values and gradients of each element of Hd at all data points.  This will be used to provide
  a way to reconstruct Hd for new expansions.
* Updated options for gfortran

## 2.7.4 (2015-04-27)

* Fixed a bug that after the QR decomposition with pivoting of projection operators in case of multi-dimensional irredicuble 
  representations the projection operators were not correctly pivoted back to the original basis.  The program is now tested
  for degenerate irreducible representations.
* Memory consumption information can now be explicitly printed out with printlvl>2.

## 2.7.3 (2015-04-18)

* Fixed a potential memory violation in EvaluateVal() in hddata.f90
* Bug fix in potlib.f90 and libutility.f90 that may cause memory leak

## 2.7.2 (2015-03-10)
* Minor bug fixes and print out improvements

## 2.7.1 (2015-02-16)

* Preliminary code to identify degenerate coordinates in the local coordinate generation process.  The rotation among them has
  not been implemented yet.
* New parameters `mng_scale_ener` and `mng_scale_grad` control if the scaled values or original values will be used by the 
  point manager.
* A new parameter `deggrdbinding` forces all gradient/couplings to be include or excluded at the same time if they are in the
  range of degeneracy.  The default is .true. which is recommended.  .false. will reproduce previous results. 

## 2.7.0 (2015-02-05)

* A procedure to identify the symmetry of each coordinate at every data point is tested to facilitate a future functionality
  that allows the fitting of points with states that carry degenerate representations.
* A point management feature is added.  The user can assign a set of points to be managed by the program by using the 
  MN keyword in `points.in` :   MN     X     I    J
  Here I and J defines the range of points to be managed.   Managed points will be excluded unless they are reasonably
  predicted by the fit.  New points that meet the criteria are added every iteration.
* Options nrmediff2 and ediffcutoff2 have been removed.  Energies are no longer scaled by energy difference, because they 
  are likely treated by a quasi-degeneracy when this is the case and the energy difference in non-diagonal representation is
  not meaningful.  This is in opposition to the gradient and coupling scaling, which only happens when energy difference is 
  large and a diagonal representation is likely used.
* Option gScaleMode has been removed.  Scaling can be turned off by setting intGradS=intGradT=0d0

## 2.6.15 (2015-01-23)

* Feature added to treat more-than-two-state quasi-degeneracy by orthogonalizing all g and h pairs, as described in JCP 141(17),174109
* A new data structure gWeight is added, which records the fitting weights of each gradient or coupling block for every data point. 
  The weight adjustment is due to high energy or large energy difference, and the relative weight assigned to w_fij and w_grad.
* Automatic assignment of phases in fixphase() now use a weighed sum of errors(using gWeight) instead of plain sum so that the important couplings
  can receive better representation.  This brings the error in line with actual fitting error.  This should make convergence smoother.

## 2.6.14 (2014-06-04)

* The program should no longer require a very large stacksize to run.

## 2.6.13 (2014-05-07)

* Fixed a problem where compilation may cause internal compiler error in earlier versions of ifort.

## 2.6.12 (2014-02-05)

* A switch `no_nad` is added to turn of couplings.

## 2.6.11 (2013-11-24)

* Minor bug fixes in findcp and readColGeom

## 2.6.10 (2013-11-4)

* Changed License to CC-BY
* Fixed a bug where gradient and energy input files are not properly closed.

## 2.6.9 (2013-10-18)

* fitinfo.csv now also contains the norm of fit gradients and couplings
* a new utility program analyse that can be used to show the internal coordinates of geometries
* Moved the internal coordinate analysis subroutine that is redundant in a number of utility programs to libutil.f90

## 2.6.8 (2013-10-12)

* Fixed torsion angle info output in findmex and findcp to correctly create list of linked torsions.
* New subroutine `getNeightbor` in potlib.f90 that allows the caller to identify the closest data point from the last evaluation.  If
distance calculation is not enabled, it returns 0.  
* tespoints utility now returns also the closest point with respect to the test point.   It now also properly handles the input argument
and use it as the geometry filename. 

## 2.6.7 (2013-10-07)

* new utilities program findmex that can be used to search for minimum energy point on N-state intersection seams on fit surface
* New subroutine OrthogonalizeGH that can be used by external programs to obtain intersection adapted representations
* adjusted bond-length threshold in findcp

## 2.6.6 (2013-10-06)

* OpenMP parallelization for evaluation of raw terms.
* New subroutine `finalizeSurfgen` in potlib.f90 that prints out the maximum deviation and evaluation count for the last trajectory

## 2.6.5 (2013-10-06)

* Updated `setsgenvars.sh` and `setsgenvars.csh` to set also the binary directory of surfgen as variable `SGENDIR`, and add surfgen directory to `PATH` variable if it is not already present.   
* Updated `findcp` utility program to correctly search for critical points.  It now use command line options to set geometry input file and initial state, so that the source code will not need ot be modified for different searches.   
* Utility programs will not be copied to binary directory upon installation.

## 2.6.4 (2013-10-05)

* When making libraries with `make libs`, recommended compilation and link flags for the program to use surfgen evaluation libraries 
will be printed out.   Two shell scripts, one for bash/sh one for csh/tcsh, will also be created in the `bin` directory that can be 
* sourced to set these environment variables and adjust stack size settings. 
PDF documentations now distributed with the repository.   
* A critical point search program `findmin` is added to the utilities.  
This program can be used to search for minima, saddle points, or conical intersections on the fit surface.
* More documentations for utility programs, along with installation scripts for each of them

## 2.6.3 (2013-10-02)

Changed format string in potlib to allow parsing of longer propagation times.

## 2.6.2 (2013-10-02)

Compatibility on Mac OS X.  Build OS X dynamic libraries with target `dylib` 

## 2.6.1 (2013-09-20)

Built-in evaluation counters added and is incremented on every EvaluateSurfgen calls.    
To retrieve the evaluation count, use GetEvalCount(count) subroutine.
The counter can be reset with the method ResetEvalCount().

## 2.6.0 (2013-09-13)

* Updated Makefile for smoother compilation on NERSC (tested for Hopper, Edison II and Carver)
* Minor modifications that allows compilation on Cray fortran compilers
* Fixed a few use-before-initialization bugs.

## 2.5.17 (2013-09-08)

* Slightly modified makefile to enable automatic detection based on $MKL_HOME variable
* Size of test decreased to more reasonable level.

## 2.5.16 (2013-09-01)

* New coordinate type anti-symmetric bend (type=1) added to describe bendings of atoms attached to a ring or chain

## 2.5.15 (2013-08-08)

* New option `cpcutoff` now allows couplings to be automatically removed above a set energy. 
works in a similar way as `gradcutoff`.  It is used to set the couplings cut off energy to be 
a different value from that of gradients.   If not set, the gradient cut off will be used instead.

## 2.5.14 (2013-08-08)

* The program also generates a fitinfo.csv file when PrintError=.true.
This file contains information needed to generate figures regarding the quality of fit.

## 2.5.13 (2013-08-05)

* Molecular information (atoms and number of states) are now written to the first record of trajectory log files 

## 2.5.12 (2013-07-21)

* Distance output now calculates the RMS average distance per non-vanishing coordinates.   A new parameter `cvanish` is
added and a coordinate is considered vanished when the norm of B-matrix elements of that coordinate is smaller than this value.
Setting it to 0 will reproduce the old distance definition.

## 2.5.11 (2013-07-15)

* Bug fix to prevent division by zero due to scaling function when distance is extremely large
* Format update for energy output.

## 2.5.10 (2013-07-04)

* A new utility program `ghplot` that creates two dimensional plots of the coupled PESs.
* Evaluation subroutine now records eigenvectors and can be retrieved with getEvec subroutine.

## 2.5.9 (2013-06-28)

* Removal of linear dependencies from exact equations
* New option gradcutoff allows gradients and coupling equations to be automatically removed when above certain energy threshold.

## 2.5.8 (2013-06-27)

* Increased the number of allowed input directories from 100 to 999.
* Minor output fixes.

## 2.5.7 (2013-06-26)

* Fixed a problem where the derivative couplings were not properly divided by energy difference.

## 2.5.6 (2013-06-26)

* Ability to plot size of couplings in evaluation libraries.

## 2.5.5 (2013-06-19)

* Updated `surfgen.in` manual pages for the options for evaluation subroutines.

## 2.5.4.1 (2013-06-18)

* Diagonal shift changed to a Morse form :  D*(1-w/w0)^2 , so that minimum will not be moved while shifting dissociation energy.

## 2.5.4 (2013-06-18)

* Diagonal shift by a list of single coordinate functions for size consistency corrections.
* During evaluation, RMS deviation per coordinate from existing point instead of total norm of deviation is displayed. 
* Partially constructed parallelization wrapper with MPI through ScaLapack.

## 2.5.3 (2013-06-17)

* Timing for normal equations and solution of Newton-Raphson equations are done separately.

## 2.5.2 (2013-06-13)

* Print out for contribution of each nascent coordinate during local internal construction.   This helps to diagnose reduced dimensionality issues.  Enabled when `printlvl>1`.

## 2.5.1 (2013-05-29)

* option `guide` allows the specification of a set of guide eigenvectors which will be used to determine state ordering at these points.
* enhanced output for detected state flipping during iterations or between loaded ckl files and automatic orderings.   
* singularities are now removed from linear equations of near degeneracy treatments

## 2.5.0 (2013-05-26)

The eigenvalue decomposition procedure has been replaced by linear equation solver `dsysv`.  This alleviates memory issues and speeds up the procedure.  
However, it can also be less reliable; specificly, it will fail when exact equations are linearly dependent.   This will be fixed in a future patch.

## 2.4.9 (2013-05-13)

* Added utilities `pauseParsing` and `resumeParsing` to allow programs to temporarily suspend parsing in evaluation subroutines.

## 2.4.8 (2013-05-12)

* Added option `printError` in group `MAKESURF`.  When true, the program will generate files `refgeom` and `error.log`, which can be used
by `libsurfgen` libraries to evaluate distances to existing points.
* Tested to integrate with ANT2012
* Rotation of Hd g and h vectors is now skipped when the corresponding rotation is forbiden by absence of data for ab initio data.
* Total fitting errors will be printed after final iteration.
* Updated manual pages for `surfgen.in`

## 2.4.7 (2013-05-10)

* Gradient ordering now also take into account all the errors of derivative couplings, in order to achieve better ordering accuracy.
* Fixed a problem with gradient ordering that can cause difficulty to obtain correct order when one of the states in the ordering group does not have gradients 

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
