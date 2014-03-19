/* This program generates geometries at (probably) some asymptote, where
 * the dissociated group (or atom) is rotated 2*pi in two planes.
 */

/* By Christopher L. Malbon (cthree-40)
 * Yarkony Group
 * The Johns Hopkins University
 * cmalbon1@jhu.edu
 */

/* Created		03-18-2014 */

/* Change Log------
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define M_PI 3.14159265358979323846

// Global Variables
int natoms; 	 	// Number of atoms
int AtomMid1;		// Atom to define midpoint
int AtomMid2;		// Atom to define midpoint
int AtomSpec;		// Atom to rotate 
int NumDisps;		// Number of displacements

// Functions
void printHeader();
void printFooter();
void readInput();
void readGeomFile(double Geometry[]);						// Reads in geometry
void printGeomFile(double Geometry[]);						// Prints geometry to stdout
void ShiftGeom(double Geometry[], double NewGeom[], double ShifVec[]);	// Shifts geometry, outputs shift vector, too
void rotationGen1( double g[], double j[] );
void RotVec1( double f[], double g[], double h[] );
void RotVec2( double f[], double g[], double h[] );
void GenRotMat( double g[], double a, double h[] );
void MatVec( int m, int n, double f[], double g[], double h[] );
void PerfRotation( double g[], double h[], int j, double i[] , int choice);
void RepCoord( double g[], double h[] );
void WriteFile( double g[], double j[], int displ, int choice );


double MidPoint(double g[], int a, int b, double h[]);


//-----------------------------------------------------------------------------------
main() // This is the main driver.
{
	printHeader(); 					// print header information
	readInput();					// read input information
	
	double geomStart[natoms * 3]; 		// Starting geometry in Row-Major Order
	double *gS;						// Pointer to geomStart
	gS = &geomStart[0];
	
	readGeomFile(gS);					// Read in geometry file
	printGeomFile(gS);				// Print original geometry file
	
	double geomOrigin[natoms * 3];		// Shifted Geometry
	double *gO;
	gO = &geomOrigin[0];
	
	double shiftVector[3];				// Shift Vector
	double *sV;
	sV = &shiftVector[0];
	
//	double RotVector1[3];				// Rotation in plane
//	double *rV;
//	rV = &RotVector1[0];
	
	ShiftGeom(gS,gO,sV);				// Shift Geometry, also returns shift vector
	rotationGen1(gO,sV);				// Perform Rotations
	
	printFooter();
}
//-----------------------------------------------------------------------------------

/*>printHeader
 */
void printHeader()
{
	printf( "\n\n\n=============================================================================\n"   );
	printf( "--------------------- SINGLE ATOM ROTATION GENERATOR-------------------------\n"   );
	printf( "-----------------------------------------------------------------------------\n\n" );
	printf( "cthree-40 \n" );
	printf( "Yarkony Group \n" );
	printf( "The Johns Hopkins University\n");
	printf( "-----------------------------------------------------------------------------\n"   );
	printf( "-----------------------------------------------------------------------------\n"   );
}
/*>printFooter
 */
void printFooter()
{
	printf( "\n\n"  );
	printf( "                    ~ DONE ~\n\n" );
	printf( "-----------------------------------------------------------------------------\n"   );
	printf( "=============================================================================\n\n\n"   );
}
/*>readInput
 */
void readInput()
{
	FILE *input; // Set file pointer
	input = fopen("input.SAR", "r");
	if ( input == NULL ){
		printf("Can't open rotation.in\n");
		exit(1);
	}
	
	printf( "\n         Reading input file...\n\n");
	
	
	while (fscanf(input, "%d\n%d\n%d\n%d%d\n", &natoms, &AtomMid1, &AtomMid2, &AtomSpec, &NumDisps ) != EOF) {
		  printf(" -- The number of atoms in the system is %d \n", natoms );
		  printf(" -- The atoms between which the origin will be is %d  and  %d \n", AtomMid1, AtomMid2 );
		  printf(" -- The atom to rotate is atom %d \n", AtomSpec );
		  printf(" -- The number of points along the rotation path will be  %d \n", NumDisps );
	}
	fclose(input);
}	

/*>readGeomFile
 * GeometryP[]	= Geometry to be read in
 */
void readGeomFile( double GeometryP[] )
{
	printf( "\n         Reading geometry file...\n\n");
	
	FILE *refgeom; //Set file pointer
	refgeom = fopen("geom.start","r");
	if ( refgeom == NULL ) {
		printf("Can't open refgeom.\n");
		exit(1);
	}
	
	int i;
	for ( i=0; i < natoms; i++ ){
		fscanf( refgeom, "%lF  %lF  %lF\n", &GeometryP[(i*3+0)], &GeometryP[(i*3+1)], &GeometryP[(i*3+2)] );
	}
	close(refgeom);
	
}

/*>printGeomFile
 * GeomP[]	= Array of natoms*3 elements containing the starting geometry
 */
void printGeomFile( double GeomP[] )
{	
	int j;
	
	printf( "\n         Starting Geometry: \n\n");
	for ( j=1; j <= natoms; j++){
		printf( " Atom%2d:    %10lF  %10lF  %10lF \n", j, GeomP[(j-1)*3+0], GeomP[(j-1)*3+1], GeomP[(j-1)*3+2] );
	}
}

/*>ShiftGeom
 * GeomP[]	= Array of natoms*3 elements containing the starting geometry
 * NewGP[]	= Array of natoms*3 elements containing the shifted geometry
 * MidVP[]	= Array of 3 elements containg the cartesian coordinates of the Midpoint
 */
void ShiftGeom( double GeomP[], double NewGP[], double MidVP[] )
{	
	int i; 
	double midPoint[3];					// Location of midpoint btwn. AtomMid1 and AtomMid2
	double *mpP;						// Pointer to midPoint
	mpP = &midPoint[0];
	
	printf( "\n Computing midpoint between atoms %d and %d ...\n\n", AtomMid1, AtomMid2 );
	MidPoint( GeomP, AtomMid1, AtomMid2, MidVP );		// Compute midpoint. Returns vector

	for (i=0; i < natoms; i++ ){
		NewGP[i*3 +0] = GeomP[i*3 +0] - MidVP[0];		// x
		NewGP[i*3 +1] = GeomP[i*3 +1] - MidVP[1];		// y
		NewGP[i*3 +2] = GeomP[i*3 +2] - MidVP[2];		// z
	}

	printf( "           Shifted Geometry: \n\n");
	for ( i=1; i<=natoms; i++ ){
		printf( " Atom%2d:    %10lF  %10lF  %10lF \n", i, NewGP[(i-1)*3+0], NewGP[(i-1)*3+1], NewGP[(i-1)*3+2] );
	}
}

/*>MidPoint
 */
double MidPoint( double gP[], int a, int b, double hP[] )
{
	int i;	// Loop variable
	int atom1;	// Location of atom1
	int atom2;	// Location of atom2
	
	atom1 = (a-1)*3;
	atom2 = (b-1)*3;
	
	for ( i=0; i<3; i++ ){
		hP[i] = ( gP[atom1+i] + gP[atom2+i] ) / 2;
	}
}
		
/*>rotationGen1
 *  Rotation of AtomSpec out of plane defined by AtomMid1 AtomMid2 and AtomSpec
 */
void rotationGen1( double gP[], double ShiftVector[] )
{
	double *At1P, *AtSP;						// Pointers to AtomMid1 and AtomSpec
	double *rV1, *rV2, *rV3;					// Pointers to Rotation vectors 1,2 and 3, respectively
	double rotVecs[9];						// Rotation vectors are stored in same array
	
	At1P = &gP[(AtomMid1-1)*3];					// Set pointers
	AtSP = &gP[(AtomSpec-1)*3];
	
	
	RotVec1( At1P, AtSP, &rotVecs[0] );							// Generate Rotation vector 1 (Orthogonal to plane)
	RotVec1( AtSP, &rotVecs[0], &rotVecs[3] );					// Generate Rotation vector 2 (In plane)
	RotVec2( &rotVecs[0], &rotVecs[3], &rotVecs[6] );				// Generate Rotation vector 3 
	
	
	double *rMP;							// Rotation Matrix Pointer
	double rotationMatrix[9];					// Rotation Matrix
	rMP = &rotationMatrix[0];
	
	double rotints;							// Intervals of rotation
	
	rotints = 2*M_PI/ NumDisps;					// Angle of rotation
	
	int i,j;
	
	for ( i=0; i<3; i++ ){
		GenRotMat( &rotVecs[i*3], rotints, rMP );					// Generate rotation matrix
		j=i+1;
		if ( i==0 ){
			printf( "\n Generating out-of-plane points...\n");
		} else if ( i==1 ){
			printf( "\n Generating in-plane points...\n");
		} else {
			printf( "\n Generating points around the R3 axis...\n");
		}
		PerfRotation( rMP, gP, NumDisps, ShiftVector, j );			// Perform rotations, geometries will be generated in files
	}

	
}

/*>PerfRotation
 * RotMatP		= Rotation Matrix [9]
 * MolGeoP		= Geometry [natoms*3]
 * disp		= Number of rotations to perform
 * ShiftVector	= Shift vector
 * choice		= which rotation to perform
 */
void PerfRotation(double RotMatP[], double MolGeoP[], int disp, double ShiftVector[], int choice )
{
	double *gSCr, *nVec;						// Scratch geometry for manipulation, new vec for special atom
	double geomScratch[natoms*3], newVec[3];			//
//	nVec = &newVec[0];			
	
	
	int i,j;
	
	for ( i=0; i<=natoms*3; i++ ){
		geomScratch[i] = MolGeoP[i];				// Set scratch geometry to original geometry
	}
		
	for ( i=1; i<=disp; i++ ){
		nVec = &newVec[0];
		MatVec( 3, 3, RotMatP, &geomScratch[(AtomSpec-1)*3], nVec );		// Perform Matrix multplication ( Rotation )
		RepCoord( geomScratch, newVec );					// Replace Coordinates of Special atom
		WriteFile( geomScratch, ShiftVector, i, choice );
	}
}

/*>WriteFile
 * Write geometry to file
 * Input:
 * Geometry		= Geometry
 * ShiftVector 	= Shift vector
 * i			= Indice to append file name with 
 * rotation 	= Curren rotation
 */
void WriteFile( double Geometry[], double ShiftVector[], int index, int rotation )
{
	int i;	
	char DirName[25];

	sprintf( DirName, "PointsR%d", rotation );				// Directory name	
	
	char FileName[25];								// File name
	FILE *gfP;										// Geometry file pointer
	
	double geomPrint[natoms*3];							// Geometry to be printed. This geometry is shifted.
	double *gPP;									// Pointer to above geometry
	gPP = &geomPrint[0];
	
	for (i=0; i<natoms; i++ ){
		geomPrint[i*3+0] = Geometry[i*3+0]+ShiftVector[0];
		geomPrint[i*3+1] = Geometry[i*3+1]+ShiftVector[1];
		geomPrint[i*3+2] = Geometry[i*3+2]+ShiftVector[2];
	}	
		
	mkdir( DirName, 0777 );								// Make directory
	
	sprintf( FileName, "PointsR%d/loopgeom.%d", rotation, index );	// Name file
	
	gfP = fopen( FileName, "w" );							// Open file for writing
		
	if ( gfP == NULL ){
		fprintf(stderr, "Cannot open loopgeom.%d", index);
		exit(1);
	}
	
	for ( i = 0; i < natoms; i++ ){
		fprintf( gfP, " %19.15lF %19.15lF %19.15lF \n", geomPrint[i*3+0], geomPrint[i*3+1], geomPrint[i*3+2] );
	}
	
	fclose(gfP);
}
	
	
	
/*>RepCoord
 * Geom	= Molecular geometry
 * NewCoord = New coordinates for AtomSpec
 */
void RepCoord( double Geometry[], double NewCoord[] )
{
	int i;
										// Replace Coordinates with rotated coords
	Geometry[(AtomSpec-1)*3+0] = NewCoord[0];
	Geometry[(AtomSpec-1)*3+1] = NewCoord[1];
	Geometry[(AtomSpec-1)*3+2] = NewCoord[2];
}

/*>MatVec
 * Rows	= Number of rows of Matrix
 * Cols	= Number of columns of Matrix
 * Matrix[] = Single dimension array containing values of matrix in row-major order ( This is C ! )
 * Vector[] = Vector to perform Mv on
 * Output[] = Output vector
 */
void MatVec( int Rows, int Cols, double Matrix[], double Vector[], double Output[] )
{
	int i, j; 				// Loop variables
	Output[0] = 0.0;
	Output[1] = 0.0;
	Output[2] = 0.0;
	for (i=0; i<Rows; i++ ){
		for (j=0; j<Cols; j++ ){
			Output[i] = Output[i] + Matrix[i*Cols + j]*Vector[j]; // Hij vj = v'i
		}
	}
}
	
		
	
	
/*>RotVec1
 */
void RotVec1(double aP[], double bP[], double hP[] )
{

	hP[0] = (aP[1]*bP[2] - aP[2]*bP[1]);
	hP[1] = (aP[2]*bP[0] - aP[0]*bP[2]);
	hP[2] = (aP[0]*bP[1] - aP[1]*bP[0]);
	
	double norm = sqrt(pow(hP[0],2) + pow(hP[1],2) + pow(hP[2],2));
	
	hP[0] = hP[0]/norm;
	hP[1] = hP[1]/norm;
	hP[2] = hP[2]/norm;
	
	norm = 0.0;
	
	norm = pow(hP[0],2) + pow(hP[1],2) + pow(hP[2],2);
	
}
/*>RotVec2
 * Generates rotation vector 3
 * vector1[]	= rotvec1
 * vector2[]	= rotvec2
 * Output[]		= rotvec3
 */
void RotVec2( double vector1[], double vector2[], double Output[] )
{
	int i;
	
	for (i=0; i<3; i++ ){
		Output[i] = vector1[i] + vector2[i];
	}
	
	double norm = sqrt(pow(Output[0],2) + pow(Output[1],2) + pow(Output[2],2));
	
	Output[0] = Output[0]/norm;
	Output[1] = Output[1]/norm;
	Output[2] = Output[2]/norm;
}
	
		
/*>GenRotMat	i	
 * aP[]	= Axis of rotation	(input)
 * interval	= Rotation angle		(input)
 * bP[]	= Rotation matrix 	(output)
 */
void GenRotMat( double aP[], double interval, double bP[] )
{
	double t = 1.0 - cos(interval);
	double S = sin(interval);
	double C = cos(interval);

	bP[0] = t*pow(aP[0],2) + C;
	bP[1] = t*aP[0]*aP[1] - S*aP[2];
	bP[2] = t*aP[0]*aP[2] + S*aP[1];
	bP[3] = t*aP[0]*aP[1] + S*aP[2];
	bP[4] = t*pow(aP[1],2) + C;
	bP[5] = t*aP[1]*aP[2] - S*aP[0];
	bP[6] = t*aP[0]*aP[2] - S*aP[1];
	bP[7] = t*aP[1]*aP[2] + S*aP[0];
	bP[8] = t*pow(aP[2],2) + C;

#ifdef DEBUG
	printf( "\n====================================================\n");
	printf( "*** DEBUG .ne. NULL ...\n");
	printf( " Debugging the Rotation Matrix generating routine, GenRotMat(...)...\n" );
	printf( " If v is the vector defining the axis of rotation, then Rv=v. \n");
	
	double TestVector[3];	// Testing vector
	double TestOutput[3];	// Output  vector from test
	int i;
	
	for (i=0;i<3;i++){
		TestVector[i] = aP[i];
	}
	
	MatVec(3, 3, bP, TestVector, TestOutput );
	for (i=0;i<3;i++){
		printf( "\n %10lF should equal %10lF \n ", TestVector[i], TestOutput[i] );
	}
	printf( "\n====================================================\n");
#endif	
}
	
		
		
		
		
	
		
