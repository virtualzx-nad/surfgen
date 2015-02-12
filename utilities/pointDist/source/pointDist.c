/* This program generates a point distribution based upon
 * bond distances given in an input file: (given at CL)
 * Written by: cthree-40
 */
#include <stdio.h>
#include <stdlib.h>
struct geom {
  double x;
  double y;
  double z;
};
struct pointGeom {
  int               index; /* index of point */
  struct geom      *gdata; /* molecular geometry */
  struct pointGeom  *next; /* pointer to next node in ll */
};
struct pointDist {
  int              index; /* index of point */
  double           *dist; /* computed distances array */
  struct pointDist *next; /* pointer to next node in ll */
};
struct coord {
  int          *atoms;      /* atoms in coordinate */
  int            type;      /* type of coordinate  */
  struct coord  *next;      /* pointer to next node in ll */
};
void     readInput1 ( char *flname, int natoms, int ngeoms,
		                    int mcoord, struct coord *monCoord );
void     readInput2 (  int  natoms, int ngeoms, struct pointGeom *data );
void    processData (  int  natoms, int ngeoms, struct pointDist *outd );
double     compDist (                           int  atom1, int  atom2 );
void      printData (               int ngeoms, struct pointDist *outd );
/* main */
main ( int argc, char *argv[] )
{
  /* 
   * natoms = number of atoms (
   * ngeoms = number of geometries
   * data   = (index; geom atom 1, geom atom 2, etc )
   * outd   = (index; dist1, dist2, etc )
   */
  int                              natoms, ngeoms;
  char                                   *inputfl;
  struct coord                          *monCoord;
  struct pointGeom                          *data;
  struct pointDist                          *outd;
  /* process command line arguments:  [input file name] */
  if ( argc != 3 ) { /* argc should be 3 */
    printf( "\n Usage; %s [input filenmame] [number of atoms]\n", argv[0] );
    printf( "\n Exiting...\n");
    exit(1); /* exit */
  } else {
    inputfl = argv[1];
    natoms  = strtol( argv[2], NULL, 10 );
    printf( "\n User input file: %s", inputfl );
    printf( "\n Number of atoms: %d", natoms  );
  };
  /* allocate first node of monCoord array */
  monCoord = malloc( sizeof( struct coord ) );
  monCoord->atoms = malloc( natoms * sizeof( int ) );
  /* read user made input file */
  readInput1( inputfl, natoms, &ngeoms, *bonds );
};
