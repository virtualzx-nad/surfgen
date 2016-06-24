! Shared Program Information
module progdata
! TYPE DEFINITIONS
! Type abpoint     
! [Purpose] 
!   Stores all information of an ab initio data point
! [Fields]
! cgeom      DOUBLE PRECISION,dimension(3*natoms)
!            Cartesian geometry of all the atoms
! igeom      DOUBLE PRECISION,dimension(ncoord)
!            Scaled global internal geometry
! energy     DOUBLE PRECISION,dimension(nstates)
!            Ab initio energies of all the states
! grads      DOUBLE PRECISIOIN,dimension(3*natoms,nstates,nstates)
!            Ab initio energy gradients and derivative couplings
! hessian    DOUBLE PRECISION,dimension(3*natoms,3*natoms)
!            Hessian matrices for one of the states
! freqs      DOUBLE PRECISION,dimension(3*natoms)
!            Harmonic frequencies of one of the states
! bmat       DOUBLE PRECISION,dimension(ncoords,3*natoms)
!            Wilson's B matrix from Cartesian to scaled global internal coords,
!            or from local linear internal coordinates to scaled global internal
!            coordinates if useIntGrad==.true. 
! cbmat      DOUBLE PRECISION,dimension(ncoords,3*natoms)
!            B matrix in the original cartesian coordinate
! eval       DOUBLE PRECISION,dimension(3*natoms)
!            eigenvalues of B^t.B
! lmat       DOUBLE PRECISION,dimension(3*natoms,3*natoms)
! scale      DOUBLE PRECISION,dimension(3*natoms)
!            Scaling factor for each of the internal motions
! nvibs      INTEGER
!            lmat Transformation to local internal coordinates, defined as the
!            eigenvectors of B^T.B
!            nvibs Number of internal degree of freedom at this geometry
!            Those eigenvectors with eigenvalue larger than threshold intGradT
!            (stored in MODULE makesurfdata) are considered internal.  
!            The non-internal modes are also stored in the last fields of lmat
! ndeggrp    INTEGER
!            Number of groups of degenerate states. Non-degenerate state counts
!            as one group by itself.
! deg_groups INTEGER,dimension(ndeggrp,2)
!            Stores the index of the first state and number of states for each
!            degenerate group            
! lb,ub      INTEGER
!            Specifies the index range between which ab initio data is available
! nsym       INTEGER
!            Number of symmetry operations under which the geometry is invariant.
! symp,symi  INTEGER,dimension(:),allocatable
!            Symmetry operations that keep geometry invariant.
!      
  TYPE abpoint
   INTEGER                                       :: id
   DOUBLE PRECISION,dimension(:),allocatable     :: cgeom
   DOUBLE PRECISION,dimension(:),allocatable     :: igeom
   DOUBLE PRECISION,dimension(:,:),allocatable   :: energy
   DOUBLE PRECISION,dimension(:),allocatable     :: eval
   DOUBLE PRECISION,dimension(:,:,:),allocatable :: grads
   DOUBLE PRECISION,dimension(:,:),allocatable   :: hessian
   DOUBLE PRECISION,dimension(:),allocatable     :: freqs
   INTEGER                                       :: lb,ub

   INTEGER                                       :: nvibs
   DOUBLE PRECISION,dimension(:,:),allocatable   :: bmat,cbmat
   DOUBLE PRECISION,dimension(:,:),allocatable   :: lmat
   DOUBLE PRECISION,dimension(:),allocatable     :: scale

   INTEGER,dimension(:,:),allocatable            :: deg_groups
   INTEGER                                       :: ndeggrp
   INTEGER                                       :: nsym
   INTEGER,dimension(:),allocatable              :: symp,symi
  end TYPE abpoint

 ! derived type for a pointer to a two dimensional integer array
 type T2DList
   integer,dimension(2)                           ::  length
   integer,dimension(:,:),allocatable             ::  List
 end type


! Coordinate Scaling Function Definitions
! types  / scaling
!    0    rij
!         0   unscaled
!         1   Morse              Exp(-c1*(r-c2))
!         2   Gaussian           Exp(-c1*(r-c2)**2)
!         3   Reciprocal         Exp(-c1*(r-c2))/r
!         4   Chasmless Reciprcl Exp(-c1*(r-c2))/(r+c2)
!   -1    tetrahedron OOP.  all 6 internuclear distances are used to scale
!         0   scaled by power of the product of reciprocal distances
!         >0  use the product of scaled distance coordinates to scale the
!             scalar triple product
!         <0  use power of the reciprocal of the SUM of scaled distances to 
!             scale the scalar triple product
!   -2    umbrella OOP, 3 distance that link one vertex to the other three are
!         used in the scaling process
!         >=0 similar to tetrahedron OOP, except 3 distances instead of 6 are used
!         <0  use the harmonic mean of the scaled dist to scale the triple product
!    1    angles
!         0   unscaled
!         1   cos(theta)
!         2   cos(theta)/(1+exp(c1*(r1^2+r2^2-c2^2)))
!    2    torsion  UNIMPLEMENTED
!         0   unscaled
!         1   sin(theta)
!         2   cos(theta)
! Order specifies the maximum order of this specific set.
! Coef(nrij,2) coefficients for type=1 or 2
 type TCoordSet
   INTEGER                             ::  Type
   INTEGER                             ::  Scaling
   INTEGER                             ::  Order
   INTEGER                             ::  nCoord
   INTEGER,DIMENSION(:),allocatable    ::  iCoord
   INTEGER,DIMENSION(4)                ::  atom
   INTEGER,DIMENSION(:,:),allocatable  ::  Coord
   DOUBLE PRECISION,DIMENSION(:),allocatable ::  Coef
 end type

! CONSTANTS
  DOUBLE PRECISION,parameter                   :: AU2CM1  = 2.194746313705D5
  DOUBLE PRECISION,parameter                   :: AU2WAVE = 2.194746313705D5*Sqrt(dble(1)/1822.88848D0)
  DOUBLE PRECISION,parameter                   :: RAD2ANG = 57.29577951308D0
  DOUBLE PRECISION,parameter                   :: aJ2Eh = 1D0/4.359748D0
  DOUBLE PRECISION,parameter                   :: ang2bohr = 5.29177249D-1
  DOUBLE PRECISION,parameter                   :: PI = 3.14159265358979323846D0
  DOUBLE PRECISION,parameter                   :: EMASS = 1822.88848D0

  INTEGER,parameter                            :: INPUTFILE=12
  INTEGER,parameter                            :: IRREPFL=14
  INTEGER,parameter                            :: PTFL=15
  INTEGER,parameter                            :: OUTFILE=31

  INTEGER,parameter                            :: MAX_ALLOWED_SYM=10

  INTEGER                                      :: nCoordSets
  type(TCoordSet),dimension(:),allocatable     :: CoordSet
  INTEGER                                      :: nCoordCond
  INTEGER,DIMENSION(:,:),ALLOCATABLE           :: CoordCond
  INTEGER,DIMENSION(:),ALLOCATABLE             :: condRHS

! GLOBAL VARIABLES
  INTEGER                                      :: printlvl
  INTEGER                                      :: natoms
  CHARACTER(255)                               :: inputfl
  CHARACTER(255)                               :: basisfl
  CHARACTER(255)                               :: inputdir

! * MOLECULE PROPERTIES *
! typeCount            :  number of types of indistinguishable atoms.
! atomList(typeCount,:):  index list of each type of atoms
! atomCount(typeCount) :  number of atom from each type
 INTEGER                  :: typeCount
 INTEGER, DIMENSION(:),allocatable            :: atoms
 INTEGER, DIMENSION(:,:),allocatable          :: atomList
 INTEGER, DIMENSION(:),allocatable            :: atomCount

! * COORDIANTE PROPERTIES *
 INTEGER         ::  NRij,NOut
!coordmap maps each scaled coordinate back to (set_index,index_in_set)
 INTEGER,DIMENSION(:,:),ALLOCATABLE           ::  coordmap

! when set to true, the pot() funciton in potlib.f90 will switch the first 2
! diabatic states on output.
 LOGICAL                                      ::  switchdiab

! Connectivities to be preserved for filtering of non-feasible permutations
 CHARACTER(255)                               ::  cntfl

CONTAINS
 !***********************************************************************
 ! take the array of atom labels and replace then with integer indices,
 ! then writes the numeric array of labels and a map between strings
 !  and numbers to the module CNPI
 SUBROUTINE genAtomList(atmgrp)
  IMPLICIT NONE
  INTEGER           :: atmgrp(natoms)
  INTEGER           :: i
  character(3)      :: tpStr
  !deallocate global arrays if already allocated
  if (allocated(atomList))deallocate(atomList)
  if (allocated(atomCount))deallocate(atomCount)
  if (allocated(atoms))deallocate(atoms)
  typeCount = maxval(atmgrp)
  allocate(atomCount(typecount))
  allocate(atomList(typecount,natoms))
  allocate(atoms(natoms))
  atoms=atmgrp
  atomCount=0
  do i= 1, natoms
    atomCount(atmGrp(i))=atomCount(atmGrp(i))+1
    atomList(atmGrp(i),atomCount(atmGrp(i)))=i
  end do!i=1,natoms
  if(printlvl>0)then
    write (tpStr,"(I3)")natoms
    print "(12X,'Atoms:     [',"//trim(tpStr)//"I3,'  ]')",atmGrp
    write (tpStr,"(I3)")typeCount
    print "(12X,'Atom Count:[',"//trim(tpStr)//"I3,'  ]')",atomCount
  end if!(printlvl>0)
  NRij=natoms*(natoms-1)/2
  NOut=nrij*(natoms-2)*(natoms-3)/12
 END SUBROUTINE !genAtomList

END MODULE progdata
