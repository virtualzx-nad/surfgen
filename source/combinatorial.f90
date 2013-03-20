!********************************************************************
! combinatorial is used to generate factorials and permutations
!********************************************************************
!Public interface
!  factl(X)   : Returns the factorial of X
!  Permutation(M): Returns a M!*M list that contains all permutations
!                  of the first M natural integers.
!********************************************************************
!Algorithms
!factl:
!  factl_data contains all the factorials of integers up to UFactl
!  When the requested integer X is higher than UFactl, list factl_data
!  be expanded and extra factorials up to the requested number will be
!  calculated. Otherwise it will read the list and return the number
!  stored at factl_data(X).
!Permutation:
!  pmtList is a list of permutations of UPmt numbers. The permutations
!  are ordered in such a way that the upper-left M!*M block of pmtList
!  contains all permutations of first M integers when M<=UPmt.
!  When M>UPmt the list will be expanded.
MODULE combinatorial
 IMPLICIT NONE
 PRIVATE
 PUBLIC   ::  factl,Permutation
 INTEGER                                    :: UFactl=0
 INTEGER,DIMENSION(:),ALLOCATABLE           :: factl_data
 INTEGER                                    :: UPmt=0
 INTEGER,DIMENSION(:,:),ALLOCATABLE         :: pmtList

CONTAINS
 !FACTL(N) returns the factorial of N
 FUNCTION factl(N)
   IMPLICIT NONE
   INTEGER,INTENT(IN)               :: N
   INTEGER                          :: factl
   integer,dimension(:),allocatable :: tmp
   integer                          :: i
   if(.not.allocated(factl_data))then !intialization
     UFactl=5
     allocate(factl_data(0:5))
     factl_data=(/1,1,2,6,24,120/)
   end if!(.not.allocated(factl_data))
   if(N<0)then
     factl=0
   elseif(N<=UFactl)then
     factl=factl_data(N)
   else !N>UFactl, the list will be expanded
     allocate(tmp(0:UFactl))
     tmp=factl_data
     deallocate(factl_data)
     allocate(factl_data(0:N))
     factl_data(0:UFactl)=tmp
     deallocate(tmp)
     do i=UFactl+1,N
       factl_data(i)=factl_data(i-1)*i
     end do!i=UFactl+1,N
     factl=factl_data(N)
     UFactl=N
   end if
 END FUNCTION factl
!********************************************************************
 !PERMUTATION(M,MFact,MList) returns the list of all permutations of
 !elements of MList
 ! M       (input) INTEGER
 !         Number of element of the list
 ! MFact   (input) INTEGER
 !         Factorial of M, also the number of permutations
 ! MList   (input) INTEGER,dimension(M)
 !         List of elements of which the permutations will be generated.
 !         Each element is supposed to be identical.  Redundant entries
 !         however will not be checked by the function.
 ! PList   (output) INTEGER,dimension(MFact,M)
 !         Each row of PList contains a different permutation of the
 !         elements in MList
 FUNCTION Permutation(M,MFact,M_List) RESULT(PList)
   IMPLICIT NONE
   INTEGER,INTENT(IN)                       :: M,MFact
   INTEGER,DIMENSION(M),INTENT(IN),OPTIONAL :: M_List
   INTEGER,DIMENSION(MFact,M)               :: PList
   integer,dimension(:,:),allocatable       :: tmp
   integer                                  :: i,j,lb,ub
   integer,dimension(M)                     :: MList

   if(present(M_List))then
     MList=M_List
   else
     do i=1,M
       MList(i)=i
     end do
   end if!(present(M_List))
   if(MFact/=factl_data(M))stop "Permutation: inconsistent M and M! input."
   if(.not.allocated(pmtList) .or. UPmt<3)then !initialization
     UPmt=3
     allocate(pmtList(6,3))
     pmtList=reshape((/1,2,1,2,3,3,2,1,3,3,1,2,3,3,2,1,2,1/),(/6,3/))
   end if!(.not.allocated(pmtList) .or. UPmt<3)
   if(M<1)RETURN
   if(M>UPmt)then
     allocate(tmp(factl_data(UPmt),UPmt))
     tmp=pmtList
     deallocate(pmtList)
     allocate(pmtList(factl_data(M),M))
     pmtList(1:factl_data(UPmt),1:UPmt)=tmp
     deallocate(tmp)
     do i=UPmt+1,M
       lb=1
       ub=factl_data(i-1)
       pmtList(lb:ub,i)=i
       do j=i-1,1,-1
         lb=lb+factl_data(i-1)
         ub=ub+factl_data(i-1)
         pmtList(lb:ub,1:j-1)=pmtList(1:factl_data(i-1),1:j-1)
         pmtList(lb:ub,j)=i
         pmtList(lb:ub,j+1:i)=pmtList(1:factl_data(i-1),j:i-1)
       end do!j=1,i-1
     end do!i=UPmt+1,M
     UPmt=M
   end if!(M>UPmt)
   do i=1,M
     PList(:,i)=MList(pmtList(1:factl_data(M),i))
   end do
 END FUNCTION Permutation
END MODULE combinatorial
