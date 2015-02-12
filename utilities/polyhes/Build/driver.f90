program driver
! Driver for optimizing an intersection with
!  polyhes in surfJay.
  implicit none
  integer              :: maxIter, icode
  integer              :: i
! Loop over polyhes
  icode = 0
  maxIter = 1
  do i=1, maxIter
     call polyhes(icode)
     print *, "icode=",icode
  end do
end program driver
