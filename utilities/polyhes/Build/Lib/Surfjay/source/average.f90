! Subroutine to compute the real average of an array
subroutine average( Array, ALength, AverageOut)
        implicit none
        integer, intent(in)   :: ALength
        real*8, dimension( ALength ), intent(in) :: Array
        real*8,  intent(out)  :: AverageOut
        integer               :: i
        real*8                :: ArrSum
        !---------------------------------------------------
        ArrSum = 0d0                      ! Zero out ArrSum
        do i=1, ALength
            ArrSum = ArrSum + Array(i)    ! Add elements
        end do
        ! Divide by ALength
        AverageOut = ArrSum / real( ALength )
        return
end subroutine
