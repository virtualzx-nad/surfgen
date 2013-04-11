#ifndef SGENVER
#define SGENVER 'Unkown'
#endif

! This subroutine is called by other part of the program to retrieve the
! version number of surfgen. The version number is reported as a 72
! character string.
SUBROUTINE getver(ver)
        IMPLICIT NONE
        CHARACTER(72),intent(OUT)       :: ver

        ver = SGENVER


END SUBROUTINE getver
