module Andersen_modul

use READ_DATA

IMPLICIT NONE

contains

    subroutine Andersen(v,temp)
    IMPLICIT NONE
    INTEGER i
    REAL*8 temp,nu,n1,n2,n3,n4,n5,n6
    REAL*8, DIMENSION(:,:) :: v
    nu=0.1/h
    sigma=sqrt(temp)
    DO i=1,n_particles
        IF(RAND().lt.nu*h)THEN
            n1=RAND();n2=RAND()
            n3=RAND();n4=RAND()
            n5=RAND();n6=RAND()
            v(n_particles,:)=(/sigma*sqrt(-2d0*log10(n1))*cos(2d0*3.1415*n2),sigma*sqrt(-2d0*log10(n1))*sin(2d0*3.1415*n2),&
                &sigma*sqrt(-2d0*log10(n3))*cos(2d0*3.1415*n4)/)
        END IF
    END DO
    RETURN
    end subroutine Andersen

end module Andersen_modul
