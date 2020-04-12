module Andersen_modul

use READ_DATA

IMPLICIT NONE

contains

    subroutine Andersen(v,temp)
    IMPLICIT NONE
    INTEGER i
    REAL*8 temp,nu,n1,n2,n3,n4
    REAL*8 temp1
    REAL*8, DIMENSION(:,:) :: v
    nu=0.1
    temp1=sqrt(temp)
    DO i=1,n_particles
        IF(RAND().lt.nu)THEN
        !Iniciem aquest bucle per tal de fer una transformació de Box - Muller i obtenir una distribució normal de les velocitats.
            n1=RAND();n2=RAND()
            n3=RAND();n4=RAND()
            v(i,1)=temp1*sqrt(-2d0*log(n1))*cos(2d0*3.1415*n2)
            v(i,2)=temp1*sqrt(-2d0*log(n1))*sin(2d0*3.1415*n2)
            v(i,3)=temp1*sqrt(-2d0*log(n3))*cos(2d0*3.1415*n4)
        END IF
    END DO
    RETURN
    end subroutine Andersen

end module Andersen_modul
