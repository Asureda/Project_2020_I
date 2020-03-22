module Andersen_modul

use READ_DATA
use matrix

IMPLICIT NONE

contains

    subroutine Andersen(v,temp)
    !Aquesta subrutina serveix per calcular les noves velocitats un cop apliquem el termostat d'Andersen 
    IMPLICIT NONE
    INTEGER i
    REAL*8 temp,nu,n1,n2,n3,n4,n5,n6
    REAL*8, DIMENSION(:,:) :: v
    nu=0.1/h
    sigma=sqrt(temp)
    !DO i=1,n_particles
        !IF(RAND().lt.nu*h)THEN
        !Iniciem aquest bucle per tal de fer una transformació de Box - Muller i obtenir una distribució normal de les velocitats.
        !taskid= identificador del processador
        IF (taskid.lt.n_working) THEN
            DO i=simple_matrix(taskid,1),simple_matrix(taskid,2)
                n1=RAND()
                n2=RAND()
                n3=RAND()
                n4=RAND()
                n5=RAND()
                n6=RAND()
                v(n_particles,:)=(/sigma*sqrt(-2d0*log10(n1))*cos(2d0*3.1415*n2),sigma*sqrt(-2d0*log10(n1))*sin(2d0*3.1415*n2),&
                &sigma*sqrt(-2d0*log10(n3))*cos(2d0*3.1415*n4)/)
            END DO
        END IF
        MPI_BARRIER(comm, ierror)
        
    RETURN
    end subroutine Andersen

end module Andersen_modul