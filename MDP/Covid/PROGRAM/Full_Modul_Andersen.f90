module Andersen_modul

use READ_DATA
use parallel_routines

IMPLICIT NONE

contains

    subroutine Andersen(v,temp)
    IMPLICIT NONE
    INTEGER i, k
    REAL*8 temp,nu,n1,n2,n3,n4 , RAND
    REAL*8 temp1
    REAL*8, DIMENSION(:,:) :: v
    nu=0.1
    temp1=sqrt(temp)
    DO i=index_matrix(taskid+1,1), index_matrix(taskid+1,2)
        IF(RAND().lt.nu)THEN
        !Iniciem aquest bucle per tal de fer una transformació de Box - Muller i obtenir una distribució normal de les velocitats.
            n1=RAND();n2=RAND()
            n3=RAND();n4=RAND()
            v(i,1)=temp1*sqrt(-2d0*log(n1))*cos(2d0*3.1415*n2)
            v(i,2)=temp1*sqrt(-2d0*log(n1))*sin(2d0*3.1415*n2)
            v(i,3)=temp1*sqrt(-2d0*log(n3))*cos(2d0*3.1415*n4)
        END IF
    END DO
    ! DO k=1,3
    ! CALL MPI_ALLGATHERV(v(index_matrix(taskid+1,1):index_matrix(taskid+1,2),k),&
    !                     & (index_matrix(taskid+1,2)-index_matrix(taskid+1,1)+1),MPI_DOUBLE_PRECISION, &
    !                     & v(:,k),num_send,desplac,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
    ! END DO
    end subroutine Andersen

end module Andersen_modul
