module Andersen_modul

use READ_DATA
use parallel_routines

IMPLICIT NONE

contains

    subroutine Andersen(v,temp)
    !Aquesta subrutina serveix per calcular les noves velocitats un cop apliquem el termostat d'Andersen 
    IMPLICIT NONE
    INTEGER i,k
    REAL*8 temp,nu,n1,n2,n3,n4,n5,n6,RAND
    REAL*8, DIMENSION(:,:) :: v
    nu=0.1/h
    sigma=sqrt(temp)
    !DO i=1,n_particles
        !IF(RAND().lt.nu*h)THEN
        !Iniciem aquest bucle per tal de fer una transformació de Box - Muller i obtenir una distribució normal de les velocitats.
        !taskid= identificador del processador
        IF (taskid.le.nworking_simple) THEN
            DO i=index_matrix(taskid,1),index_matrix(taskid,2)
                CALL RANDOM_SEED()
                CALL RANDOM_NUMBER(n1)
                CALL RANDOM_NUMBER(n2)
                CALL RANDOM_NUMBER(n3)
                CALL RANDOM_NUMBER(n4)
                CALL RANDOM_NUMBER(n5)
                CALL RANDOM_NUMBER(n6)
                v(n_particles,:)=(/sigma*sqrt(-2d0*log10(n1))*cos(2d0*3.1415*n2),sigma*sqrt(-2d0*log10(n1))*sin(2d0*3.1415*n2),&
                &sigma*sqrt(-2d0*log10(n3))*cos(2d0*3.1415*n4)/)
            END DO
        END IF
        !print*,'andersen before algather from proc',taskid,'force', v(1,:)
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        DO k=1,3
            CALL MPI_ALLGATHERV(v(index_matrix(taskid,1):index_matrix(taskid,2),k), &
                                &(index_matrix(taskid,2)-index_matrix(taskid,1)+1), &
                                &MPI_DOUBLE_PRECISION,&
                                & v(:,k),num_send, desplac,&
                                & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierror )
        END DO
        !print*,'andersen after algather',taskid,'force', v(1,:)
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        
    RETURN
    end subroutine Andersen

end module Andersen_modul