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
    DO i=index_matrix(taskid,1),index_matrix(taskid,2)
        CALL RANDOM_NUMBER(n5)
        IF(n5.lt.nu)THEN
        !CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(n1)
            CALL RANDOM_NUMBER(n2)
            CALL RANDOM_NUMBER(n3)
            CALL RANDOM_NUMBER(n4)
            
            v(i,:)=(/sigma*sqrt(-2d0*log(n1))*cos(2d0*3.1415*n2),sigma*sqrt(-2d0*log(n1))*sin(2d0*3.1415*n2),&
                  &sigma*sqrt(-2d0*log(n3))*cos(2d0*3.1415*n4)/)
         END IF
    END DO
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    end subroutine Andersen

end module Andersen_modul
