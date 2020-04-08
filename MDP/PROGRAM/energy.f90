! Modul del cutoff
MODULE Interaction_Cutoff_Modul
use READ_DATA
use Lennard_Jones
use PBC
implicit none
contains

SUBROUTINE INTERACTION_CUTOFF(r,F,cutoff)
! Calcul de l'energia d'interaccio total donades una subrutina externa pel potencial
! Lennard-Jones i per les condicions periodiques de contorn. Tambe calculem la pressio
    IMPLICIT NONE
    INTEGER :: i,j,k
    REAL*8 :: cutoff,pot
    REAL*8 :: dx,dy,dz,d,ff
    REAL*8, DIMENSION(:,:) :: r, F
    !print*,'in inter'
    F=0d0
    potential=0d0
    !print*,'hola-1'
    pressure=0.0
    !print*,'hola'
    IF (paral_double.eqv..TRUE.)THEN
        IF (taskid.le.nworking_simple) THEN
            DO i=index_matrix(taskid,1),index_matrix(taskid,2)
                DO j=1,n_particles
                    IF (j.ne.i) THEN
                        dx=PBC1(r(i,1)-r(j,1),L)
                        dy=PBC1(r(i,2)-r(j,2),L)
                        dz=PBC1(r(i,3)-r(j,3),L)
                        d=(dx**2d0+dy**2d0+dz**2d0)**0.5

                        CALL L_J(d,ff,pot,cutoff)
                        F(i,1)=F(i,1)+ff*dx
                        F(i,2)=F(i,2)+ff*dy
                        F(i,3)=F(i,3)+ff*dz
                                !F(j,1)=F(j,1)-ff*dx
                                !F(j,2)=F(j,2)-ff*dy
                                !F(j,3)=F(j,3)-ff*dz
                        !WE ARE DOUBLE COUNTING POTENTIAL AND PRESSURE
                        !MUST DIVIDE BY TWO LATER !!!!
                        potential=potential+pot
                        pressure=pressure+(ff*dx**2d0+ff*dy**2d0+ff*dz**2d0)
                    END IF
                END DO
            END DO
        END IF
        !SINCE WE COMPUTE A SMALL PART OF THE FORCE MATRIX WE MUST
        !SEND EACH PART TO ALL OTHER PROCESSATORS
        print*,'before algather from proc',taskid,'force', f(1,:)
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        DO k=1,3
            CALL MPI_ALLGATHERV(f(index_matrix(taskid,1):index_matrix(taskid,2),k), &
                                &(index_matrix(taskid,2)-index_matrix(taskid,1)+1), &
                                &MPI_DOUBLE_PRECISION,&
                                & f(:,k),num_send, desplac,&
                                & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierror ) ! OPCIONAL crec, amb posicions es suficient
        END DO
        print*,'after algather',taskid,'force', f(1,:)
        call MPI_REDUCE( potential, potential,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
        call MPI_REDUCE(pressure,pressure,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    ELSE if(paral_double.eqv..false.)then
        if (taskid.eq.1) then
            DO i=1,n_particles
              DO j=i+1,n_particles
                dx=PBC1(r(i,1)-r(j,1),L)
                dy=PBC1(r(i,2)-r(j,2),L)
                dz=PBC1(r(i,3)-r(j,3),L)
                d=(dx**2d0+dy**2d0+dz**2d0)**0.5

                        CALL L_J(d,ff,pot,cutoff)
                        F(i,1)=F(i,1)+ff*dx
                        F(i,2)=F(i,2)+ff*dy
                        F(i,3)=F(i,3)+ff*dz
                        F(j,1)=F(j,1)-ff*dx
                        F(j,2)=F(j,2)-ff*dy
                        F(j,3)=F(j,3)-ff*dz
                        potential=potential+pot
                        pressure=pressure+(ff*dx**2d0+ff*dy**2d0+ff*dz**2d0)
                    END DO
                END DO
        end if
        !HERE WE COMPUTE THE FORCE, POTENTIAL AND PRESSURE IN ONE PROCESSATOR
        !WE MUST SEND THE FORCE TO ALL OTHER PROCESSATORS
        print*,'before algather'
        DO k=1,3
            CALL MPI_BCAST(F(1:n_particles,k),n_particles,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierror)
        END DO
        print*,'after algather'
    END IF

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    !call MPI_REDUCE(pressure,pressure,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
    !call MPI_REDUCE(potential,potential,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
END SUBROUTINE INTERACTION_CUTOFF

END MODULE Interaction_Cutoff_Modul
