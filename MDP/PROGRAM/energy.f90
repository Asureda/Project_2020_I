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
    REAL*8 :: cutoff,pot,potential1
    REAL*8 :: dx,dy,dz,d,ff
    REAL*8, DIMENSION(:,:) :: r, F
    F=0d0
    potential=0d0
    pressure=0d0

    DO i=index_matrix(taskid,1),index_matrix(taskid,2)
        DO j=1,n_particles
            IF (j.ne.i) THEN
                dx=PBC1(r(i,1)-r(j,1),L)
                dy=PBC1(r(i,2)-r(j,2),L)
                dz=PBC1(r(i,3)-r(j,3),L)
                d=(dx**2d0+dy**2d0+dz**2d0)**0.5

                CALL L_J(d,ff,pot,cutoff)
                F(i,1)=F(i,1)+ff*dx/2.
                F(i,2)=F(i,2)+ff*dy/2.
                F(i,3)=F(i,3)+ff*dz/2.
                F(j,1)=F(j,1)-ff*dx/2.
                F(j,2)=F(j,2)-ff*dy/2.
                F(j,3)=F(j,3)-ff*dz/2.
                potential=potential+pot
                pressure=pressure+(ff*dx**2d0+ff*dy**2d0+ff*dz**2d0)
            END IF
        END DO
    END DO

    !call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(potential,potential,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(pressure,pressure,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    !call MPI_BARRIER(MPI_COMM_WORLD,ierror)
END SUBROUTINE INTERACTION_CUTOFF

END MODULE Interaction_Cutoff_Modul
