MODULE Verlet_Algorithm
use READ_DATA
use Interaction_Cutoff_Modul
use PBC
use parallel_routines
implicit none
contains
SUBROUTINE VELO_VERLET(r,v,F)
    INTEGER i, k
    REAL*8 r(:,:),v(:,:),r0(n_particles,3),v0(n_particles,3),f0(n_particles,3)
    REAL*8 F(:,:),cutoff
    kinetic=0d0
    cutoff=0.99*L*5d-1
    CALL INTERACTION_CUTOFF(r,F,cutoff)

    DO i=index_matrix(taskid+1,1), index_matrix(taskid+1,2)
        r(i,:)=r(i,:)+v(i,:)*h+5d-1*F(i,:)*h*h
        v(i,:)=v(i,:)+5d-1*(F(i,:))*h
        r(i,1)=PBC2(r(i,1),L)
        r(i,2)=PBC2(r(i,2),L)
        r(i,3)=PBC2(r(i,3),L)
    END DO
    DO k=1,3
    CALL MPI_ALLGATHERV(r(index_matrix(taskid+1,1):index_matrix(taskid+1,2),k),&
                        & (index_matrix(taskid+1,2)-index_matrix(taskid+1,1)+1),MPI_DOUBLE_PRECISION, &
                        & r(:,k),num_send,desplac,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
    END DO
    CALL INTERACTION_CUTOFF(r,F,cutoff)

    DO i=index_matrix(taskid+1,1), index_matrix(taskid+1,2)
        v(i,:)=v(i,:)+5d-1*(F(i,:))*h
    END DO
   ! DO k=1,3
   !  ! CALL MPI_ALLGATHERV(v(index_matrix(taskid+1,1):index_matrix(taskid+1,2),k),&
   !  !                     & (index_matrix(taskid+1,2)-index_matrix(taskid+1,1)+1),MPI_DOUBLE_PRECISION, &
   !  !                     & v(:,k),num_send,desplac,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
   !  CALL MPI_GATHERV(v(index_matrix(taskid+1,1):index_matrix(taskid+1,2),k),&
   !                       & (index_matrix(taskid+1,2)-index_matrix(taskid+1,1)+1),MPI_DOUBLE_PRECISION, &
   !                       & v(:,k),num_send,desplac,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
   !  END DO

      !  IF(taskid==0) then
      !    DO i = 1,n_particles
      !      kinetic=kinetic+5d-1*(v(i,1)**2d0+v(i,2)**2d0+v(i,3)**2d0)
      !   end do
      ! end if
     ! DO i = index_matrix(taskid+1,1), index_matrix(taskid+1,2)
     !       kinetic=kinetic+5d-1*(v(i,1)**2d0+v(i,2)**2d0+v(i,3)**2d0)
     !    end do
     ! call MPI_REDUCE(kinetic,kinetic,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror) !(No funciona be)


END SUBROUTINE
END MODULE Verlet_Algorithm
