! Modul de l'algoritme d'integració Verlet
! Integració de les equacions de Newton mitjançant l'integració de Verlet
! El modul consta de dos parts, en primer lloc el càlcul de les velocitats i posicions a partir de les configuracions inicials
! En segon lloc calcular les noves forces a partir de les posicions anteriors i tot seguit fer el càlcul de les noves velocitats v(t+h) a partir de les
! noves forces.
MODULE Verlet_Algorithm
use READ_DATA
use Interaction_Cutoff_Modul
use PBC
use parallel_routines
implicit none
contains
SUBROUTINE VELO_VERLET(r,v,F,iter)
    INTEGER i,k,iter,num,master
    REAL*8 r(:,:),v(:,:),r0(n_particles,3),v0(n_particles,3),f0(n_particles,3)
    REAL*8 F(:,:),cutoff,kinetic1
    cutoff=0.99*L*5d-1
    r0=r
    v0=v
    f0=f
    !print*,'initial verlet',taskid,KINETIC_ENERGY(v)
    CALL INTERACTION_CUTOFF(r,F0,cutoff)
    DO i=index_matrix(taskid,1),index_matrix(taskid,2)
        r(i,:)=r(i,:)+v(i,:)*h+5d-1*F0(i,:)*h*h
        r(i,1)=PBC2(r(i,1),L)
        r(i,2)=PBC2(r(i,2),L)
        r(i,3)=PBC2(r(i,3),L)
    END DO
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    DO k=1,3
    CALL MPI_ALLGATHERV(r(index_matrix(taskid,1):index_matrix(taskid,2),k),&
                        & (index_matrix(taskid,2)-index_matrix(taskid,1)+1),MPI_DOUBLE_PRECISION, &
                        & r(:,k),num_send,desplac,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
    END DO
    
    !print*,'r',r(100,:)

    CALL INTERACTION_CUTOFF(r,F,cutoff)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    kinetic=0d0
    kinetic1=0d0

    DO i=index_matrix(taskid,1),index_matrix(taskid,2)
        v(i,:)=v(i,:)+5d-1*(F(i,:)+F0(i,:))*h
        !kinetic(taskid)=kinetic(taskid)+5d-1*(v(i,1)**2d0+v(i,2)**2d0+v(i,3)**2d0)
        kinetic1=kinetic1+5d-1*(v(i,1)**2d0+v(i,2)**2d0+v(i,3)**2d0)
    END DO
    DO k=1,3
    CALL MPI_ALLGATHERV(v(index_matrix(taskid,1):index_matrix(taskid,2),k),&
                        & (index_matrix(taskid,2)-index_matrix(taskid,1)+1),MPI_DOUBLE_PRECISION, &
                        & v(:,k),num_send,desplac,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
    END DO
    DO k=1,3
    CALL MPI_ALLGATHERV(f(index_matrix(taskid,1):index_matrix(taskid,2),k),&
                        & (index_matrix(taskid,2)-index_matrix(taskid,1)+1),MPI_DOUBLE_PRECISION, &
                        & f(:,k),num_send,desplac,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
    END DO
    num=1
    master=0
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(kinetic1,kinetic,num,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
    !print*,kinetic,taskid,iter,'reduce'
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    
    
    !print*,'final verlet',taskid,KINETIC_ENERGY(v),kinetic
END SUBROUTINE
END MODULE Verlet_Algorithm
