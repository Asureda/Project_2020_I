module parallel_routines
    use read_data
    IMPLICIT NONE

    contains
    SUBROUTINE simple_loop_matrix()
        IMPLICIT NONE
        integer :: i,j,a

        a=INT(REAL(n_particles)/REAL(numproc))
        !print*,'particules',n_particles,'CPUs',numproc,'part/CPU',a
        allocate(index_matrix(numproc,2),desplac(numproc),num_send(numproc))
        index_matrix=0

        !IF (paral_simple.eqv..TRUE.) THEN
          if(n_particles>numproc)then
            nworking_simple=numproc
            DO i=1,numproc
              index_matrix(i,1)=a*(i-1)+1
              index_matrix(i,2)=i*a
            ENDDO
            index_matrix(numproc,2)=n_particles
          ! else
          !   nworking_simple=n_particles
          !   DO i=1,n_particles
          !     index_matrix(i,1)=i
          !     index_matrix(i,2)=i
          !   ENDDO
          endif
        !else
          !nworking_simple=1
          !index_matrix(1,1)=1
          !index_matrix(1,2)=n_particles
        !END IF
        
        DO i=2,numproc
          desplac(i)=index_matrix(i-1,2)
          num_send(i)=index_matrix(i,2)-index_matrix(i,1)+1
        END DO
        desplac(1)=0
        num_send(1)=index_matrix(1,2)-index_matrix(1,1)+1
    END SUBROUTINE
END module
