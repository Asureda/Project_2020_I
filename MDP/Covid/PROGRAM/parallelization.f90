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


            DO i=1,numproc
              index_matrix(i,1)=a*(i-1)+1
              if ( i /= numproc) then
              index_matrix(i,2)=i*a
              else
              index_matrix(numproc,2)=n_particles
            end if
            END DO

        DO i=1,numproc
          if(i.eq.1)then
            desplac(1)=0
          else
            desplac(i)=index_matrix(i-1,2)
            end if
            num_send(i)=index_matrix(i,2)-index_matrix(i,1)+1
        end do
    END SUBROUTINE
END module
