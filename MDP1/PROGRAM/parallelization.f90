module parallel_routines
    use read_data
    IMPLICIT NONE

    contains
    SUBROUTINE simple_loop_matrix()
        IMPLICIT NONE
        integer :: i,j,a
        a=INT(REAL(n_particles)/REAL(numproc))
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
                desplac(i)=0
                else
                desplac(i)=index_matrix(i-1,2)
                end if
                num_send(i)=index_matrix(i,2)-index_matrix(i,1)+1
            end do

        !print*,'CPU',taskid+1,'nº particles', index_matrix(taskid+1,2)-index_matrix(taskid+1,1)+1
        !print*,'CPU',taskid+1,'relative displacement', desplac(taskid+1)
    END SUBROUTINE
END module
