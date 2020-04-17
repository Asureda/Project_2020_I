module parallel_routines
    use read_data
    IMPLICIT NONE

    contains
    SUBROUTINE simple_loop_matrix()
        IMPLICIT NONE
        integer :: i,j,a, partition
        integer, dimension(:,:), allocatable :: index_matrix2
        integer, allocatable :: pairindex(:,:)
        a=INT(REAL(n_particles)/REAL(numproc))
        allocate(index_matrix(numproc,2),desplac(numproc),num_send(numproc))
        allocate(index_matrix2(numproc,2))
        integer,intent(in):: Pair((N*(N-1))/2,2)

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

            partition = INT(REAL(n_particles*(n_particles-1)/2d0)/REAL(numproc))

            DO i=1,numproc
              index_matrix2(i,1)=partition*(i-1)+1
              if ( i /= numproc) then
              index_matrix2(i,2)=i*partition
              else
              index_matrix2(numproc,2)=partition
            end if
            END DO


            !Define vector of pairs of particles
              n = 1
              do i = 1,npar-1
                  do j = i+1,npar
                      pairindex(n,:) = (/i,j/)
                      n = n + 1
                  enddo
              enddo

        !print*,'CPU',taskid+1,'nº particles', index_matrix(taskid+1,2)-index_matrix(taskid+1,1)+1
        !print*,'CPU',taskid+1,'relative displacement', desplac(taskid+1)
    END SUBROUTINE
END module
