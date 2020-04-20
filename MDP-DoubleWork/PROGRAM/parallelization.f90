module parallel_routines
    use read_data
    IMPLICIT NONE

    contains
    SUBROUTINE simple_loop_matrix()
        IMPLICIT NONE
        integer :: i,j,a,res
        a=INT(REAL(n_particles)/REAL(numproc))
        res=mod(n_particles,numproc)
        allocate(index_matrix(numproc,2),desplac(numproc),num_send(numproc))

            DO i=1,numproc
              IF((i-1)<res)then
                index_matrix(i,1)=(a+1)*(i-1)+1
                index_matrix(i,2)=index_matrix(i,1)+a
              else
              index_matrix(i,1)=(i-1)*a+res+1
              index_matrix(i,2)=index_matrix(i,1)+a-1
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
        if(taskid==0) then
          DO i =1,numproc
            print*,index_matrix(i,:)
          END DO
          print*,desplac(:)
          print*,num_send(:)
        end if
        ! print*,'CPU',taskid+1,'nº particles', index_matrix(taskid+1,2)-index_matrix(taskid+1,1)+1
        ! print*,'CPU',taskid+1,'relative displacement', desplac(taskid+1)
    END SUBROUTINE
END module
