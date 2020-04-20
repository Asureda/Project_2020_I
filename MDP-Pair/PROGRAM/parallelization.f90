module parallel_routines
    use read_data
    IMPLICIT NONE

    contains
    SUBROUTINE simple_loop_matrix()
        IMPLICIT NONE
        integer :: i,j,a,res,b,b_res,ii
        a=INT(REAL(n_particles)/REAL(numproc))
        res=mod(n_particles,numproc)
        allocate(index_matrix(numproc,2),desplac(numproc),num_send(numproc))
        allocate(index_matrix2(numproc,2))
        allocate(pairindex((n_particles*(n_particles-1))/2,2))

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
        ! if(taskid==0) then
        !   DO i =1,numproc
        !     print*,index_matrix(i,:)
        !   END DO
        !   print*,desplac(:)
        !   print*,num_send(:)
        ! end if

        !Define vector of pairs of particles
        ii = 1
        do i = 1,n_particles-1
            do j = i+1,n_particles
                pairindex(ii,:) = (/i,j/)
                ii = ii + 1
            enddo
        enddo
      b =  INT(REAL((n_particles*(n_particles-1))/2)/REAL(numproc))
      b_res=mod((n_particles*(n_particles-1))/2,numproc)
      DO i=1,numproc
        IF((i-1)<b_res)then
          index_matrix2(i,1)=(b+1)*(i-1)+1
          index_matrix2(i,2)=index_matrix2(i,1)+b
        else
        index_matrix2(i,1)=(i-1)*b+b_res+1
        index_matrix2(i,2)=index_matrix2(i,1)+b-1
        end if
      END DO


    END SUBROUTINE
END module
