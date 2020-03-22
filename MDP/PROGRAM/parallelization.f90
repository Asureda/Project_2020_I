module parallel_routines
    use read_vars
    IMPLICIT NONE

    contains
    SUBROUTINE simple_loop_matrix()
        IMPLICIT NONE
        integer :: i,j,a

        a=INT(REAL(numproc)/REAL(n_particles))
        allocate(index_matrix(numproc,2))
        index_matrix=0

        if(n_particles>numproc)then
          DO i=1,numproc
            index_matrix(i,1)=a*(i-1)+1
            index_matrix(i,2)=i*a
          ENDDO
          if (mod(numproc,n_particles)/=0)then
            index_matrix(numproc)=n_particles-numproc*a
          endif
        else
          DO i=1,n_particles
            index_matrix(i,1)=i
            index_matrix(i,2)=i+1
          ENDDO
        endif

    END SUBROUTINE

    SUBROUTINE doble_loop_matrix()
        IMPLICIT NONE
        integer :: i,j, aa, ites
        real*8 :: fact

        ites=int(FACT(numproc)/FACT(numproc-2)/2)
        aa=INT(REAL(ites)/REAL(numproc))
        allocate(double_matrix(numproc,4))
        DO i=1,n_particules
          DO j=i+1,n_particules
            double_matrix(i,1)=
            double_matrix(i,2)=
            double_matrix(i,3)=
            double_matrix(i,4)=
          ENDDO
        ENDDO
    END SUBROUTINE

    real*8 function FACT(n)
      implicit NONE
      integer :: n,i

      fact = 1.0
      DO i=2,n
        facte=fact*i
      ENDDO
    END FUNCTION FACT

END module
