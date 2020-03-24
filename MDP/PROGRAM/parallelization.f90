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
          nworking_simple=numproc
          DO i=1,numproc
            index_matrix(i,1)=a*(i-1)+1
            index_matrix(i,2)=i*a
          ENDDO
          index_matrix(numproc)=n_particles
        else
          nworking_simple=n_particles
          DO i=1,n_particles
            index_matrix(i,1)=i
            index_matrix(i,2)=i
          ENDDO
        endif

    END SUBROUTINE

    SUBROUTINE double_loop_matrix()
        IMPLICIT NONE
        integer :: i,j, aa, ites, suma
        real*8 :: fact

        ites=int(FACT(numproc)/FACT(numproc-2)/2)
        aa=INT(REAL(ites)/REAL(numproc))
        allocate(double_matrix(numproc,4))
        if(ites>numproc)then
          suma=0
          DO i=1,n_particules
            DO j=i+1,n_particules
              suma=suma+1
              if(mod(suma,aa)==0)then
                double_matrix(i,1)=i
                double_matrix(i,3)=j
                if(j==n_particles)then
                  double_matrix(i,2)=i+1
                  double_matrix(i,4)=i+2
                else
                  double_matrix(i,2)=i
                  double_matrix(i,4)=j+1
                endif
              endif
            ENDDO
          ENDDO
          double_matrix(numproc,2)=n_particles
          double_matrix(numproc,4)=n_particles
        else
          suma=0
          DO i=1,n_particules
            DO j=i+1,n_particules
              double_matrix(i,1)=i
              double_matrix(i,2)=i
              double_matrix(i,3)=j
              double_matrix(i,4)=j
            ENDDO
          ENDDO
          double_matrix(numproc,2)=n_particles
          double_matrix(numproc,4)=n_particles
        endif
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
