module parallel_routines
    use read_data
    IMPLICIT NONE

    contains
    SUBROUTINE simple_loop_matrix()
        IMPLICIT NONE
        integer :: i,j,a

        a=INT(REAL(n_particles)/REAL(numproc))
        print*,'particules',n_particles,'CPUs',numproc,'part/CPU',a
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

    SUBROUTINE double_loop_matrix()
        IMPLICIT NONE
        integer :: i,j, aa, ites, suma,suma2

        ites=0
        DO i=1,n_particles
          DO j=i+1,n_particles
            ites=ites+1
          END DO
        END DO
        aa=INT(REAL(ites)/REAL(numproc))
        print*,'doble bucle',ites,'CPUs',numproc,'ites/CPU',aa, 'test', aa*numproc
        allocate(double_matrix(numproc,4))
        double_matrix=0
        IF (paral_double.eqv..TRUE.) THEN
        if(ites.gt.numproc)then
          nworking_double=numproc
          suma=0
          suma2=0
          DO i=1,n_particles
            DO j=i+1,n_particles
              suma=suma+1
              if((mod(suma,aa).eq.0).and.(suma2.lt.numproc-1)) then
                !print*,'test',suma,aa,i,j
                suma2=suma2+1
                double_matrix(suma2,2)=i
                double_matrix(suma2,4)=j
                if(j==n_particles)then
                  double_matrix(suma2+1,1)=i+1
                  double_matrix(suma2+1,3)=i+2
                else
                  double_matrix(suma2+1,1)=i
                  double_matrix(suma2+1,3)=j+1
                endif
              endif
            ENDDO
          ENDDO
          !print*,'suma............................',suma,suma2
          double_matrix(1,1)=1
          double_matrix(1,3)=2
          double_matrix(numproc,2)=n_particles
          double_matrix(numproc,4)=n_particles
        else
          nworking_double=ites
          suma=0
          DO i=1,n_particles
            DO j=i+1,n_particles
              suma=suma+1
              double_matrix(suma,1)=i
              double_matrix(suma,2)=i
              double_matrix(suma,3)=j
              double_matrix(suma,4)=j
            ENDDO
          ENDDO
          double_matrix(numproc,2)=n_particles
          double_matrix(numproc,4)=n_particles
        endif
       END IF
    END SUBROUTINE
END module
