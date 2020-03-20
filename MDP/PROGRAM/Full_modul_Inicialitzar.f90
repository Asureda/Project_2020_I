module Inicialitzar

use READ_DATA

IMPLICIT NONE

contains

    subroutine FCC_Initialize(r)
    !En aquesta primera subrutina, definirem les posicions inicials de la nostra configuració. 
    INTEGER :: n,i,j,k
    REAL*8 :: r(:,:)
    n=1
    !print*,n
    !Definim aquesta configuració en les tres dimensions de l'espai
    DO i=0,M-1
        DO j=0,M-1
            DO k=0,M-1
                !print*,size(positions)
                r(n,:)=a*(/i,j,k/)
                !print*,positions(n,:)
                r(n+1,:)=r(n,:)+a*(/0.5,0.5,0.0/)
                r(n+2,:)=r(n,:)+a*(/0.5,0.0,0.5/)
                r(n+3,:)=r(n,:)+a*(/0.0,0.5,0.5/)
                n=n+4
            END DO
        END DO
    END DO
    PRINT*, 'particles positioned', n-1, 'of a total imput', n_particles
    RETURN
    end subroutine FCC_Initialize

    subroutine Uniform_velocity(v,T)
    !Aquesta subrutina, ens permetrà seguir una distribució normal al crear la matriu de les velocitats de la simulació.
    !Sempre es farà a una temperatura concreta
    INTEGER :: i,j,seed
    REAL*8 :: v(:,:),vi,vtot,T
    seed=13
    CALL SRAND(seed)
    DO i=1,3
        vtot=0
        DO j=1,n_particles-1
            vi=2*RAND()-1
            v(j,i)=vi
            vtot=vtot+vi
        END DO
        v(n_particles,i)=-vtot
    END DO
    CALL VELO_RESCALING_MOD(v,T)
    RETURN
    end subroutine Uniform_velocity
    subroutine Velo_Rescaling_mod(v,T)
    !Amb aquesta subrutina trobarem la velocitat, un cop hem reescalat, a una temperatura i un número de partícules concret
    IMPLICIT NONE
    REAL*8 v(:,:),T,alpha,KINETIC_LOC
    KINETIC_LOC=KINETIC_ENERGY(v)
    alpha=sqrt(3d0*n_particles*T/(2d0*KINETIC_LOC))
    v=alpha*v
    end subroutine Velo_Rescaling_mod

end module Inicialitzar