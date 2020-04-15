!GRUP I: Àlex, Oriol, Laia, Sílvia i Elena

MODULE Distribucio_Uniforme_vel

use READ_DATA
use Reescala_velocitats

implicit none

contains

SUBROUTINE UNIFORM_VELO(v,T)

    !OBJECTIU: Aconseguir una distribució uniforme de les velocitats
    !          Després, a partir de la subrutina del reescalat de velocitats (VELO_RESCALING), fer - li el reescalat a una T.
    
    !INPUTS: Matriu de velocitats(v), temperatura(T)
    
    !OUTPUTS: Matriu de velocitats(v)

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
     
     !Fem el reescalat a la temperatura marcada
     CALL VELO_RESCALING(v,T)
     RETURN
END SUBROUTINE UNIFORM_VELO

END MODULE Distribucio_Uniforme_vel
