! Modul per la distribucio uniforme de les velocitats
MODULE Distribucio_Uniforme_vel
use READ_DATA
use Reescala_velocitats
implicit none
contains
SUBROUTINE UNIFORM_VELO(v,T)
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
     
     !Resacling the velocities to the temperature
     CALL VELO_RESCALING(v,T)
     RETURN
END SUBROUTINE UNIFORM_VELO

END MODULE Distribucio_Uniforme_vel