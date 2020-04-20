!GRUP I: Àlex, Oriol, Laia, Sílvia i Elena
module Reescala_velocitats

use READ_DATA

IMPLICIT NONE

contains

!Cridem la funció definida a l'arxiu input, per tal de definir l'energia cinètica:

    subroutine Velo_Rescaling(v,T)
    !OBJECTIU: Amb aquesta subrutina trobarem la velocitat, un cop hem reescalat amb el valor alpha,
    ! a una temperatura i un número de partícules concret.
    
    !INPUT: velocitats i temperatura
    
    !OUTPUTS: velocitats
    
    IMPLICIT NONE
    REAL*8 v(:,:),T,alpha
    alpha=sqrt(3d0*n_particles*T/(2d0*KINETIC))
    v=alpha*v
    end subroutine Velo_Rescaling

end module Reescala_velocitats
