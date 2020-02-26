module Reescala_velocitats

use READ_DATA

IMPLICIT NONE

contains

    subroutine Velo_Rescaling(v,T)
    IMPLICIT NONE
    REAL*8 v(:,:),T,alpha
    alpha=sqrt(3d0*n_particles*T/(2d0*KINETIK))
    v=alpha*v
    end subroutine Velo_Rescaling

end module Reescala_velocitats