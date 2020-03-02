module Distribucio_Radial

use READ_DATA
!Cridem el mòdul PBC definit  a l'arxiu PBC
use PBC

IMPLICIT NONE

contains

    subroutine RAD_DIST_INTER(r,vec)
    !Aquesta subrutina, ens permetrà calcular a partir dels diferencials de distància, el valor de la distribució radial
    IMPLICIT NONE
    INTEGER i,coef,j
    REAL*8 dist,vec(:,:), r(:,:),dx,dy,dz
    DO i=1,n_particles
        DO j=1,n_particles
            IF (i.ne.j) THEN
                !Calculem el diferencial en l'espai de les tres dimensions sobre les que treballem
                dx=PBC1(r(j,1)-r(i,1),L)
                dy=PBC1(r(j,2)-r(i,2),L)
                dz=PBC1(r(j,3)-r(i,3),L)
                dist=(dx**2d0+dy**2d0+dz**2d0)**0.5
                coef=int(0.5+dist/dx_radial)
                IF ((coef.gt.0).and.(coef.le.n_radial)) THEN
                    vec(coef)=vec(coef)+1.0
                END IF
            END IF
        END DO
    END DO
    RETURN
    end subroutine RAD_DIST_INTER

    subroutine RAD_DIST_FINAL(vec,n_gr_meas)
    !Aquesta subrutina ens permet 
    IMPLICIT NONE
    INTEGER i,n_gr_meas
    REAL*8 vec(0:n_radial+1),result(0:n_radial+1),aux
    vec=vec/(1d0*n_gr_meas)
    DO i=2,n_radial
        !Fem aquest bucle, des de 2 fins a n_radial per analitzar cada 
        aux=(density*4d0*3.1415*((((i)*dx_radial)**3d0)-(((i-1)*dx_radial)**3d0)))/3d0
        result(i)=vec(i)/aux
    END DO
    vec=result 
    RETURN
    end subroutine RAD_DIST_FINAL

end module Distribucio_Radial