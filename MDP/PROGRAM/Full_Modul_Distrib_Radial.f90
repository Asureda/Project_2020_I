module Distribucio_Radial

use READ_DATA
!Cridem el mòdul PBC definit  a l'arxiu PBC
use PBC
use parallel_routines

IMPLICIT NONE

contains

    subroutine RAD_DIST_INTER(r,vec)
    !Aquesta subrutina, ens permetrà calcular a partir dels diferencials de distància, el valor de la distribució radial
    IMPLICIT NONE
    INTEGER i,coef,j
    REAL*8 dist,vec(:), r(:,:),dx,dy,dz
    !DO i=1,n_particles
        !DO j=1,n_particles
        !taskid= identificador del processador
        !IF (paral_double.eqv..TRUE.) THEN
        IF (taskid.le.nworking) THEN
            DO i=simple_matrix(taskid,1),simple_matrix(taskid,2)
                 DO j = i+1, n_particles
                    IF (i.ne.j) THEN
                        !Calculem el diferencial en l'espai de les tres dimensions sobre les que treballem
                        dx=PBC1(r(j,1)-r(i,1),L)
                        dy=PBC1(r(j,2)-r(i,2),L)
                        dz=PBC1(r(j,3)-r(i,3),L)
                        dist=(dx**2d0+dy**2d0+dz**2d0)**0.5
                        coef=int(0.5+dist/dx_radial)
                        IF ((coef.gt.0).and.(coef.le.n_radial)) THEN
                            vec(coef)=vec(coef)+2.0
                        END IF
                    END IF
                END DO
            END DO
        END IF
        ELSE
            DO i=1,n_particles
                DO j=i+1,n_particles
                    IF (i.ne.j) THEN
                        !Calculem el diferencial en l'espai de les tres dimensions sobre les que treballem
                        dx=PBC1(r(j,1)-r(i,1),L)
                        dy=PBC1(r(j,2)-r(i,2),L)
                        dz=PBC1(r(j,3)-r(i,3),L)
                        dist=(dx**2d0+dy**2d0+dz**2d0)**0.5
                        coef=int(0.5+dist/dx_radial)
                        IF ((coef.gt.0).and.(coef.le.n_radial)) THEN
                            vec(coef)=vec(coef)+2.0
                        END IF
                    END IF
                END DO
            END DO
        END IF
        !MPI_BARRIER(comm , ierror)

    end subroutine RAD_DIST_INTER

    subroutine RAD_DIST_FINAL(vec,n_gr_meas)
    !Per tal de fer un càlcul correcte, s'han de fer un cúmul de g(r)
    !Aquesta subrutina ens permet calcular el promig dels valors anteriors
    IMPLICIT NONE
    INTEGER i,n_gr_meas
    REAL*8 vec(:),result(0:n_radial+1),aux
    vec=vec/(1d0*n_gr_meas*n_particles)
    result=0d0
    DO i=2,n_radial
        aux=(density*4d0*3.1415*((((i)*dx_radial)**3d0)-(((i-1)*dx_radial)**3d0)))/3d0
        result(i)=vec(i)/aux
    END DO
    vec=0d0
    vec=result

    end subroutine RAD_DIST_FINAL

end module Distribucio_Radial