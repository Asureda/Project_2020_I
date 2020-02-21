!!Creem els mòduls de les subrutines, posicions i velocitats inicials, Andersen i funcions de distribució radial:

module Iniciatitzar

use READ_DATA

IMPLICIT NONE

contains

    subroutine FCC_Initialize(r)
    INTEGER :: n,i,j,k
    REAL*8 :: positions(:,:)
    n=1
    !print*,n
    DO i=0,M-1
        DO j=0,M-1
            DO k=0,M-1
                !print*,size(positions)
                positions(n,:)=a*(/i,j,k/)
                !print*,positions(n,:)
                positions(n+1,:)=positions(n,:)+a*(/0.5,0.5,0.0/)
                positions(n+2,:)=positions(n,:)+a*(/0.5,0.0,0.5/)
                positions(n+3,:)=positions(n,:)+a*(/0.0,0.5,0.5/)
                n=n+4
            END DO
        END DO
    END DO
    PRINT*, 'particles positioned', n-1, 'of a total imput', n_particles
    RETURN
    end subroutine FCC_Initialize

    subroutine Uniform_velocity(v,T)
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
    end subroutine Uniform_velocity

end module Inicialitzar

module Andersen_modul

use READ_DATA

IMPLICIT NONE

contains

    subroutine Andersen(v,temp)
    IMPLICIT NONE
    INTEGER i
    REAL*8 temp,nu,n1,n2,n3,n4,n5,n6
    REAL*8, DIMENSION(:,:) :: v
    nu=0.1/h
    sigma=sqrt(temp)
    DO i=1,n_particles
        IF(RAND().lt.nu*dt)THEN
            n1=RAND();n2=RAND()
            n3=RAND();n4=RAND()
            n5=RAND();n6=RAND()
            v(n_particles,:)=(/sigma*sqrt(-2d0*log10(n1))*cos(2d0*3.1415*n2),sigma*sqrt(-2d0*log10(n1))*sin(2d0*3.1415*n2),&
                &sigma*sqrt(-2d0*log10(n3))*cos(2d0*3.1415*n4)/)
        END IF
    END DO
    RETURN
    end subroutine Andersen

end module Andersen_modul

module Distribucio_Radial

use READ_DATA
use LJ_PBC

IMPLICIT NONE

contains

    subroutine RAD_DIST_INTER(r,vec)
    IMPLICIT NONE
    INTEGER i,coef,j
    REAL*8 dist,vec(:,:), r(:,:),dx,dy,dz
    DO i=1,n_particles
        DO j=1,n_particles
            IF (i.ne.j) THEN
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
    IMPLICIT NONE
    INTEGER i,n_gr_meas
    REAL*8 vec(0:n_radial+1),result(0:n_radial+1),aux
    vec=vec/(1d0*n_gr_meas)
    DO i=2,n_radial
        aux=(density*4d0*3.1415*((((i)*dx_radial)**3d0)-(((i-1)*dx_radial)**3d0)))/3d0
        result(i)=vec(i)/aux
    END DO
    vec=result 
    RETURN
    end subroutine RAD_DIST_FINAL

end module Distribucio_Radial

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
