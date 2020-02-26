! Moduls de Lennard Jones, Periodic Bounday Conditions, Interaction Cutoff, velocitats distribuïdes uniformement, Velocity Verlet


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


! Modul del cutoff
MODULE Interaction_Cutoff_Modul
use READ_DATA
use Lennard_Jones
use PBC
implicit none
contains
SUBROUTINE INTERACTION_CUTOFF(r,F,cutoff)
!COMPUTING THE TOTAL INTERACTION ENERGY GIVEN AN EXTENAL FUNCTION FOR POTENCIAL
!AND BOUNDARY CONDITIONS
    IMPLICIT NONE
    INTEGER :: i,j
    REAL*8 :: cutoff,pot
    REAL*8 :: dx,dy,dz,d,ff
    REAL*8, DIMENSION(:,:) :: r, F
    !print*,'in inter'
    F=0d0
    potential=0d0
    !print*,'hola-1'
    pressure=0.0
    !print*,'hola'
    DO i=1,n_particles
        DO j=i+1,n_particles
            dx=PBC1(r(i,1)-r(j,1),L)
            dy=PBC1(r(i,2)-r(j,2),L)
            dz=PBC1(r(i,3)-r(j,3),L)
            d=(dx**2d0+dy**2d0+dz**2d0)**0.5
            
            CALL L_J(d,ff,pot,cutoff)
            F(i,1)=F(i,1)+ff*dx
            F(i,2)=F(i,2)+ff*dy
            F(i,3)=F(i,3)+ff*dz
            F(j,1)=F(j,1)-ff*dx
            F(j,2)=F(j,2)-ff*dy
            F(j,3)=F(j,3)-ff*dz
            potential=potential+pot
            !print*,'hola1'
            pressure=pressure+(ff*dx**2d0+ff*dy**2d0+ff*dz**2d0)
            !print*,'hola2'
        END DO
    END DO
    !print*,'en inter'
    RETURN
END SUBROUTINE INTERACTION_CUTOFF

END MODULE Interaction_Cutoff_Modul


! Modul del Lennard-Jones
MODULE Lennard_Jones
use READ_DATA
implicit none
contains
SUBROUTINE L_J(d,ff,pot,cutoff)
    IMPLICIT NONE
    REAL*8 d,cutoff,ff,pot
    ff=(48d0/d**14d0)-(24d0/d**8d0)
    pot=0d0
    IF (d<cutoff) THEN
        pot=4d0*((1d0/d**12d0)-(1d0/d**6d0))
    END IF
    RETURN
END SUBROUTINE L_J

END MODULE Lennard_Jones


! Modul de les condicions periodiques de contorn
MODULE PBC
use READ_DATA
implicit none
contains
FUNCTION PBC1(x,L)
!FUNCTION THAT RETURNS THE DISTANCE GIVEN PBC AND L
    IMPLICIT NONE
    REAL*8 x,L,PBC1
    PBC1=x-int(2d0*x/L)*L
    RETURN
END FUNCTION
FUNCTION PBC2(x,L)
!FUNCTION THAT RETURNS THE DISTANCE GIVEN PBC AND L
    IMPLICIT NONE
    REAL*8 x,L,PBC2
    IF(x.lt.0)THEN
        x=L+mod(x,L)
    ELSE IF(x.gt.L) THEN
        x=mod(x,L)
    END IF
    PBC2=x
    RETURN
END FUNCTION
END MODULE PBC

! Modul de l'algoritme d'integració Verlet
MODULE Verlet_Algorithm
use READ_DATA
use Interaction_Cutoff_Modul
use PBC
implicit none
contains
SUBROUTINE VELO_VERLET(r,v,F)
    INTEGER i
    REAL*8 r(:,:),v(:,:)
    REAL*8 F(:,:),cutoff
    cutoff=0.99*L*5d-1
    DO i=1,n_particles
        v(i,:)=v(i,:)+5d-1*F(i,:)*h
    END DO
    DO i=1,n_particles
        r(i,:)=r(i,:)+v(i,:)*h+5d-1*F(i,:)*h*h
        r(i,1)=PBC2(r(i,1),L)
        r(i,2)=PBC2(r(i,2),L)
        r(i,3)=PBC2(r(i,3),L)
    END DO
    CALL INTERACTION_CUTOFF(r,F,cutoff)
    kinetic=0d0
    DO i=1,n_particles
        v(i,:)=v(i,:)+5d-1*F(i,:)*h
        kinetic=kinetic+5d-1*(v(i,1)**2d0+v(i,2)**2d0+v(i,3)**2d0)
    END DO
    !print*,'out verlet'
    RETURN
END SUBROUTINE
END MODULE Verlet_Algorithm

! MODULE Kinetic_Energy
!     use READ_DATA
!     FUNCTION KINETIK(velocity)
!     IMPLICIT NONE
!     INTEGER n_particles,M,i
!     REAL*8 :: density,L,a
!     REAL*8 :: velocity(n_particles,3),KINETIK
!     COMMON/PARAMETERS/n_particles,M,density,L,a
!     KINETIK=0d0
!     DO i=1,n_particles
!         KINETIK=KINETIK
!     END DO
!     RETURN
! END FUNCTION KINETIK
! END MODULE Kinetic_Energy