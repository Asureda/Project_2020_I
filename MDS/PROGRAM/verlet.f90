!GRUP I: Àlex, Oriol, Laia, Sílvia i Elena

MODULE Verlet_Algorithm

use READ_DATA
use Interaction_Cutoff_Modul
use PBC

implicit none

contains

SUBROUTINE VELO_VERLET(r,v,F)

    ! Primer, farem una petita descripció sobre els passos que seguirem per tractar aquest algorisme:
    ! Integració de les equacions de Newton mitjançant l'integració de Verlet
    ! El modul consta de dos parts, en primer lloc el càlcul de les velocitats i posicions a partir de les configuracions inicials
    ! En segon lloc calcular les noves forces a partir de les posicions anteriors i tot seguit fer el càlcul de les noves velocitats v(t+h) a partir de les
    ! noves forces.

    !OBJECTIU: Apliquem l'algorisme de Verlet, explicat anteriorment.
    
    !INPUTS: matriu de posicions(r), matriu de velocitats(v) i matriu de forces(F).
    !        Les matrius inicials de posicions, velocitats i forces també les utilitzarem en aquesta subrutina. 
    !        També, cridem a la subrutina que inlcourà el cutoff i ens definirà la PBC.
    
    !OUTPUTS: matriu de posicions(r), matriu de velocitats(v) i l'energia cinètica.

    INTEGER i
    REAL*8 r(:,:),v(:,:),r0(n_particles,3),v0(n_particles,3),f0(n_particles,3)
    REAL*8 F(:,:),cutoff
    cutoff=0.99*L*5d-1
    r0=r
    v0=v
    f0=f
    CALL INTERACTION_CUTOFF(r,F0,cutoff)
    DO i=1,n_particles
        r(i,:)=r0(i,:)+v0(i,:)*h+5d-1*F0(i,:)*h*h
        r(i,1)=PBC2(r(i,1),L)
        r(i,2)=PBC2(r(i,2),L)
        r(i,3)=PBC2(r(i,3),L)
    END DO
    CALL INTERACTION_CUTOFF(r,F,cutoff)
    kinetic=0d0
    DO i=1,n_particles
        v(i,:)=v0(i,:)+5d-1*(F(i,:)+F0(i,:))*h
        kinetic=kinetic+5d-1*(v(i,1)**2d0+v(i,2)**2d0+v(i,3)**2d0)
    END DO
    RETURN
END SUBROUTINE
END MODULE Verlet_Algorithm
