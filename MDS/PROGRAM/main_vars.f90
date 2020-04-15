!GRUP I: Àlex, Oriol, Laia, Sílvia i Elena

MODULE ALLOCATE_VARS

    USE READ_DATA
    
    IMPLICIT NONE
    
    INTEGER seed,n_verlet,n_gr_meas,i
    REAL*8 time_re,energy_re,dist_re,temp_re,press_re,n_mols,total_mass,rho_re,t,temp_instant,cutoff_aux
    REAL*8,DIMENSION(:,:),ALLOCATABLE :: r, v, f
    REAL*8,DIMENSION(:),ALLOCATABLE :: g_r

    CONTAINS
    
    SUBROUTINE INITIALIZE_VARS()
    
        !OBJECTIU: Alocatar i donar valors a les variables necessàries per fer la simulació
        !          També, fixarem alguns paràmetres, per tal de poder fer els gràfics amb valors reals.
        
        !INPUTS: matriu de posicions(r), matriu de velocitats(v), matriu de forces(F), vector de
        !        la distribució radial(g_r) i els factors de reescalat.
        
        !OUTPUTS: matriu de posicions(r), matriu de velocitats(v), matriu de forces(F), vector de
        !         la distribució radial(g_r) i els factors de reescalat.
    
        IMPLICIT NONE

        ALLOCATE(r(n_particles,3),v(n_particles,3),F(n_particles,3),g_r(0:n_radial+1))
        seed=1996

        
        !Factors de Reescalat
        temp_re=epsilon/k_B   !Kelvin
        energy_re=epsilon     !kJ/mol
        dist_re=sigma         !Angstroms
        time_re=0.1*sqrt(mass*sigma**2d0/epsilon)  !Picosegons
        press_re=1d33*epsilon/(n_avog*sigma**3d0)  !Pascals
        n_mols=n_particles/n_avog
        total_mass=n_mols*mass
        rho_re=mass*1d24/(sigma**3d0*n_avog)

        !Valors d'integració de l'algorisme de Verlet
        n_verlet=int((t_b-t_a)/h)  !Time Steps
        t=t_a
        g_r=0d0
        n_gr_meas=0

    END SUBROUTINE
END MODULE
