MODULE READ_DATA
    IMPLICIT NONE
    !Variables from parameters.dat
    INTEGER :: n_particles
    REAL*8 :: density,t_b,h,sigma,epsilon,mass,T_ini, T_therm, dx_radial
    LOGICAL :: is_thermostat

    !Variables from config.dat
    INTEGER :: n_meas,n_meas_gr,n_meas_time_ev,n_melting
    REAL*8 :: T_therm_prov
    LOGICAL :: is_print_thermo,is_compute_gr,is_time_evol

    !Variables from constants.dat
    REAL*8 :: k_b,n_avog

    !Additional variables
    INTEGER :: M,n_radial
    REAL*8 :: L,a,T_a,kinetic,potential,pressure
    !LOGICAL :: is_thermostat
    CONTAINS
    SUBROUTINE READ_ALL_DATA()
        IMPLICIT NONE
        !---------------------------------------------------
        !      READING PARAMETERS FILE
        !---------------------------------------------------
        OPEN(11,FILE='parameters.dat',status='OLD')
        READ(11,*)n_particles
        READ(11,*)density
        READ(11,*)t_b
        READ(11,*)h
        READ(11,*)sigma
        READ(11,*)epsilon
        READ(11,*)mass
        READ(11,*)T_ini
        READ(11,*)is_thermostat
        READ(11,*)T_therm
        READ(11,*)dx_radial
        CLOSE(11)
        !--------------------------------------------------
        ! READING CONFIGURATION FILE
        !--------------------------------------------------
        OPEN(12,FILE='config.dat',status='OLD')
        READ(12,*)T_therm_prov
        READ(12,*)n_melting
        READ(12,*)is_print_thermo
        READ(12,*)n_meas
        READ(12,*)is_compute_gr
        READ(12,*)n_meas_gr
        READ(12,*)is_time_evol
        READ(12,*)n_meas_time_ev
        CLOSE(12)

        OPEN(13,FILE='constants.dat',status='OLD')
        READ(13,*)k_b
        READ(13,*)n_avog
        CLOSE(13)
    END SUBROUTINE
    SUBROUTINE OTHER_GLOBAL_VARS()
        IMPLICIT NONE
        L=((n_particles*1d0)/density)**(1d0/3d0)     !Computing the magnitude of the system
        M=nint(((n_particles*1d0)/4d0)**(1d0/3d0))   !Number of cells
        a=L/(M*1d0)                                  !dimension of the cells
        n_radial=int(L/dx_radial)
        ta=0d0
        kinetik=0d0
        potential=0d0
        pressure=0d0
    END SUBROUTINE
    FUNCTION KINETIC_ENERGY(v)
        IMPLICIT NONE
        INTEGER i
        REAL*8 :: v(:,:),KINETIC_ENERGY
        KINETIC_ENERGY=0d0
        DO i=1,n_particles
            KINETIC_ENERGY=KINETIC_ENERGY+5d-1*(v(i,1)**2d0+v(i,2)**2d0+v(i,3)**2d0)
        END DO
        RETURN
    END FUNCTION KINETIK

END MODULE

MODULE ALLOCATE_VARS
    USE READ_DATA
    IMPLICIT NONE
    INTEGER seed
    REAL*8 temp_re,energy_re,dist_re,temp_re,press_re,n_mols,total_mass,rho_re
    REAL*8,DIMENSION(:,:),ALLOCATABLE :: positions, velocities, force
    REAL*8,DIMENSION(:),ALLOCATABLE :: g_r

    CONTAINS
    SUBROUTINE INITIALIZE_VARS()
        IMPLICIT NONE

        ALLOCATE(r(n_particles,3),v(n_particles,3),F(n_particles,3),g_r(0:n_radial+1))
        seed=1996

        
        !dimensional factors
        temp_re=epsilon/k_B   !Kelvin
        energy_re=epsilon     !kJ/mol
        dist_re=sigma         !Angstroms
        time_re=1d-13*sqrt(mass*sigma**2d0/epsilon)  !Seconds
        press_re=1d27*epsilon/(n_avog*sigma**3d0)          !Pascals
        n_mols=n_particles/n_avog
        total_mass=n_mols*mass
        rho_re=mass*1d24/(sigma**3d0*n_avog)

    END SUBROUTINE
END MODULE

