MODULE READ_DATA
    IMPLICIT NONE
    !Variables from parameters.dat
    INTEGER :: n_particles
    REAL*8 :: density,t_b,h,sigma,epsilon,mass,T_ini, T_therm, dx_radial
    LOGICAL :: is_thermostat


    CONTAINS
    SUBROUTINE READ_ALL_DATA()
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
    END SUBROUTINE
END MODULE
