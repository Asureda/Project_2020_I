PROGRAM SEQUENTIAL_MD
  use READ_DATA
  use ALLOCATE_VARS
  use Inicialitzar
  use Distribucio_Uniforme_vel
  use Interaction_Cutoff_Modul
  use Verlet_Algorithm
  use Andersen_modul
  use Distribucio_Radial
  use Reescala_velocitats
  use SAMPLE
  IMPLICIT NONE
  call srand(seed)
  call read_all_data()
  call other_global_vars()
  call INITIALIZE_VARS()


  call FCC_Initialize(r)
  call Uniform_velocity(v)
  call VELO_RESCALING_MOD(v,T_therm_prov)

   DO i=1,n_melting
    call velo_verlet(r,v,F)
     call andersen(v,T_therm_prov)
   END DO
  !call VELO_RESCALING_MOD(v,T_ini) -> HO DEIXO AQUI PER SI A CAS
  open(51,file='thermodynamics_reduced.dat')
  open(52,file='thermodynamics_real.dat')
  open(53,file='distrib_funct.dat')
  open(54,file='positions.xyz')

  DO i=1,n_verlet
    t=t_a+i*h
    call VELO_VERLET(r,v,F)
    if(is_thermostat.eqv..true.)THEN
      call andersen(v,T_therm)
    end if
  CALL SAMPLES()
  end do
  CALL gdr()
  print*,'PROGRAM END'
END PROGRAM SEQUENTIAL_MD
