PROGRAM SEQUENTIAL_MD
  use READ_DATA
  use ALLOCATE_VARS
  use Inicialitzar
  use Interaction_Cutoff_Modul
  use Verlet_Algorithm
  use Andersen_modul
  use Distribucio_Radial
  use Reescala_velocitats
  use parallel_routines
  use SAMPLE
  IMPLICIT NONE
  integer master, k
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
  master = 0

  call srand(seed)
  call read_all_data()
  call other_global_vars()
  call INITIALIZE_VARS()
  call simple_loop_matrix()


  call FCC_Initialize(r)
  call Uniform_velocity(v)
  DO k= 1,3
    call MPI_BCAST(r(:,k), n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
    call MPI_BCAST(v(:,k), n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
  END DO
  !call VELO_RESCALING_MOD(v,T_therm_prov)
  call VELO_RESCALING_MOD(v,T_therm)
  DO k= 1,3
    call MPI_BCAST(v(:,k), n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
  END DO

   DO i=1,10
    call velo_verlet(r,v,F)
     !call andersen(v,T_therm_prov)
   END DO
  !call VELO_RESCALING_MOD(v,T_ini) -> HO DEIXO AQUI PER SI A CAS
  if (taskid==0) then
  open(51,file='thermodynamics_reduced.dat')
  open(52,file='thermodynamics_real.dat')
  open(53,file='distrib_funct.dat')
  open(54,file='positions.xyz')
  end if

  DO i=1,n_verlet
    t=t_a+i*h
    call VELO_VERLET(r,v,F)
    !if(is_thermostat.eqv..true.)THEN
    !call andersen(v,T_therm)
    !end if

   if (taskid==0) then
   CALL SAMPLES()
   end if
  end do
   if (taskid==0) then
   CALL gdr()
   end if
  print*,'PROGRAM END'
  call MPI_FINALIZE(ierror)
END PROGRAM SEQUENTIAL_MD
