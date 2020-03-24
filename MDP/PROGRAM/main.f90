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
  use parallel_routines

  IMPLICIT NONE
  INTEGER k, master

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

  master=0
  if(taskid==master)then
    call srand(seed)
    !LLEGIM EL FITXER INPUT AMB LES SEGÜENTS DADES:
    !PARAMETRES DE DENSITAT, MASSA, TEMPERATURA DE REFERÈNCIA, TEMPERATURA DEL BANYS ETC.
    call read_all_data()
    !CALCULEM ELS PARÀMETRES GLOBALS DEL LATTICE: NÚMERO DE PARTÍCULES, LONGITUD DE LA CAIXA,
    !DISTANCIA ENTRE PARTICULES ETC.
    call other_global_vars()
    !PARALLELIZATION MPI SUBROUTINES IN ORDER TO DISTRIBUTE THE PARTICLES AMONG THE PROCESSORS
    call simple_loop_matrix()
    DO k=1,numproc
      print*,index_matrix(k,:)
    ENDDO
    call double_loop_matrix()
    DO k=1,numproc
      print*,double_matrix(k,:)
    ENDDO
  endif
  stop

  !INICITALITZEM LES VARIABLES D'ESTAT EN UNITATS REDUÏDES
  call INITIALIZE_VARS()
  !DEFINIM LA CONFIGURACIÓ INICIAL DE LES PARTICULES COM UNA XARXA FCC
  call FCC_Initialize(r)
  !LI DONEM UNA VELOCITAT INICIAL A LES PARTICULES (VELOCITATS INICIALS RANDOM)
  call Uniform_velocity(v,T_ini)
  !FEM UN REESCALATGE DE LES VELOCITATS A LA TEMPERATURA INICIAL
  !LI DONEM UNA TEMPERATURA INICIAL SUFICIENTMENT GRAN COM PER DESFER LA ESTRUCTURA
  !CRISTALINA (MELTING)
  call VELO_RESCALING_MOD(v,T_therm_prov)
  !UN COP CALCULAT EL NÚMERO DE ITERACIONS NECESSARIES PER FONDRE EL SÒLID INICIAL
  !APLIQUEM EL TERMOSTAT DE ANDERSEN TANTS COPS COM SIGUIN NECESSARIS
  !cutoff_aux=0.99*L*5d-1
  !CALL INTERACTION_CUTOFF(r,F,cutoff_aux)
  DO i=1,n_melting
    call velo_verlet(r,v,F) !EN UNA REGIÓ LxL AMB UNES CONDICIONS DE CONTORN PERIODIQUES
                            ! EN FUNCIO DE LES FORCES D'INTERACCIÓ S'ACTUALITZEN LES VELOCITATS
                            ! I LES POSICIONS DE LES PARTÍCULES
    call andersen(v,T_therm_prov) !AMB EL TERMOSTAT RECALCULEM LES VELOCITATS ARA EN FUNCIO
                                  ! DE LES TEMPERATURES
  end do
  print*,'FINAL MELTING'
  !AMB EL SÒLID FOS I LES PARTICULES MOVENT-SE COM UN FLUID LES VELOCITATS ES REESCALEN CALCULANT
  !L'ENERGIA CINÈTICA DEGUDA A LA TEMPERATURA DE LES PARTÍCULES
  !COPIEM ELS PRIMERS RESULTATS DE LES PARTICULES COM A FLUID, VELOCITAT, POSICIONS, TEMPERATURES I
  !PRESSIÓ, EN UNITATS REDUÏDES I NO REDUÏDES I LES POSICIONS DE LES PARTÍCULES
  !I LES ESCRIBIIM EN UN FITXER OUTPUT
  call Velo_Rescaling(v,T_ini)
  open(51,file='thermodynamics_reduced.dat')
  open(52,file='thermodynamics_real.dat')
  open(53,file='distrib_funct.dat')
  open(54,file='positions.xyz')
  !APLIQUEM L'ALGORITME DE VERLET I EL TERMOSTAT D'ANDERSEN PER OBTENIR
  !VELOCITAT, POSICIONS, TEMPERATURES I
  !PRESSIÓ, EN UNITATS REDUÏDES I NO REDUÏDES, I LES POSICIONS DE LES PARTÍCULES
  !I LES ESCRIBIM EN UN FITXER OUTPUT, PER N TIME STEPS D'UN INTERVAL DE TEMPS
  !cutoff_aux=0.99*L*5d-1
  !CALL INTERACTION_CUTOFF(r,F,cutoff_aux)
  pressure=(density*temp_instant+pressure/(3d0*L**3d0))
  print*,'pres',press_re,pressure,pressure*press_re
  DO i=1,n_verlet
    t=t_a+i*h
    call VELO_VERLET(r,v,F)
    if(is_thermostat.eqv..true.)THEN
      call andersen(v,T_therm)
    end if
  !PER OBTENIR LA DISTRIBUCIÓ RADIAL DE LES PARTÍCULES A CADA TIME STEP
  !DE LES PARTÍCULES DE LA REGIÓ DE LA CAIXA LI APLIQUEM LA FUNCIÍ G EN FUNCIÓ
  !DEL RADI
    if((mod(i,n_meas).eq.0).and.(is_print_thermo.eqv..true.))then
      temp_instant=2d0*kinetic/(3d0*n_particles)
      pressure=(density*temp_instant+pressure/(3d0*L**3d0))
      write(51,*)t,kinetic,potential,(kinetic+potential),temp_instant,pressure
      write(52,*)t*time_re,kinetic*energy_re,potential*energy_re,(kinetic+potential)*energy_re,temp_instant*&
                                                                                    &temp_re,pressure*press_re
    endif
    if((mod(i,n_meas_gr).eq.0).and.(is_compute_gr.eqv..true.))then
      call RAD_DIST_INTER(r,g_r) !càlcul g(r) a cada pas
      n_gr_meas=n_gr_meas+1
    endif
    IF(is_time_evol.eqv..TRUE.)THEN
      WRITE(54,*)n_particles
      WRITE(54,*)
      DO k=1,n_particles
        WRITE(54,*)'X',r(k,:)
      END DO
    END IF
  enddo
  if((is_compute_gr.eqv..true.))then
    call RAD_DIST_INTER(r,g_r)
    n_gr_meas=n_gr_meas+1
    call RAD_DIST_FINAL(g_r,n_gr_meas) !càlcul g(r) com a cúmul
    do k=1,n_radial
      write(53,*)dx_radial*k,g_r(k)
    enddo
  endif
  print*,'PROGRAM END'
END PROGRAM SEQUENTIAL_MD
