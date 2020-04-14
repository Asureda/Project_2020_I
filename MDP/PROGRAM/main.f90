PROGRAM SEQUENTIAL_MD
  !use MPI


  use READ_DATA
  use ALLOCATE_VARS
  use Inicialitzar
  use Interaction_Cutoff_Modul
  use Verlet_Algorithm
  use Andersen_modul
  use Distribucio_Radial
  use Reescala_velocitats
  use parallel_routines

  IMPLICIT NONE
  INTEGER k, master_task
  REAL*8 starttime, endtime
  !#############################################################!
  !            INITIALIZING MPI                                 !
  !#############################################################!
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
  starttime = MPI_WTIME()
  taskid=taskid+1
  master_task=1

!#############################################################!
!            READING AND SETTING VARIABLES                    !
!#############################################################!

  !call srand(seed)
  !LLEGIM EL FITXER INPUT AMB LES SEGÜENTS DADES:
  !PARAMETRES DE DENSITAT, MASSA, TEMPERATURA DE REFERÈNCIA, TEMPERATURA DEL BANYS ETC.
  call read_all_data()
  !CALCULEM ELS PARÀMETRES GLOBALS DEL LATTICE: NÚMERO DE PARTÍCULES, LONGITUD DE LA CAIXA,
  !DISTANCIA ENTRE PARTICULES ETC.
  call other_global_vars()
  !PARALLELIZATION MPI SUBROUTINES IN ORDER TO DISTRIBUTE THE PARTICLES AMONG THE PROCESSORS
  call simple_loop_matrix()
  !INICITALITZEM LES VARIABLES D'ESTAT EN UNITATS REDUÏDES
  call INITIALIZE_VARS()
  !DEFINIM LA CONFIGURACIÓ INICIAL DE LES PARTICULES COM UNA XARXA FCC
!#############################################################!
!            iNITIAL CONFIGURATION                            !
!#############################################################!
  if(taskid.eq.master_task)then
    call FCC_Initialize(r)
    !LI DONEM UNA VELOCITAT INICIAL A LES PARTICULES (VELOCITATS INICIALS RANDOM)
    call Uniform_velocity(v,T_therm_prov)
    !FEM UN REESCALATGE DE LES VELOCITATS A LA TEMPERATURA INICIAL
    !LI DONEM UNA TEMPERATURA INICIAL SUFICIENTMENT GRAN COM PER DESFER LA ESTRUCTURA
    !CRISTALINA (MELTING)
    !call VELO_RESCALING_MOD(v,T_therm_prov)
  endif
  DO k=1,3
        CALL MPI_BCAST(v(1:n_particles,k),n_particles,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierror)
        CALL MPI_BCAST(r(1:n_particles,k),n_particles,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierror)
    END DO
    f=0

  !#############################################################!
  !            MELTING THE FCC STRUCTURE                        !
  !#############################################################!
  DO i=1,n_melting
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call velo_verlet(r,v,F,i)
    IF(taskid.eq.1)THEN
    temp_instant=2d0*kinetic/(3d0*n_particles)
    pressure=(density*temp_instant+pressure/(2d0*3d0*L**3d0))
    !print*,'1',temp_instant,kinetic/(1d0*n_particles)
    !print*,'2',potential/(2d0*n_particles),pressure
    endif
    !print*,'kinetic from main',kinetic,taskid
                            !EN UNA REGIÓ LxL AMB UNES CONDICIONS DE CONTORN PERIODIQUES
                            ! EN FUNCIO DE LES FORCES D'INTERACCIÓ S'ACTUALITZEN LES VELOCITATS
                            ! I LES POSICIONS DE LES PARTÍCULES
    call andersen(v,T_therm_prov) !AMB EL TERMOSTAT RECALCULEM LES VELOCITATS ARA EN FUNCIO
                                  ! DE LES TEMPERATURES
  end do
  !AMB EL SÒLID FOS I LES PARTICULES MOVENT-SE COM UN FLUID LES VELOCITATS ES REESCALEN CALCULANT
  !L'ENERGIA CINÈTICA DEGUDA A LA TEMPERATURA DE LES PARTÍCULES
  !COPIEM ELS PRIMERS RESULTATS DE LES PARTICULES COM A FLUID, VELOCITAT, POSICIONS, TEMPERATURES I
  !PRESSIÓ, EN UNITATS REDUÏDES I NO REDUÏDES I LES POSICIONS DE LES PARTÍCULES
  !I LES ESCRIBIIM EN UN FITXER OUTPUT
  IF(taskid.eq.master_task)THEN
    open(41,file='kin_pot_red.dat')
    open(42,file='total_red.dat')
    open(43,file='temp_pres_red.dat')
    open(44,file='kin_pot_real.dat')
    open(45,file='total_real.dat')
    open(46,file='temp_pres_real.dat')
    open(51,file='thermodynamics_red.dat')
    open(52,file='thermodynamics_real.dat')
    open(53,file='distrib_funct.dat')
    open(54,file='positions.xyz')
    print*,'-----------------------------------------------'
  END IF
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  !   call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  ! print*,taskid,r(10,:)
  ! call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  ! print*,'----------------------'
  ! call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  ! print*,taskid,v(10,:)
  ! call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  ! print*,'----------------------'
  ! call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  ! print*,taskid,f(10,:)
  ! call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  ! print*,'----------------------'
  ! call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  ! stop
  !APLIQUEM L'ALGORITME DE VERLET I EL TERMOSTAT D'ANDERSEN PER OBTENIR
  !VELOCITAT, POSICIONS, TEMPERATURES I
  !PRESSIÓ, EN UNITATS REDUÏDES I NO REDUÏDES, I LES POSICIONS DE LES PARTÍCULES
  !I LES ESCRIBIM EN UN FITXER OUTPUT, PER N TIME STEPS D'UN INTERVAL DE TEMPS
  DO k=1,3
        CALL MPI_BCAST(v(1:n_particles,k),n_particles,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierror)
        CALL MPI_BCAST(r(1:n_particles,k),n_particles,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierror)
        CALL MPI_BCAST(F(1:n_particles,k),n_particles,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierror)
  END DO
  print*,n_verlet

!#############################################################!
!      SYMULATION : VELOCITY VERLET INTEGRATION               !
!#############################################################!
  DO i=1,n_verlet
    IF(taskid.eq.master_task) THEN
    t=t_a+i*h
    !print*,'verlet',i
    END IF
    CALL MPI_BARRIER(comm,ierror)

    !VELOCITY VERLET STEP
    call VELO_VERLET(r,v,F,i)
    !ANDERSEN THERMOSTAT
    if(is_thermostat.eqv..true.)THEN
      call andersen(v,T_therm)
    end if
  !PER OBTENIR LA DISTRIBUCIÓ RADIAL DE LES PARTÍCULES A CADA TIME STEP
  !DE LES PARTÍCULES DE LA REGIÓ DE LA CAIXA LI APLIQUEM LA FUNCIÍ G EN FUNCIÓ
  !DEL RADI
  !#############################################################!
  !        SAVING DATA TO FILES                                 !
  !#############################################################!
      IF (taskid.eq.master_task) THEN
        if((mod(i,n_meas).eq.0).and.(is_print_thermo.eqv..true.))then
          temp_instant=2d0*kinetic/(3d0*n_particles)
          pressure=(density*temp_instant+pressure/(2d0*3d0*L**3d0))
          !print*,'1',temp_instant,kinetic/(1d0*n_particles)
          !print*,'2',potential/(2d0*n_particles),pressure
          !print*,'pressió',pressure/2.,temp_instant
          !print*,'Kinetic, Potential',kinetic/(1d0*n_particles),potential/(2d0*n_particles),temp_instant
          !print*,'r',r(150,1),v(150,1),f(150,1)
          !print*,'Energia',kinetic/(1d0*n_particles)+potential/(2d0*n_particles)
!write(51,*)t,kinetic/(1d0*n_particles),potential/(2d0*n_particles),(kinetic/(1d0*n_particles)+potential/(2d0*n_particles)),temp_instant,pressure
!write(52,*)t*time_re,kinetic*energy_re,potential*energy_re,(kinetic+potential)*energy_re,temp_instant*temp_re,pressure*press_re
          write(41,*) t,kinetic/(1d0*n_particles),potential/(2d0*n_particles)
          write(42,*) t,kinetic/(1d0*n_particles)+potential/(2d0*n_particles)
          write(43,*) t,temp_instant,pressure
          write(44,*) t,kinetic/(1d0*n_particles)*energy_re,potential/(2d0*n_particles)*energy_re
          write(45,*) t,(kinetic/(1d0*n_particles)+potential/(2d0*n_particles))*energy_re
          write(46,*) t,temp_instant*temp_re,pressure*press_re
!print*,t,kinetic,potential,kinetic+potential,temp_instant,pressure
!print*,t,kinetic,potential,(kinetic+potential),temp_instant,pressure
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
      END IF
  !------------------------------------
  enddo   ! FINAL VERLET
  !------------------------------------
  
  !COPMUTING THE FINAL RADIAL DISTRIBUTION FUNCTION
  IF (taskid.eq.master_task)then
    if((is_compute_gr.eqv..true.))then
      call RAD_DIST_INTER(r,g_r)
      n_gr_meas=n_gr_meas+1
      call RAD_DIST_FINAL(g_r,n_gr_meas) !càlcul g(r) com a cúmul
      do k=1,n_radial
        write(53,*)dx_radial*k,g_r(k)
      enddo
    endif    
  END IF

  !#############################################################!
  !        CLOSING MPI AND PROGRAM                              !
  !#############################################################!
  endtime = MPI_WTIME()
  IF (taskid.eq.master_task) THEN
    print*,'time = ',endtime-starttime
    print*,'PROGRAM END'
    !print*,f
  END IF
  CALL MPI_FINALIZE(ierror)
END PROGRAM SEQUENTIAL_MD
