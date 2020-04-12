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
  IMPLICIT NONE
  INTEGER k
  call srand(seed)
  call read_all_data()
  call other_global_vars()
  call INITIALIZE_VARS()


  call FCC_Initialize(r)
  call Uniform_velocity(v)
  call VELO_RESCALING_MOD(v,T_therm_prov)

   DO i=1,n_melting
    call velo_verlet(r,v,F)
    !if(is_thermostat.eqv..true.)then
     call andersen(v,T_therm_prov)
    !end if
   END DO


  !call VELO_RESCALING_MOD(v,T_ini) -> HO DEIXO AQUI PER SI A CAS
  open(51,file='thermodynamics_reduced.dat')
  open(52,file='thermodynamics_real.dat')
  open(53,file='distrib_funct.dat')
  open(54,file='positions.xyz')

  pressure=(density*temp_instant+pressure/(3d0*L**3d0))

  DO i=1,n_verlet
    t=t_a+i*h
    call VELO_VERLET(r,v,F)
    if(is_thermostat.eqv..true.)THEN
      call andersen(v,T_therm)
    end if

    if((mod(i,n_meas).eq.0).and.(is_print_thermo.eqv..true.))then
      temp_instant=2d0*kinetic/(3d0*n_particles)
      pressure=(density*temp_instant+pressure/(3d0*L**3d0))
      kinetic = kinetic/n_particles
      potential = potential/n_particles
      write(51,*)t,kinetic,potential,(kinetic+potential),temp_instant,pressure
      write(52,*)t*time_re,kinetic*energy_re,potential*energy_re,(kinetic+potential)*energy_re,temp_instant*&
                                                                                       &temp_re,pressure*press_re
      !print*,'potential',potential
      !print*,'kinetic',kinetic
      !print*,'energy',potential + kinetic
      !print*,'Temperatura red, real', temp_instant, temp_instant*temp_re
    endif

    !-----------------------------------------------------------------------------
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
