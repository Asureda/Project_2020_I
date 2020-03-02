PROGRAM SEQUENTIAL_MD
  use READ_DATA
  use ALLOCATE_VARS
  use Inicialitzar
  use Distribucio_Uniforme_vel
  use Verlet_Algorithm
  use Andersen_modul
  use Distribucio_Radial
  use Reescala_velocitats

  IMPLICIT NONE

  call srand(seed)
  call read_all_data()
  call other_global_vars()
  call INITIALIZE_VARS()
  call FCC_Initialize(r)
  call Uniform_velocity(v,T_ini)
  call VELO_RESCALING_MOD(v,T_therm_prov)
  print*,'inicial melting'
  DO i=1,n_melting
    call velo_verlet(r,v,F)
    call andersen(v,T_therm_prov)
  end do
  print*,'FINAL MELTING'
  call Velo_Rescaling(v,T_ini)
  open(51,file='thermodynamics_reduced.dat')
  open(52,file='thermodynamics_real.dat')
  open(53,file='distriv_funct.dat')
  open(54,file='positions.xyz')
  DO i=1,n_verlet
    t=t_a+i*h
    call VELO_VERLET(r,v,F)
    if(is_thermostat.eqv..true.)THEN
      call andersen(v,T_therm)
    end if
    if((mod(i,n_meas).eq.0).and.(is_print_thermo.eqv..true.))then
      temp_instant=2d0*kinetic/(3d0*n_particles)
      pressure=(density*temp_instant*pressure/(3d0*L**3d0))
      write(51,*)t,kinetic,potential,(kinetic+potential),temp_instant,pressure
      write(52,*)t*time_re,kinetic*energy_re,potential*energy_re,(kinetic+potential)*energy_re,temp_instant*&
                                                                                    &temp_re,pressure*press_re
    endif
    ! el seguent if sera x fer un fitxer xyz pel vmd :)
    !if()THEN
    !endif
    if((mod(i,n_meas_gr).eq.0).and.(is_compute_gr.eqv..true.))then
      call RAD_DIST_INTER(r,g_r)
      n_gr_meas=n_gr_meas+1
    endif
  enddo
  if((is_compute_gr.eqv..true.))then
    call RAD_DIST_INTER(r,g_r)
    n_gr_meas=n_gr_meas+1
    call RAD_DIST_FINAL(g_r,n_gr_meas)
    do i=1,n_radial
      write(53,*)dx_radial*i,g_r(i)
    enddo
  endif
  !if()THEN
    !DO i=1,n_particules
      !write(54,*)'C',r(i,:)
    !enddo
  !endif
  print*,'PROGRAM END'









END PROGRAM SEQUENTIAL_MD
