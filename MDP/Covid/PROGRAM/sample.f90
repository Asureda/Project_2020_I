MODULE SAMPLE
use READ_DATA
use ALLOCATE_VARS
use Distribucio_Radial
IMPLICIT NONE
INTEGER kk
contains
SUBROUTINE SAMPLES()
  if((mod(i,n_meas).eq.0).and.(is_print_thermo.eqv..true.))then
    temp_instant=2d0*kinetic/(3d0*n_particles)
    pressure=(density*temp_instant+pressure/(3d0*L**3d0))
    kinetic = kinetic/n_particles
    potential = potential/n_particles
    !print*, potential
    write(51,*)t,kinetic,potential,(kinetic+potential),temp_instant,pressure
    write(52,*)t*time_re,kinetic*energy_re,potential*energy_re,(kinetic+potential)*energy_re,temp_instant*&
                                                                                       &temp_re,pressure*press_re
  endif

  if((mod(i,n_meas_gr).eq.0).and.(is_compute_gr.eqv..true.))then
    call RAD_DIST_INTER(r,g_r) !càlcul g(r) a cada pas
    n_gr_meas=n_gr_meas+1
  endif
  IF(is_time_evol.eqv..TRUE.)THEN
  !  WRITE(54,*)n_particles
  !  WRITE(54,*)
    DO kk=1,n_particles
  !    WRITE(54,*)'X',r(kk,:)
    END DO
  END IF
END SUBROUTINE SAMPLES

SUBROUTINE gdr()
  if((is_compute_gr.eqv..true.))then
  call RAD_DIST_INTER(r,g_r)
  n_gr_meas=n_gr_meas+1
  call RAD_DIST_FINAL(g_r,n_gr_meas) !càlcul g(r) com a cúmul
  do kk=1,n_radial
  !  write(53,*)dx_radial*kk,g_r(kk)
  enddo
  endif
END SUBROUTINE gdr

END MODULE SAMPLE
