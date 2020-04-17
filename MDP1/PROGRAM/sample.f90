MODULE SAMPLE
use READ_DATA
use ALLOCATE_VARS
use Distribucio_Radial
IMPLICIT NONE
INTEGER kk
contains
SUBROUTINE SAMPLES()
  if(taskid==0) then
  if((mod(nstep,n_meas).eq.0).and.(is_print_thermo.eqv..true.))then
    temp_instant=2d0*kinetic/(3d0*n_particles)
    pressure=(density*temp_instant+pressure/(2d0*3d0*L**3d0))
    !kinetic = kinetic/n_particles
    !potential = potential/n_particles
          write(41,*) t,kinetic/(1d0*n_particles),potential/(2d0*n_particles)
          write(42,*) t,kinetic/(1d0*n_particles)+potential/(2d0*n_particles)
          write(43,*) t,temp_instant,pressure
          write(44,*) t*time_re,kinetic/(1d0*n_particles)*energy_re,potential/(2d0*n_particles)*energy_re
          write(45,*) t*time_re,(kinetic/(1d0*n_particles)+potential/(2d0*n_particles))*energy_re
          write(46,*) t*time_re,temp_instant*temp_re,pressure*press_re
  endif
end if

if (taskid==0) then
  if((mod(nstep,n_meas_gr).eq.0).and.(is_compute_gr.eqv..true.))then
    call RAD_DIST_INTER(r,g_r) !càlcul g(r) a cada pas
    n_gr_meas=n_gr_meas+1
  endif
  IF(is_time_evol.eqv..TRUE.)THEN
    WRITE(54,*)n_particles
    WRITE(54,*)
    DO kk=1,n_particles
      WRITE(54,*)'X',r(kk,:)
    END DO
  END IF
end if
END SUBROUTINE SAMPLES

SUBROUTINE gdr()
  if (taskid==0) then
  if((is_compute_gr.eqv..true.))then
  call RAD_DIST_INTER(r,g_r)
  n_gr_meas=n_gr_meas+1
  call RAD_DIST_FINAL(g_r,n_gr_meas) !càlcul g(r) com a cúmul
  do kk=1,n_radial
    write(53,*)dx_radial*kk,g_r(kk)
  enddo
  endif
end if
END SUBROUTINE gdr

END MODULE SAMPLE
