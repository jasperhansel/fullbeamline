!+
! Subroutine time_runge_kutta_periodic_kick_hook (orbit, ele, param, stop_time, init_needed)
!
! Custom routine to add a kick to a particle at periodic times.
!
! For example, this routine could be used to add the kick due to a passing beam
! on a residual gas ion that is being tracked.
!
! The time_runge_kutta routine will initially call this routine at the beginning
! of tracking through an element with init_needed = True_int$
! This routine should then set stop_time to the time at which time_runge_kutta
! should again call this routine. 
! When this routine is called with init_needed /= True_int$, this routine
! should give the appropriate kick to the orbit and set stop_time to the next time when
! time_runge_kutta should call this routine.
!
! Note: If this routine sets orbit%state to anything but alive$, this will cause time_runge_kutta
! to stop tracking and return to the calling routine.
! 
! Input:
!   orbit       -- coord_struct: Particle orbit.
!     %vec(6)       -- "time" units: See convert_particle_coordinates_s_to_t for more details.
!   z_phase     -- real(rp): z phase space coordinate.
!   ele         -- ele_struct: Element to propagate the geometry through.
!   param       -- lat_param_struct: Branch parameters.
!   init_needed -- integer: Initialization needed? See above for details.
!
! Output:
!   orbit       -- coord_struct: Possibly modified particle orbit.
!   z_phase     -- real(rp): Possibly modified z phase space coordinate.
!   stop_time   -- real(rp): Set to time when time_runge_kutta should next call this routine.
!                    Set to real_garbage$ to prevent time_runge_kutta from stopping and calling this routine.
!-

subroutine time_runge_kutta_periodic_kick_hook (orbit, z_phase, ele, param, stop_time, init_needed)

use bmad, except_dummy => time_runge_kutta_periodic_kick_hook

implicit none

type (coord_struct) orbit
type (ele_struct) ele
type (lat_param_struct) param

real(rp) z_phase, stop_time
integer init_needed

! Init needed?

if (init_needed == true_int$) then
  stop_time = real_garbage$        ! Set stop_time
  return
endif

! Give particle a kick

! And set stop_time to the next time time_runge_kutta should next call this routine.

stop_time = real_garbage$

end subroutine time_runge_kutta_periodic_kick_hook
