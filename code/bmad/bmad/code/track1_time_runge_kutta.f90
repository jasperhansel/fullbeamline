!+ 
! Subroutine track1_time_runge_kutta(start_orb, ele, param, end_orb, err_flag, track)
!
! Routine to track a particle through an element using 
! Runge-Kutta time-based tracking. Converts to and from element
! coordinates before and after tracking.
!
! Modules Needed:
!   use time_tracker_mod
!
! Input:
!   start_orb   -- coord_struct: starting position, t-based global
!   ele         -- ele_struct: element
!   param       -- lat_param_struct: lattice parameters
!
! Output:
!   end_orb     -- coord_struct: end position, t-based global
!   err_flag    -- Logical: Set True if there is an error. False otherwise
!   track       -- track_struct (optional): Contains array of the step-by-step particle
!                    trajectory along with the field at these positions.
!                    When tracking through multiple elements, the trajectory in an element
!                    is appended to the existing trajectory. To reset: Set track%n_pt = -1.
!-

subroutine track1_time_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)

use time_tracker_mod, dummy => track1_time_runge_kutta
use track1_mod, dummy2 => track1_time_runge_kutta

implicit none

type (coord_struct) :: start_orb, end_orb, start_orb_saved
type (coord_struct) :: ele_origin
type (lat_param_struct), target :: param
type (ele_struct), target :: ele
type (ele_struct), pointer :: hard_ele
type (track_struct), optional :: track
type (em_field_struct) :: saved_field

real(rp) vec(6), d_radius
real(rp) s_lab, s0, s1, s2, ds_ref, del_s, p0c_save, s_body, z_phase
real(rp) s_edge_track, s_edge_hard, rf_time, beta_ref, r, dref_time

integer :: i, hard_end, t_dir

logical :: abs_time, err_flag, err, set_spin

character(*), parameter :: r_name = 'track1_time_runge_kutta'

!---------------------------------

if (ele%key /= patch$ .and. ele%value(l$) == 0) then
  call track_a_zero_length_element (start_orb, ele, param, end_orb, err_flag, track)
  return
endif

!

err_flag = .true.

t_dir = 1
start_orb_saved = start_orb
end_orb = start_orb
rf_time = particle_rf_time (end_orb, ele, .true., end_orb%s - ele%s_start)
beta_ref = ele%value(p0c$) / ele%value(E_tot$)
set_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$ .and. &
            (ele%field_calc == bmad_standard$ .or. ele%field_calc == fieldmap$) .and. &
            is_true(ele%value(spin_fringe_on$)))

! s_lab is longitudinal position relative to entrance end (element body coords).
! Remember that if an element has reversed orientation, +s direction in body coords is opposite to
! the reference trajectory +s direction.
! Adjust to match element edges if close enough


s_lab =  end_orb%s - ele%s_start

if (abs(s_lab) < bmad_com%significant_length) then
  s_lab = 0
else if (abs(s_lab - ele%value(l$))  < bmad_com%significant_length) then
  s_lab = ele%value(l$)
endif

!------
! Check wall

call check_aperture_limit (end_orb, ele, in_between$, param)
if (end_orb%state /= alive$) then
  call out_io (s_info$, r_name, "PARTICLE STARTED IN REGION OUTSIDE OF WALL: "//trim(ele%name))
  end_orb%state = lost$
  return
endif

!------
! Convert particle to element coordinates
! Kicks and multipoles should be turned off in offset_particle

! Interior start, reference momentum is at the end.
if (end_orb%location == inside$) then
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., s_pos =s_lab, s_out = s_body, set_spin = set_spin)
  if (ele%value(l$) < 0) t_dir = -1

elseif (ele%key == patch$) then
  if (start_orb_saved%location == inside$) then
    call out_io (s_error$, r_name, 'TIME-RUNGE-KUTTA TRACKING STARTING WITH A INSIDE A PATCH IS NOT PERMITED: ' // ele%name)
    return
  endif
  call track_a_patch (ele, end_orb, .false., s0, ds_ref)
  end_orb%vec(5) = end_orb%vec(5) + (ds_ref + s0 * end_orb%direction * ele%orientation) * end_orb%beta / beta_ref 
  if (s0*ele%orientation > 0) t_dir = -1
  s_body = s0

! Particle is at an end.
else
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., s_out = s_body, set_spin = set_spin)
  if (ele%value(l$) < 0) t_dir = -1
endif

! Convert orbit coords to time based.

call convert_particle_coordinates_s_to_t (end_orb, s_body, ele%orientation, z_phase)

!

if (present(track)) then
  call save_a_step (track, ele, param, .false., start_orb_saved, s_lab, .true., rf_time = rf_time)
endif

! Track through element

if ((ele%key == lcavity$ .or. ele%key == rfcavity$) .and. bmad_com%use_hard_edge_drifts .and. &
                  ele%field_calc == bmad_standard$ .and. ele%value(l$) < ele%value(l_hard_edge$)) then
  call out_io (s_error$, r_name, 'TIME-RUNGE-KUTTA TRACKING THROUGH RF CAVITY: ' // ele%name, &
                          'WILL NOT BE ACCURATE SINCE THE LENGTH IS LESS THAN THE HARD EDGE MODEL LENGTH.')
endif

call odeint_bmad_time(end_orb, z_phase, ele, param, t_dir, rf_time, err, track)

if (err) return

!------
! Convert back to s-based coordinates.
! The particle is either dead or is alive and at an element end.

if (end_orb%location == upstream_end$) then
  end_orb%p0c = ele%value(p0c_start$)
  call convert_particle_coordinates_t_to_s(end_orb, z_phase, ele, s_body)
  end_orb%direction = -1  ! In case t_to_s conversion confused by roundoff error.
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., s_pos = s_body, set_spin = set_spin)

elseif (end_orb%location == downstream_end$) then
  end_orb%p0c = ele%value(p0c$)
  call convert_particle_coordinates_t_to_s(end_orb, z_phase, ele, s_body)
  end_orb%direction = 1  ! In case t_to_s conversion confused by roundoff error
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., s_pos = s_body, set_spin = set_spin)

elseif (end_orb%state /= alive$) then
  ! Particle is lost in the interior of the element.
  ! The reference energy in the interior is equal to the ref energy at the exit end of the element.
  if (end_orb%vec(6) < 0) then
    end_orb%p0c = ele%value(p0c$)
    end_orb%direction = -1
  else 
    end_orb%p0c = ele%value(p0c$)
    end_orb%direction = 1
  end if

  call convert_particle_coordinates_t_to_s(end_orb, z_phase, ele, s_body)
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., s_pos = s_body, set_spin = set_spin)

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE LEAVING ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! Set relativistic beta

call convert_pc_to (end_orb%p0c * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)

! The z value computed in odeint_bmad_time is off for elements where the particle changes energy is not 
! constant and when the ref time is not calcuated via the length of the reference orbit (as with a wiggler). 
! In this case, make the needed correction. Odeint_bmad_time uses a reference time
! assuming that the reference velocity is constant and equal to the velocity at the final energy.
! If the particle has started inside the element then just assume that the reference time is linear in s.
! This is not a great approximation but absolute time tracking should be used here which makes the
! tracking independent of z.

if (ele%key /= patch$) then
  r = abs(end_orb%s - start_orb_saved%s) / ele%value(l$)
  dref_time = ele%value(l$) / (beta_ref * c_light)
  end_orb%vec(5) = end_orb%vec(5) + r * (ele%value(delta_ref_time$) - dref_time) * end_orb%beta * c_light
endif

err_flag = .false.

end subroutine

