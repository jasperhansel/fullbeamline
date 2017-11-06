!+
! Subroutine track1 (start_orb, ele, param, end_orb, track, err_flag, ignore_radiation, mat6, make_matrix)
!
! Particle tracking through a single element. 
! This includes synchrotron radiation and space charge kicks if enabled by the appropriate
! bmad_com%radiation_damping_on, etc. components. See the code for more details. 
!
! Note: the ignore_radiation argument is meant as a "temporary" override to turn off
! all radiation and space-charge effects independent of the settings of the bmad_com%xxx
! switches. This is used by routines that are tracking the reference particle.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb   -- Coord_struct: Starting position.
!   ele         -- Ele_struct: Element to track through.
!   param       -- lat_param_struct: Reference particle info.
!   track       -- track_struct, optional: Structure holding existing track.
!   ignore_radiation
!               -- Logical, optional: If present and True then do not include radiation
!                  effects along with space charge effects. 
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before the element. Only used with bmad_standard tracking.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   end_orb     -- Coord_struct: End position.
!   track       -- track_struct, optional: Structure holding the track information if the 
!                    tracking method does tracking step-by-step.
!                    When tracking through multiple elements, the trajectory in an element
!                    is appended to the existing trajectory. To reset: Set track%n_pt = -1.
!   err_flag    -- Logical, optional: Set true if there is an error. False otherwise.
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix including the element. Only used with bmad_standard tracking.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!-

recursive subroutine track1 (start_orb, ele, param, end_orb, track, err_flag, ignore_radiation, mat6, make_matrix)

use bmad, except_dummy1 => track1
use mad_mod, only: track1_mad
use boris_mod, only: track1_boris
use space_charge_mod, except_dummy2 => track1
use spin_mod, except_dummy3 => track1

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb
type (ele_struct)   :: ele
type (ele_struct), pointer :: lord
type (lat_param_struct) :: param
type (track_struct), optional :: track

real(rp), optional :: mat6(6,6)
real(rp) p0c_start, z_start
integer tracking_method, stm

character(8), parameter :: r_name = 'track1'

logical, optional :: make_matrix
logical, optional :: err_flag, ignore_radiation
logical err, do_extra, finished, radiation_included, do_spin_tracking

! Use start2_orb since start_orb is strictly input.

if (present(err_flag)) err_flag = .true.
start2_orb = start_orb

do_extra = .not. logic_option(.false., ignore_radiation)

! symp_lie_bmad tracking does include radiation effects

radiation_included = (ele%tracking_method == symp_lie_bmad$) 

! For historical reasons, the calling routine may not have correctly 
! initialized the starting orbit. If so, we do an init here.
! Time Runge-Kutta is tricky so do not attempt to do a set.

if (start2_orb%state == not_set$) then
  if (ele%tracking_method == time_runge_kutta$) then
    call out_io (s_error$, r_name, 'STARTING ORBIT NOT PROPERLY INITIALIZED! [NEEDS AN INIT_COORD CALL?]')
    return
  endif
  call init_coord(start2_orb, start2_orb, ele, upstream_end$, particle = default_tracking_species(param)) 
endif

! Set start2_orb%location appropriate for tracking through the element.
! For time runge-kutta the particle may be starting inside the element.
! In this case, do not set the location.

if (ele%tracking_method /= time_runge_kutta$ .or. start2_orb%location /= inside$) then
  if (start2_orb%direction == 1) then
    start2_orb%location = upstream_end$
    start2_orb%s = ele%s_start
  else
    start2_orb%location = downstream_end$
    start2_orb%s = ele%s
  endif
endif

if (start2_orb%species /= photon$ .and. start2_orb%state == alive$) then
  err = .false.

  if (ele%key == beginning_ele$) then
    if (abs(ele%value(p0c$) - start2_orb%p0c) > 1d-15*start2_orb%p0c) err = .true. ! For e_gun case
  else if (start2_orb%location == upstream_end$) then
    if (abs(ele%value(p0c_start$) - start2_orb%p0c) > 1d-15*start2_orb%p0c) err = .true.
  else if (start2_orb%location == downstream_end$) then
    if (abs(ele%value(p0c$) - start2_orb%p0c) > 1d-15*start2_orb%p0c) err = .true. ! For e_gun case
  endif

  if (err) then
    call out_io (s_error$, r_name, 'STARTING ORBIT NOT PROPERLY INITIALIZED! [NEEDS AN INIT_COORD CALL?]')
    return
  endif
endif

! Photons get the z-coordinate reset to zero.

if (start2_orb%species == photon$ .and. start2_orb%location == downstream_end$) then
  start2_orb%vec(5) = 0
endif

! Check for particle lost even before we begin

if (start2_orb%state /= alive$) then
  end_orb = start_orb
  end_orb%vec = 0
  if (present(err_flag)) err_flag = .false.
  return
endif

! If a particle is inside the element then only time_runge_kutta
! can handle this situation.

if (start2_orb%state == inside$ .and. ele%tracking_method /= time_runge_kutta$) then
  call out_io (s_error$, r_name, 'PARTICLE''S STARTING POSITION IS INSIDE THE ELEMENT! ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! Custom tracking if the custom routine is to do everything.

finished = .false.
call track1_preprocess (start2_orb, ele, param, err, finished, radiation_included, track)
if (err) return
if (finished) then
  end_orb = start2_orb
  if (present(err_flag)) err_flag = err
  return
endif

! Init

if (bmad_com%auto_bookkeeper) call attribute_bookkeeper (ele, param)

! check for particles outside aperture.

call check_aperture_limit (start2_orb, ele, first_track_edge$, param)

if (start2_orb%state /= alive$) then
  end_orb = start2_orb
  if (present(err_flag)) err_flag = .false.
  return
endif

! Radiation damping and/or fluctuations for the 1st half of the element.

if ((bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) .and. &
                                           .not. radiation_included .and. ele%is_on .and. do_extra) then
  call track1_radiation (start2_orb, ele, param, start2_orb, start_edge$, z_start) 
endif

! bmad_standard handles the case when the element is turned off.

do_spin_tracking = bmad_com%spin_tracking_on

tracking_method = ele%tracking_method
if (.not. ele%is_on) tracking_method = bmad_standard$

select case (tracking_method)

case (bmad_standard$)
  if (start2_orb%species == photon$) then
    call track1_bmad_photon (start2_orb, ele, param, end_orb, err)
  else
    call track1_bmad (start2_orb, ele, param, end_orb, err, mat6, make_matrix)
  endif

  if (present(track)) call add_to_track()
  if (err) return

case (runge_kutta$) 
  call track1_runge_kutta (start2_orb, ele, param, end_orb, err, track)
  if (err) return

case (linear$) 
  call track1_linear (start2_orb, ele, param, end_orb)
  if (present(track)) call add_to_track()

case (taylor$) 
  call track1_taylor (start2_orb, ele, param, end_orb)
  if (present(track)) call add_to_track()

case (symp_map$) 
  call track1_symp_map (start2_orb, ele, param, end_orb)
  if (present(track)) call add_to_track()

case (symp_lie_bmad$) 
  call symp_lie_bmad (ele, param, start2_orb, end_orb, .false., track)

case (symp_lie_ptc$)
  call track1_symp_lie_ptc (start2_orb, ele, param, end_orb, track)
  stm = ele%spin_tracking_method
  if (stm == tracking$ .or. stm == symp_lie_ptc$) do_spin_tracking = .false.

case (boris$) 
  call track1_boris (start2_orb, ele, param, end_orb, err, track)
  if (err) return

case (mad$)
  call track1_mad (start2_orb, ele, param, end_orb)
  if (present(track)) call add_to_track()

case (custom$)
  call track1_custom (start2_orb, ele, param, end_orb, err, finished, track)
  if (err) return
  if (finished) then
    if (present(err_flag)) err_flag = err
    return
  endif

case (time_runge_kutta$)
  call track1_time_runge_kutta (start2_orb, ele, param, end_orb, err, track)
  if (err) return

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN TRACKING_METHOD: \i0\ ', ele%tracking_method)
  if (global_com%exit_on_error) call err_exit
  return

end select

! Check

if (orbit_too_large (end_orb, param)) return

! spin tracking. Must do after regular tracking in the case of spin_tracking_method = bmad_standard
 
if (do_spin_tracking) call track1_spin (start2_orb, ele, param, end_orb)

! Set ix_ele. If the element is a slice_slave then the appropriate ix_ele is given by the lord.

if (ele%slave_status == slice_slave$) then
  lord => pointer_to_lord (ele, 1)
  end_orb%ix_ele = lord%ix_ele
else
  end_orb%ix_ele = ele%ix_ele
endif

if (tracking_method /= time_runge_kutta$) then
  if (end_orb%state /= alive$) then
    end_orb%location = inside$
  elseif (start2_orb%direction == 1) then
    end_orb%location = downstream_end$
  else
    end_orb%location = upstream_end$
  endif
endif

! Radiation damping and/or fluctuations for the last half of the element

if ((bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) .and. &
                                            .not. radiation_included .and. ele%is_on .and. do_extra) then
  call track1_radiation (end_orb, ele, param, end_orb, end_edge$, z_start) 
endif

! space charge

if (bmad_com%space_charge_on .and. do_extra) call track1_ultra_rel_space_charge (ele, param, end_orb)

! check for particles outside aperture

call check_aperture_limit (end_orb, ele, second_track_edge$, param)

call track1_postprocess (start2_orb, ele, param, end_orb)

if (present(err_flag)) err_flag = .false.

!--------------------------------------------------------------------
contains

! Add the the track. 

subroutine add_to_track()
if (start_orb%direction == 1) then
  call save_a_step(track, ele, param, .false., start_orb, 0.0_rp)
  call save_a_step(track, ele, param, .false., end_orb, ele%value(l$))
else
  call save_a_step(track, ele, param, .false., start_orb, ele%value(l$))
  call save_a_step(track, ele, param, .false., end_orb, 0.0_rp)
endif
end subroutine add_to_track

end subroutine
