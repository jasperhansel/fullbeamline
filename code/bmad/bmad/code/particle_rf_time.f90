!+
! Function particle_rf_time (orbit, ele, apply_hard_edge_offset, s_rel) result (time)
!
! Routine to return the reference time used to calculate the phase of
! time-dependent EM fields.
!
! Essentially: The oscillations of EM fields are synched relative to the absolute clock if 
! absolute time tracking is used or are synched relative to the reference particle
! if relative time tracking is used. See the Bmad manual for more details.
!
! Also see set_particle_from_rf_time which is the inverse of this routine.
!
! Input:
!   orbit   -- Coord_struct: Particle coordinates
!   ele     -- ele_struct: Element being tracked through.
!   apply_hard_edge_offset -- logical: If True then account for offset due to a hard edge model.
!   s_rel   -- real(rp), optional: Longitudinal position relative to the upstream edge of the element.
!               Needed for relative time tracking when the particle is inside the element. Default is 0.
!
! Ouput:
!   time    -- Real(rp): Current time.
!-

function particle_rf_time (orbit, ele, apply_hard_edge_offset, s_rel) result (time)

use multipass_mod, dummy => particle_rf_time

implicit none

type (coord_struct) orbit
type (ele_struct), target :: ele
type (ele_struct), pointer :: ref_ele
type (ele_pointer_struct), allocatable :: chain(:)

real(rp) time, s_hard_offset, beta0
real(rp), optional :: s_rel
integer ix_pass, n_links
logical apply_hard_edge_offset, abs_time

character(*), parameter :: r_name = 'particle_rf_time'

!

ref_ele => ele
if (ref_ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) ref_ele => pointer_to_lord (ref_ele, 1)

call multipass_chain(ref_ele, ix_pass, n_links, chain)
if (ix_pass > 1) ref_ele => chain(1)%ele

! With absolute time tracking the reference time is relative to the reference time of the element.
! This way the phase does not have to be adjusted when switching between absolute and relative time tracking.
! Note: With a multipass_slave, use the reference time of the pass element.
! Note: e_gun uses absolute time tracking to get around the problem when orbit%beta = 0.

if (absolute_time_tracking(ele)) then
  time = orbit%t - ref_ele%value(ref_time_start$)

else
  if (orbit%beta == 0) then
    call out_io (s_fatal$, r_name, 'PARTICLE IN NON E-GUN ELEMENT HAS VELOCITY = 0!')
    if (global_com%exit_on_error) call err_exit
    time = orbit%t  ! Just to keep on going
    return
  endif
  time = -orbit%vec(5) / (orbit%beta * c_light)

  if (present(s_rel)) then
    beta0 = ele%value(p0c$) / ele%value(E_tot$)
    time = time + (s_rel + ele%s_start - ref_ele%s_start) / (c_light * beta0)
  endif
endif

!

if (apply_hard_edge_offset .and. has_attribute(ref_ele, 'L_HARD_EDGE')) then
  s_hard_offset = (ref_ele%value(l$) - hard_edge_model_length(ref_ele)) / 2  ! Relative to entrance end of the cavity
  beta0 = ele%value(p0c_start$) / ele%value(E_tot_start$)
  time = time - s_hard_offset / (c_light * beta0)
endif

end function particle_rf_time

