!+
! Subroutine set_particle_from_rf_time (rf_time, ele, apply_hard_edge_offset, orbit)
!
! Routine to set the z (vec(5)) and t components of a particle orbit given a time related to the "RF clock".
! Also see particle_rf_time which is the inverse of this routine.
!
! Input:
!   rf_time   -- rel(rp): Time related to the RF clock.
!   ele       -- rel(rp): Element that uses this clock.
!   apply_hard_edge_offset -- logical: If True then account for offset due to a hard edge model.
!
! Ouput:
!   orbit     -- coord_struct: Particle orbit.
!-

subroutine set_particle_from_rf_time (rf_time, ele, apply_hard_edge_offset, orbit)

use multipass_mod, dummy => set_particle_from_rf_time

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (ele_struct), pointer :: ref_ele
type (ele_pointer_struct), allocatable :: chain(:)

real(rp) rf_time, s_hard_offset
integer ix_pass, n_links
logical apply_hard_edge_offset

!

ref_ele => ele
if (ref_ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) ref_ele => pointer_to_lord (ref_ele, 1)

call multipass_chain(ref_ele, ix_pass, n_links, chain)
if (ix_pass > 1) ref_ele => chain(1)%ele

orbit%t = rf_time + ref_ele%value(ref_time_start$)

!

orbit%vec(5) = -rf_time * orbit%beta * c_light

!

if (apply_hard_edge_offset .and. has_attribute(ele, 'L_HARD_EDGE')) then
  s_hard_offset = (ele%value(l$) - hard_edge_model_length(ele)) / 2  ! Relative to entrance end of the cavity
  orbit%t = orbit%t + s_hard_offset / (c_light * orbit%beta)
  orbit%vec(5) = orbit%vec(5) + c_light * orbit%beta * s_hard_offset
endif



end subroutine set_particle_from_rf_time
