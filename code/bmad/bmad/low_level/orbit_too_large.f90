!+
! Function orbit_too_large (orbit, param) result (is_too_large)
!
! Routine to check if an orbit is too large.
! This routine is used to prevent floating point overflow.
! Too large is defined by:
!   |x|, |y| > bmad_com%max_aperture_limit or 
!   |pz| > 1 (photons only) or 
!   px^2 + py^2 > (1 + pz)^2 (non-photons only)
!
! Also see:
!   check_aperture_limit
!
! Input:
!   orbit         -- coord_struct: Particle orbit.
!
! Output:
!   orbit         -- coord_struct: Particle orbit.
!     %state          -- Particle status. Not set if orbit is OK.
!   is_too_large  -- logical: True if orbit is too large. False otherwise.
!   param         -- lat_param_struct, optional: 
!     %unstable_factor  -- Set if orbit is too large. Otherwise not set
!-

function orbit_too_large (orbit, param) result (is_too_large)

use bmad_interface, dummy => orbit_too_large

implicit none

type (coord_struct) orbit
type (lat_param_struct), optional :: param

logical is_too_large
real(rp) rel_p

! Assume the worst

is_too_large = .true.

! Test aperture

if (orbit%vec(1) > bmad_com%max_aperture_limit) then
  orbit%state = lost_pos_x_aperture$
  if (present(param)) param%unstable_factor = (abs(orbit%vec(1)) - bmad_com%max_aperture_limit) / bmad_com%max_aperture_limit
  return
elseif (-orbit%vec(1) > bmad_com%max_aperture_limit) then
  orbit%state = lost_neg_x_aperture$
  if (present(param)) param%unstable_factor = (abs(orbit%vec(1)) - bmad_com%max_aperture_limit) / bmad_com%max_aperture_limit
  return
endif

if (orbit%vec(3) > bmad_com%max_aperture_limit) then
  orbit%state = lost_pos_y_aperture$
  if (present(param)) param%unstable_factor = (abs(orbit%vec(3)) - bmad_com%max_aperture_limit) / bmad_com%max_aperture_limit
  return
elseif (-orbit%vec(3) > bmad_com%max_aperture_limit) then
  orbit%state = lost_neg_y_aperture$
  if (present(param)) param%unstable_factor = (abs(orbit%vec(3)) - bmad_com%max_aperture_limit) / bmad_com%max_aperture_limit
  return
endif

! Test photons

if (orbit%species == photon$) then
  if (abs(orbit%vec(6)) > 1) then
    orbit%state = lost_pz_aperture$
    if (present(param)) param%unstable_factor = abs(orbit%vec(6)) - 1
    return
  endif

! charged particle test

else
  rel_p = 1 + orbit%vec(6)
  if (rel_p < 0) return

  if (orbit%vec(2)**2 + orbit%vec(4)**2 > rel_p**2) then
    if (present(param)) param%unstable_factor = sqrt(orbit%vec(2)**2 + orbit%vec(4)**2 - rel_p**2)

    if (abs(orbit%vec(2)) > abs(orbit%vec(4))) then
      if (orbit%vec(2) > 0) then
        orbit%state = lost_pos_x_aperture$
      else
        orbit%state = lost_neg_x_aperture$
      endif

    else
      if (orbit%vec(4) > 0) then
        orbit%state = lost_pos_y_aperture$
      else
        orbit%state = lost_neg_y_aperture$
      endif
    endif

    return
  endif

endif

! Passed tests.

is_too_large = .false.

end function orbit_too_large

