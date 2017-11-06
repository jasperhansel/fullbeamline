!+
! Subroutine calc_z_tune (lat, ix_branch)
!
! Routine to calculate the synchrotron tune from the full 6X6 1-turn matrix.
!
! Note: The tune will be negative above transition which corresponds to
! counter-clockwise rotation in z-pz space.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat       -- lat_struct: Lat
!   ix_branch -- integer, optional: Branch of lattice to analyze. Default is 0.
!
! Output:
!   lat -- lat_struct
!     branch(ix_branch)%z%tune            -- Synchrotron tune (radians)
!     branch(ix_branch)%param%t1_with_RF  -- 6x6 1-turn matrix.
!-

subroutine calc_z_tune (lat, ix_branch)

use bmad_interface, except_dummy => calc_z_tune
use nrtype
use nr

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

real(rp) a(6,6), wr(6), wi(6), cos_z, denom

integer, optional :: ix_branch
integer i, sgn

!

branch => lat%branch(integer_option(0, ix_branch))

call transfer_matrix_calc (lat, a, ix_branch = ix_branch)
branch%param%t1_with_RF = a

denom = 2 * (a(5,5)*a(6,6) - a(5,6)*a(6,5))
if (denom == 0) then
  branch%z%tune = 0
  return
endif

cos_z = (a(5,5) + a(6,6)) / denom
sgn = sign_of(a(5,6))

if (cos_z > 0.9999999) then
  branch%z%tune = 0
  return
endif

call balanc(a)
call elmhes(a)
call hqr(a,wr,wi)

! we need to find which eigen-value is closest to the z_tune

i = minloc(abs(wr-cos_z), 1)
branch%z%tune = sgn * abs(atan2(wi(i),wr(i)))

end subroutine
