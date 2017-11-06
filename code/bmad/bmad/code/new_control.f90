!+
! Subroutine new_control (lat, ix_ele)
!
! Routine to create a new control element.
!
! Note: This routine may reallocate the lat%branch arrays so
! any pointers to elements will need to be repointed.
!
! Modules Needed:
!   use bmad
!
! Input:
!     lat -- lat_struct: Lat used
!
! Output
!     ix_ele -- Integer: Index of the new control element
!-

subroutine new_control (lat, ix_ele)

use bmad_interface, except_dummy => new_control

implicit none

type (lat_struct)  lat
integer ix_ele

!

lat%n_ele_max = lat%n_ele_max + 1
ix_ele = lat%n_ele_max
if (ix_ele > ubound(lat%ele, 1))  call allocate_lat_ele_array (lat)

end subroutine
