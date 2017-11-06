!+
! Subroutine track1_postprocess (start_orb, ele, param, end_orb)
!
! Dummy routine for post processing after the track1 routine is done.
!
! Also see:
!   track1_preprocess
!   track1_custom
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- Coord_struct: Starting position.
!   ele    -- Ele_struct: Element.
!   param  -- lat_param_struct: Lattice parameters.
!   end_orb   -- Coord_struct: End position.
!
! Output:
!   end_orb   -- Coord_struct: End position.
!-

subroutine track1_postprocess (start_orb, ele, param, end_orb)

use bmad, except_dummy => track1_postprocess

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param

character(*), parameter :: r_name = 'track1_postprocess'

!

end subroutine
