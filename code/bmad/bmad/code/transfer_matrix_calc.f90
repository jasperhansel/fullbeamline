!+         
! Subroutine transfer_matrix_calc (lat, xfer_mat, xfer_vec, ix1, ix2, ix_branch, one_turn)
!
! Subroutine to calculate the 6X6 transfer matrix between two elements.
! To calculate transfer maps see the routine: transfer_map_calc.
!
! The transfer matrix is from the end of element ix1 to the end of element ix2.
! If ix1 and ix2 are not present, the full 1-turn matrix is calculated.
!
! If ix2 < ix1 and lat%param%geometry is closed$ then the
! calculation will "wrap around" the lattice end.
! For example, if ix1 = 900 and ix2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If ix2 < ix1 and lat%param%geometry is open$ then the backwards
! transfer matrix is computed.
!
! If ix1 = ix2 then xfr_mat is the unit matrix except if one_turn = True and 
! the lattice is circular.
!
! To get the one-turn matrix for a circular lattice, omit ix2 from the argument list.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat   -- lat_struct: Lattice used in the calculation.
!     %ele(:)%mat6  -- Transfer matrices used in the calculation.
!   ix1   -- Integer, optional: Element start index for the calculation.
!              Default is 0.
!   ix2   -- Integer, optional: Element end index for the calculation.
!              Defaults:
!                If ix1 is not present: ix2 = lat%n_ele_track
!                If ix1 is present and lattice is circular: Calculate the 
!                  one-turn matrix from ix1 back to ix1.
!   ix_branch -- Integer, optional: Branch index. Default is 0.
!   one_turn   -- Logical, optional: If present and True, and ix1 = ix2, and the lattice
!                   is circular: Construct the one-turn matrix from ix1 back to ix1.
!                   Otherwise mat6 is unchanged or the unit map if unit_start = T.
!                   Default = False.
!
! Output:
!    xfr_mat(6,6) -- Real(rp): Transfer matrix.
!    xfr_vec(6)   -- Real(rp), optional: 0th order part of the transfer map.
!-

subroutine transfer_matrix_calc (lat, xfer_mat, xfer_vec, ix1, ix2, ix_branch, one_turn)

use bmad_interface, except_dummy => transfer_matrix_calc
use sim_utils, only: integer_option

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

real(rp) :: xfer_mat(:,:)
real(rp), optional :: xfer_vec(:)
real(rp) rf_mat(6,6)

integer, optional :: ix1, ix2, ix_branch
integer i, i1, i2

logical vec_present, one_turn_this
logical, optional :: one_turn

character(*), parameter :: r_name = 'transfer_matrix_calc'

!

vec_present = present(xfer_vec)
if (vec_present) xfer_vec = 0

call mat_make_unit (xfer_mat)

branch => lat%branch(integer_option(0, ix_branch))  
i1 = integer_option(0, ix1) 
i2 = integer_option(branch%n_ele_track, ix2) 
one_turn_this = logic_option (.false., one_turn)

if (i1 < 0 .or. i1 > branch%n_ele_track .or. i2 < 0 .or. i2 > branch%n_ele_track) then
  call out_io (s_fatal$, r_name, 'ELEMENT INDEXES OUT OF BOUNDS: \i0\, \i0\ ', &
                                 i_array = [i1, i2])
  if (global_com%exit_on_error) call err_exit
  xfer_mat = 0
  return
endif 

if (branch%param%geometry == closed$ .and. present(ix1) .and. .not. present(ix2)) then
  i2 = i1
endif

if (i1 == i2 .and. .not. one_turn_this) return

! Normal case with i1 < i2

if (i1 < i2) then  
  do i = i1+1, i2
    call add_on_to_xfer_mat
  enddo

! Here when i2 < i1 ...
! For a circular lattice push through the origin.

elseif (branch%param%geometry == closed$ .or. one_turn_this) then
  do i = i1+1, branch%n_ele_track
    call add_on_to_xfer_mat
  enddo
  do i = 1, i2
    call add_on_to_xfer_mat
  enddo

! For a linear lattice compute the backwards matrix

else
  do i = i2+1, i1 
    call add_on_to_xfer_mat
  enddo
  call mat_inverse (xfer_mat, xfer_mat)
  xfer_vec = -matmul(xfer_mat, xfer_vec)

endif

!--------------------------------------------------------
contains

subroutine add_on_to_xfer_mat

xfer_mat = matmul (branch%ele(i)%mat6, xfer_mat)
if (vec_present) xfer_vec = matmul(branch%ele(i)%mat6, xfer_vec) + branch%ele(i)%vec0

end subroutine

end subroutine
