!+
! Subroutine mat_make_unit (mat)
!
! routine to create a unit matrix.
!
! Output:
!   mat(:,:) -- Real(rp): Square matrix
!-

subroutine mat_make_unit (mat)

  use precision_def

  implicit none

  integer i
  real(rp) mat(:,:)

!

  mat = 0
  forall (i = 1:size(mat, 1)) mat(i,i) = 1

end subroutine

