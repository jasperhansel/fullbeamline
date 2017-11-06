SUBROUTINE cplx_lubksb(a,indx,b)

  USE nrtype; USE nrutil, ONLY : assert_eq

  IMPLICIT NONE

  COMPLEX(dp), DIMENSION(:,:), INTENT(IN) :: a
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
  COMPLEX(dp), DIMENSION(:), INTENT(INOUT) :: b
  INTEGER(I4B) :: i,n,ii,ll
  COMPLEX(dp) :: summ

  n=assert_eq(size(a,1),size(a,2),size(indx),'cplx_lubksb')

  !- indx is ignored because the pivoting approach used in nr does not work
  !- when the columns are complex conjugate pairs

  !- Diagonal elements of L are assumed to be 1.0
  !- First use forward substitution to solve Ly=b
  do i=2,n
    ! the conjg below un-does the conjg done by fortran's dot_product.
    b(i) = b(i) - DOT_PRODUCT(CONJG(a(i,1:i-1)),b(1:i-1))
  enddo

  b(n) = b(n)/a(n,n)
  do i=n-1,1,-1
    b(i) = ( b(i) - DOT_PRODUCT(CONJG(a(i,i+1:n)),b(i+1:n)) ) / a(i,i)
  enddo

END SUBROUTINE cplx_lubksb
