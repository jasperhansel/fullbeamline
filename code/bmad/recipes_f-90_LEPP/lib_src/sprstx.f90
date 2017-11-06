	SUBROUTINE sprstx_dp(sa,x,b)
	USE nrtype; USE nrutil, ONLY : assert_eq,scatter_add
	IMPLICIT NONE
	TYPE(sprs2_dp), INTENT(IN) :: sa
	REAL(dp), DIMENSION (:), INTENT(IN) :: x
	REAL(dp), DIMENSION (:), INTENT(OUT) :: b
	INTEGER(I4B) :: ndum
	ndum=assert_eq(sa%n,size(x),size(b),'sprstx_dp')
	b=0.0_dp
	call scatter_add(b,sa%val*x(sa%irow),sa%jcol)
	END SUBROUTINE sprstx_dp

