	FUNCTION fpoly(x,n)
	USE nrtype; USE nrutil, ONLY : geop
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: n
	REAL(dp), DIMENSION(n) :: fpoly
	fpoly=geop(1.0_dp,x,n)
	END FUNCTION fpoly
