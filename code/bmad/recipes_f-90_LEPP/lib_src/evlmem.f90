	FUNCTION evlmem(fdt,d,xms)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: fdt,xms
	REAL(dp), DIMENSION(:), INTENT(IN) :: d
	REAL(dp) :: evlmem
	COMPLEX(dpc) :: z,zz
	REAL(DP) :: theta
	theta=TWOPI_D*fdt
	z=cmplx(cos(theta),sin(theta),kind=dpc)
	zz=1.0_dp-z*poly(z,d)
	evlmem=xms/abs(zz)**2
	END FUNCTION evlmem
