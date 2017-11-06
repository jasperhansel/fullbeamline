	FUNCTION erfc_s(x)
	USE nrtype
	USE nr, ONLY : gammp,gammq
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: erfc_s
	erfc_s=merge(1.0_dp+gammp(0.5_dp,x**2),gammq(0.5_dp,x**2), x < 0.0)
	END FUNCTION erfc_s


	FUNCTION erfc_v(x)
	USE nrtype
	USE nr, ONLY : gammp,gammq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x
	REAL(dp), DIMENSION(size(x)) :: erfc_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	mask = (x < 0.0)
	erfc_v=merge(1.0_dp+gammp(spread(0.5_dp,1,size(x)), &
		merge(x,0.0_dp,mask)**2),gammq(spread(0.5_dp,1,size(x)), &
		merge(x,0.0_dp,.not. mask)**2),mask)
	END FUNCTION erfc_v
