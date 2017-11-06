	FUNCTION erfcc_s(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: erfcc_s
	REAL(dp) :: t,z
	REAL(dp), DIMENSION(10) :: coef = (/-1.26551223_dp,1.00002368_dp,&
		0.37409196_dp,0.09678418_dp,-0.18628806_dp,0.27886807_dp,&
		-1.13520398_dp,1.48851587_dp,-0.82215223_dp,0.17087277_dp/)
	z=abs(x)
	t=1.0_dp/(1.0_dp+0.5_dp*z)
	erfcc_s=t*exp(-z*z+poly(t,coef))
	if (x < 0.0) erfcc_s=2.0_dp-erfcc_s
	END FUNCTION erfcc_s


	FUNCTION erfcc_v(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x
	REAL(dp), DIMENSION(size(x)) :: erfcc_v,t,z
	REAL(dp), DIMENSION(10) :: coef = (/-1.26551223_dp,1.00002368_dp,&
		0.37409196_dp,0.09678418_dp,-0.18628806_dp,0.27886807_dp,&
		-1.13520398_dp,1.48851587_dp,-0.82215223_dp,0.17087277_dp/)
	z=abs(x)
	t=1.0_dp/(1.0_dp+0.5_dp*z)
	erfcc_v=t*exp(-z*z+poly(t,coef))
	where (x < 0.0) erfcc_v=2.0_dp-erfcc_v
	END FUNCTION erfcc_v
