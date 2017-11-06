	FUNCTION pythag_dp(a,b)
	USE nrtype
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp) :: pythag_dp
	REAL(dp) :: absa,absb
	absa=abs(a)
	absb=abs(b)
	if (absa > absb) then
		pythag_dp=absa*sqrt(1.0_dp+(absb/absa)**2)
	else
		if (absb == 0.0) then
			pythag_dp=0.0
		else
			pythag_dp=absb*sqrt(1.0_dp+(absa/absb)**2)
		end if
	end if
	END FUNCTION pythag_dp

