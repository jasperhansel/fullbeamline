	FUNCTION gammp_s(a,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,x
	REAL(dp) :: gammp_s
	call assert( x >= 0.0,  a > 0.0, 'gammp_s args')
	if (x<a+1.0_dp) then
		gammp_s=gser(a,x)
	else
		gammp_s=1.0_dp-gcf(a,x)
	end if
	END FUNCTION gammp_s


	FUNCTION gammp_v(a,x)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: a,x
	REAL(dp), DIMENSION(size(x)) :: gammp_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(x),'gammp_v')
	call assert( all(x >= 0.0),  all(a > 0.0), 'gammp_v args')
	mask = (x<a+1.0_dp)
	gammp_v=merge(gser(a,merge(x,0.0_dp,mask)), &
		1.0_dp-gcf(a,merge(x,0.0_dp,.not. mask)),mask)
	END FUNCTION gammp_v
