	FUNCTION gammq_s(a,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,x
	REAL(dp) :: gammq_s
	call assert( x >= 0.0,  a > 0.0, 'gammq_s args')
	if (x<a+1.0_dp) then
		gammq_s=1.0_dp-gser(a,x)
	else
		gammq_s=gcf(a,x)
	end if
	END FUNCTION gammq_s


	FUNCTION gammq_v(a,x)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: a,x
	REAL(dp), DIMENSION(size(a)) :: gammq_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(x),'gammq_v')
	call assert( all(x >= 0.0),  all(a > 0.0), 'gammq_v args')
	mask = (x<a+1.0_dp)
	gammq_v=merge(1.0_dp-gser(a,merge(x,0.0_dp,mask)), &
		gcf(a,merge(x,0.0_dp,.not. mask)),mask)
	END FUNCTION gammq_v
