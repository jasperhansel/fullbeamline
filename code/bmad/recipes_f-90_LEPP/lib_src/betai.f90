	FUNCTION betai_s(a,b,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : betacf,gammln
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b,x
	REAL(dp) :: betai_s
	REAL(dp) :: bt
	call assert(x >= 0.0, x <= 1.0, 'betai_s arg')
	if (x == 0.0 .or. x == 1.0) then
		bt=0.0
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)&
			+a*log(x)+b*log(1.0_dp-x))
	end if
	if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
		betai_s=bt*betacf(a,b,x)/a
	else
		betai_s=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
	end if
	END FUNCTION betai_s


	FUNCTION betai_v(a,b,x)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : betacf,gammln
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: a,b,x
	REAL(dp), DIMENSION(size(a)) :: betai_v
	REAL(dp), DIMENSION(size(a)) :: bt
	LOGICAL(LGT), DIMENSION(size(a)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(b),size(x),'betai_v')
	call assert(all(x >= 0.0), all(x <= 1.0), 'betai_v arg')
	where (x == 0.0 .or. x == 1.0)
		bt=0.0
	elsewhere
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)&
			+a*log(x)+b*log(1.0_dp-x))
	end where
	mask=(x < (a+1.0_dp)/(a+b+2.0_dp))
	betai_v=bt*betacf(merge(a,b,mask),merge(b,a,mask),&
		merge(x,1.0_dp-x,mask))/merge(a,b,mask)
	where (.not. mask) betai_v=1.0_dp-betai_v
	END FUNCTION betai_v
