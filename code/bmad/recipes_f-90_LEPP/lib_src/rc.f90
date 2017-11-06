	FUNCTION rc_s(x,y)
	USE nrtype; USE nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x,y
	REAL(dp) :: rc_s
	REAL(dp), PARAMETER :: ERRTOL=0.04_dp,TINY=1.69e-38_dp,&
		SQRTNY=1.3e-19_dp,BIG=3.0e37_dp,TNBG=TINY*BIG,&
		COMP1=2.236_dp/SQRTNY,COMP2=TNBG*TNBG/25.0_dp,&
		THIRD=1.0_dp/3.0_dp,&
		C1=0.3_dp,C2=1.0_dp/7.0_dp,C3=0.375_dp,C4=9.0_dp/22.0_dp
	REAL(dp) :: alamb,ave,s,w,xt,yt
	call assert( (/x >= 0.0,y /= 0.0,x+abs(y) >= TINY,x+abs(y) <= BIG, &
		y >= -COMP1 .or. x <= 0.0 .or. x >= COMP2/),'rc_s')
	if (y > 0.0) then
		xt=x
		yt=y
		w=1.0
	else
		xt=x-y
		yt=-y
		w=sqrt(x)/sqrt(xt)
	end if
	do
		alamb=2.0_dp*sqrt(xt)*sqrt(yt)+yt
		xt=0.25_dp*(xt+alamb)
		yt=0.25_dp*(yt+alamb)
		ave=THIRD*(xt+yt+yt)
		s=(yt-ave)/ave
		if (abs(s) <= ERRTOL) exit
	end do
	rc_s=w*(1.0_dp+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
	END FUNCTION rc_s


	FUNCTION rc_v(x,y)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x,y
	REAL(dp), DIMENSION(size(x)) :: rc_v
	REAL(dp), PARAMETER :: ERRTOL=0.04_dp,TINY=1.69e-38_dp,&
		SQRTNY=1.3e-19_dp,BIG=3.0e37_dp,TNBG=TINY*BIG,&
		COMP1=2.236_dp/SQRTNY,COMP2=TNBG*TNBG/25.0_dp,&
		THIRD=1.0_dp/3.0_dp,&
		C1=0.3_dp,C2=1.0_dp/7.0_dp,C3=0.375_dp,C4=9.0_dp/22.0_dp
	REAL(dp), DIMENSION(size(x)) :: alamb,ave,s,w,xt,yt
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(x),size(y),'rc_v')
	call assert( (/all(x >= 0.0),all(y /= 0.0),all(x+abs(y) >= TINY), &
		all(x+abs(y) <= BIG),all(y >= -COMP1 .or. x <= 0.0  &
		.or. x >= COMP2) /),'rc_v')
	where (y > 0.0)
		xt=x
		yt=y
		w=1.0
	elsewhere
		xt=x-y
		yt=-y
		w=sqrt(x)/sqrt(xt)
	end where
	converged=.false.
	do
		where (.not. converged)
			alamb=2.0_dp*sqrt(xt)*sqrt(yt)+yt
			xt=0.25_dp*(xt+alamb)
			yt=0.25_dp*(yt+alamb)
			ave=THIRD*(xt+yt+yt)
			s=(yt-ave)/ave
			converged = (abs(s) <= ERRTOL)
		end where
		if (all(converged)) exit
	end do
	rc_v=w*(1.0_dp+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
	END FUNCTION rc_v
