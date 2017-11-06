	SUBROUTINE flmoon(n,nph,jd,frac)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n,nph
	INTEGER(I4B), INTENT(OUT) :: jd
	REAL(dp), INTENT(OUT) :: frac
	REAL(dp), PARAMETER :: RAD=PI/180.0_dp
	INTEGER(I4B) :: i
	REAL(dp) :: am,as,c,t,t2,xtra
	c=n+nph/4.0_dp
	t=c/1236.85_dp
	t2=t**2
	as=359.2242_dp+29.105356_dp*c
	am=306.0253_dp+385.816918_dp*c+0.010730_dp*t2
	jd=2415020+28*n+7*nph
	xtra=0.75933_dp+1.53058868_dp*c+(1.178e-4_dp-1.55e-7_dp*t)*t2
	select case(nph)
		case(0,2)
			xtra=xtra+(0.1734_dp-3.93e-4_dp*t)*sin(RAD*as)-0.4068_dp*sin(RAD*am)
		case(1,3)
			xtra=xtra+(0.1721_dp-4.0e-4_dp*t)*sin(RAD*as)-0.6280_dp*sin(RAD*am)
		case default
			call nrerror('flmoon: nph is unknown')
	end select
	i=int(merge(xtra,xtra-1.0_dp, xtra >= 0.0))
	jd=jd+i
	frac=xtra-i
	END SUBROUTINE flmoon
