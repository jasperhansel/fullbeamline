	FUNCTION rd_s(x,y,z)
	USE nrtype; USE nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x,y,z
	REAL(dp) :: rd_s
	REAL(dp), PARAMETER :: ERRTOL=0.05_dp,TINY=1.0e-25_dp,BIG=4.5e21_dp,&
		C1=3.0_dp/14.0_dp,C2=1.0_dp/6.0_dp,C3=9.0_dp/22.0_dp,&
		C4=3.0_dp/26.0_dp,C5=0.25_dp*C3,C6=1.5_dp*C4
	REAL(dp) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
		ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
	call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, &
		'rd_s args')
	xt=x
	yt=y
	zt=z
	sum=0.0
	fac=1.0
	do
		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		sum=sum+fac/(sqrtz*(zt+alamb))
		fac=0.25_dp*fac
		xt=0.25_dp*(xt+alamb)
		yt=0.25_dp*(yt+alamb)
		zt=0.25_dp*(zt+alamb)
		ave=0.2_dp*(xt+yt+3.0_dp*zt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave
		if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
	end do
	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.0_dp*eb
	ee=ed+ec+ec
	rd_s=3.0_dp*sum+fac*(1.0_dp+ed*(-C1+C5*ed-C6*delz*ee)&
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
	END FUNCTION rd_s


	FUNCTION rd_v(x,y,z)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x,y,z
	REAL(dp), DIMENSION(size(x)) :: rd_v
	REAL(dp), PARAMETER :: ERRTOL=0.05_dp,TINY=1.0e-25_dp,BIG=4.5e21_dp,&
		C1=3.0_dp/14.0_dp,C2=1.0_dp/6.0_dp,C3=9.0_dp/22.0_dp,&
		C4=3.0_dp/26.0_dp,C5=0.25_dp*C3,C6=1.5_dp*C4
	REAL(dp), DIMENSION(size(x)) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
		ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(x),size(y),size(z),'rd_v')
	call assert(all(min(x,y) >= 0.0), all(min(x+y,z) >= TINY), &
		all(max(x,y,z) <= BIG), 'rd_v args')
	xt=x
	yt=y
	zt=z
	sum=0.0
	fac=1.0
	converged=.false.
	do
		where (.not. converged)
			sqrtx=sqrt(xt)
			sqrty=sqrt(yt)
			sqrtz=sqrt(zt)
			alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
			sum=sum+fac/(sqrtz*(zt+alamb))
			fac=0.25_dp*fac
			xt=0.25_dp*(xt+alamb)
			yt=0.25_dp*(yt+alamb)
			zt=0.25_dp*(zt+alamb)
			ave=0.2_dp*(xt+yt+3.0_dp*zt)
			delx=(ave-xt)/ave
			dely=(ave-yt)/ave
			delz=(ave-zt)/ave
			converged = (all(max(abs(delx),abs(dely),abs(delz)) <= ERRTOL))
		end where
		if (all(converged)) exit
	end do
	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.0_dp*eb
	ee=ed+ec+ec
	rd_v=3.0_dp*sum+fac*(1.0_dp+ed*(-C1+C5*ed-C6*delz*ee)&
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
	END FUNCTION rd_v
