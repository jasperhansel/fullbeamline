	FUNCTION betacf_s(a,b,x)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b,x
	REAL(dp) :: betacf_s
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(dp), PARAMETER :: EPS=epsilon(x), FPMIN=tiny(x)/EPS
	REAL(dp) :: aa,c,d,del,h,qab,qam,qap
	INTEGER(I4B) :: m,m2
	qab=a+b
	qap=a+1.0_dp
	qam=a-1.0_dp
	c=1.0
	d=1.0_dp-qab*x/qap
	if (abs(d) < FPMIN) d=FPMIN
	d=1.0_dp/d
	h=d
	do m=1,MAXIT
		m2=2*m
		aa=m*(b-m)*x/((qam+m2)*(a+m2))
		d=1.0_dp+aa*d
		if (abs(d) < FPMIN) d=FPMIN
		c=1.0_dp+aa/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_dp/d
		h=h*d*c
		aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
		d=1.0_dp+aa*d
		if (abs(d) < FPMIN) d=FPMIN
		c=1.0_dp+aa/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_dp/d
		del=d*c
		h=h*del
		if (abs(del-1.0_dp) <= EPS) exit
	end do
	if (m > MAXIT)&
		call nrerror('a or b too big, or MAXIT too small in betacf_s')
	betacf_s=h
	END FUNCTION betacf_s


	FUNCTION betacf_v(a,b,x)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: a,b,x
	REAL(dp), DIMENSION(size(x)) :: betacf_v
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(dp), PARAMETER :: EPS=epsilon(x), FPMIN=tiny(x)/EPS
	REAL(dp), DIMENSION(size(x)) :: aa,c,d,del,h,qab,qam,qap
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	INTEGER(I4B) :: m
	INTEGER(I4B), DIMENSION(size(x)) :: m2
	m=assert_eq(size(a),size(b),size(x),'betacf_v')
	qab=a+b
	qap=a+1.0_dp
	qam=a-1.0_dp
	c=1.0
	d=1.0_dp-qab*x/qap
	where (abs(d) < FPMIN) d=FPMIN
	d=1.0_dp/d
	h=d
	converged=.false.
	do m=1,MAXIT
		where (.not. converged)
			m2=2*m
			aa=m*(b-m)*x/((qam+m2)*(a+m2))
			d=1.0_dp+aa*d
			d=merge(FPMIN,d, abs(d)<FPMIN )
			c=1.0_dp+aa/c
			c=merge(FPMIN,c, abs(c)<FPMIN )
			d=1.0_dp/d
			h=h*d*c
			aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
			d=1.0_dp+aa*d
			d=merge(FPMIN,d, abs(d)<FPMIN )
			c=1.0_dp+aa/c
			c=merge(FPMIN,c, abs(c)<FPMIN )
			d=1.0_dp/d
			del=d*c
			h=h*del
			converged = (abs(del-1.0_dp) <= EPS)
		end where
		if (all(converged)) exit
	end do
	if (m > MAXIT)&
		call nrerror('a or b too big, or MAXIT too small in betacf_v')
	betacf_v=h
	END FUNCTION betacf_v
