	SUBROUTINE laguer(a,x,its)
	USE nrtype; USE nrutil, ONLY : nrerror,poly,poly_term
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: its
	COMPLEX(dpc), INTENT(INOUT) :: x
	COMPLEX(dpc), DIMENSION(:), INTENT(IN) :: a
	REAL(dp), PARAMETER :: EPS=epsilon(1.0_dp)
	INTEGER(I4B), PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
	INTEGER(I4B) :: iter,m
	REAL(dp) :: abx,abp,abm,err
	COMPLEX(dpc) :: dx,x1,f,g,h,sq,gp,gm,g2
	COMPLEX(dpc), DIMENSION(size(a)) :: b,d
	REAL(dp), DIMENSION(MR) :: frac = &
		(/ 0.5_dp,0.25_dp,0.75_dp,0.13_dp,0.38_dp,0.62_dp,0.88_dp,1.0_dp /)
	m=size(a)-1
	do iter=1,MAXIT
		its=iter
		abx=abs(x)
		b(m+1:1:-1)=poly_term(a(m+1:1:-1),x)
		d(m:1:-1)=poly_term(b(m+1:2:-1),x)
		f=poly(x,d(2:m))
		err=EPS*poly(abx,abs(b(1:m+1)))
		if (abs(b(1)) <= err) RETURN
		g=d(1)/b(1)
		g2=g*g
		h=g2-2.0_dp*f/b(1)
		sq=sqrt((m-1)*(m*h-g2))
		gp=g+sq
		gm=g-sq
		abp=abs(gp)
		abm=abs(gm)
		if (abp < abm) gp=gm
		if (max(abp,abm) > 0.0) then
			dx=m/gp
		else
			dx=exp(cmplx(log(1.0_dp+abx),iter,kind=dpc))
		end if
		x1=x-dx
		if (x == x1) RETURN
		if (mod(iter,MT) /= 0) then
			x=x1
		else
			x=x-dx*frac(iter/MT)
		end if
	end do
	call nrerror('laguer: too many iterations')
	END SUBROUTINE laguer
