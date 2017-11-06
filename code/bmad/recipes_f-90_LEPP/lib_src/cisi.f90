	SUBROUTINE cisi(x,ci,si)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp), INTENT(OUT) :: ci,si
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(dp), PARAMETER :: EPS=epsilon(x),FPMIN=4.0_dp*tiny(x),&
		BIG=huge(x)*EPS,TMIN=2.0
	INTEGER(I4B) :: i,k
	REAL(dp) :: a,err,fact,sign,sum,sumc,sums,t,term
	COMPLEX(dpc) :: h,b,c,d,del
	LOGICAL(LGT) :: odd
	t=abs(x)
	if (t == 0.0) then
		si=0.0
		ci=-BIG
		RETURN
	end if
	if (t > TMIN) then
		b=cmplx(1.0_dp,t,kind=dpc)
		c=BIG
		d=1.0_dp/b
		h=d
		do i=2,MAXIT
			a=-(i-1)**2
			b=b+2.0_dp
			d=1.0_dp/(a*d+b)
			c=b+a/c
			del=c*d
			h=h*del
			if (absc(del-1.0_dp) <= EPS) exit
		end do
		if (i > MAXIT) call nrerror('continued fraction failed in cisi')
		h=cmplx(cos(t),-sin(t),kind=dpc)*h
		ci=-real(h)
		si=PIO2+aimag(h)
	else
		if (t < sqrt(FPMIN)) then
			sumc=0.0
			sums=t
		else
			sum=0.0
			sums=0.0
			sumc=0.0
			sign=1.0
			fact=1.0
			odd=.true.
			do k=1,MAXIT
				fact=fact*t/k
				term=fact/k
				sum=sum+sign*term
				err=term/abs(sum)
				if (odd) then
					sign=-sign
					sums=sum
					sum=sumc
				else
					sumc=sum
					sum=sums
				end if
				if (err < EPS) exit
				odd=.not. odd
			end do
			if (k > MAXIT) call nrerror('MAXIT exceeded in cisi')
		end if
		si=sums
		ci=sumc+log(t)+EULER
	end if
	if (x < 0.0) si=-si
	CONTAINS
!BL
	FUNCTION absc(z)
	IMPLICIT NONE
	COMPLEX(dpc), INTENT(IN) :: z
	REAL(dp) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
	END SUBROUTINE cisi
