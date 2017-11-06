	SUBROUTINE frenel(x,s,c)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp), INTENT(OUT) :: s,c
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(dp), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x),BIG=huge(x)*EPS,&
		XMIN=1.5
	INTEGER(I4B) :: k,n
	REAL(dp) :: a,ax,fact,pix2,sign,sum,sumc,sums,term,test
	COMPLEX(dpc) :: b,cc,d,h,del,cs
	LOGICAL(LGT) :: odd
	ax=abs(x)
	if (ax < sqrt(FPMIN)) then
		s=0.0
		c=ax
	else if (ax <= XMIN) then
		sum=0.0
		sums=0.0
		sumc=ax
		sign=1.0
		fact=PIO2*ax*ax
		odd=.true.
		term=ax
		n=3
		do k=1,MAXIT
			term=term*fact/k
			sum=sum+sign*term/n
			test=abs(sum)*EPS
			if (odd) then
				sign=-sign
				sums=sum
				sum=sumc
			else
				sumc=sum
				sum=sums
			end if
			if (term < test) exit
			odd=.not. odd
			n=n+2
		end do
		if (k > MAXIT) call nrerror('frenel: series failed')
		s=sums
		c=sumc
	else
		pix2=PI*ax*ax
		b=cmplx(1.0_dp,-pix2,kind=dpc)
		cc=BIG
		d=1.0_dp/b
		h=d
		n=-1
		do k=2,MAXIT
			n=n+2
			a=-n*(n+1)
			b=b+4.0_dp
			d=1.0_dp/(a*d+b)
			cc=b+a/cc
			del=cc*d
			h=h*del
			if (absc(del-1.0_dp) <= EPS) exit
		end do
		if (k > MAXIT) call nrerror('cf failed in frenel')
		h=h*cmplx(ax,-ax,kind=dpc)
		cs=cmplx(0.5_dp,0.5_dp,kind=dpc)*(1.0_dp-&
			cmplx(cos(0.5_dp*pix2),sin(0.5_dp*pix2),kind=dpc)*h)
		c=real(cs)
		s=aimag(cs)
	end if
	if (x < 0.0) then
		c=-c
		s=-s
	end if
	CONTAINS
!BL
	FUNCTION absc(z)
	IMPLICIT NONE
	COMPLEX(dpc), INTENT(IN) :: z
	REAL(dp) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
	END SUBROUTINE frenel
