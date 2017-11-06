	SUBROUTINE gaulag(x,w,alf)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: alf
	REAL(dp), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-13_dp
	INTEGER(I4B) :: its,j,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(dp) :: anu
	REAL(dp), PARAMETER :: C1=9.084064e-01_dp,C2=5.214976e-02_dp,&
		C3=2.579930e-03_dp,C4=3.986126e-03_dp
	REAL(dp), DIMENSION(size(x)) :: rhs,r2,r3,theta
	REAL(DP), DIMENSION(size(x)) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION(size(x)) :: unfinished
	n=assert_eq(size(x),size(w),'gaulag')
	anu=4.0_dp*n+2.0_dp*alf+2.0_dp
	rhs=arth(4*n-1,-4,n)*PI/anu
	r3=rhs**(1.0_dp/3.0_dp)
	r2=r3**2
	theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
	z=anu*cos(theta)**2
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=1.0
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=((2.0_dp*j-1.0_dp+alf-z)*p2-(j-1.0_dp+alf)*p3)/j
			end where
		end do
		where (unfinished)
			pp=(n*p1-(n+alf)*p2)/z
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS*z)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gaulag')
	x=z
	w=-exp(gammln(alf+n)-gammln(real(n,dp)))/(pp*n*p2)
	END SUBROUTINE gaulag
