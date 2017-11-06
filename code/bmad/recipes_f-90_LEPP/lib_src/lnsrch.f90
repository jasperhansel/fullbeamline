	SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,vabs
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: xold,g
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
	REAL(dp), INTENT(IN) :: fold,stpmax
	REAL(dp), DIMENSION(:), INTENT(OUT) :: x
	REAL(dp), INTENT(OUT) :: f
	LOGICAL(LGT), INTENT(OUT) :: check
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp) :: func
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
	END INTERFACE
	REAL(dp), PARAMETER :: ALF=1.0e-4_dp,TOLX=epsilon(x)
	INTEGER(I4B) :: ndum
	REAL(dp) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
	ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
	check=.false.
	pabs=vabs(p(:))
	if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
	slope=dot_product(g,p)
	if (slope >= 0.0) call nrerror('roundoff problem in lnsrch')
	alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_dp))
	alam=1.0
	do
		x(:)=xold(:)+alam*p(:)
		f=func(x)
		if (alam < alamin) then
			x(:)=xold(:)
			check=.true.
			RETURN
		else if (f <= fold+ALF*alam*slope) then
			RETURN
		else
			if (alam == 1.0) then
				tmplam=-slope/(2.0_dp*(f-fold-slope))
			else
				rhs1=f-fold-alam*slope
				rhs2=f2-fold-alam2*slope
				a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
				b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
					(alam-alam2)
				if (a == 0.0) then
					tmplam=-slope/(2.0_dp*b)
				else
					disc=b*b-3.0_dp*a*slope
					if (disc < 0.0) then
						tmplam=0.5_dp*alam
					else if (b <= 0.0) then
						tmplam=(-b+sqrt(disc))/(3.0_dp*a)
					else
						tmplam=-slope/(b+sqrt(disc))
					end if
				end if
				if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
			end if
		end if
		alam2=alam
		f2=f
		alam=max(tmplam,0.1_dp*alam)
	end do
	END SUBROUTINE lnsrch
