	FUNCTION golden(ax,bx,cx,func,tol,xmin)
	USE nrtype
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: ax,bx,cx,tol
	REAL(dp), INTENT(OUT) :: xmin
	REAL(dp) :: golden
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: func
		END FUNCTION func
	END INTERFACE
	REAL(dp), PARAMETER :: R=0.61803399_dp,C=1.0_dp-R
	REAL(dp) :: f1,f2,x0,x1,x2,x3
	x0=ax
	x3=cx
	if (abs(cx-bx) > abs(bx-ax)) then
		x1=bx
		x2=bx+C*(cx-bx)
	else
		x2=bx
		x1=bx-C*(bx-ax)
	end if
	f1=func(x1)
	f2=func(x2)
	do
		if (abs(x3-x0) <= tol*(abs(x1)+abs(x2))) exit
		if (f2 < f1) then
			call shft3(x0,x1,x2,R*x2+C*x3)
			call shft2(f1,f2,func(x2))
		else
			call shft3(x3,x2,x1,R*x1+C*x0)
			call shft2(f2,f1,func(x1))
		end if
	end do
	if (f1 < f2) then
		golden=f1
		xmin=x1
	else
		golden=f2
		xmin=x2
	end if
	CONTAINS
!BL
	SUBROUTINE shft2(a,b,c)
	REAL(dp), INTENT(OUT) :: a
	REAL(dp), INTENT(INOUT) :: b
	REAL(dp), INTENT(IN) :: c
	a=b
	b=c
	END SUBROUTINE shft2
!BL
	SUBROUTINE shft3(a,b,c,d)
	REAL(dp), INTENT(OUT) :: a
	REAL(dp), INTENT(INOUT) :: b,c
	REAL(dp), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft3
	END FUNCTION golden
