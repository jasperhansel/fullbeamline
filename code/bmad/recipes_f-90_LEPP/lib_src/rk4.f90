	SUBROUTINE rk4(y,dydx,x,h,yout,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: y,dydx
	REAL(dp), INTENT(IN) :: x,h
	REAL(dp), DIMENSION(:), INTENT(OUT) :: yout
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: ndum
	REAL(dp) :: h6,hh,xh
	REAL(dp), DIMENSION(size(y)) :: dym,dyt,yt
	ndum=assert_eq(size(y),size(dydx),size(yout),'rk4')
	hh=h*0.5_dp
	h6=h/6.0_dp
	xh=x+hh
	yt=y+hh*dydx
	call derivs(xh,yt,dyt)
	yt=y+hh*dyt
	call derivs(xh,yt,dym)
	yt=y+h*dym
	dym=dyt+dym
	call derivs(x+h,yt,dyt)
	yout=y+h6*(dydx+dyt+2.0_dp*dym)
	END SUBROUTINE rk4
