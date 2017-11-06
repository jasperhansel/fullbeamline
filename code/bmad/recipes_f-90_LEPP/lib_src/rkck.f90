	SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: y,dydx
	REAL(dp), INTENT(IN) :: x,h
	REAL(dp), DIMENSION(:), INTENT(OUT) :: yout,yerr
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
	REAL(dp), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
	REAL(dp), PARAMETER :: A2=0.2_dp,A3=0.3_dp,A4=0.6_dp,A5=1.0_dp,&
		A6=0.875_dp,B21=0.2_dp,B31=3.0_dp/40.0_dp,B32=9.0_dp/40.0_dp,&
		B41=0.3_dp,B42=-0.9_dp,B43=1.2_dp,B51=-11.0_dp/54.0_dp,&
		B52=2.5_dp,B53=-70.0_dp/27.0_dp,B54=35.0_dp/27.0_dp,&
		B61=1631.0_dp/55296.0_dp,B62=175.0_dp/512.0_dp,&
		B63=575.0_dp/13824.0_dp,B64=44275.0_dp/110592.0_dp,&
		B65=253.0_dp/4096.0_dp,C1=37.0_dp/378.0_dp,&
		C3=250.0_dp/621.0_dp,C4=125.0_dp/594.0_dp,&
		C6=512.0_dp/1771.0_dp,DC1=C1-2825.0_dp/27648.0_dp,&
		DC3=C3-18575.0_dp/48384.0_dp,DC4=C4-13525.0_dp/55296.0_dp,&
		DC5=-277.0_dp/14336.0_dp,DC6=C6-0.25_dp
	ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
	ytemp=y+B21*h*dydx
	call derivs(x+A2*h,ytemp,ak2)
	ytemp=y+h*(B31*dydx+B32*ak2)
	call derivs(x+A3*h,ytemp,ak3)
	ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
	call derivs(x+A4*h,ytemp,ak4)
	ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
	call derivs(x+A5*h,ytemp,ak5)
	ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
	call derivs(x+A6*h,ytemp,ak6)
	yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
	yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
	END SUBROUTINE rkck
