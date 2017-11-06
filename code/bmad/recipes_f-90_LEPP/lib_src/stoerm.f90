	SUBROUTINE stoerm(y,d2y,xs,htot,nstep,yout,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: y,d2y
	REAL(dp), INTENT(IN) :: xs,htot
	INTEGER(I4B), INTENT(IN) :: nstep
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
	INTEGER(I4B) :: neqn,neqn1,nn,nv
	REAL(dp) :: h,h2,halfh,x
	REAL(dp), DIMENSION(size(y)) :: ytemp
	nv=assert_eq(size(y),size(d2y),size(yout),'stoerm')
	neqn=nv/2
	neqn1=neqn+1
	h=htot/nstep
	halfh=0.5_dp*h
	ytemp(neqn1:nv)=h*(y(neqn1:nv)+halfh*d2y(1:neqn))
	ytemp(1:neqn)=y(1:neqn)+ytemp(neqn1:nv)
	x=xs+h
	call derivs(x,ytemp,yout)
	h2=h*h
	do nn=2,nstep
		ytemp(neqn1:nv)=ytemp(neqn1:nv)+h2*yout(1:neqn)
		ytemp(1:neqn)=ytemp(1:neqn)+ytemp(neqn1:nv)
		x=x+h
		call derivs(x,ytemp,yout)
	end do
	yout(neqn1:nv)=ytemp(neqn1:nv)/h+halfh*yout(1:neqn)
	yout(1:neqn)=ytemp(1:neqn)
	END SUBROUTINE stoerm
