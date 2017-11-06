	FUNCTION rtnewt(funcd,x1,x2,xacc)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x1,x2,xacc
	REAL(dp) :: rtnewt
	INTERFACE
		SUBROUTINE funcd(x,fval,fderiv)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), INTENT(OUT) :: fval,fderiv
		END SUBROUTINE funcd
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=20
	INTEGER(I4B) :: j
	REAL(dp) :: df,dx,f
	rtnewt=0.5_dp*(x1+x2)
	do j=1,MAXIT
		call funcd(rtnewt,f,df)
		dx=f/df
		rtnewt=rtnewt-dx
		if ((x1-rtnewt)*(rtnewt-x2) < 0.0)&
			call nrerror('rtnewt: values jumped out of brackets')
		if (abs(dx) < xacc) RETURN
	end do
	call nrerror('rtnewt exceeded maximum iterations')
	END FUNCTION rtnewt
