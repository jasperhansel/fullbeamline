	FUNCTION rtbis(func,x1,x2,xacc)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x1,x2,xacc
	REAL(dp) :: rtbis
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=40
	INTEGER(I4B) :: j
	REAL(dp) :: dx,f,fmid,xmid
	fmid=func(x2)
	f=func(x1)
	if (f*fmid >= 0.0) call nrerror('rtbis: root must be bracketed')
	if (f < 0.0) then
		rtbis=x1
		dx=x2-x1
	else
		rtbis=x2
		dx=x1-x2
	end if
	do j=1,MAXIT
		dx=dx*0.5_dp
		xmid=rtbis+dx
		fmid=func(xmid)
		if (fmid <= 0.0) rtbis=xmid
		if (abs(dx) < xacc .or. fmid == 0.0) RETURN
	end do
	call nrerror('rtbis: too many bisections')
	END FUNCTION rtbis
