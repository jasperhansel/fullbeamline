	FUNCTION qromb(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint,trapzd
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp) :: qromb
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(dp), PARAMETER :: EPS=1.0e-6_dp
	REAL(dp), DIMENSION(JMAXP) :: h,s
	REAL(dp) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromb,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_dp*h(j)
	end do
	call nrerror('qromb: too many steps')
	END FUNCTION qromb
