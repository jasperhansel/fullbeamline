	FUNCTION qromo(func,a,b,choose)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp) :: qromo
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		END FUNCTION func
!BL
		SUBROUTINE choose(funk,aa,bb,s,n)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: aa,bb
		REAL(dp), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			USE nrtype
			IMPLICIT NONE
			REAL(dp), DIMENSION(:), INTENT(IN) :: x
			REAL(dp), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE choose
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=14,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(dp), PARAMETER :: EPS=1.0e-6
	REAL(dp), DIMENSION(JMAXP) :: h,s
	REAL(dp) :: dqromo
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call choose(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromo,dqromo)
			if (abs(dqromo) <= EPS*abs(qromo)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=h(j)/9.0_dp
	end do
	call nrerror('qromo: too many steps')
	END FUNCTION qromo
