	FUNCTION qsimp(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : trapzd
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp) :: qsimp
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20
	REAL(dp), PARAMETER :: EPS=1.0e-6_dp
	INTEGER(I4B) :: j
	REAL(dp) :: os,ost,st
	ost=0.0
	os= 0.0
	do j=1,JMAX
		call trapzd(func,a,b,st,j)
		qsimp=(4.0_dp*st-ost)/3.0_dp
		if (j > 5) then
			if (abs(qsimp-os) < EPS*abs(os) .or. &
				(qsimp == 0.0 .and. os == 0.0)) RETURN
		end if
		os=qsimp
		ost=st
	end do
	call nrerror('qsimp: too many steps')
	END FUNCTION qsimp
