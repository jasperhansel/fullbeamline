	SUBROUTINE pcshft(a,b,d)
	USE nrtype; USE nrutil, ONLY : geop
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: d
	INTEGER(I4B) :: j,n
	REAL(dp), DIMENSION(size(d)) :: dd
	REAL(dp) :: x
	n=size(d)
	dd=d*geop(1.0_dp,2.0_dp/(b-a),n)
	x=-0.5_dp*(a+b)
	d(1)=dd(n)
	d(2:n)=0.0
	do j=n-1,1,-1
		d(2:n+1-j)=d(2:n+1-j)*x+d(1:n-j)
		d(1)=d(1)*x+dd(j)
	end do
	END SUBROUTINE pcshft
