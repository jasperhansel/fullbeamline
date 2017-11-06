	FUNCTION chder(a,b,c)
	USE nrtype; USE nrutil, ONLY : arth,cumsum
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp), DIMENSION(:), INTENT(IN) :: c
	REAL(dp), DIMENSION(size(c)) :: chder
	INTEGER(I4B) :: n
	REAL(dp) :: con
	REAL(dp), DIMENSION(size(c)) :: temp
	n=size(c)
	temp(1)=0.0
	temp(2:n)=2.0_dp*arth(n-1,-1,n-1)*c(n:2:-1)
	chder(n:1:-2)=cumsum(temp(1:n:2))
	chder(n-1:1:-2)=cumsum(temp(2:n:2))
	con=2.0_dp/(b-a)
	chder=chder*con
	END FUNCTION chder
