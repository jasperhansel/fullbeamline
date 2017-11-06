	SUBROUTINE midsql(funk,aa,bb,s,n)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: aa,bb
	REAL(dp), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION funk(x)
		USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: funk
		END FUNCTION funk
	END INTERFACE
	REAL(dp) :: a,b,del
	INTEGER(I4B) :: it
	REAL(dp), DIMENSION(2*3**(n-2)) :: x
	b=sqrt(bb-aa)
	a=0.0
	if (n == 1) then
		s=(b-a)*sum(func( (/0.5_dp*(a+b)/) ))
	else
		it=3**(n-2)
		del=(b-a)/(3.0_dp*it)
		x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it)
		x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
		s=s/3.0_dp+del*sum(func(x))
	end if
	CONTAINS
!BL
		FUNCTION func(x)
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		func=2.0_dp*x*funk(aa+x**2)
		END FUNCTION func
	END SUBROUTINE midsql
