	FUNCTION fredin(x,a,b,t,f,w,g,ak)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp), DIMENSION(:), INTENT(IN) :: x,t,f,w
	REAL(dp), DIMENSION(size(x)) :: fredin
	INTERFACE
		FUNCTION g(t)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: t
		REAL(dp), DIMENSION(size(t)) :: g
		END FUNCTION g
!BL
		FUNCTION ak(t,s)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: t,s
		REAL(dp), DIMENSION(size(t),size(s)) :: ak
		END FUNCTION ak
	END INTERFACE
	INTEGER(I4B) :: n
	n=assert_eq(size(f),size(t),size(w),'fredin')
	fredin=g(x)+matmul(ak(x,t),w*f)
	END FUNCTION fredin
