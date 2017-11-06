	SUBROUTINE gaucof(a,b,amu0,x,w)
	USE nrtype; USE nrutil, ONLY : assert_eq,unit_matrix
	USE nr, ONLY : eigsrt,tqli
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: a,b
	REAL(dp), INTENT(IN) :: amu0
	REAL(dp), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(dp), DIMENSION(size(a),size(a)) :: z
	INTEGER(I4B) :: n
	n=assert_eq(size(a),size(b),size(x),size(w),'gaucof')
	b(2:n)=sqrt(b(2:n))
	call unit_matrix(z)
	call tqli(a,b,z)
	call eigsrt(a,z)
	x=a
	w=amu0*z(1,:)**2
	END SUBROUTINE gaucof
