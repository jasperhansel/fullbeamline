	SUBROUTINE fixrts(d)
	USE nrtype
	USE nr, ONLY : zroots
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: d
	INTEGER(I4B) :: i,m
	LOGICAL(LGT) :: polish
	COMPLEX(dpc), DIMENSION(size(d)+1) :: a
	COMPLEX(dpc), DIMENSION(size(d)) :: roots
	m=size(d)
	a(m+1)=cmplx(1.0_dp,kind=dpc)
	a(m:1:-1)=cmplx(-d(1:m),kind=dpc)
	polish=.true.
	call zroots(a(1:m+1),roots,polish)
	where (abs(roots) > 1.0) roots=1.0_dp/conjg(roots)
	a(1)=-roots(1)
	a(2:m+1)=cmplx(1.0_dp,kind=dpc)
	do i=2,m
		a(2:i)=a(1:i-1)-roots(i)*a(2:i)
		a(1)=-roots(i)*a(1)
	end do
	d(m:1:-1)=-real(a(1:m))
	END SUBROUTINE fixrts
