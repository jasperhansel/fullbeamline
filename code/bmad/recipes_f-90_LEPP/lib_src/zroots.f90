	SUBROUTINE zroots(a,roots,polish)
	USE nrtype; USE nrutil, ONLY : assert_eq,poly_term
	USE nr, ONLY : laguer,indexx
	IMPLICIT NONE
	COMPLEX(dpc), DIMENSION(:), INTENT(IN) :: a
	COMPLEX(dpc), DIMENSION(:), INTENT(OUT) :: roots
	LOGICAL(LGT), INTENT(IN) :: polish
	REAL(dp), PARAMETER :: EPS=1.0e-6_dp
	INTEGER(I4B) :: j,its,m
	INTEGER(I4B), DIMENSION(size(roots)) :: indx
	COMPLEX(dpc) :: x
	COMPLEX(dpc), DIMENSION(size(a)) :: ad
	m=assert_eq(size(roots),size(a)-1,'zroots')
	ad(:)=a(:)
	do j=m,1,-1
		x=cmplx(0.0_dp,kind=dpc)
		call laguer(ad(1:j+1),x,its)
		if (abs(aimag(x)) <= 2.0_dp*EPS**2*abs(real(x))) &
			x=cmplx(real(x),kind=dpc)
		roots(j)=x
		ad(j:1:-1)=poly_term(ad(j+1:2:-1),x)
	end do
	if (polish) then
		do j=1,m
			call laguer(a(:),roots(j),its)
		end do
	end if
	call indexx(real(roots),indx)
	roots=roots(indx)
	END SUBROUTINE zroots
