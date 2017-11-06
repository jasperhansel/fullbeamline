	SUBROUTINE fdjac(x,fvec,df)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: fvec
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: x
	REAL(dp), DIMENSION(:,:), INTENT(OUT) :: df
	INTERFACE
		FUNCTION funcv(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: funcv
		END FUNCTION funcv
	END INTERFACE
	REAL(dp), PARAMETER :: EPS=1.0e-4_dp
	INTEGER(I4B) :: j,n
	REAL(dp), DIMENSION(size(x)) :: xsav,xph,h
	n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
	xsav=x
	h=EPS*abs(xsav)
	where (h == 0.0) h=EPS
	xph=xsav+h
	h=xph-xsav
	do j=1,n
		x(j)=xph(j)
		df(:,j)=(funcv(x)-fvec(:))/h(j)
		x(j)=xsav(j)
	end do
	END SUBROUTINE fdjac
