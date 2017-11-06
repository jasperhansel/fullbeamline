	SUBROUTINE quadct(x,y,xx,yy,fa,fb,fc,fd)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x,y
	REAL(dp), DIMENSION(:), INTENT(IN) :: xx,yy
	REAL(dp), INTENT(OUT) :: fa,fb,fc,fd
	INTEGER(I4B) :: na,nb,nc,nd,nn
	REAL(dp) :: ff
	nn=assert_eq(size(xx),size(yy),'quadct')
	na=count(yy(:) > y .and. xx(:) > x)
	nb=count(yy(:) > y .and. xx(:) <= x)
	nc=count(yy(:) <= y .and. xx(:) <= x)
	nd=nn-na-nb-nc
	ff=1.0_dp/nn
	fa=ff*na
	fb=ff*nb
	fc=ff*nc
	fd=ff*nd
	END SUBROUTINE quadct
