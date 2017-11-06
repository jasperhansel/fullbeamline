	SUBROUTINE sphbes_s(n,x,sj,sy,sjp,syp)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessjy
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(dp), INTENT(IN) :: x
	REAL(dp), INTENT(OUT) :: sj,sy,sjp,syp
	REAL(dp), PARAMETER :: RTPIO2=1.253314137315500_dp
	REAL(dp) :: factor,order,rj,rjp,ry,ryp
	call assert(n >= 0, x > 0.0, 'sphbes_s args')
	order=n+0.5_dp
	call bessjy(x,order,rj,ry,rjp,ryp)
	factor=RTPIO2/sqrt(x)
	sj=factor*rj
	sy=factor*ry
	sjp=factor*rjp-sj/(2.0_dp*x)
	syp=factor*ryp-sy/(2.0_dp*x)
	END SUBROUTINE sphbes_s


	SUBROUTINE sphbes_v(n,x,sj,sy,sjp,syp)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessjy
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(dp), DIMENSION(:), INTENT(IN) :: x
	REAL(dp), DIMENSION(:), INTENT(OUT) :: sj,sy,sjp,syp
	REAL(dp), PARAMETER :: RTPIO2=1.253314137315500_dp
	REAL(dp) :: order
	REAL(dp), DIMENSION(size(x)) :: factor,rj,rjp,ry,ryp
	call assert(n >= 0,  all(x > 0.0), 'sphbes_v args')
	order=n+0.5_dp
	call bessjy(x,order,rj,ry,rjp,ryp)
	factor=RTPIO2/sqrt(x)
	sj=factor*rj
	sy=factor*ry
	sjp=factor*rjp-sj/(2.0_dp*x)
	syp=factor*ryp-sy/(2.0_dp*x)
	END SUBROUTINE sphbes_v
