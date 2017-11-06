	FUNCTION chebev_s(a,b,c,x)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b,x
	REAL(dp), DIMENSION(:), INTENT(IN) :: c
	REAL(dp) :: chebev_s
	INTEGER(I4B) :: j,m
	REAL(dp) :: d,dd,sv,y,y2
	if ((x-a)*(x-b) > 0.0) call nrerror('x not in range in chebev_s')
	m=size(c)
	d=0.0
	dd=0.0
	y=(2.0_dp*x-a-b)/(b-a)
	y2=2.0_dp*y
	do j=m,2,-1
		sv=d
		d=y2*d-dd+c(j)
		dd=sv
	end do
	chebev_s=y*d-dd+0.5_dp*c(1)
	END FUNCTION chebev_s


	FUNCTION chebev_v(a,b,c,x)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp), DIMENSION(:), INTENT(IN) :: c,x
	REAL(dp), DIMENSION(size(x)) :: chebev_v
	INTEGER(I4B) :: j,m
	REAL(dp), DIMENSION(size(x)) :: d,dd,sv,y,y2
	if (any((x-a)*(x-b) > 0.0)) call nrerror('x not in range in chebev_v')
	m=size(c)
	d=0.0
	dd=0.0
	y=(2.0_dp*x-a-b)/(b-a)
	y2=2.0_dp*y
	do j=m,2,-1
		sv=d
		d=y2*d-dd+c(j)
		dd=sv
	end do
	chebev_v=y*d-dd+0.5_dp*c(1)
	END FUNCTION chebev_v
