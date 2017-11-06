	FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : spline,splint
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x1a,x2a
	REAL(dp), DIMENSION(:,:), INTENT(IN) :: ya,y2a
	REAL(dp), INTENT(IN) :: x1,x2
	REAL(dp) :: splin2
	INTEGER(I4B) :: j,m,ndum
	REAL(dp), DIMENSION(size(x1a)) :: yytmp,y2tmp2
	m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m')
	ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum')
	do j=1,m
		yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
	end do
	call spline(x1a,yytmp,1.0e30_dp,1.0e30_dp,y2tmp2)
	splin2=splint(x1a,yytmp,y2tmp2,x1)
	END FUNCTION splin2
