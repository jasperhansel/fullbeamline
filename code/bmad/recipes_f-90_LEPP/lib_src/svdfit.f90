	SUBROUTINE svdfit(x,y,sig,a,v,w,chisq,funcs)
	USE nrtype; USE nrutil, ONLY : assert_eq,vabs
	USE nr, ONLY : svbksb,svdcmp
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x,y,sig
	REAL(dp), DIMENSION(:), INTENT(OUT) :: a,w
	REAL(dp), DIMENSION(:,:), INTENT(OUT) :: v
	REAL(dp), INTENT(OUT) :: chisq
	INTERFACE
		FUNCTION funcs(x,n)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: n
		REAL(dp), DIMENSION(n) :: funcs
		END FUNCTION funcs
	END INTERFACE
	REAL(dp), PARAMETER :: TOL=1.0e-5_dp
	INTEGER(I4B) :: i,ma,n
	REAL(dp), DIMENSION(size(x)) :: b,sigi
	REAL(dp), DIMENSION(size(x),size(a)) :: u,usav
	n=assert_eq(size(x),size(y),size(sig),'svdfit: n')
	ma=assert_eq(size(a),size(v,1),size(v,2),size(w),'svdfit: ma')
	sigi=1.0_dp/sig
	b=y*sigi
	do i=1,n
		usav(i,:)=funcs(x(i),ma)
	end do
	u=usav*spread(sigi,dim=2,ncopies=ma)
	usav=u
	call svdcmp(u,w,v)
	where (w < TOL*maxval(w)) w=0.0
	call svbksb(u,w,v,b,a)
	chisq=vabs(matmul(usav,a)-b)**2
	END SUBROUTINE svdfit
