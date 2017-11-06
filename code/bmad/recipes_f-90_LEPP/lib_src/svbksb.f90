	SUBROUTINE svbksb_dp(u,w,v,b,x)
	USE nrtype; USE nrutil, ONLY : assert_eq
	REAL(dp), DIMENSION(:,:), INTENT(IN) :: u,v
	REAL(dp), DIMENSION(:), INTENT(IN) :: w,b
	REAL(dp), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B) :: mdum,ndum
	REAL(dp), DIMENSION(size(x)) :: tmp
	mdum=assert_eq(size(u,1),size(b),'svbksb_dp: mdum')
	ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),&
		'svbksb_dp: ndum')
	where (w /= 0.0)
		tmp=matmul(b,u)/w
	elsewhere
		tmp=0.0
	end where
	x=matmul(v,tmp)
	END SUBROUTINE svbksb_dp

