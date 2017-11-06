	SUBROUTINE memcof(data,xms,d)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(OUT) :: xms
	REAL(dp), DIMENSION(:), INTENT(IN) :: data
	REAL(dp), DIMENSION(:), INTENT(OUT) :: d
	INTEGER(I4B) :: k,m,n
	REAL(dp) :: denom,pneum
	REAL(dp), DIMENSION(size(data)) :: wk1,wk2,wktmp
	REAL(dp), DIMENSION(size(d)) :: wkm
	m=size(d)
	n=size(data)
	xms=dot_product(data,data)/n
	wk1(1:n-1)=data(1:n-1)
	wk2(1:n-1)=data(2:n)
	do k=1,m
		pneum=dot_product(wk1(1:n-k),wk2(1:n-k))
		denom=dot_product(wk1(1:n-k),wk1(1:n-k))+ &
			dot_product(wk2(1:n-k),wk2(1:n-k))
		d(k)=2.0_dp*pneum/denom
		xms=xms*(1.0_dp-d(k)**2)
		d(1:k-1)=wkm(1:k-1)-d(k)*wkm(k-1:1:-1)
		if (k == m) RETURN
		wkm(1:k)=d(1:k)
		wktmp(2:n-k)=wk1(2:n-k)
		wk1(1:n-k-1)=wk1(1:n-k-1)-wkm(k)*wk2(1:n-k-1)
		wk2(1:n-k-1)=wk2(2:n-k)-wkm(k)*wktmp(2:n-k)
	end do
	call nrerror('never get here in memcof')
	END SUBROUTINE memcof
