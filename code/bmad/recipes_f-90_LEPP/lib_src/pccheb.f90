	FUNCTION pccheb(d)
	USE nrtype; USE nrutil, ONLY : arth,cumprod,geop
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: d
	REAL(dp), DIMENSION(size(d)) :: pccheb
	INTEGER(I4B) :: k,n
	REAL(dp), DIMENSION(size(d)) :: denom,numer,pow
	n=size(d)
	pccheb(1)=2.0_dp*d(1)
	pow=geop(1.0_dp,2.0_dp,n)
	numer(1)=1.0
	denom(1)=1.0
	denom(2:(n+3)/2)=cumprod(arth(1.0_dp,1.0_dp,(n+1)/2))
	pccheb(2:n)=0.0
	do k=2,n
		numer(2:(k+3)/2)=cumprod(arth(k-1.0_dp,-1.0_dp,(k+1)/2))
		pccheb(k:1:-2)=pccheb(k:1:-2)+&
			d(k)/pow(k-1)*numer(1:(k+1)/2)/denom(1:(k+1)/2)
	end do
	END FUNCTION pccheb
