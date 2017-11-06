	FUNCTION predic(data,d,nfut)
	USE nrtype
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: data,d
	INTEGER(I4B), INTENT(IN) :: nfut
	REAL(dp), DIMENSION(nfut) :: predic
	INTEGER(I4B) :: j,ndata,m
	REAL(dp) :: discrp,sm
	REAL(dp), DIMENSION(size(d)) :: reg
	m=size(d)
	ndata=size(data)
	reg(1:m)=data(ndata:ndata+1-m:-1)
	do j=1,nfut
		discrp=0.0
		sm=discrp+dot_product(d,reg)
		reg=eoshift(reg,-1,sm)
		predic(j)=sm
	end do
	END FUNCTION predic
