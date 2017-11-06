	FUNCTION bessi_s(n,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessi0
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: bessi_s
	INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4B) :: j,m
	REAL(dp) :: bi,bim,bip,tox
	call assert(n >= 2, 'bessi_s args')
	bessi_s=0.0
	if (x*x <= 8.0_dp*tiny(x)) RETURN
	tox=2.0_dp/abs(x)
	bip=0.0
	bi=1.0
	m=2*((n+int(sqrt(real(IACC*n,dp)))))
	do j=m,1,-1
		bim=bip+j*tox*bi
		bip=bi
		bi=bim
		if (exponent(bi) > IEXP) then
			bessi_s=scale(bessi_s,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end if
		if (j == n) bessi_s=bip
	end do
	bessi_s=bessi_s*bessi0(x)/bi
	if (x < 0.0 .and. mod(n,2) == 1) bessi_s=-bessi_s
	END FUNCTION bessi_s


	FUNCTION bessi_v(n,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessi0
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(dp), DIMENSION(:), INTENT(IN) :: x
	REAL(dp), DIMENSION(size(x)) :: bessi_v
	INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4B) :: j,m
	REAL(dp), DIMENSION(size(x)) :: bi,bim,bip,tox
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	call assert(n >= 2, 'bessi_v args')
	bessi_v=0.0
	mask = (x <= 8.0_dp*tiny(x))
	tox=2.0_dp/merge(2.0_dp,abs(x),mask)
	bip=0.0
	bi=1.0_dp
	m=2*((n+int(sqrt(real(IACC*n,dp)))))
	do j=m,1,-1
		bim=bip+j*tox*bi
		bip=bi
		bi=bim
		where (exponent(bi) > IEXP)
			bessi_v=scale(bessi_v,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end where
		if (j == n) bessi_v=bip
	end do
	bessi_v=bessi_v*bessi0(x)/bi
	where (mask) bessi_v=0.0_dp
	where (x < 0.0 .and. mod(n,2) == 1) bessi_v=-bessi_v
	END FUNCTION bessi_v
