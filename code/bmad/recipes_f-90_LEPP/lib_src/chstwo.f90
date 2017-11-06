	SUBROUTINE chstwo(bins1,bins2,knstrn,df,chsq,prob)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : gammq
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: knstrn
	REAL(dp), INTENT(OUT) :: df,chsq,prob
	REAL(dp), DIMENSION(:), INTENT(IN) :: bins1,bins2
	INTEGER(I4B) :: ndum
	LOGICAL(LGT), DIMENSION(size(bins1)) :: nzeromask
	ndum=assert_eq(size(bins1),size(bins2),'chstwo')
	nzeromask = bins1(:) /= 0.0 .or. bins2(:) /= 0.0
	chsq=sum((bins1(:)-bins2(:))**2/(bins1(:)+bins2(:)),mask=nzeromask)
	df=count(nzeromask)-knstrn
	prob=gammq(0.5_dp*df,0.5_dp*chsq)
	END SUBROUTINE chstwo
