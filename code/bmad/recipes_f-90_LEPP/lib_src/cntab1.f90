	SUBROUTINE cntab1(nn,chisq,df,prob,cramrv,ccc)
	USE nrtype; USE nrutil, ONLY : outerprod
	USE nr, ONLY : gammq
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
	REAL(dp), INTENT(OUT) :: chisq,df,prob,cramrv,ccc
	REAL(dp), PARAMETER :: TINY=1.0e-30_dp
	INTEGER(I4B) :: nni,nnj
	REAL(dp) :: sumn
	REAL(dp), DIMENSION(size(nn,1)) :: sumi
	REAL(dp), DIMENSION(size(nn,2)) :: sumj
	REAL(dp), DIMENSION(size(nn,1),size(nn,2)) :: expctd
	sumi(:)=sum(nn(:,:),dim=2)
	sumj(:)=sum(nn(:,:),dim=1)
	sumn=sum(sumi(:))
	nni=size(sumi)-count(sumi(:) == 0.0)
	nnj=size(sumj)-count(sumj(:) == 0.0)
	df=nni*nnj-nni-nnj+1
	expctd(:,:)=outerprod(sumi(:),sumj(:))/sumn
	chisq=sum((nn(:,:)-expctd(:,:))**2/(expctd(:,:)+TINY))
	prob=gammq(0.5_dp*df,0.5_dp*chisq)
	cramrv=sqrt(chisq/(sumn*min(nni-1,nnj-1)))
	ccc=sqrt(chisq/(chisq+sumn))
	END SUBROUTINE cntab1
