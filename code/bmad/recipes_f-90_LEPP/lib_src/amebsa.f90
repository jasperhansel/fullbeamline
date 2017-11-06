	SUBROUTINE amebsa(p,y,pb,yb,ftol,func,iter,temptr)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,swap
	USE nr, ONLY : ran1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: iter
	REAL(dp), INTENT(INOUT) :: yb
	REAL(dp), INTENT(IN) :: ftol,temptr
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: y,pb
	REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: NMAX=200
	INTEGER(I4B) :: ihi,ndim
	REAL(dp) :: yhi
	REAL(dp), DIMENSION(size(p,2)) :: psum
	call amebsa_private
	CONTAINS
!BL
	SUBROUTINE amebsa_private
	INTEGER(I4B) :: i,ilo,inhi
	REAL(dp) :: rtol,ylo,ynhi,ysave,ytry
	REAL(dp), DIMENSION(size(y)) :: yt,harvest
	ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,size(pb),'amebsa')
	psum(:)=sum(p(:,:),dim=1)
	do
		call ran1(harvest)
		yt(:)=y(:)-temptr*log(harvest)
		ilo=iminloc(yt(:))
		ylo=yt(ilo)
		ihi=imaxloc(yt(:))
		yhi=yt(ihi)
		yt(ihi)=ylo
		inhi=imaxloc(yt(:))
		ynhi=yt(inhi)
		rtol=2.0_dp*abs(yhi-ylo)/(abs(yhi)+abs(ylo))
		if (rtol < ftol .or. iter < 0) then
			call swap(y(1),y(ilo))
			call swap(p(1,:),p(ilo,:))
			RETURN
		end if
		ytry=amotsa(-1.0_dp)
		iter=iter-1
		if (ytry <= ylo) then
			ytry=amotsa(2.0_dp)
			iter=iter-1
		else if (ytry >= ynhi) then
			ysave=yhi
			ytry=amotsa(0.5_dp)
			iter=iter-1
			if (ytry >= ysave) then
				p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
				do i=1,ndim+1
					if (i /= ilo) y(i)=func(p(i,:))
				end do
				iter=iter-ndim
				psum(:)=sum(p(:,:),dim=1)
			end if
		end if
	end do
	END SUBROUTINE amebsa_private
!BL
	FUNCTION amotsa(fac)
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: fac
	REAL(dp) :: amotsa
	REAL(dp) :: fac1,fac2,yflu,ytry,harv
	REAL(dp), DIMENSION(size(p,2)) :: ptry
	fac1=(1.0_dp-fac)/ndim
	fac2=fac1-fac
	ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
	ytry=func(ptry)
	if (ytry <= yb) then
		pb(:)=ptry(:)
		yb=ytry
	end if
	call ran1(harv)
	yflu=ytry+temptr*log(harv)
	if (yflu < yhi) then
		y(ihi)=ytry
		yhi=yflu
		psum(:)=psum(:)-p(ihi,:)+ptry(:)
		p(ihi,:)=ptry(:)
	end if
	amotsa=yflu
	END FUNCTION amotsa
	END SUBROUTINE amebsa
