	SUBROUTINE amoeba(p,y,ftol,func,iter)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(dp), INTENT(IN) :: ftol
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: y
	REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=5000
	REAL(dp), PARAMETER :: TINY=1.0e-10
	INTEGER(I4B) :: ihi,ndim
	REAL(dp), DIMENSION(size(p,2)) :: psum
	call amoeba_private
	CONTAINS
!BL
	SUBROUTINE amoeba_private
	IMPLICIT NONE
	INTEGER(I4B) :: i,ilo,inhi
	REAL(dp) :: rtol,ysave,ytry,ytmp
	ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
	iter=0
	psum(:)=sum(p(:,:),dim=1)
	do
		ilo=iminloc(y(:))
		ihi=imaxloc(y(:))
		ytmp=y(ihi)
		y(ihi)=y(ilo)
		inhi=imaxloc(y(:))
		y(ihi)=ytmp
		rtol=2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
		if (rtol < ftol) then
			call swap(y(1),y(ilo))
			call swap(p(1,:),p(ilo,:))
			RETURN
		end if
		if (iter >= ITMAX) call nrerror('ITMAX exceeded in amoeba')
		ytry=amotry(-1.0_dp)
		iter=iter+1
		if (ytry <= y(ilo)) then
			ytry=amotry(2.0_dp)
			iter=iter+1
		else if (ytry >= y(inhi)) then
			ysave=y(ihi)
			ytry=amotry(0.5_dp)
			iter=iter+1
			if (ytry >= ysave) then
				p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
				do i=1,ndim+1
					if (i /= ilo) y(i)=func(p(i,:))
				end do
				iter=iter+ndim
				psum(:)=sum(p(:,:),dim=1)
			end if
		end if
	end do
	END SUBROUTINE amoeba_private
!BL
	FUNCTION amotry(fac)
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: fac
	REAL(dp) :: amotry
	REAL(dp) :: fac1,fac2,ytry
	REAL(dp), DIMENSION(size(p,2)) :: ptry
	fac1=(1.0_dp-fac)/ndim
	fac2=fac1-fac
	ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
	ytry=func(ptry)
	if (ytry < y(ihi)) then
		y(ihi)=ytry
		psum(:)=psum(:)-p(ihi,:)+ptry(:)
		p(ihi,:)=ptry(:)
	end if
	amotry=ytry
	END FUNCTION amotry
	END SUBROUTINE amoeba
