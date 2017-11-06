MODULE rkdumb_path
	USE nrtype
	REAL(dp), DIMENSION(:), ALLOCATABLE:: xx
	REAL(dp), DIMENSION(:,:), ALLOCATABLE :: y
END MODULE rkdumb_path

	SUBROUTINE rkdumb(vstart,x1,x2,nstep,derivs)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : rk4
	USE rkdumb_path
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: vstart
	REAL(dp), INTENT(IN) :: x1,x2
	INTEGER(I4B), INTENT(IN) :: nstep
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: k
	REAL(dp) :: h,x
	REAL(dp), DIMENSION(size(vstart)) :: dv,v
	v(:)=vstart(:)
	if (allocated(xx)) deallocate(xx)
	if (allocated(y)) deallocate(y)
	allocate(xx(nstep+1))
	allocate(y(size(vstart),nstep+1))
	y(:,1)=v(:)
	xx(1)=x1
	x=x1
	h=(x2-x1)/nstep
	do k=1,nstep
		call derivs(x,v,dv)
		call rk4(v,dv,x,h,v,derivs)
		if (x+h == x) call nrerror('stepsize not significant in rkdumb')
		x=x+h
		xx(k+1)=x
		y(:,k+1)=v(:)
	end do
	END SUBROUTINE rkdumb
