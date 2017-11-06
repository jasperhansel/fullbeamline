	FUNCTION dawson_s(x)
	USE nrtype; USE nrutil, ONLY : arth,geop
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: dawson_s
	INTEGER(I4B), PARAMETER :: NMAX=6
	REAL(dp), PARAMETER :: H=0.4_dp,A1=2.0_dp/3.0_dp,A2=0.4_dp,&
		A3=2.0_dp/7.0_dp
	INTEGER(I4B) :: i,n0
	REAL(dp) :: ec,x2,xp,xx
	REAL(dp), DIMENSION(NMAX) :: d1,d2,e1
	REAL(dp), DIMENSION(NMAX), SAVE :: c=(/ (0.0_dp,i=1,NMAX) /)
	if (c(1) == 0.0) c(1:NMAX)=exp(-(arth(1,2,NMAX)*H)**2)
	if (abs(x) < 0.2_dp) then
		x2=x**2
		dawson_s=x*(1.0_dp-A1*x2*(1.0_dp-A2*x2*(1.0_dp-A3*x2)))
	else
		xx=abs(x)
		n0=2*nint(0.5_dp*xx/H)
		xp=xx-real(n0,dp)*H
		ec=exp(2.0_dp*xp*H)
		d1=arth(n0+1,2,NMAX)
		d2=arth(n0-1,-2,NMAX)
		e1=geop(ec,ec**2,NMAX)
		dawson_s=0.5641895835477563_dp*sign(exp(-xp**2),x)*&
			sum(c*(e1/d1+1.0_dp/(d2*e1)))
	end if
	END FUNCTION dawson_s

	FUNCTION dawson_v(x)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x
	REAL(dp), DIMENSION(size(x)) :: dawson_v
	INTEGER(I4B), PARAMETER :: NMAX=6
	REAL(dp), PARAMETER :: H=0.4_dp,A1=2.0_dp/3.0_dp,A2=0.4_dp,&
		A3=2.0_dp/7.0_dp
	INTEGER(I4B) :: i,n
	REAL(dp), DIMENSION(size(x)) :: x2
	REAL(dp), DIMENSION(NMAX), SAVE :: c=(/ (0.0_dp,i=1,NMAX) /)
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	if (c(1) == 0.0) c(1:NMAX)=exp(-(arth(1,2,NMAX)*H)**2)
	mask = (abs(x) >= 0.2_dp)
	dawson_v=dawsonseries_v(x,mask)
	where (.not. mask)
		x2=x**2
		dawson_v=x*(1.0_dp-A1*x2*(1.0_dp-A2*x2*(1.0_dp-A3*x2)))
	end where
	CONTAINS
!BL
	FUNCTION dawsonseries_v(xin,mask)
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: xin
	LOGICAL(LGT), DIMENSION(size(xin)), INTENT(IN) :: mask
	REAL(dp), DIMENSION(size(xin)) :: dawsonseries_v
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: n0
	REAL(dp), DIMENSION(:), ALLOCATABLE :: d1,d2,e1,e2,sm,xp,xx,x
	n=count(mask)
	if (n == 0) RETURN
	allocate(n0(n),d1(n),d2(n),e1(n),e2(n),sm(n),xp(n),xx(n),x(n))
	x=pack(xin,mask)
	xx=abs(x)
	n0=2*nint(0.5_dp*xx/H)
	xp=xx-real(n0,dp)*H
	e1=exp(2.0_dp*xp*H)
	e2=e1**2
	d1=n0+1.0_dp
	d2=d1-2.0_dp
	sm=0.0
	do i=1,NMAX
		sm=sm+c(i)*(e1/d1+1.0_dp/(d2*e1))
		d1=d1+2.0_dp
		d2=d2-2.0_dp
		e1=e2*e1
	end do
	sm=0.5641895835477563_dp*sign(exp(-xp**2),x)*sm
	dawsonseries_v=unpack(sm,mask,0.0_dp)
	deallocate(n0,d1,d2,e1,e2,sm,xp,xx)
	END FUNCTION dawsonseries_v
	END FUNCTION dawson_v
