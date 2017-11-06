	MODULE quad3d_qgaus_mod
	USE nrtype
	PRIVATE
	PUBLIC quad3d_qgaus
	REAL(dp) :: xsav,ysav
	INTERFACE
		FUNCTION func(x,y,z)
		USE nrtype
		REAL(dp), INTENT(IN) :: x,y
		REAL(dp), DIMENSION(:), INTENT(IN) :: z
		REAL(dp), DIMENSION(size(z)) :: func
		END FUNCTION func
!BL
		FUNCTION y1(x)
		USE nrtype
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: y1
		END FUNCTION y1
!BL
		FUNCTION y2(x)
		USE nrtype
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: y2
		END FUNCTION y2
!BL
		FUNCTION z1(x,y)
		USE nrtype
		REAL(dp), INTENT(IN) :: x,y
		REAL(dp) :: z1
		END FUNCTION z1
!BL
		FUNCTION z2(x,y)
		USE nrtype
		REAL(dp), INTENT(IN) :: x,y
		REAL(dp) :: z2
		END FUNCTION z2
	END INTERFACE
	CONTAINS
!BL
	FUNCTION h(x)
	REAL(dp), DIMENSION(:), INTENT(IN) :: x
	REAL(dp), DIMENSION(size(x)) :: h
	INTEGER(I4B) :: i
	do i=1,size(x)
		xsav=x(i)
		h(i)=qgaus(g,y1(xsav),y2(xsav))
	end do
	END FUNCTION h
!BL
	FUNCTION g(y)
	REAL(dp), DIMENSION(:), INTENT(IN) :: y
	REAL(dp), DIMENSION(size(y)) :: g
	INTEGER(I4B) :: j
	do j=1,size(y)
		ysav=y(j)
		g(j)=qgaus(f,z1(xsav,ysav),z2(xsav,ysav))
	end do
	END FUNCTION g
!BL
	FUNCTION f(z)
	REAL(dp), DIMENSION(:), INTENT(IN) :: z
	REAL(dp), DIMENSION(size(z)) :: f
	f=func(xsav,ysav,z)
	END FUNCTION f
!BL
	RECURSIVE FUNCTION qgaus(func,a,b)
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp) :: qgaus
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(dp) :: xm,xr
	REAL(dp), DIMENSION(5) :: dx, w = (/ 0.2955242247_dp,0.2692667193_dp,&
		0.2190863625_dp,0.1494513491_dp,0.0666713443_dp /),&
		x = (/ 0.1488743389_dp,0.4333953941_dp,0.6794095682_dp,&
		0.8650633666_dp,0.9739065285_dp /)
	xm=0.5_dp*(b+a)
	xr=0.5_dp*(b-a)
	dx(:)=xr*x(:)
	qgaus=xr*sum(w(:)*(func(xm+dx)+func(xm-dx)))
	END FUNCTION qgaus
!BL
	SUBROUTINE quad3d_qgaus(x1,x2,ss)
	REAL(dp), INTENT(IN) :: x1,x2
	REAL(dp), INTENT(OUT) :: ss
	ss=qgaus(h,x1,x2)
	END SUBROUTINE quad3d_qgaus
	END MODULE quad3d_qgaus_mod

	MODULE quad3d_qromb_mod
	USE nrtype
	PRIVATE
	PUBLIC quad3d_qromb
	REAL(dp) :: xsav,ysav
	INTERFACE
		FUNCTION func(x,y,z)
		USE nrtype
		REAL(dp), INTENT(IN) :: x,y
		REAL(dp), DIMENSION(:), INTENT(IN) :: z
		REAL(dp), DIMENSION(size(z)) :: func
		END FUNCTION func
!BL
		FUNCTION y1(x)
		USE nrtype
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: y1
		END FUNCTION y1
!BL
		FUNCTION y2(x)
		USE nrtype
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: y2
		END FUNCTION y2
!BL
		FUNCTION z1(x,y)
		USE nrtype
		REAL(dp), INTENT(IN) :: x,y
		REAL(dp) :: z1
		END FUNCTION z1
!BL
		FUNCTION z2(x,y)
		USE nrtype
		REAL(dp), INTENT(IN) :: x,y
		REAL(dp) :: z2
		END FUNCTION z2
	END INTERFACE
	CONTAINS
!BL
	FUNCTION h(x)
	REAL(dp), DIMENSION(:), INTENT(IN) :: x
	REAL(dp), DIMENSION(size(x)) :: h
	INTEGER(I4B) :: i
	do i=1,size(x)
		xsav=x(i)
		h(i)=qromb(g,y1(xsav),y2(xsav))
	end do
	END FUNCTION h
!BL
	FUNCTION g(y)
	REAL(dp), DIMENSION(:), INTENT(IN) :: y
	REAL(dp), DIMENSION(size(y)) :: g
	INTEGER(I4B) :: j
	do j=1,size(y)
		ysav=y(j)
		g(j)=qromb(f,z1(xsav,ysav),z2(xsav,ysav))
	end do
	END FUNCTION g
!BL
	FUNCTION f(z)
	REAL(dp), DIMENSION(:), INTENT(IN) :: z
	REAL(dp), DIMENSION(size(z)) :: f
	f=func(xsav,ysav,z)
	END FUNCTION f
!BL
	RECURSIVE FUNCTION qromb(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp) :: qromb
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(dp), PARAMETER :: EPS=3.0e-6_dp
	REAL(dp), DIMENSION(JMAXP) :: h,s
	REAL(dp) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromb,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_dp*h(j)
	end do
	call nrerror('qromb: too many steps')
	END FUNCTION qromb
!BL
	RECURSIVE SUBROUTINE trapzd(func,a,b,s,n)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(dp) :: del,fsum
	INTEGER(I4B) :: it
	if (n == 1) then
		s=0.5_dp*(b-a)*sum(func( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func(arth(a+0.5_dp*del,del,it)))
		s=0.5_dp*(s+del*fsum)
	end if
	END SUBROUTINE trapzd
!BL
	SUBROUTINE quad3d_qromb(x1,x2,ss)
	REAL(dp), INTENT(IN) :: x1,x2
	REAL(dp), INTENT(OUT) :: ss
	ss=qromb(h,x1,x2)
	END SUBROUTINE quad3d_qromb
	END MODULE quad3d_qromb_mod
