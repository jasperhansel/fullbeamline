	SUBROUTINE stiff(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq,diagadd,nrerror
	USE nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: y
	REAL(dp), DIMENSION(:), INTENT(IN) :: dydx,yscal
	REAL(dp), INTENT(INOUT) :: x
	REAL(dp), INTENT(IN) :: htry,eps
	REAL(dp), INTENT(OUT) :: hdid,hnext
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE jacobn(x,y,dfdx,dfdy)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dfdx
		REAL(dp), DIMENSION(:,:), INTENT(OUT) :: dfdy
		END SUBROUTINE jacobn
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXTRY=40
	REAL(dp), PARAMETER :: SAFETY=0.9_dp,GROW=1.5_dp,PGROW=-0.25_dp,&
		SHRNK=0.5_dp,PSHRNK=-1.0_dp/3.0_dp,ERRCON=0.1296_dp,&
		GAM=1.0_dp/2.0_dp,&
		A21=2.0_dp,A31=48.0_dp/25.0_dp,A32=6.0_dp/25.0_dp,C21=-8.0_dp,&
		C31=372.0_dp/25.0_dp,C32=12.0_dp/5.0_dp,&
		C41=-112.0_dp/125.0_dp,C42=-54.0_dp/125.0_dp,&
		C43=-2.0_dp/5.0_dp,B1=19.0_dp/9.0_dp,B2=1.0_dp/2.0_dp,&
		B3=25.0_dp/108.0_dp,B4=125.0_dp/108.0_dp,E1=17.0_dp/54.0_dp,&
		E2=7.0_dp/36.0_dp,E3=0.0_dp,E4=125.0_dp/108.0_dp,&
		C1X=1.0_dp/2.0_dp,C2X=-3.0_dp/2.0_dp,C3X=121.0_dp/50.0_dp,&
		C4X=29.0_dp/250.0_dp,A2X=1.0_dp,A3X=3.0_dp/5.0_dp
	INTEGER(I4B) :: jtry,ndum
	INTEGER(I4B), DIMENSION(size(y)) :: indx
	REAL(dp), DIMENSION(size(y)) :: dfdx,dytmp,err,g1,g2,g3,g4,ysav
	REAL(dp), DIMENSION(size(y),size(y)) :: a,dfdy
	REAL(dp) :: d,errmax,h,xsav
	ndum=assert_eq(size(y),size(dydx),size(yscal),'stiff')
	xsav=x
	ysav(:)=y(:)
	call jacobn(xsav,ysav,dfdx,dfdy)
	h=htry
	do jtry=1,MAXTRY
		a(:,:)=-dfdy(:,:)
		call diagadd(a,1.0_dp/(GAM*h))
		call ludcmp(a,indx,d)
		g1=dydx+h*C1X*dfdx
		call lubksb(a,indx,g1)
		y=ysav+A21*g1
		x=xsav+A2X*h
		call derivs(x,y,dytmp)
		g2=dytmp+h*C2X*dfdx+C21*g1/h
		call lubksb(a,indx,g2)
		y=ysav+A31*g1+A32*g2
		x=xsav+A3X*h
		call derivs(x,y,dytmp)
		g3=dytmp+h*C3X*dfdx+(C31*g1+C32*g2)/h
		call lubksb(a,indx,g3)
		g4=dytmp+h*C4X*dfdx+(C41*g1+C42*g2+C43*g3)/h
		call lubksb(a,indx,g4)
		y=ysav+B1*g1+B2*g2+B3*g3+B4*g4
		err=E1*g1+E2*g2+E3*g3+E4*g4
		x=xsav+h
		if (x == xsav) call &
			nrerror('stepsize not significant in stiff')
		errmax=maxval(abs(err/yscal))/eps
		if (errmax <= 1.0) then
			hdid=h
			hnext=merge(SAFETY*h*errmax**PGROW, GROW*h, &
				errmax > ERRCON)
			RETURN
		else
			hnext=SAFETY*h*errmax**PSHRNK
			h=sign(max(abs(hnext),SHRNK*abs(h)),h)
		end if
	end do
	call nrerror('exceeded MAXTRY in stiff')
	END SUBROUTINE stiff
