	SUBROUTINE stifbs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,cumsum,iminloc,nrerror,&
		outerdiff,outerprod,upper_triangle
	USE nr, ONLY : simpr,pzextr
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: y
	REAL(dp), DIMENSION(:), INTENT(IN) :: dydx,yscal
	REAL(dp), INTENT(IN) :: htry,eps
	REAL(dp), INTENT(INOUT) :: x
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
	INTEGER(I4B), PARAMETER :: IMAX=8, KMAXX=IMAX-1
	REAL(dp), PARAMETER :: SAFE1=0.25_dp,SAFE2=0.7_dp,REDMAX=1.0e-5_dp,&
		REDMIN=0.7_dp,TINY=1.0e-30_dp,SCALMX=0.1_dp
	INTEGER(I4B) :: k,km,ndum
	INTEGER(I4B), DIMENSION(IMAX) :: nseq = (/ 2,6,10,14,22,34,50,70 /)
	INTEGER(I4B), SAVE :: kopt,kmax,nvold=-1
	REAL(dp), DIMENSION(KMAXX,KMAXX), SAVE :: alf
	REAL(dp), DIMENSION(KMAXX) :: err
	REAL(dp), DIMENSION(IMAX), SAVE :: a
	REAL(dp), SAVE :: epsold = -1.0
	REAL(dp) :: eps1,errmax,fact,h,red,scale,wrkmin,xest
	REAL(dp), SAVE :: xnew
	REAL(dp), DIMENSION(size(y)) :: dfdx,yerr,ysav,yseq
	REAL(dp), DIMENSION(size(y),size(y)) :: dfdy
	LOGICAL(LGT) :: reduct
	LOGICAL(LGT), SAVE :: first=.true.
	ndum=assert_eq(size(y),size(dydx),size(yscal),'stifbs')
	if (eps /= epsold .or. nvold /= size(y)) then
		hnext=-1.0e29_dp
		xnew=-1.0e29_dp
		eps1=SAFE1*eps
		a(:)=cumsum(nseq,1)
		where (upper_triangle(KMAXX,KMAXX)) alf=eps1** &
			(outerdiff(a(2:),a(2:))/outerprod(arth( &
			3.0_dp,2.0_dp,KMAXX),(a(2:)-a(1)+1.0_dp)))
		epsold=eps
		nvold=size(y)
		a(:)=cumsum(nseq,1+nvold)
		do kopt=2,KMAXX-1
			if (a(kopt+1) > a(kopt)*alf(kopt-1,kopt)) exit
		end do
		kmax=kopt
	end if
	h=htry
	ysav(:)=y(:)
	call jacobn(x,y,dfdx,dfdy)
	if (h /= hnext .or. x /= xnew) then
		first=.true.
		kopt=kmax
	end if
	reduct=.false.
	main_loop: do
		do k=1,kmax
			xnew=x+h
			if (xnew == x) call nrerror('step size underflow in stifbs')
			call simpr(ysav,dydx,dfdx,dfdy,x,h,nseq(k),yseq,derivs)
			xest=(h/nseq(k))**2
			call pzextr(k,xest,yseq,y,yerr)
			if (k /= 1) then
				errmax=maxval(abs(yerr(:)/yscal(:)))
				errmax=max(TINY,errmax)/eps
				km=k-1
				err(km)=(errmax/SAFE1)**(1.0_dp/(2*km+1))
			end if
			if (k /= 1 .and. (k >= kopt-1 .or. first)) then
				if (errmax < 1.0) exit main_loop
				if (k == kmax .or. k == kopt+1) then
					red=SAFE2/err(km)
					exit
				else if (k == kopt) then
					if (alf(kopt-1,kopt) < err(km)) then
						red=1.0_dp/err(km)
						exit
					end if
				else if (kopt == kmax) then
					if (alf(km,kmax-1) < err(km)) then
						red=alf(km,kmax-1)*SAFE2/err(km)
						exit
					end if
				else if (alf(km,kopt) < err(km)) then
					red=alf(km,kopt-1)/err(km)
					exit
				end if
			end if
		end do
		red=max(min(red,REDMIN),REDMAX)
		h=h*red
		reduct=.true.
	end do main_loop
	x=xnew
	hdid=h
	first=.false.
	kopt=1+iminloc(a(2:km+1)*max(err(1:km),SCALMX))
	scale=max(err(kopt-1),SCALMX)
	wrkmin=scale*a(kopt)
	hnext=h/scale
	if (kopt >= k .and. kopt /= kmax .and. .not. reduct) then
		fact=max(scale/alf(kopt-1,kopt),SCALMX)
		if (a(kopt+1)*fact <= wrkmin) then
			hnext=h/fact
			kopt=kopt+1
		end if
	end if
	END SUBROUTINE stifbs
