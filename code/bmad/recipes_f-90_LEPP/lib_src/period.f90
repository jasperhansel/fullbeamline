	SUBROUTINE period(x,y,ofac,hifac,px,py,jmax,prob)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror
	USE nr, ONLY : avevar
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: jmax
	REAL(dp), INTENT(IN) :: ofac,hifac
	REAL(dp), INTENT(OUT) :: prob
	REAL(dp), DIMENSION(:), INTENT(IN) :: x,y
	REAL(dp), DIMENSION(:), POINTER :: px,py
	INTEGER(I4B) :: i,n,nout
	REAL(dp) :: ave,cwtau,effm,expy,pnow,sumc,sumcy,&
		sums,sumsh,sumsy,swtau,var,wtau,xave,xdif,xmax,xmin
	REAL(DP), DIMENSION(size(x)) :: tmp1,tmp2,wi,wpi,wpr,wr
	LOGICAL(LGT), SAVE :: init=.true.
	n=assert_eq(size(x),size(y),'period')
	if (init) then
		init=.false.
		nullify(px,py)
	else
		if (associated(px)) deallocate(px)
		if (associated(py)) deallocate(py)
	end if
	nout=0.5_dp*ofac*hifac*n
	allocate(px(nout),py(nout))
	call avevar(y(:),ave,var)
	if (var == 0.0) call nrerror('zero variance in period')
	xmax=maxval(x(:))
	xmin=minval(x(:))
	xdif=xmax-xmin
	xave=0.5_dp*(xmax+xmin)
	pnow=1.0_dp/(xdif*ofac)
	tmp1(:)=TWOPI_D*((x(:)-xave)*pnow)
	wpr(:)=-2.0_dp*sin(0.5_dp*tmp1)**2
	wpi(:)=sin(tmp1(:))
	wr(:)=cos(tmp1(:))
	wi(:)=wpi(:)
	do i=1,nout
		px(i)=pnow
		sumsh=dot_product(wi,wr)
		sumc=dot_product(wr(:)-wi(:),wr(:)+wi(:))
		wtau=0.5_dp*atan2(2.0_dp*sumsh,sumc)
		swtau=sin(wtau)
		cwtau=cos(wtau)
		tmp1(:)=wi(:)*cwtau-wr(:)*swtau
		tmp2(:)=wr(:)*cwtau+wi(:)*swtau
		sums=dot_product(tmp1,tmp1)
		sumc=dot_product(tmp2,tmp2)
		sumsy=dot_product(y(:)-ave,tmp1)
		sumcy=dot_product(y(:)-ave,tmp2)
		tmp1(:)=wr(:)
		wr(:)=(wr(:)*wpr(:)-wi(:)*wpi(:))+wr(:)
		wi(:)=(wi(:)*wpr(:)+tmp1(:)*wpi(:))+wi(:)
		py(i)=0.5_dp*(sumcy**2/sumc+sumsy**2/sums)/var
		pnow=pnow+1.0_dp/(ofac*xdif)
	end do
	jmax=imaxloc(py(1:nout))
	expy=exp(-py(jmax))
	effm=2.0_dp*nout/ofac
	prob=effm*expy
	if (prob > 0.01_dp) prob=1.0_dp-(1.0_dp-expy)**effm
	END SUBROUTINE period
