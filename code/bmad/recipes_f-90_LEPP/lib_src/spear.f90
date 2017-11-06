	SUBROUTINE spear(data1,data2,d,zd,probd,rs,probrs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : betai,erfcc,sort2
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(dp), INTENT(OUT) :: d,zd,probd,rs,probrs
	INTEGER(I4B) :: n
	REAL(dp) :: aved,df,en,en3n,fac,sf,sg,t,vard
	REAL(dp), DIMENSION(size(data1)) :: wksp1,wksp2
	n=assert_eq(size(data1),size(data2),'spear')
	wksp1(:)=data1(:)
	wksp2(:)=data2(:)
	call sort2(wksp1,wksp2)
	call crank(wksp1,sf)
	call sort2(wksp2,wksp1)
	call crank(wksp2,sg)
	wksp1(:)=wksp1(:)-wksp2(:)
	d=dot_product(wksp1,wksp1)
	en=n
	en3n=en**3-en
	aved=en3n/6.0_dp-(sf+sg)/12.0_dp
	fac=(1.0_dp-sf/en3n)*(1.0_dp-sg/en3n)
	vard=((en-1.0_dp)*en**2*(en+1.0_dp)**2/36.0_dp)*fac
	zd=(d-aved)/sqrt(vard)
	probd=erfcc(abs(zd)/SQRT2)
	rs=(1.0_dp-(6.0_dp/en3n)*(d+(sf+sg)/12.0_dp))/sqrt(fac)
	fac=(1.0_dp+rs)*(1.0_dp-rs)
	if (fac > 0.0) then
		t=rs*sqrt((en-2.0_dp)/fac)
		df=en-2.0_dp
		probrs=betai(0.5_dp*df,0.5_dp,df/(df+t**2))
	else
		probrs=0.0
	end if
	CONTAINS
!BL
	SUBROUTINE crank(w,s)
	USE nrtype; USE nrutil, ONLY : arth,array_copy
	IMPLICIT NONE
	REAL(dp), INTENT(OUT) :: s
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: w
	INTEGER(I4B) :: i,n,ndum,nties
	INTEGER(I4B), DIMENSION(size(w)) :: tstart,tend,tie,idx
	n=size(w)
	idx(:)=arth(1,1,n)
	tie(:)=merge(1,0,w == eoshift(w,-1))
	tie(1)=0
	w(:)=idx(:)
	if (all(tie == 0)) then
		s=0.0
		RETURN
	end if
	call array_copy(pack(idx(:),tie(:)<eoshift(tie(:),1)),tstart,nties,ndum)
	tend(1:nties)=pack(idx(:),tie(:)>eoshift(tie(:),1))
	do i=1,nties
		w(tstart(i):tend(i))=(tstart(i)+tend(i))/2.0_dp
	end do
	tend(1:nties)=tend(1:nties)-tstart(1:nties)+1
	s=sum(tend(1:nties)**3-tend(1:nties))
	END SUBROUTINE crank
	END SUBROUTINE spear
