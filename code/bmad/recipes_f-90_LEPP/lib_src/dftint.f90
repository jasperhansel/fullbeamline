	SUBROUTINE dftint(func,a,b,w,cosint,sinint)
	USE nrtype; USE nrutil, ONLY : arth
	USE nr, ONLY : dftcor,polint,realft
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b,w
	REAL(dp), INTENT(OUT) :: cosint,sinint
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: M=64,NDFT=1024,MPOL=6
	INTEGER(I4B) :: nn
	INTEGER(I4B), SAVE :: init=0
	INTEGER(I4B), DIMENSION(MPOL) :: nnmpol
	REAL(dp) :: c,cdft,cerr,corfac,corim,corre,en,s,sdft,serr
	REAL(dp), SAVE :: delta
	REAL(dp), DIMENSION(MPOL) :: cpol,spol,xpol
	REAL(dp), DIMENSION(NDFT), SAVE :: data
	REAL(dp), DIMENSION(8), SAVE :: endpts
	REAL(dp), SAVE :: aold=-1.0e30_dp,bold=-1.0e30_dp
	if (init /= 1 .or. a /= aold .or. b /= bold) then
		init=1
		aold=a
		bold=b
		delta=(b-a)/M
		data(1:M+1)=func(a+arth(0,1,M+1)*delta)
		data(M+2:NDFT)=0.0
		endpts(1:4)=data(1:4)
		endpts(5:8)=data(M-2:M+1)
		call realft(data(1:NDFT),1)
		data(2)=0.0
	end if
	en=w*delta*NDFT/TWOPI+1.0_dp
	nn=min(max(int(en-0.5_dp*MPOL+1.0_dp),1),NDFT/2-MPOL+1)
	nnmpol=arth(nn,1,MPOL)
	cpol(1:MPOL)=data(2*nnmpol(:)-1)
	spol(1:MPOL)=data(2*nnmpol(:))
	xpol(1:MPOL)=nnmpol(:)
	call polint(xpol,cpol,en,cdft,cerr)
	call polint(xpol,spol,en,sdft,serr)
	call dftcor(w,delta,a,b,endpts,corre,corim,corfac)
	cdft=cdft*corfac+corre
	sdft=sdft*corfac+corim
	c=delta*cos(w*a)
	s=delta*sin(w*a)
	cosint=c*cdft-s*sdft
	sinint=s*cdft+c*sdft
	END SUBROUTINE dftint
