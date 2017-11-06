	SUBROUTINE newt(x,check)
	USE nrtype; USE nrutil, ONLY : nrerror,vabs
	USE nr, ONLY : fdjac,lnsrch,lubksb,ludcmp
	USE fminln
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: x
	LOGICAL(LGT), INTENT(OUT) :: check
	INTEGER(I4B), PARAMETER :: MAXITS=200
	REAL(dp), PARAMETER :: TOLF=1.0e-4_dp,TOLMIN=1.0e-6_dp,TOLX=epsilon(x),&
		STPMX=100.0
	INTEGER(I4B) :: its
	INTEGER(I4B), DIMENSION(size(x)) :: indx
	REAL(dp) :: d,f,fold,stpmax
	REAL(dp), DIMENSION(size(x)) :: g,p,xold
	REAL(dp), DIMENSION(size(x)), TARGET :: fvec
	REAL(dp), DIMENSION(size(x),size(x)) :: fjac
	fmin_fvecp=>fvec
	f=fmin(x)
	if (maxval(abs(fvec(:))) < 0.01_dp*TOLF) then
		check=.false.
		RETURN
	end if
	stpmax=STPMX*max(vabs(x(:)),real(size(x),dp))
	do its=1,MAXITS
		call fdjac(x,fvec,fjac)
		g(:)=matmul(fvec(:),fjac(:,:))
		xold(:)=x(:)
		fold=f
		p(:)=-fvec(:)
		call ludcmp(fjac,indx,d)
		call lubksb(fjac,indx,p)
		call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin)
		if (maxval(abs(fvec(:))) < TOLF) then
			check=.false.
			RETURN
		end if
		if (check) then
			check=(maxval(abs(g(:))*max(abs(x(:)),1.0_dp) / &
				max(f,0.5_dp*size(x))) < TOLMIN)
			RETURN
		end if
		if (maxval(abs(x(:)-xold(:))/max(abs(x(:)),1.0_dp)) < TOLX) &
			RETURN
	end do
	call nrerror('MAXITS exceeded in newt')
	END SUBROUTINE newt
