	FUNCTION fleg(x,nl)
	USE nrtype
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: nl
	REAL(dp), DIMENSION(nl) :: fleg
	INTEGER(I4B) :: j
	REAL(dp) :: d,f1,f2,twox
	fleg(1)=1.0
	fleg(2)=x
	if (nl > 2) then
		twox=2.0_dp*x
		f2=x
		d=1.0
		do j=3,nl
			f1=d
			f2=f2+twox
			d=d+1.0_dp
			fleg(j)=(f2*fleg(j-1)-f1*fleg(j-2))/d
		end do
	end if
	END FUNCTION fleg
