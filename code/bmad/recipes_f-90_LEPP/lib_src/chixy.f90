MODULE chixyfit
	USE nrtype; USE nrutil, ONLY : nrerror
	REAL(dp), DIMENSION(:), POINTER :: xxp,yyp,sxp,syp,wwp
	REAL(dp) :: aa,offs
CONTAINS
!BL
	FUNCTION chixy(bang)
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: bang
	REAL(dp) :: chixy
	REAL(dp), PARAMETER :: BIG=1.0e30_dp
	REAL(dp) :: avex,avey,sumw,b
	if (.not. associated(wwp)) call nrerror("chixy: bad pointers")
	b=tan(bang)
	wwp(:)=(b*sxp(:))**2+syp(:)**2
	where (wwp(:) < 1.0/BIG)
		wwp(:)=BIG
	elsewhere
		wwp(:)=1.0_dp/wwp(:)
	end where
	sumw=sum(wwp)
	avex=dot_product(wwp,xxp)/sumw
	avey=dot_product(wwp,yyp)/sumw
	aa=avey-b*avex
	chixy=sum(wwp(:)*(yyp(:)-aa-b*xxp(:))**2)-offs
	END FUNCTION chixy
END MODULE chixyfit
