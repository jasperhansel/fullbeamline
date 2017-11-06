	SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
	USE nrtype
	USE nr, ONLY : chebev
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP), INTENT(OUT) :: gam1,gam2,gampl,gammi
	INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
	REAL(dp) :: xx
	REAL(dp), DIMENSION(7) :: c1=(/-1.142022680371168_dp,&
		6.5165112670737e-3_dp,3.087090173086e-4_dp,-3.4706269649e-6_dp,&
		6.9437664e-9_dp,3.67795e-11_dp,-1.356e-13_dp/)
	REAL(dp), DIMENSION(8) :: c2=(/1.843740587300905_dp,&
		-7.68528408447867e-2_dp,1.2719271366546e-3_dp,&
		-4.9717367042e-6_dp, -3.31261198e-8_dp,2.423096e-10_dp,&
		-1.702e-13_dp,-1.49e-15_dp/)
	xx=8.0_dp*x*x-1.0_dp
	gam1=chebev(-1.0_dp,1.0_dp,c1(1:NUSE1),xx)
	gam2=chebev(-1.0_dp,1.0_dp,c2(1:NUSE2),xx)
	gampl=gam2-x*gam1
	gammi=gam2+x*gam1
	END SUBROUTINE beschb_s


	SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
	USE nrtype
	USE nr, ONLY : chebev
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
	INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
	REAL(dp), DIMENSION(size(x)) :: xx
	REAL(dp), DIMENSION(7) :: c1=(/-1.142022680371168_dp,&
		6.5165112670737e-3_dp,3.087090173086e-4_dp,-3.4706269649e-6_dp,&
		6.9437664e-9_dp,3.67795e-11_dp,-1.356e-13_dp/)
	REAL(dp), DIMENSION(8) :: c2=(/1.843740587300905_dp,&
		-7.68528408447867e-2_dp,1.2719271366546e-3_dp,&
		-4.9717367042e-6_dp, -3.31261198e-8_dp,2.423096e-10_dp,&
		-1.702e-13_dp,-1.49e-15_dp/)
	xx=8.0_dp*x*x-1.0_dp
	gam1=chebev(-1.0_dp,1.0_dp,c1(1:NUSE1),xx)
	gam2=chebev(-1.0_dp,1.0_dp,c2(1:NUSE2),xx)
	gampl=gam2-x*gam1
	gammi=gam2+x*gam1
	END SUBROUTINE beschb_v
