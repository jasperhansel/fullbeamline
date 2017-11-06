	FUNCTION elle_s(phi,ak)
	USE nrtype
	USE nr, ONLY : rd,rf
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: phi,ak
	REAL(dp) :: elle_s
	REAL(dp) :: cc,q,s
	s=sin(phi)
	cc=cos(phi)**2
	q=(1.0_dp-s*ak)*(1.0_dp+s*ak)
	elle_s=s*(rf(cc,q,1.0_dp)-((s*ak)**2)*rd(cc,q,1.0_dp)/3.0_dp)
	END FUNCTION elle_s


	FUNCTION elle_v(phi,ak)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : rd,rf
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: phi,ak
	REAL(dp), DIMENSION(size(phi)) :: elle_v
	REAL(dp), DIMENSION(size(phi)) :: cc,q,s
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(phi),size(ak),'elle_v')
	s=sin(phi)
	cc=cos(phi)**2
	q=(1.0_dp-s*ak)*(1.0_dp+s*ak)
	elle_v=s*(rf(cc,q,spread(1.0_dp,1,size(phi)))-((s*ak)**2)*&
		rd(cc,q,spread(1.0_dp,1,size(phi)))/3.0_dp)
	END FUNCTION elle_v
