	FUNCTION ellf_s(phi,ak)
	USE nrtype
	USE nr, ONLY : rf
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: phi,ak
	REAL(dp) :: ellf_s
	REAL(dp) :: s
	s=sin(phi)
	ellf_s=s*rf(cos(phi)**2,(1.0_dp-s*ak)*(1.0_dp+s*ak),1.0_dp)
	END FUNCTION ellf_s


	FUNCTION ellf_v(phi,ak)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : rf
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: phi,ak
	REAL(dp), DIMENSION(size(phi)) :: ellf_v
	REAL(dp), DIMENSION(size(phi)) :: s
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(phi),size(ak),'ellf_v')
	s=sin(phi)
	ellf_v=s*rf(cos(phi)**2,(1.0_dp-s*ak)*(1.0_dp+s*ak),&
		spread(1.0_dp,1,size(phi)))
	END FUNCTION ellf_v
