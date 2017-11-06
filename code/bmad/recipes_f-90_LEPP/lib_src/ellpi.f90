	FUNCTION ellpi_s(phi,en,ak)
	USE nrtype
	USE nr, ONLY : rf,rj
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: phi,en,ak
	REAL(dp) :: ellpi_s
	REAL(dp) :: cc,enss,q,s
	s=sin(phi)
	enss=en*s*s
	cc=cos(phi)**2
	q=(1.0_dp-s*ak)*(1.0_dp+s*ak)
	ellpi_s=s*(rf(cc,q,1.0_dp)-enss*rj(cc,q,1.0_dp,1.0_dp+enss)/3.0_dp)
	END FUNCTION ellpi_s


	FUNCTION ellpi_v(phi,en,ak)
	USE nrtype
	USE nr, ONLY : rf,rj
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: phi,en,ak
	REAL(dp), DIMENSION(size(phi)) :: ellpi_v
	REAL(dp), DIMENSION(size(phi)) :: cc,enss,q,s
	s=sin(phi)
	enss=en*s*s
	cc=cos(phi)**2
	q=(1.0_dp-s*ak)*(1.0_dp+s*ak)
	ellpi_v=s*(rf(cc,q,spread(1.0_dp,1,size(phi)))-enss*&
		rj(cc,q,spread(1.0_dp,1,size(phi)),1.0_dp+enss)/3.0_dp)
	END FUNCTION ellpi_v
