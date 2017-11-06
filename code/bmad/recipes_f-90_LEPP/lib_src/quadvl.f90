	SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
	USE nrtype
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x,y
	REAL(dp), INTENT(OUT) :: fa,fb,fc,fd
	REAL(dp) :: qa,qb,qc,qd
	qa=min(2.0_dp,max(0.0_dp,1.0_dp-x))
	qb=min(2.0_dp,max(0.0_dp,1.0_dp-y))
	qc=min(2.0_dp,max(0.0_dp,x+1.0_dp))
	qd=min(2.0_dp,max(0.0_dp,y+1.0_dp))
	fa=0.25_dp*qa*qb
	fb=0.25_dp*qb*qc
	fc=0.25_dp*qc*qd
	fd=0.25_dp*qd*qa
	END SUBROUTINE quadvl
