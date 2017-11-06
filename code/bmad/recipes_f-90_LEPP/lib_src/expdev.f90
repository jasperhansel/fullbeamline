	SUBROUTINE expdev_s(harvest)
	USE nrtype
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(dp), INTENT(OUT) :: harvest
	REAL(dp) :: dum
	call ran1(dum)
	harvest=-log(dum)
	END SUBROUTINE expdev_s

	SUBROUTINE expdev_v(harvest)
	USE nrtype
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(OUT) :: harvest
	REAL(dp), DIMENSION(size(harvest)) :: dum
	call ran1(dum)
	harvest=-log(dum)
	END SUBROUTINE expdev_v
