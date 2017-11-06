	FUNCTION qtrap(func,a,b, eps_in)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : trapzd
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp) :: qtrap, eps
  real(dp), optional :: eps_in
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20
	REAL(dp), PARAMETER :: EPS_std=1.0e-6_dp
	REAL(dp) :: olds
	INTEGER(I4B) :: j
  if (present(eps_in)) then
    eps = eps_in
  else
    eps = eps_std
  endif
	olds=0.0
	do j=1,JMAX
		call trapzd(func,a,b,qtrap,j)
		if (j > 5) then
			if (abs(qtrap-olds) <= EPS*abs(olds) .or. &
				(qtrap == 0.0 .and. olds == 0.0)) RETURN
		end if
		olds=qtrap
	end do
	call nrerror('qtrap: too many steps')
	END FUNCTION qtrap
