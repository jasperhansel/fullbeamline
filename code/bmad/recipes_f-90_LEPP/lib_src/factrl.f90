	FUNCTION factrl_s(n)
	USE nrtype; USE nrutil, ONLY : arth,assert,cumprod
	USE nr, ONLY : gammln
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(dp) :: factrl_s
	INTEGER(I4B), SAVE :: ntop=0
	INTEGER(I4B), PARAMETER :: NMAX=32
	REAL(dp), DIMENSION(NMAX), SAVE :: a
	call assert(n >= 0, 'factrl_s arg')
	if (n < ntop) then
		factrl_s=a(n+1)
	else if (n < NMAX) then
		ntop=NMAX
		a(1)=1.0
		a(2:NMAX)=cumprod(arth(1.0_dp,1.0_dp,NMAX-1))
		factrl_s=a(n+1)
	else
		factrl_s=exp(gammln(n+1.0_dp))
	end if
	END FUNCTION factrl_s

	FUNCTION factrl_v(n)
	USE nrtype; USE nrutil, ONLY : arth,assert,cumprod
	USE nr, ONLY : gammln
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
	REAL(dp), DIMENSION(size(n)) :: factrl_v
	LOGICAL(LGT), DIMENSION(size(n)) :: mask
	INTEGER(I4B), SAVE :: ntop=0
	INTEGER(I4B), PARAMETER :: NMAX=32
	REAL(dp), DIMENSION(NMAX), SAVE :: a
	call assert(all(n >= 0), 'factrl_v arg')
	if (ntop == 0) then
		ntop=NMAX
		a(1)=1.0
		a(2:NMAX)=cumprod(arth(1.0_dp,1.0_dp,NMAX-1))
	end if
	mask = (n >= NMAX)
	factrl_v=unpack(exp(gammln(pack(n,mask)+1.0_dp)),mask,0.0_dp)
	where (.not. mask) factrl_v=a(n+1)
	END FUNCTION factrl_v
