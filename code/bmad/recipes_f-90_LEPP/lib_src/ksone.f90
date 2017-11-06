	SUBROUTINE ksone(data,func,d,prob)
	USE nrtype; USE nrutil, ONLY : arth
	USE nr, ONLY : probks,sort
	IMPLICIT NONE
	REAL(dp), INTENT(OUT) :: d,prob
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: data
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B) :: n
	REAL(dp) :: en
	REAL(dp), DIMENSION(size(data)) :: fvals
	REAL(dp), DIMENSION(size(data)+1) :: temp
	call sort(data)
	n=size(data)
	en=n
	fvals(:)=func(data(:))
	temp=arth(0,1,n+1)/en
	d=maxval(max(abs(temp(1:n)-fvals(:)), &
		abs(temp(2:n+1)-fvals(:))))
	en=sqrt(en)
	prob=probks((en+0.12_dp+0.11_dp/en)*d)
	END SUBROUTINE ksone
