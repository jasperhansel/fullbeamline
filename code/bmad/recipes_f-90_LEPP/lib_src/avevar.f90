	SUBROUTINE avevar(data,ave,var)
	USE nrtype
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: data
	REAL(dp), INTENT(OUT) :: ave,var
	INTEGER(I4B) :: n
	REAL(dp), DIMENSION(size(data)) :: s
	n=size(data)
	ave=sum(data(:))/n
	s(:)=data(:)-ave
	var=dot_product(s,s)
	var=(var-sum(s)**2/n)/(n-1)
	END SUBROUTINE avevar
