	SUBROUTINE ttest(data1,data2,t,prob)
	USE nrtype
	USE nr, ONLY : avevar,betai
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(dp), INTENT(OUT) :: t,prob
	INTEGER(I4B) :: n1,n2
	REAL(dp) :: ave1,ave2,df,var,var1,var2
	n1=size(data1)
	n2=size(data2)
	call avevar(data1,ave1,var1)
	call avevar(data2,ave2,var2)
	df=n1+n2-2
	var=((n1-1)*var1+(n2-1)*var2)/df
	t=(ave1-ave2)/sqrt(var*(1.0_dp/n1+1.0_dp/n2))
	prob=betai(0.5_dp*df,0.5_dp,df/(df+t**2))
	END SUBROUTINE ttest
