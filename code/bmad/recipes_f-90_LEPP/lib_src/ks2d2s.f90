	SUBROUTINE ks2d2s(x1,y1,x2,y2,d,prob)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : pearsn,probks,quadct
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
	REAL(dp), INTENT(OUT) :: d,prob
	INTEGER(I4B) :: j,n1,n2
	REAL(dp) :: d1,d2,dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,r2,rr,sqen
	n1=assert_eq(size(x1),size(y1),'ks2d2s: n1')
	n2=assert_eq(size(x2),size(y2),'ks2d2s: n2')
	d1=0.0
	do j=1,n1
		call quadct(x1(j),y1(j),x1,y1,fa,fb,fc,fd)
		call quadct(x1(j),y1(j),x2,y2,ga,gb,gc,gd)
		d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
	end do
	d2=0.0
	do j=1,n2
		call quadct(x2(j),y2(j),x1,y1,fa,fb,fc,fd)
		call quadct(x2(j),y2(j),x2,y2,ga,gb,gc,gd)
		d2=max(d2,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
	end do
	d=0.5_dp*(d1+d2)
	sqen=sqrt(real(n1,dp)*real(n2,dp)/real(n1+n2,dp))
	call pearsn(x1,y1,r1,dum,dumm)
	call pearsn(x2,y2,r2,dum,dumm)
	rr=sqrt(1.0_dp-0.5_dp*(r1**2+r2**2))
	prob=probks(d*sqen/(1.0_dp+rr*(0.25_dp-0.75_dp/sqen)))
	END SUBROUTINE ks2d2s
