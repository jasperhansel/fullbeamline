	PROGRAM sfroid
	USE nrtype; USE nrutil, ONLY : arth
	USE nr, ONLY : plgndr,solvde
	USE sfroid_data
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NE=3,NB=1
	INTEGER(I4B) :: itmax
	INTEGER(I4B), DIMENSION(NE) :: indexv
	REAL(dp) :: conv,slowc
	REAL(dp), DIMENSION(M) :: deriv,fac1,fac2
	REAL(dp), DIMENSION(NE) :: scalv
	REAL(dp), DIMENSION(NE,M) :: y
	itmax=100
	conv=5.0e-6_dp
	slowc=1.0
	h=1.0_dp/(M-1)
	c2=0.0
	write(*,*) 'ENTER M,N'
	read(*,*) mm,n
	indexv(1:3)=merge( (/ 1, 2, 3 /), (/ 2, 1, 3 /), (mod(n+mm,2) == 1) )
	anorm=1.0
	if (mm /= 0) then
		anorm=(-0.5_dp)**mm*product(&
			arth(n+1,1,mm)*arth(real(n,dp),-1.0_dp,mm)/arth(1,1,mm))
	end if
	x(1:M-1)=arth(0,1,M-1)*h
	fac1(1:M-1)=1.0_dp-x(1:M-1)**2
	fac2(1:M-1)=fac1(1:M-1)**(-mm/2.0_dp)
	y(1,1:M-1)=plgndr(n,mm,x(1:M-1))*fac2(1:M-1)
	deriv(1:M-1)=-((n-mm+1)*plgndr(n+1,mm,x(1:M-1))-(n+1)*&
		x(1:M-1)*plgndr(n,mm,x(1:M-1)))/fac1(1:M-1)
	y(2,1:M-1)=mm*x(1:M-1)*y(1,1:M-1)/fac1(1:M-1)+deriv(1:M-1)*fac2(1:M-1)
	y(3,1:M-1)=n*(n+1)-mm*(mm+1)
	x(M)=1.0
	y(1,M)=anorm
	y(3,M)=n*(n+1)-mm*(mm+1)
	y(2,M)=(y(3,M)-c2)*y(1,M)/(2.0_dp*(mm+1.0_dp))
	scalv(1:3)=(/ abs(anorm), max(abs(anorm),y(2,M)), max(1.0_dp,y(3,M)) /)
	do
		write (*,*) 'ENTER C**2 OR 999 TO END'
		read (*,*) c2
		if (c2 == 999.0) exit
		call solvde(itmax,conv,slowc,scalv,indexv,NB,y)
		write (*,*) ' M = ',mm,'  N = ',n,&
			'  C**2 = ',c2,'  LAMBDA = ',y(3,1)+mm*(mm+1)
	end do
	END PROGRAM sfroid
