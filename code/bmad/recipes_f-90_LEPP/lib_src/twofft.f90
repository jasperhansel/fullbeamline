	SUBROUTINE twofft(data1,data2,fft1,fft2)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : four1
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: data1,data2
	COMPLEX(dpc), DIMENSION(:), INTENT(OUT) :: fft1,fft2
	INTEGER(I4B) :: n,n2
	COMPLEX(dpc), PARAMETER :: C1=(0.5_dp,0.0_dp), C2=(0.0_dp,-0.5_dp)
	COMPLEX, DIMENSION(size(data1)/2+1) :: h1,h2
	n=assert_eq(size(data1),size(data2),size(fft1),size(fft2),'twofft')
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in twofft')
	fft1=cmplx(data1,data2,kind=dpc)
	call four1(fft1,1)
	fft2(1)=cmplx(aimag(fft1(1)),0.0_dp,kind=dpc)
	fft1(1)=cmplx(real(fft1(1)),0.0_dp,kind=dpc)
	n2=n/2+1
	h1(2:n2)=C1*(fft1(2:n2)+conjg(fft1(n:n2:-1)))
	h2(2:n2)=C2*(fft1(2:n2)-conjg(fft1(n:n2:-1)))
	fft1(2:n2)=h1(2:n2)
	fft1(n:n2:-1)=conjg(h1(2:n2))
	fft2(2:n2)=h2(2:n2)
	fft2(n:n2:-1)=conjg(h2(2:n2))
	END SUBROUTINE twofft
