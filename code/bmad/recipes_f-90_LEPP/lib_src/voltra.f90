	SUBROUTINE voltra(t0,h,t,f,g,ak)
	USE nrtype; USE nrutil, ONLY : array_copy,assert_eq,unit_matrix
	USE nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: t0,h
	REAL(dp), DIMENSION(:), INTENT(OUT) :: t
	REAL(dp), DIMENSION(:,:), INTENT(OUT) :: f
	INTERFACE
		FUNCTION g(t)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: t
		REAL(dp), DIMENSION(:), POINTER :: g
		END FUNCTION g
!BL
		FUNCTION ak(t,s)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: t,s
		REAL(dp), DIMENSION(:,:), POINTER :: ak
		END FUNCTION ak
	END INTERFACE
	INTEGER(I4B) :: i,j,n,ncop,nerr,m
	INTEGER(I4B), DIMENSION(size(f,1)) :: indx
	REAL(dp) :: d
	REAL(dp), DIMENSION(size(f,1)) :: b
	REAL(dp), DIMENSION(size(f,1),size(f,1)) :: a
	n=assert_eq(size(f,2),size(t),'voltra: n')
	t(1)=t0
	call array_copy(g(t(1)),f(:,1),ncop,nerr)
	m=assert_eq(size(f,1),ncop,ncop+nerr,'voltra: m')
	do i=2,n
		t(i)=t(i-1)+h
		b=g(t(i))+0.5_dp*h*matmul(ak(t(i),t(1)),f(:,1))
		do j=2,i-1
			b=b+h*matmul(ak(t(i),t(j)),f(:,j))
		end do
		call unit_matrix(a)
		a=a-0.5_dp*h*ak(t(i),t(i))
		call ludcmp(a,indx,d)
		call lubksb(a,indx,b)
		f(:,i)=b(:)
	end do
	END SUBROUTINE voltra
