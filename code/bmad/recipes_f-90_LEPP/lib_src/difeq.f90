	SUBROUTINE difeq(k,k1,k2,jsf,is1,isf,indexv,s,y)
	USE nrtype
	USE sfroid_data
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: is1,isf,jsf,k,k1,k2
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
	REAL(dp), DIMENSION(:,:), INTENT(OUT) :: s
	REAL(dp), DIMENSION(:,:), INTENT(IN) :: y
	REAL(dp) :: temp,temp2
	INTEGER(I4B), DIMENSION(3) :: indexv3
	indexv3(1:3)=3+indexv(1:3)
	if (k == k1) then
		if (mod(n+mm,2) == 1) then
			s(3,indexv3(1:3))= (/ 1.0_dp, 0.0_dp, 0.0_dp /)
			s(3,jsf)=y(1,1)
		else
			s(3,indexv3(1:3))= (/ 0.0_dp, 1.0_dp, 0.0_dp /)
			s(3,jsf)=y(2,1)
		end if
	else if (k > k2) then
		s(1,indexv3(1:3))= (/ -(y(3,M)-c2)/(2.0_dp*(mm+1.0_dp)),&
			1.0_dp, -y(1,M)/(2.0_dp*(mm+1.0_dp)) /)
		s(1,jsf)=y(2,M)-(y(3,M)-c2)*y(1,M)/(2.0_dp*(mm+1.0_dp))
		s(2,indexv3(1:3))=(/ 1.0_dp, 0.0_dp, 0.0_dp /)
		s(2,jsf)=y(1,M)-anorm
	else
		s(1,indexv(1:3))=(/ -1.0_dp, -0.5_dp*h, 0.0_dp /)
		s(1,indexv3(1:3))=(/ 1.0_dp, -0.5_dp*h, 0.0_dp /)
		temp=h/(1.0_dp-(x(k)+x(k-1))**2*0.25_dp)
		temp2=0.5_dp*(y(3,k)+y(3,k-1))-c2*0.25_dp*(x(k)+x(k-1))**2
		s(2,indexv(1:3))=(/ temp*temp2*0.5_dp,&
			-1.0_dp-0.5_dp*temp*(mm+1.0_dp)*(x(k)+x(k-1)),&
			0.25_dp*temp*(y(1,k)+y(1,k-1)) /)
		s(2,indexv3(1:3))=s(2,indexv(1:3))
		s(2,indexv3(2))=s(2,indexv3(2))+2.0_dp
		s(3,indexv(1:3))=(/ 0.0_dp, 0.0_dp, -1.0_dp /)
		s(3,indexv3(1:3))=(/ 0.0_dp, 0.0_dp, 1.0_dp /)
		s(1,jsf)=y(1,k)-y(1,k-1)-0.5_dp*h*(y(2,k)+y(2,k-1))
		s(2,jsf)=y(2,k)-y(2,k-1)-temp*((x(k)+x(k-1))*&
			0.5_dp*(mm+1.0_dp)*(y(2,k)+y(2,k-1))-temp2*&
			0.5_dp*(y(1,k)+y(1,k-1)))
		s(3,jsf)=y(3,k)-y(3,k-1)
	end if
	END SUBROUTINE difeq
