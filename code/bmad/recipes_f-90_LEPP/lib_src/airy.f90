	SUBROUTINE airy(x,ai,bi,aip,bip)
	USE nrtype
	USE nr, ONLY : bessik,bessjy
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp), INTENT(OUT) :: ai,bi,aip,bip
	REAL(dp) :: absx,ri,rip,rj,rjp,rk,rkp,rootx,ry,ryp,z
	REAL(dp), PARAMETER :: THIRD=1.0_dp/3.0_dp,TWOTHR=2.0_dp/3.0_dp, &
		ONOVRT=0.5773502691896258_dp
	absx=abs(x)
	rootx=sqrt(absx)
	z=TWOTHR*absx*rootx
	if (x > 0.0) then
		call bessik(z,THIRD,ri,rk,rip,rkp)
		ai=rootx*ONOVRT*rk/PI
		bi=rootx*(rk/PI+2.0_dp*ONOVRT*ri)
		call bessik(z,TWOTHR,ri,rk,rip,rkp)
		aip=-x*ONOVRT*rk/PI
		bip=x*(rk/PI+2.0_dp*ONOVRT*ri)
	else if (x < 0.0) then
		call bessjy(z,THIRD,rj,ry,rjp,ryp)
		ai=0.5_dp*rootx*(rj-ONOVRT*ry)
		bi=-0.5_dp*rootx*(ry+ONOVRT*rj)
		call bessjy(z,TWOTHR,rj,ry,rjp,ryp)
		aip=0.5_dp*absx*(ONOVRT*ry+rj)
		bip=0.5_dp*absx*(ONOVRT*rj-ry)
	else
		ai=0.3550280538878172_dp
		bi=ai/ONOVRT
		aip=-0.2588194037928068_dp
		bip=-aip/ONOVRT
	end if
	END SUBROUTINE airy
