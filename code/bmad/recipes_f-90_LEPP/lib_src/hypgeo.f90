MODULE hypgeo_info
	USE nrtype
	COMPLEX(dpc) :: hypgeo_aa,hypgeo_bb,hypgeo_cc,hypgeo_dz,hypgeo_z0
END MODULE hypgeo_info


	FUNCTION hypgeo(a,b,c,z)
	USE nrtype
	USE hypgeo_info
	USE nr, ONLY : bsstep,hypdrv,hypser,odeint
	IMPLICIT NONE
	COMPLEX(dpc), INTENT(IN) :: a,b,c,z
	COMPLEX(dpc) :: hypgeo
	REAL(dp), PARAMETER :: EPS=1.0e-6_dp
	COMPLEX(dpc), DIMENSION(2) :: y
	REAL(dp), DIMENSION(4) :: ry
	if (real(z)**2+aimag(z)**2 <= 0.25) then
		call hypser(a,b,c,z,hypgeo,y(2))
		RETURN
	else if (real(z) < 0.0) then
		hypgeo_z0=cmplx(-0.5_dp,0.0_dp,kind=dpc)
	else if (real(z) <= 1.0) then
		hypgeo_z0=cmplx(0.5_dp,0.0_dp,kind=dpc)
	else
		hypgeo_z0=cmplx(0.0_dp,sign(0.5_dp,aimag(z)),kind=dpc)
	end if
	hypgeo_aa=a
	hypgeo_bb=b
	hypgeo_cc=c
	hypgeo_dz=z-hypgeo_z0
	call hypser(hypgeo_aa,hypgeo_bb,hypgeo_cc,hypgeo_z0,y(1),y(2))
	ry(1:4:2)=real(y)
	ry(2:4:2)=aimag(y)
	call odeint(ry,0.0_dp,1.0_dp,EPS,0.1_dp,0.0001_dp,hypdrv,bsstep)
	y=cmplx(ry(1:4:2),ry(2:4:2),kind=dpc)
	hypgeo=y(1)
	END FUNCTION hypgeo
