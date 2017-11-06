MODULE pwtcom
	USE nrtype
	INTEGER(I4B), SAVE :: ncof=0,ioff,joff
	REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE :: cc,cr
END MODULE pwtcom

	SUBROUTINE pwtset(n)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE pwtcom
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(dp) :: sig
	REAL(dp), PARAMETER :: &
		c4(4)=(/&
		0.4829629131445341_dp, 0.8365163037378079_dp, &
		0.2241438680420134_dp,-0.1294095225512604_dp /), &
		c12(12)=(/&
		0.111540743350_dp, 0.494623890398_dp, 0.751133908021_dp, &
		0.315250351709_dp,-0.226264693965_dp,-0.129766867567_dp, &
		0.097501605587_dp, 0.027522865530_dp,-0.031582039318_dp, &
		0.000553842201_dp, 0.004777257511_dp,-0.001077301085_dp /), &
		c20(20)=(/&
		0.026670057901_dp, 0.188176800078_dp, 0.527201188932_dp, &
		0.688459039454_dp, 0.281172343661_dp,-0.249846424327_dp, &
		-0.195946274377_dp, 0.127369340336_dp, 0.093057364604_dp, &
		-0.071394147166_dp,-0.029457536822_dp, 0.033212674059_dp, &
		0.003606553567_dp,-0.010733175483_dp, 0.001395351747_dp, &
		0.001992405295_dp,-0.000685856695_dp,-0.000116466855_dp, &
		0.000093588670_dp,-0.000013264203_dp /)
	if (allocated(cc)) deallocate(cc)
	if (allocated(cr)) deallocate(cr)
	allocate(cc(n),cr(n))
	ncof=n
	ioff=-n/2
	joff=-n/2
	sig=-1.0
	select case(n)
		case(4)
			cc=c4
		case(12)
			cc=c12
		case(20)
			cc=c20
		case default
			call nrerror('unimplemented value n in pwtset')
	end select
	cr(n:1:-1) = cc
	cr(n:1:-2) = -cr(n:1:-2)
	END SUBROUTINE pwtset
