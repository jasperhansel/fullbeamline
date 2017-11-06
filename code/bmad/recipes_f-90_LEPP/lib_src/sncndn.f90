	SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: uu,emmc
	REAL(dp), INTENT(OUT) :: sn,cn,dn
	REAL(dp), PARAMETER :: CA=0.0003_dp
	INTEGER(I4B), PARAMETER :: MAXIT=13
	INTEGER(I4B) :: i,ii,l
	REAL(dp) :: a,b,c,d,emc,u
	REAL(dp), DIMENSION(MAXIT) :: em,en
	LOGICAL(LGT) :: bo
	emc=emmc
	u=uu
	if (emc /= 0.0) then
		bo=(emc < 0.0)
		if (bo) then
			d=1.0_dp-emc
			emc=-emc/d
			d=sqrt(d)
			u=d*u
		end if
		a=1.0
		dn=1.0
		do i=1,MAXIT
			l=i
			em(i)=a
			emc=sqrt(emc)
			en(i)=emc
			c=0.5_dp*(a+emc)
			if (abs(a-emc) <= CA*a) exit
			emc=a*emc
			a=c
		end do
		if (i > MAXIT) call nrerror('sncndn: convergence failed')
		u=c*u
		sn=sin(u)
		cn=cos(u)
		if (sn /= 0.0) then
			a=cn/sn
			c=a*c
			do ii=l,1,-1
				b=em(ii)
				a=c*a
				c=dn*c
				dn=(en(ii)+a)/(b+a)
				a=c/b
			end do
			a=1.0_dp/sqrt(c**2+1.0_dp)
			sn=sign(a,sn)
			cn=c*sn
		end if
		if (bo) then
			a=dn
			dn=cn
			cn=a
			sn=sn/d
		end if
	else
		cn=1.0_dp/cosh(u)
		dn=cn
		sn=tanh(u)
	end if
	END SUBROUTINE sncndn
