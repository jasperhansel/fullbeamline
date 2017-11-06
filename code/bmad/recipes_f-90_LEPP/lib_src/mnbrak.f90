	SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
	USE nrtype; USE nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(dp), INTENT(INOUT) :: ax,bx
	REAL(dp), INTENT(OUT) :: cx,fa,fb,fc
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: func
		END FUNCTION func
	END INTERFACE
	REAL(dp), PARAMETER :: GOLD=1.618034_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp
	REAL(dp) :: fu,q,r,u,ulim
	fa=func(ax)
	fb=func(bx)
	if (fb > fa) then
		call swap(ax,bx)
		call swap(fa,fb)
	end if
	cx=bx+GOLD*(bx-ax)
	fc=func(cx)
	do
		if (fb < fc) RETURN
		r=(bx-ax)*(fb-fc)
		q=(bx-cx)*(fb-fa)
		u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*sign(max(abs(q-r),TINY),q-r))
		ulim=bx+GLIMIT*(cx-bx)
		if ((bx-u)*(u-cx) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				ax=bx
				fa=fb
				bx=u
				fb=fu
				RETURN
			else if (fu > fb) then
				cx=u
				fc=fu
				RETURN
			end if
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		else if ((cx-u)*(u-ulim) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				bx=cx
				cx=u
				u=cx+GOLD*(cx-bx)
				call shft(fb,fc,fu,func(u))
			end if
		else if ((u-ulim)*(ulim-cx) >= 0.0) then
			u=ulim
			fu=func(u)
		else
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		end if
		call shft(ax,bx,cx,u)
		call shft(fa,fb,fc,fu)
	end do
	CONTAINS
!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(dp), INTENT(OUT) :: a
	REAL(dp), INTENT(INOUT) :: b,c
	REAL(dp), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END SUBROUTINE mnbrak
