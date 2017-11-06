	FUNCTION dbrent(ax,bx,cx,func,dfunc,tol,xmin)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: ax,bx,cx,tol
	REAL(dp), INTENT(OUT) :: xmin
	REAL(dp) :: dbrent
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: func
		END FUNCTION func
!BL
		FUNCTION dfunc(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(dp), PARAMETER :: ZEPS=1.0e-3_dp*epsilon(ax)
	INTEGER(I4B) :: iter
	REAL(dp) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,&
		u,u1,u2,v,w,x,xm
	LOGICAL :: ok1,ok2
	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	fx=func(x)
	fv=fx
	fw=fx
	dx=dfunc(x)
	dv=dx
	dw=dx
	do iter=1,ITMAX
		xm=0.5_dp*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0_dp*tol1
		if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) exit
		if (abs(e) > tol1) then
			d1=2.0_dp*(b-a)
			d2=d1
			if (dw /= dx) d1=(w-x)*dx/(dx-dw)
			if (dv /= dx) d2=(v-x)*dx/(dx-dv)
			u1=x+d1
			u2=x+d2
			ok1=((a-u1)*(u1-b) > 0.0) .and. (dx*d1 <= 0.0)
			ok2=((a-u2)*(u2-b) > 0.0) .and. (dx*d2 <= 0.0)
			olde=e
			e=d
			if (ok1 .or. ok2) then
				if (ok1 .and. ok2) then
					d=merge(d1,d2, abs(d1) < abs(d2))
				else
					d=merge(d1,d2,ok1)
				end if
				if (abs(d) <= abs(0.5_dp*olde)) then
					u=x+d
					if (u-a < tol2 .or. b-u < tol2) &
						d=sign(tol1,xm-x)
				else
					e=merge(a,b, dx >= 0.0)-x
					d=0.5_dp*e
				end if
			else
				e=merge(a,b, dx >= 0.0)-x
				d=0.5_dp*e
			end if
		else
			e=merge(a,b, dx >= 0.0)-x
			d=0.5_dp*e
		end if
		if (abs(d) >= tol1) then
			u=x+d
			fu=func(u)
		else
			u=x+sign(tol1,d)
			fu=func(u)
			if (fu > fx) exit
		end if
		du=dfunc(u)
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			call mov3(v,fv,dv,w,fw,dw)
			call mov3(w,fw,dw,x,fx,dx)
			call mov3(x,fx,dx,u,fu,du)
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu <= fw .or. w == x) then
				call mov3(v,fv,dv,w,fw,dw)
				call mov3(w,fw,dw,u,fu,du)
			else if (fu <= fv .or. v == x .or. v == w) then
				call mov3(v,fv,dv,u,fu,du)
			end if
		end if
	end do
	if (iter > ITMAX) call nrerror('dbrent: exceeded maximum iterations')
	xmin=x
	dbrent=fx
	CONTAINS
!BL
	SUBROUTINE mov3(a,b,c,d,e,f)
	REAL(dp), INTENT(IN) :: d,e,f
	REAL(dp), INTENT(OUT) :: a,b,c
	a=d
	b=e
	c=f
	END SUBROUTINE mov3
	END FUNCTION dbrent
