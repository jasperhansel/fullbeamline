	SUBROUTINE dftcor(w,delta,a,b,endpts,corre,corim,corfac)
	USE nrtype; USE nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: w,delta,a,b
	REAL(dp), INTENT(OUT) :: corre,corim,corfac
	REAL(dp), DIMENSION(:), INTENT(IN) :: endpts
	REAL(dp) :: a0i,a0r,a1i,a1r,a2i,a2r,a3i,a3r,arg,c,cl,cr,s,sl,sr,t,&
		t2,t4,t6
	REAL(DP) :: cth,ctth,spth2,sth,sth4i,stth,th,th2,th4,&
		tmth2,tth4i
	th=w*delta
	call assert(a < b, th >= 0.0, th <= PI_D, 'dftcor args')
	if (abs(th) < 5.0e-2_dp) then
		t=th
		t2=t*t
		t4=t2*t2
		t6=t4*t2
		corfac=1.0_dp-(11.0_dp/720.0_dp)*t4+(23.0_dp/15120.0_dp)*t6
		a0r=(-2.0_dp/3.0_dp)+t2/45.0_dp+(103.0_dp/15120.0_dp)*t4-&
			(169.0_dp/226800.0_dp)*t6
		a1r=(7.0_dp/24.0_dp)-(7.0_dp/180.0_dp)*t2+(5.0_dp/3456.0_dp)*t4&
			-(7.0_dp/259200.0_dp)*t6
		a2r=(-1.0_dp/6.0_dp)+t2/45.0_dp-(5.0_dp/6048.0_dp)*t4+t6/64800.0_dp
		a3r=(1.0_dp/24.0_dp)-t2/180.0_dp+(5.0_dp/24192.0_dp)*t4-t6/259200.0_dp
		a0i=t*(2.0_dp/45.0_dp+(2.0_dp/105.0_dp)*t2-&
			(8.0_dp/2835.0_dp)*t4+(86.0_dp/467775.0_dp)*t6)
		a1i=t*(7.0_dp/72.0_dp-t2/168.0_dp+(11.0_dp/72576.0_dp)*t4-&
			(13.0_dp/5987520.0_dp)*t6)
		a2i=t*(-7.0_dp/90.0_dp+t2/210.0_dp-(11.0_dp/90720.0_dp)*t4+&
			(13.0_dp/7484400.0_dp)*t6)
		a3i=t*(7.0_dp/360.0_dp-t2/840.0_dp+(11.0_dp/362880.0_dp)*t4-&
			(13.0_dp/29937600.0_dp)*t6)
	else
		cth=cos(th)
		sth=sin(th)
		ctth=cth**2-sth**2
		stth=2.0_dp*sth*cth
		th2=th*th
		th4=th2*th2
		tmth2=3.0_dp-th2
		spth2=6.0_dp+th2
		sth4i=1.0_dp/(6.0_dp*th4)
		tth4i=2.0_dp*sth4i
		corfac=tth4i*spth2*(3.0_dp-4.0_dp*cth+ctth)
		a0r=sth4i*(-42.0_dp+5.0_dp*th2+spth2*(8.0_dp*cth-ctth))
		a0i=sth4i*(th*(-12.0_dp+6.0_dp*th2)+spth2*stth)
		a1r=sth4i*(14.0_dp*tmth2-7.0_dp*spth2*cth)
		a1i=sth4i*(30.0_dp*th-5.0_dp*spth2*sth)
		a2r=tth4i*(-4.0_dp*tmth2+2.0_dp*spth2*cth)
		a2i=tth4i*(-12.0_dp*th+2.0_dp*spth2*sth)
		a3r=sth4i*(2.0_dp*tmth2-spth2*cth)
		a3i=sth4i*(6.0_dp*th-spth2*sth)
	end if
	cl=a0r*endpts(1)+a1r*endpts(2)+a2r*endpts(3)+a3r*endpts(4)
	sl=a0i*endpts(1)+a1i*endpts(2)+a2i*endpts(3)+a3i*endpts(4)
	cr=a0r*endpts(8)+a1r*endpts(7)+a2r*endpts(6)+a3r*endpts(5)
	sr=-a0i*endpts(8)-a1i*endpts(7)-a2i*endpts(6)-a3i*endpts(5)
	arg=w*(b-a)
	c=cos(arg)
	s=sin(arg)
	corre=cl+c*cr-s*sr
	corim=sl+s*cr+c*sr
	END SUBROUTINE dftcor
