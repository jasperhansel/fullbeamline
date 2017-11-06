!+
! This module was originally in obj_src/sphfpt.f90.
! This was a problem since this forced compiliation of the sphfpt program and
! then the cmake system would put the program in the library.
! The solution was to separate the module from the program.
!-

MODULE sphfpt_data
	USE nrtype
	INTEGER(I4B) :: m,n
	REAL(dp) :: c2,dx,gamma
END MODULE sphfpt_data

MODULE sphfpt_caller
	USE nrtype
	INTEGER(I4B) :: nn2
	REAL(dp) :: x1,x2,xf
END MODULE sphfpt_caller

