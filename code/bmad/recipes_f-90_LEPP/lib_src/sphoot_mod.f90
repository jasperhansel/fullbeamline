!+
! This module was originally in obj_src/sphoot.f90.
! This was a problem since this forced compiliation of the sphoot program and
! then the cmake system would put the program in the library.
! The solution was to separate the module from the program.
!-

MODULE sphoot_data
	USE nrtype
	INTEGER(I4B) :: m,n
	REAL(dp) :: c2,dx,gamma
END MODULE sphoot_data

MODULE sphoot_caller
	USE nrtype
	INTEGER(I4B) :: nvar
	REAL(dp) :: x1,x2
END MODULE sphoot_caller

