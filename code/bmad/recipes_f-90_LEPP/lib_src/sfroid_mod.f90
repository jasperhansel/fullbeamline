!+
! This module was originally in obj_src/sfroid.f90.
! This was a problem since this forced compiliation of the sfroid program and
! then the cmake system would put the program in the library.
! The solution was to separate the module from the program.
!-

MODULE sfroid_data
	USE nrtype
	INTEGER(I4B), PARAMETER :: M=41
	INTEGER(I4B) :: mm,n
	REAL(dp) :: anorm,c2,h
	REAL(dp), DIMENSION(M) :: x
END MODULE sfroid_data

