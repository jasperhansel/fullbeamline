!+
! Below are dummy routines for Numerical Recipes. 
! This is only needed when using dynamic linking since, in this case, 
! all external references must be satisfied.
!-

! See recipes/lib_src/quad3d.f90

FUNCTION y1(x)
use precision_def
REAL(dp), INTENT(IN) :: x
REAL(dp) :: y1
y1 = 0
END FUNCTION y1

FUNCTION y2(x)
use precision_def
REAL(dp), INTENT(IN) :: x
REAL(dp) :: y2
y2 = 0
END FUNCTION y2

FUNCTION z1(x,y)
use precision_def
REAL(dp), INTENT(IN) :: x,y
REAL(dp) :: z1
z1 = 0
END FUNCTION z1

FUNCTION z2(x,y)
use precision_def
REAL(dp), INTENT(IN) :: x,y
REAL(dp) :: z2
z2 = 0
END FUNCTION z2
