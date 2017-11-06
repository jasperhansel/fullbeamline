module modulo2_mod

use precision_def

!+
! Function modulo2 (x, amp) result (mod2)
!
! Function to return
!     mod2 = x + 2 * n * amp
! where n is an integer chosen such that
!    -amp <= mod2 < amp
!
! Modules needed:
!   sim_utils_interface
! 
! Input:
!   x    -- Real(sp), Real(dp), or Integer
!   amp  -- Real(sp), Real(dp), or Integer: Must be positive.
!
! Output:
!   mod2 -- Real(sp), Real(dp), or Integer: Result
!-

interface modulo2
  module procedure modulo2_sp 
  module procedure modulo2_dp
  module procedure modulo2_int
end interface

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function modulo2_sp (x, amp) result (mod2)
!
! Function to return
!     mod2 = x + 2 * n * amp
! where n is an integer chosen such that
!    -amp <= mod2 < amp
!
! Modules needed:
!   sim_utils_interface
! 
! Input:
!   x    -- Real(sp): 
!   amp  -- Real(sp): Must be positive.
!
! Output:
!   mod2 -- Real(sp): Result
!-

elemental function modulo2_sp (x, amp) result (mod2)


  use precision_def

  implicit none

  real(sp), intent(in) :: x, amp
  real(sp) mod2

!

  mod2 = modulo (x, 2*amp)
  if (mod2 >= amp) mod2 = mod2 - 2*amp

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function modulo2_dp (x, amp) result (mod2)
!
! Function to return
!     mod2 = x + 2 * n * amp
! where n is an integer chosen such that
!    -amp <= mod2 < amp
!
! Modules needed:
!   sim_utils_interface
! 
! Input:
!   x    -- Real(dp): 
!   amp  -- Real(dp): Must be positive.
!
! Output:
!   mod2 -- Real(dp): Result
!-

elemental function modulo2_dp (x, amp) result (mod2)


  use precision_def

  implicit none

  real(dp), intent(in) :: x, amp
  real(dp) mod2

!

  mod2 = modulo (x, 2*amp)
  if (mod2 >= amp) mod2 = mod2 - 2*amp

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function modulo2_int (x, amp) result (mod2)
!
! Function to return
!     mod2 = x + 2 * n * amp
! where n is an integer chosen such that
!    -amp <= mod2 < amp
!
! Modules needed:
!   sim_utils_interface
! 
! Input:
!   x    -- Integer: 
!   amp  -- Integer: Must be positive.
!
! Output:
!   mod2 -- Integer: Result
!-

elemental function modulo2_int (x, amp) result (mod2)


  use precision_def

  implicit none

  integer, intent(in) :: x, amp
  integer mod2

!

  mod2 = modulo (x, 2*amp)
  if (mod2 >= amp) mod2 = mod2 - 2*amp

end function


end module
