!+
! Subroutine transfer_ac_kick (ac_in, ac_out)
!
! Routine to transfer the ac_kick infor from ac_in to ac_out.
!
! Input:
!   ac_in    -- ac_kicker_struct, pointer: Input 
!
! Output:
!   ac_out   -- ac_kicker_struct, pointer: Gets set equal to ac_in
!-

subroutine transfer_ac_kick (ac_in, ac_out)

use bmad_struct

implicit none

type (ac_kicker_struct), pointer :: ac_in, ac_out
integer n

!

if (associated(ac_in)) then
  if (.not. associated(ac_out)) allocate(ac_out)

  if (allocated(ac_in%amp_vs_time)) then
    n = size(ac_in%amp_vs_time)
    if (allocated(ac_out%amp_vs_time)) then
      if (size(ac_out%amp_vs_time) /= n) deallocate(ac_out%amp_vs_time)
    endif
    if (.not. allocated(ac_out%amp_vs_time)) allocate(ac_out%amp_vs_time(n))
  else
    if (allocated(ac_out%amp_vs_time)) deallocate(ac_out%amp_vs_time)
  endif
    
  if (allocated(ac_in%frequencies)) then
    n = size(ac_in%frequencies)
    if (allocated(ac_out%frequencies)) then
      if (size(ac_out%frequencies) /= n) deallocate(ac_out%frequencies)
    endif
    if (.not. allocated(ac_out%frequencies)) allocate(ac_out%frequencies(n))
  else
    if (allocated(ac_out%frequencies)) deallocate(ac_out%frequencies)
  endif

  ac_out = ac_in

else
  if (associated(ac_out)) deallocate(ac_out)
endif

end subroutine transfer_ac_kick
