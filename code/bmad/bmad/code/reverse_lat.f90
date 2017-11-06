!+
! Subroutine reverse_lat (lat_in lat_rev, track_antiparticle)
!
! Routine to create a reversed lattice.
! That is a lattice with the tracking elements in reversed order and with all the elements having a reversed orientation.
! For lord elements the element positions are unshifted and for tracking elements:
!   lat_in%branch(ib)%ele(ie) => lat_rev%branch(ib)%ele(nt+1-ie)    ! where nt = lat_in%branch(ib)%n_ele_track
! Note: The beginning element %ele(0) is unshifted.
!
! Input:
!   lat_in             -- lat_struct: Input lattice to reverse.
!   track_antiparticle -- logical, optional: Set the particle species of the reversed lat to 
!                            the anti-particle of lat_in? Default is True.
!
! Output:
!   lat_rev            -- lat_struct: Reversed lattice. IMPORTANT! The actual lat_rev 
!                            argument must be different from the actual arg for lat_in.
!-

subroutine reverse_lat (lat_in, lat_rev, track_antiparticle)

use bmad, dummy => reverse_lat

implicit none

type (lat_struct), target :: lat_in, lat_rev
type (ele_struct), pointer :: ele_in, ele_rev, lord
type (branch_struct), pointer :: branch_in, branch_rev
type (control_struct), pointer :: con

integer i, n, ie, ib, nt, n_con, i1, i2
integer :: ix_con(size(lat_in%control))

logical, optional :: track_antiparticle
logical err_flag

!

call kill_ptc_layouts(lat_rev)    ! Cleanup lat_rev
call allocate_branch_array(lat_rev, ubound(lat_in%branch, 1))
call transfer_lat_parameters (lat_in, lat_rev)

!

do ib = 0, ubound(lat_in%branch, 1)
  branch_in => lat_in%branch(ib)
  call allocate_lat_ele_array (lat_rev, ubound(branch_in%ele, 1), ib)
  branch_rev => lat_rev%branch(ib)
  branch_rev%lat => lat_rev
  call transfer_branch_parameters(branch_in, branch_rev)
  if (logic_option(.true., track_antiparticle)) branch_rev%param%default_tracking_species = &
                                                  antiparticle(branch_in%param%default_tracking_species)
  branch_rev%param%t1_with_RF = 0  ! Init
  branch_rev%param%t1_no_RF = 0    ! Init

  nt = branch_in%n_ele_track

  ! Beginning element

  ele_rev => branch_rev%ele(0)
  ele_rev = branch_in%ele(0)
  ele_rev%orientation = -ele_rev%orientation
  ele_rev%s = 0
  call transfer_twiss (branch_in%ele(nt), ele_rev)

  ! All other elements

  do ie = 1, nt
    ele_in => branch_in%ele(nt+1-ie)
    ele_rev => branch_rev%ele(ie)
    ele_rev = ele_in
    ele_rev%orientation = -ele_rev%orientation
    ele_rev%ix_ele = ie
    if (associated(ele_rev%taylor(1)%term)) call kill_taylor(ele_rev%taylor)
    if (associated(ele_rev%ptc_genfield%field)) call kill_ptc_genfield(ele_rev%ptc_genfield%field)
  enddo
enddo

! Transfer lords

do ie = lat_in%n_ele_track+1, lat_in%n_ele_max
  ele_rev => lat_rev%ele(ie)
  ele_rev = lat_in%ele(ie)
enddo

! Correct control information

call reallocate_control(lat_rev, size(lat_in%control))
lat_rev%control = lat_in%control
lat_rev%ic = lat_in%ic

do i = 1, lat_rev%n_control_max
  con => lat_rev%control(i)
  if (con%slave%ix_ele <= lat_rev%n_ele_track .and. con%slave%ix_ele /= 0) con%slave%ix_ele = lat_rev%n_ele_track+1-con%slave%ix_ele
  if (con%lord%ix_ele <= lat_rev%n_ele_track .and. con%lord%ix_ele /= 0)  con%lord%ix_ele  = lat_rev%n_ele_track+1-con%lord%ix_ele
enddo

! Slaves of a super lord must be in assending sequence.
! ix_con keeps track of the switching.
! Also: adjust s-position of lords.

n_con = size(lat_in%control)
forall (i = 1:n_con) ix_con(i) = i 

do i = lat_rev%n_ele_track+1, lat_rev%n_ele_max
  lord => lat_rev%ele(i)
  if (lord%lord_status /= super_lord$) cycle
  i1 = lord%ix1_slave
  i2 = lord%ix1_slave+lord%n_slave-1
  lat_rev%control(i1:i2) = lat_rev%control(i2:i1:-1)
  ix_con(i1:i2) = ix_con(i2:i1:-1)
enddo

n = lat_in%n_ic_max
lat_rev%ic(1:n) = ix_con(lat_rev%ic(1:n))

! Finish

call set_flags_for_changed_attribute(lat_rev)
call s_calc(lat_rev)
call lat_sanity_check (lat_rev, err_flag)
call lattice_bookkeeper (lat_rev)

end subroutine 
