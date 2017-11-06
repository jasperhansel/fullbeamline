module superimpose_mod

use bmad_struct
use bmad_interface
use lat_ele_loc_mod
use bookkeeper_mod

private delete_underscore, adjust_super_slave_names, adjust_drift_names, split_this_lat

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!+
! Subroutine add_superimpose (lat, super_ele_in, ix_branch, err_flag, super_ele_out,
!                                         save_null_drift, create_jumbo_slave, ix_insert)
!
! Routine to superimpose an element. If the element can be inserted
! into the lat without making a super_lord element then this will be done.
!
! Note: This routine, since it handles only one superposition, is not sufficient for
!   superposition in a multipass region. For historical reasons, the extra code needed 
!   is buried in the parser_add_superimpose code. If you need to do multipass superpositions 
!   please contact David Sagan and this situation will be rectified.
!
! Note: Bookkeeping like recalculating reference energies and recalculating transfer matrices 
!   is *not* done by this routine.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat              -- lat_struct: Lat to modify.
!   super_ele_in     -- Ele_struct: Element to superimpose.
!         %s            -- Position of end of element.
!                          Negative distances mean distance from the end.
!   ix_branch        -- Integer: Branch index to put element.
!   save_null_drift  -- Logical, optional: Save a copy of a drift to be split as a null_ele?
!                         This is useful if further superpositions might use this drift as a 
!                         reference element. After all superpositions are done, 
!                         remove_eles_from_lat can be called to remove all null_eles.
!                         Default is False.
!   create_jumbo_slave
!                   -- Logical, optional: Default is False. If True then super_slaves
!                       that are created that have super_ele_in as their super_lord are
!                       em_field elements.
!   ix_insert       -- integer, optional: If present and positive, and super_ele_in has zero length,
!                       use ix_insert as the index to insert super_ele_in at. ix_insert is useful 
!                       when superposing next to another zero length element and you want to make sure 
!                       that the superimposed element is on the correct side of the existing element.
!                       Note: As a sanity check, it is an error if super_ele_in%s does not match the 
!                       s-position at ix_insert.
!
! Output:
!   lat           -- lat_struct: Modified lat.
!   err_flag      -- Logical :: Set True if there is an error. False otherwise
!   super_ele_out -- Ele_struct, pointer, optional :: Pointer to the super element in the lattice.
!-

subroutine add_superimpose (lat, super_ele_in, ix_branch, err_flag, super_ele_out, &
                                              save_null_drift, create_jumbo_slave, ix_insert)

implicit none

type (lat_struct), target :: lat
type (ele_struct)  super_ele_in
type (ele_struct), pointer, optional ::  super_ele_out
type (ele_struct) super_saved, slave_saved, drift, null_ele
type (ele_struct), pointer :: slave, lord, slave2, ele0, ele
type (control_struct)  sup_con(100)
type (control_struct), pointer ::  cntl
type (branch_struct), pointer :: branch
type (lat_ele_loc_struct), pointer :: loc

real(rp) s1, s2, s1_lat, s2_lat, s1_lat_fudge, s2_lat_fudge, s1_in, s2_in
real(rp) ds_small, l_super

integer, optional :: ix_insert
integer i, j, jj, k, ix, n, i2, ic, n_con, ixs, ix_branch, ii
integer ix1_split, ix2_split, ix_super, ix_super_con, ix_ic
integer ix_slave, ixn, ixc, ix_1lord, ix_lord_max_old

logical, optional :: save_null_drift, create_jumbo_slave
logical err_flag, setup_lord, split1_done, split2_done, all_drift, err, zero_length_lord

character(100) name
character(20) fmt
character(20) :: r_name = "add_superimpose"

!-------------------------------------------------------------------------
! Check for negative length

err_flag = .true.
l_super = super_ele_in%value(l$)

if (l_super < 0) then
  call out_io (s_abort$, r_name, &
                  'Superposition of element with negative length not allowed!', &
                  'Element: ' // super_ele_in%name, &
                  'Length: \es10.2\ ', r_array = [l_super] )
  if (global_com%exit_on_error) call err_exit
  return
endif

! We need a copy of super_ele_in since the actual argument may be in the lat
! and split_lat can then overwrite it.

drift%key = drift$

call transfer_ele (super_ele_in, super_saved)
super_saved%slave_status = free$
super_saved%n_lord = 0
super_saved%ic1_lord = 0

branch => lat%branch(ix_branch)

! s1 is the entrance edge of the superimpose.
! s2 is the exit edge of the superimpose.
! For a lat a superimpose can wrap around the ends of the lattice.

ix_lord_max_old = lat%n_ele_max

s1_lat = branch%ele(0)%s                 ! normally this is zero.
s1_lat_fudge = s1_lat - bmad_com%significant_length

s2_lat = branch%ele(branch%n_ele_track)%s
s2_lat_fudge = s2_lat + bmad_com%significant_length

s1_in = super_saved%s - l_super
s2_in = super_saved%s                 

s1 = s1_in
s2 = s2_in

if (s1 < s1_lat_fudge) then
  if (branch%param%geometry == open$) call out_io (s_warn$, &
         r_name, 'Superimpose is being wrapped around linear lattice for: ' // super_saved%name)
  s1 = s1 + branch%param%total_length
endif

if (s2 > s2_lat_fudge) then
  if (branch%param%geometry == open$) call out_io (s_warn$, &
         r_name, 'Superimpose is being wrapped around linear lattice for: ' // super_saved%name)
  s2 = s2 - branch%param%total_length
endif

if (s1 < s1_lat_fudge .or. s2 < s1_lat_fudge .or. s1 > s2_lat_fudge .or. s2 > s2_lat_fudge) then
  call out_io (s_abort$, r_name, &
    'SUPERIMPOSE POSITION BEYOUND END OF LATTICE FOR ELEMENT: ' // super_saved%name, &
    'ELEMENT WANTS TO BE PLACED AT: [\F10.1\, \F10.1\] ', &
    'lATTICE EXTENT:                [\f10.1\, \F10.1\]', r_array = [s1_in, s2_in, s1_lat, s2_lat])
  if (global_com%exit_on_error) call err_exit
  return
endif

!-------------------------------------------------------------------------
! If the element has zero length, just insert it in the tracking part of the lattice list.

! Note: Important to set super_saved%lord_status before calling insert_element since
! this affects the status flag setting.

if (super_saved%value(l$) == 0) then
  super_saved%lord_status  = not_a_lord$ 
  ix = integer_option(-1, ix_insert)
  if (ix > 0) then
    if (abs(branch%ele(ix-1)%s - s1) > bmad_com%significant_length / 10) then
      call out_io (s_abort$, r_name, 'SUPERIMPOSE CONFUSION WITH ZERO LENGTH SUPERPOSITION.')
      if (global_com%exit_on_error) call err_exit
      return
    endif
    ix1_split = ix - 1
  else
    call split_lat (lat, s1, ix_branch, ix1_split, split1_done, &
                                check_sanity = .false., save_null_drift = save_null_drift, err_flag = err)
    if (err) return
  endif
  call insert_element (lat, super_saved, ix1_split+1, ix_branch)
  ix_super = ix1_split + 1
  if (present(super_ele_out)) super_ele_out => branch%ele(ix_super)
  call adjust_super_slave_names (lat, lat%n_ele_track+1, lat%n_ele_max)
  call adjust_drift_names (lat, branch%ele(ix1_split))
  err_flag = .false.
  return
endif

!-------------------------------------------------------------------------
! Split lat at begining and end of the superimpose.
! The correct order of splitting is important since we are creating elements
! so that the numbering of the elments after the split changes.

! if superimpose wraps around 0 ...
if (s2 < s1) then
  if (split_this_lat(2, branch, s2, ix2_split, split2_done, save_null_drift, create_jumbo_slave)) return
  if (split_this_lat(1, branch, s1, ix1_split, split1_done, save_null_drift, create_jumbo_slave)) return

! no wrap case...
else                  
  if (split_this_lat(1, branch, s1, ix1_split, split1_done, save_null_drift, create_jumbo_slave)) return
  if (split_this_lat(2, branch, s2, ix2_split, split2_done, save_null_drift, create_jumbo_slave)) return
endif

! If the element has zero length then need to insert a zero length element

zero_length_lord = .false.
if (abs(l_super) < 10*bmad_com%significant_length .and. .not. logic_option(.false., create_jumbo_slave)) then
  zero_length_lord = .true.
  ds_small = 10 * bmad_com%significant_length

  if (branch%ele(ix1_split)%value(l$) > ds_small) then
    call split_lat (branch%lat, s1-ds_small, branch%ix_branch, ix1_split, split2_done, &
                                              .false., .false., save_null_drift, err, choose_max = .true.)
    ix2_split = ix2_split + 1

  elseif (branch%ele(ix1_split+1)%value(l$) > ds_small) then
    call split_lat (branch%lat, s1+ds_small, branch%ix_branch, ix2_split, split2_done, &
                                              .false., .false., save_null_drift, err, choose_max = .false.)

  else
    call out_io (s_fatal$, r_name, 'CONFUSED SUPERPOSITION WITH ELEMENT OF SMALL LENGTH!')
    if (global_com%exit_on_error) call err_exit
  endif

  if (err) return
  branch%ele(ix1_split)%value(l$) = branch%ele(ix1_split)%value(l$) + ds_small - l_super ! Reset to original size - l_super
  branch%ele(ix1_split)%s = branch%ele(ix1_split)%s + ds_small - l_super
  branch%ele(ix2_split)%value(l$) = l_super                                                
  branch%ele(ix2_split)%s_start = branch%ele(ix2_split)%s - branch%ele(ix2_split)%value(l$)
endif

! The splits may not be done exactly at s1 and s2 since split_lat avoids
! making extremely small "runt" elements. We thus adjust the lord length slightly to
! keep everything consistant. 

if (split1_done) super_saved%value(l$) = super_saved%value(l$) - (branch%ele(ix1_split)%s - s1) 
if (split2_done) super_saved%value(l$) = super_saved%value(l$) + (branch%ele(ix2_split)%s - s2) 

! Get rid of duplicate saved drifts if present

n = lat%n_ele_max
if (lat%ele(n)%key == null_ele$ .and. lat%ele(n-1)%key == null_ele$ .and. &
        lat%ele(n)%s == lat%ele(n-1)%s .and. lat%ele(n)%name == lat%ele(n-1)%name) then
  lat%ele(n)%key = -1
  call remove_eles_from_lat (lat, .false.)
endif

! zero length elements at the edges of the superimpose region can be excluded
! from the region

if (.not. zero_length_lord) then
  do 
    if (branch%ele(ix1_split+1)%value(l$) /= 0) exit
    ix1_split = ix1_split + 1
    if (ix1_split > branch%n_ele_track) ix1_split = 0
  enddo

  do
    if (branch%ele(ix2_split)%value(l$) /= 0) exit
    ix2_split = ix2_split - 1
    if (ix2_split == -1) ix2_split = branch%n_ele_track
  enddo
endif

! If there are null_ele elements in the superimpose region then just move them
! out of the way to the lord section of the branch. This prevents unnecessary
! splitting.

i = ix1_split
do
  i = i + 1
  if (i > branch%n_ele_track) i = 0
  if (branch%ele(i)%key == null_ele$) then

    branch%n_ele_max = branch%n_ele_max + 1
    ix = branch%n_ele_max
    if (ix > ubound(branch%ele, 1))  call allocate_lat_ele_array (lat, ix_branch = ix_branch)

    branch%ele(ix) = branch%ele(i)       ! copy null_ele
    do ic = branch%ele(i)%ic1_lord, branch%ele(i)%ic1_lord+branch%ele(i)%n_lord-1
      j = lat%ic(ic)
      lat%control(j)%slave%ix_ele = ix ! point to new null_ele.
    enddo
    branch%ele(i)%key = -1  ! Mark old null_ele for deletion
    call remove_eles_from_lat (lat, .false.)
    i = i - 1
    if (ix2_split > i) ix2_split = ix2_split - 1
    if (ix1_split > i) ix1_split = ix1_split - 1
  endif
  if (i == ix2_split) exit
enddo

! If element overlays only drifts then just 
! insert it in the tracking part of the lat list.

all_drift = (ix2_split > ix1_split)
do i = ix1_split+1, ix2_split
  if (branch%ele(i)%key /= drift$) all_drift = .false.
  if (branch%ele(i)%slave_status == super_slave$ .or. branch%ele(i)%slave_status == multipass_slave$) all_drift = .false.
  if (.not. all_drift) exit
enddo

if (all_drift) then  
  do i = ix1_split+2, ix2_split    ! remove all drifts but one
    branch%ele(i)%key = -1    ! mark for deletion
  enddo
  call remove_eles_from_lat(lat, .false.)    ! And delete
  ix_super = ix1_split + 1
  branch%ele(ix_super) = super_saved
  branch%ele(ix_super)%lord_status  = not_a_lord$
  call set_flags_for_changed_attribute (branch%ele(ix_super))
  if (present(super_ele_out)) super_ele_out => branch%ele(ix_super)

  if (split1_done .and. branch%ele(ix_super-1)%key == drift$) call adjust_drift_names (lat, branch%ele(ix_super-1))
  if (split2_done .and. branch%ele(ix_super+1)%key == drift$) call adjust_drift_names (lat, branch%ele(ix_super+1))

  err_flag = .false.
  return
endif

! Only possibility left means we have to set up a super_lord element 
! representing the superimposed element for the superposition...

! First: It is not legal for an element to be simultaneously a multipass_slave and a super_slave.
! Thus if the elements to be superimposed upon are multipass_slaves,
! we need to make them super_slaves and create a corresponding super_lord.

do i = ix1_split+1, ix2_split
  slave => branch%ele(i)
  if (slave%slave_status /= multipass_slave$) cycle
  ! Create a lord for this multipass_slave
  call new_control(lat, ixs)
  slave => branch%ele(i) ! need this if branch%ele was reallocated
  lord => lat%ele(ixs)
  lord = slave
  lord%lord_status = super_lord$
  ! Point control info to this new lord
  do j = 1, lat%n_control_max
    if (lat%control(j)%slave%ix_ele == i) lat%control(j)%slave%ix_ele = ixs
  enddo
  ! Now put in the info to make the original element a super_slave
  call add_lattice_control_structs (lord, n_add_slave = 1)
  ix = lord%ix1_slave
  lat%control(ix)%slave%ix_ele = i
  lat%control(ix)%slave%ix_branch = ix_branch
  slave%slave_status = super_slave$
  slave%name = trim(slave%name) // '#1'
  slave%ic1_lord = 0   ! So add_lattice_control_structs does the right thing
  call add_lattice_control_structs (slave, n_add_lord = 1)
  ic = slave%ic1_lord
  lat%ic(ic) = ix
  !
  if (allocated(ele_loc_com%branch)) then
    do ii = lbound(ele_loc_com%branch, 1), ubound(ele_loc_com%branch, 1)
      do jj = lbound(ele_loc_com%branch(ii)%ele, 1), ubound(ele_loc_com%branch(ii)%ele, 1)
        loc => ele_loc_com%branch(ii)%ele(jj)
        if (loc%ix_branch /= branch%ix_branch .or. loc%ix_ele /= slave%ix_ele) cycle
        loc%ix_branch = lord%ix_branch
        loc%ix_ele = lord%ix_ele
      enddo
    enddo
  endif
enddo

! Now to create the superimposed element super_lord.

ix_super = lat%n_ele_max + 1
lat%n_ele_max = ix_super
if (lat%n_ele_max > ubound(lat%ele, 1)) call allocate_lat_ele_array(lat)
lat%ele(ix_super) = super_saved
lat%ele(ix_super)%lord_status = super_lord$
call set_flags_for_changed_attribute (lat%ele(ix_super))
if (present(super_ele_out)) super_ele_out => lat%ele(ix_super)

ix_super_con = 0

!-------------------------------------------------------------------------
! If create_jumbo_slave = T:
! If any existing super_lords extend past the region to be superimposed upon then
! extend this region to include the existing super_lord. 

if (logic_option(.false., create_jumbo_slave)) then

  do
    slave => branch%ele(ix1_split+1)
    if (slave%slave_status /= super_slave$) exit
    do i = 1, slave%n_lord
      lord => pointer_to_lord(slave, i, ix_slave = ix)
      if (ix == 1) cycle
      slave2 => pointer_to_slave(lord, 1)
      ix1_split = slave2%ix_ele - 1
      exit
    enddo
    if (i == slave%n_lord + 1) exit
  enddo

  do
    slave => branch%ele(ix2_split)
    if (slave%slave_status /= super_slave$) exit
    do i = 1, slave%n_lord
      lord => pointer_to_lord(slave, i, ix_slave = ix)
      if (ix == lord%n_slave) cycle
      slave2 => pointer_to_slave(lord, lord%n_slave)
      ix2_split = slave2%ix_ele
      exit
    enddo
    if (i == slave%n_lord + 1) exit
  enddo

endif

!-------------------------------------------------------------------------
! Go through the list of elements being superimposed upon.
! Zero length elements (markers, etc.) do not get involved here.

ix_slave = ix1_split

do 

  ix_slave = ix_slave + 1
  if (ix_slave == ix2_split + 1) exit
  if (ix_slave == branch%n_ele_track + 1) ix_slave = 1

  slave => branch%ele(ix_slave)
  call transfer_ele(slave, slave_saved)
  if (slave_saved%value(l$) == 0) cycle

  select case (slave_saved%key)
  case (multipole$, ab_multipole$)
   if (slave_saved%value(l$) == 0) cycle
  end select

  ! Do we need to set up a super lord to control this slave element?
  ! if yes then create the super lord element

  setup_lord = (slave%slave_status /= super_slave$ .and. (slave%key /= drift$ .or. slave%n_lord /= 0))

  if (setup_lord) then
    call new_control (lat, ixn)
    slave => branch%ele(ix_slave)   ! need this if branch%ele was reallocated
    lat%ele(ixn) = slave_saved
    lat%ele(ixn)%lord_status = super_lord$
    ixc = lat%n_control_max + 1
    if (ixc > size(lat%control)) call reallocate_control(lat, ixc+100)
    lat%ele(ixn)%ix1_slave = ixc
    lat%ele(ixn)%n_slave = 1
    lat%control(ixc)%lord%ix_ele = ixn
    lat%control(ixc)%slave = lat_ele_loc_struct(ix_slave, ix_branch)
    lat%n_control_max = ixc

    do j = lat%ele(ixn)%ic1_lord, lat%ele(ixn)%ic1_lord+lat%ele(ixn)%n_lord-1
      jj = lat%ic(j)
      lat%control(jj)%slave%ix_ele = ixn
    enddo

    ic = lat%n_ic_max + 1
    slave%ic1_lord = ic
    slave%n_lord = 2
    lat%n_ic_max = ic + 1
    lat%ic(ic) = ixc 

  else
    call add_lattice_control_structs (slave, n_add_lord = 1)
  endif

  ! Components like %wall3d, %em_field are contained in the lord so
  ! deallocate any of these components in the future slave.
  ! Note that %a_pole and %b_pole components are the exception

  call deallocate_ele_pointers (slave, nullify_branch = .false., dealloc_poles = .false.)
  slave%slave_status = super_slave$

  ! add control info for main super lord to list

  ix_super_con = ix_super_con + 1
  sup_con(ix_super_con)%slave = lat_ele_loc_struct(ix_slave, ix_branch)
  sup_con(ix_super_con)%lord%ix_ele = ix_super
  sup_con(ix_super_con)%ix_attrib = 0

  ! change the element key

  if (logic_option(.false., create_jumbo_slave)) then
    slave%key = em_field$
  else
    call calc_super_slave_key(slave_saved, super_saved, slave)
    if (slave%key <= 0) then
      call out_io (s_abort$, r_name, &
              'I DO NOT KNOW HOW TO SUPERIMPOSE ELEMENT: "' // trim(super_saved%name) // &
                                                 '" OF TYPE: ' // key_name(super_saved%key), &
              'UPON: "' // trim(slave_saved%name) // '" OF TYPE: ' // key_name(slave_saved%key))
      if (global_com%exit_on_error) call err_exit 
      return                   
    endif
  endif

  call set_flags_for_changed_attribute (lat%ele(ix_super))

enddo

! Special case where elements on either side of the superimpose have the same
! name

if (split1_done .and. split2_done .and. &
              branch%ele(ix1_split)%name == branch%ele(ix2_split+1)%name) then
  branch%ele(ix1_split)%name = trim(branch%ele(ix1_split)%name) // '#1'
  branch%ele(ix2_split+1)%name = trim(branch%ele(ix2_split+1)%name) // '#2'
  call delete_underscore (branch%ele(ix1_split))
  call delete_underscore (branch%ele(ix2_split+1))
endif

! transfer control info from sup_con array

ixc = lat%n_control_max + 1
n_con = ixc + ix_super_con - 1
if (n_con > size(lat%control)) call reallocate_control(lat, n_con+500) 
lat%ele(ix_super)%ix1_slave = ixc
lat%ele(ix_super)%n_slave = ix_super_con
if (present(super_ele_out)) super_ele_out => lat%ele(ix_super)

do k = 1, ix_super_con
  lat%control(k+ixc-1) = sup_con(k)
  ix_slave = lat%control(k+ixc-1)%slave%ix_ele
  i2 = branch%ele(ix_slave)%ic1_lord+branch%ele(ix_slave)%n_lord-1
  lat%ic(i2) = k+ixc-1
enddo

lat%n_control_max = n_con

! order slave elements in the super_lord list to be in the correct order

call s_calc (lat)  ! just in case superimpose extended before beginning of lattice.
call order_super_lord_slaves (lat, ix_super)

!-------------------------------------------------------------------------
! If create_jumbo_slave = T:

if (logic_option(.false., create_jumbo_slave)) then

  lat%ele(ix_super)%s = super_saved%s ! correct s-shift from call to s_calc 
  lat%ele(ix_super)%s_start = lat%ele(ix_super)%s - lat%ele(ix_super)%value(l$)

  ! All slave elements that do not have a zero length element in between 
  ! them have to be combined into one. 

  nullify (ele0)
  i = ix1_split
  do
    i = i + 1
    if (i == ix2_split + 1) exit
    if (i == branch%n_ele_track + 1) i = 1

    if (branch%ele(i)%value(l$) == 0) cycle

    ! Mark first element in a "string".
    if (.not. associated (ele0)) then
      if (branch%ele(i)%value(l$) == 0) cycle ! Skip zero length elements
      ele0 => branch%ele(i)
    endif


    ! If at end of string then mark all elements but ele0 for deletion.

    if (i == ix2_split .or. i == branch%n_ele_track .or. branch%ele(i+1)%value(l$) == 0) then
      ele0%value(l$) = branch%ele(i)%s - branch%ele(ele0%ix_ele-1)%s
      ele0%s = branch%ele(i)%s
      ele0%s_start = ele0%s - ele0%value(l$)

      do j = ele0%ix_ele+1, i
        ele => branch%ele(j)
        ele%key = -1

        ! Check lords and make sure that ele0 is a super_slave of the lord.
        ! If not, make ele0 a super_slave.
        do k = 1, ele%n_lord
          lord => pointer_to_lord(ele, k)
          if (lord%lord_status /= super_lord$) cycle
          slave => pointer_to_slave(lord, 1)
          if (slave%ix_ele <= ele0%ix_ele) cycle

          ! Need to make ele0 a super_slave of this lord.
          call add_lattice_control_structs(ele0, n_add_lord = 1)
          call add_lattice_control_structs(lord, n_add_slave = 1, add_at_end = .false.)
          
          cntl => lat%control(lord%ix1_slave)
          cntl%lord%ix_ele = lord%ix_ele
          cntl%slave%ix_ele = ele0%ix_ele
          cntl%slave%ix_branch = ele0%ix_branch
        enddo
      enddo
      nullify(ele0)
    endif

  enddo

  ! Compute lord_pad2 and lord_pad1 values.

  i = ix1_split
  do
    i = i + 1
    if (i == ix2_split + 1) exit
    if (i == branch%n_ele_track + 1) i = 1
    ele => branch%ele(i)

    do k = 1, ele%n_lord
      lord => pointer_to_lord(ele, k)

      slave => pointer_to_slave(lord, 1) 
      lord%value(lord_pad1$) = lord%s_start - slave%s_start
      if (abs(lord%value(lord_pad1$)) < bmad_com%significant_length) lord%value(lord_pad1$) = 0

      slave => pointer_to_slave(lord, lord%n_slave) 
      lord%value(lord_pad2$) = slave%s - lord%s
      if (abs(lord%value(lord_pad2$)) < bmad_com%significant_length) lord%value(lord_pad2$) = 0

    enddo
  enddo

  ! And remove unwanted super_slave elements

  call remove_eles_from_lat (lat, .false.)

endif

!-------------------------------------------------------------------------
! Bookkeeping and adjust names
! Do name adjust first in case there is a bookkeeping error

call adjust_super_slave_names (lat, ix_lord_max_old+1, lat%n_ele_max)
call adjust_drift_names (lat, branch%ele(ix1_split))
call adjust_drift_names (lat, branch%ele(ix2_split+1))

call control_bookkeeper(lat, lat%ele(ix_super))

err_flag = .false.

end subroutine add_superimpose

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

function split_this_lat (which, branch, s_here, ix_split, split_done, &
                                  save_null_drift, create_jumbo_slave) result (err)

implicit none

type (branch_struct) branch

real(rp) s_here

integer ix_split, which

logical split_done, err, choose_max
logical, optional :: save_null_drift, create_jumbo_slave

! Try to split so that the minimum number of elements are to be superimposed upon.

if (which == 1) then
  choose_max = .true.
else
  choose_max = .false.
endif

! If creating a jumbo slave then only split at drift elements

if (logic_option(.false., create_jumbo_slave)) then
  split_done = .false.
  ix_split = element_at_s(branch%lat, s_here, choose_max, branch%ix_branch, err)
  if (err) return
  if (branch%ele(ix_split)%key /= drift$) then
    if (which == 1) ix_split = ix_split - 1
    return
  endif
endif

call split_lat (branch%lat, s_here, branch%ix_branch, ix_split, split_done, &
                                                            .false., .false., save_null_drift, err, choose_max = choose_max)

end function split_this_lat 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Modify: "#\" -> "\"
!         "##" -> "#"

subroutine delete_underscore(ele)

use bmad_struct

implicit none

type (ele_struct) ele
integer ix

!

ix = index(ele%name, '#\')  ! '
if (ix /= 0) ele%name = ele%name(1:ix-1) // ele%name(ix+1:)

ix = index(ele%name, '##')
if (ix /= 0) ele%name = ele%name(1:ix-1) // ele%name(ix+1:)

end subroutine delete_underscore

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine adjust_super_slave_names (lat, ix1_lord, ix2_lord, first_time)
!
! Routine to adjust the names of the slaves.
! This routine is used by add_superimpose and is not meant for general use.
!-

recursive subroutine adjust_super_slave_names (lat, ix1_lord, ix2_lord, first_time)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: lord, slave, lord2
integer ix1_lord, ix2_lord
integer i, j, k, ix, n_unique

character(40), allocatable :: slave_names(:)
character(100) name
integer, allocatable :: n_slave_names(:), ix_slave_names(:)

logical, optional :: first_time

! First time through...
! The idea is to prevent an infinite circle by making sure that an element is not
! touched more than once. So initially mark all elements as not touched.

if (.not. present(first_time)) then
  do i = 0, ubound(lat%branch, 1)
    lat%branch(i)%ele%bmad_logic = .false.  ! Have adjusted?
  enddo
endif

!

do i = ix1_lord, ix2_lord
  lord => lat%ele(i)
  if (lord%lord_status /= super_lord$) cycle
  if (lord%bmad_logic) cycle
  lord%bmad_logic = .true.

  allocate (slave_names(lord%n_slave), n_slave_names(lord%n_slave), ix_slave_names(lord%n_slave))
  n_unique = 0
  n_slave_names = 0

  slave_loop: do j = 1, lord%n_slave
    slave => pointer_to_slave(lord, j)
    if (slave%bmad_logic) cycle
    slave%bmad_logic = .true.

    name = ''
    do k = 1, slave%n_lord
      lord2 => pointer_to_lord(slave, k)
      name = trim(name) //  '\' // lord2%name     !'
    enddo
    slave%name = name(2:len(slave%name)+1)

    do k = 1, n_unique
      if (slave%name == slave_names(k)) exit
    enddo

    n_unique = max (n_unique, k)
    slave_names(k) = slave%name
    n_slave_names(k) = n_slave_names(k) + 1

  enddo slave_loop

  !

  ix_slave_names = 0
  do j = 1, lord%n_slave
    slave => pointer_to_slave(lord, j)

    do k = 1, n_unique
      if (slave%name == slave_names(k)) then
        if (index(slave%name, '\') /= 0 .and. n_slave_names(k) == 1) exit       !'
        ix = min(len_trim(slave%name), len(slave%name) - 2)
        ix_slave_names(k) = ix_slave_names(k) + 1
        write (slave%name, '(2a, i0)') slave%name(1:ix), '#', ix_slave_names(k)
      endif
    enddo
    
    do k = 1, slave%n_lord
      lord2 => pointer_to_lord(slave, k)
      if (.not. lord2%bmad_logic) call adjust_super_slave_names (lat, lord2%ix_ele, lord2%ix_ele, .false.)
    enddo

  enddo

  deallocate (slave_names, n_slave_names, ix_slave_names)

enddo

end subroutine adjust_super_slave_names

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine adjust_drift_names (lat, drift_ele)
!
! If drifts have been split due to superimpose then the split drifts
! can have some very convoluted names. Also multipass lord drifts. Here we do some cleanup.
! Collect all free elements with the same name before the '#' and renumber.
!-

subroutine adjust_drift_names (lat, drift_ele)

implicit none

type (lat_struct), target :: lat
type (ele_struct) drift_ele
type (ele_struct), pointer :: ele, slave
type (branch_struct), pointer :: branch

real(rp) drift_id
integer i, ib, k, ixx

character(40) d_name

!

if (drift_ele%key /= drift$) return
drift_id = drift_ele%value(drift_id$)
if (drift_id == 0) return

d_name = drift_ele%name
i = index(d_name, '#')
if (i /= 0) d_name = d_name(:i-1)

!

ixx = 0
do ib = 0, ubound(lat%branch, 1)

  branch => lat%branch(ib)
  do i = 1, branch%n_ele_max

    ele => branch%ele(i)
    if ((ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) .and. ele%lord_status /= multipass_lord$) cycle
    if (ele%value(drift_id$) /= drift_id) cycle

    ixx = ixx + 1
    write (ele%name, '(2a, i0)') trim(d_name), '#', ixx

    if (ele%lord_status == multipass_lord$) then
      do k = 1, ele%n_slave
        slave => pointer_to_slave(ele, k)
        write (slave%name, '(2a, i0)') trim(ele%name), '\', k  ! '
      enddo
    endif

  enddo
enddo 

end subroutine adjust_drift_names

end module
