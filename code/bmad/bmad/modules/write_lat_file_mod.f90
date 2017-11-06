module write_lat_file_mod

use multipole_mod
use multipass_mod
use element_modeling_mod
use binary_parser_mod

private re_str, rchomp, cmplx_re_str, write_line_element, array_re_str
private write_lat_in_sad_format

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_bmad_lattice_file (bmad_file, lat, err, output_form)
!
! Subroutine to write a Bmad lattice file using the information in
! a lat_struct. Optionally only part of the lattice can be generated.
! Also see: write_lattice_in_foreign_format
!
! Modules needed:
!   use write_lat_file_mod
!
! Input:
!   bmad_file     -- Character(*): Name of the output lattice file.
!   lat           -- lat_struct: Holds the lattice information.
!   ix_start      -- Integer, optional: Starting index of lat%ele(i)
!                       used for output.
!   ix_end        -- Integer, optional: Ending index of lat%ele(i)
!                       used for output.
!   output_form   -- integer, optional: 
!                       binary$ -> Write field info in binary form. All else ascii. Default.
!                       ascii$  -> All ascii.
!
! Output:
!   err    -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_bmad_lattice_file (bmad_file, lat, err, output_form)

implicit none

type multipass_region_ele_struct
  integer ix_region
  logical region_start_pt
  logical region_stop_pt
end type

type multipass_region_branch_struct
  type (multipass_region_ele_struct), allocatable :: ele(:)
end type

type multipass_region_lat_struct
  type (multipass_region_branch_struct), allocatable :: branch(:)
end type

type (multipass_region_lat_struct), target :: m_region
type (multipass_region_ele_struct), pointer :: e_region(:)

type (ele_attribute_struct) attrib
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch, branch2
type (ele_struct), pointer :: ele, super, slave, lord, lord2, s1, s2, multi_lord, slave2, ele2, ele_dflt, ele0
type (ele_struct), target :: ele_default(n_key$)
type (ele_pointer_struct), allocatable :: named_eles(:)
type (control_struct), pointer :: ctl, ctl2
type (taylor_term_struct) tm
type (multipass_all_info_struct), target :: m_info
type (multipass_ele_info_struct), pointer :: e_info
type (wake_lr_mode_struct), pointer :: lr
type (wake_sr_mode_struct), parameter :: sr0 = wake_sr_mode_struct()
type (wake_sr_mode_struct), pointer :: sr
type (ele_pointer_struct), pointer :: ss1(:), ss2(:)
type (cylindrical_map_struct), pointer :: cl_map
type (cartesian_map_struct), pointer :: ct_map
type (cartesian_map_term1_struct), pointer :: ct_term
type (grid_field_struct), pointer :: g_field
type (taylor_field_struct), pointer :: t_field
type (em_taylor_term_struct), pointer :: t_term
type (wall3d_section_struct), pointer :: section
type (wall3d_vertex_struct), pointer :: v
type (bmad_common_struct), parameter :: bmad_com_default = bmad_common_struct()
type (ac_kicker_struct), pointer :: ac
type (expression_atom_struct), pointer :: stack(:)
type (str_indexx_struct) str_index

real(rp) s0, x_lim, y_lim, val

character(*) bmad_file
character(4000) line
character(2000) line2
character(200) wake_name, file_name, path, basename
character(200), allocatable :: sr_wake_name(:), lr_wake_name(:)
character(100) string
character(60) alias
character(40) name, look_for, attrib_name
character(40), allocatable :: names(:)
character(16) polar, dependence
character(40) angle
character(4) end_str, last
character(1), parameter :: xyz(3) = ['x', 'y', 'z']
character(*), parameter :: r_name = 'write_bmad_lattice_file'

integer, optional :: output_form
integer i, j, k, n, ix, iu, im, ix_ptr, iu2, iuw, ios, ixs, n_sr, n_lr, ix1, ie, ib, ic
integer unit(6), n_names, ix_match, ie2, id1, id2, id3, j1, j2, ip, it
integer ix_slave, ix_ss, ix_l, ix_r, ix_pass
integer ix_lord, ix_super, default_val, imax, ibr
integer, allocatable :: an_indexx(:)

logical, optional :: err
logical unit_found, write_term, found, in_multi_region, expand_branch_out
logical x_lim_good, y_lim_good, is_default, need_new_region, err_flag, has_been_added

! Init...
! Init default parameters

do i = 1, size(ele_default)
  call init_ele (ele_default(i), i)
  call set_ele_defaults (ele_default(i))
enddo

! Count the number of foreign wake files

if (present(err)) err = .true.

n_sr = 0
n_lr = 0
do ie = 1, lat%n_ele_max
  ele => lat%ele(ie)
  if (.not. associated(ele%wake)) cycle
  if (ele%wake%sr_file(1:6) == 'xsif::') n_sr = n_sr + 1 
  if (ele%wake%lr_file(1:6) == 'xsif::') n_lr = n_lr + 1  
enddo
call re_allocate(sr_wake_name, n_sr)
call re_allocate(lr_wake_name, n_lr)

n_sr = 0
n_lr = 0

! Open the file

iu = lunget()
call fullfilename (bmad_file, file_name)
open (iu, file = file_name, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(bmad_file))
  return
endif

ix = splitfilename(file_name, path, basename)
if (path == '') path = '.'

! Attribute aliases

if (allocated(lat%attribute_alias)) then
  do i = 1, size(lat%attribute_alias)
    alias = lat%attribute_alias(i)
    ix = index(alias, '=')
    write (iu, '(4a)') 'parameter[', alias(1:ix-1), '] = ', alias(ix+1:) 
  enddo
endif

! Non-elemental stuff

if (lat%title /= ' ')            write (iu, '(4a)')    'title, "', trim(lat%title), '"'
if (lat%lattice /= ' ')          write (iu, '(4a)')    'parameter[lattice]     = "', trim(lat%lattice), '"'

write (iu, '(4a)') 'parameter[geometry] = ', geometry_name(lat%param%geometry)

if (lat%input_taylor_order /= 0) write (iu, '(a, i0)') 'parameter[taylor_order] = ', lat%input_taylor_order

write (iu, '(a)')
write (iu, '(4a)')    'parameter[p0c]                    = ', trim(re_str(lat%ele(0)%value(p0c_start$)))
write (iu, '(4a)')    'parameter[particle]               = ', trim(species_name(lat%param%particle))

if (lat%param%n_part /= 0)             write (iu, '(2a)') 'parameter[n_part]                 = ', trim(re_str(lat%param%n_part))

write (iu, '(a, l1)') 'parameter[absolute_time_tracking]    = ', lat%absolute_time_tracking
ele => lat%ele(lat%n_ele_track)
if (ele%name /= 'END' .or. ele%key /= marker$) then
  write (iu, '(a)') 'parameter[no_end_marker]          =  T'
endif

if (lat%photon_type /= incoherent$) then
  write (iu, '(3a)') 'parameter[photon_type] = ', photon_type_name(lat%photon_type)
endif

if (bmad_com%aperture_limit_on .neqv. bmad_com_default%aperture_limit_on) &
            write (iu, '(a, l1)') 'bmad_com[aperture_limit_on] = ', bmad_com%aperture_limit_on

if (bmad_com%use_hard_edge_drifts .neqv. bmad_com_default%use_hard_edge_drifts) &
            write (iu, '(a, l1)') 'bmad_com[use_hard_edge_drifts] = ', bmad_com%use_hard_edge_drifts

if (abs(bmad_com%ptc_cut_factor - bmad_com_default%ptc_cut_factor) >  bmad_com_default%ptc_cut_factor/1d6) &
            write (iu, '(a, l1)') 'bmad_com[ptc_cut_factor] = ', bmad_com%ptc_cut_factor

if (bmad_com%electric_dipole_moment /= 0) &
            write (iu, '(a, l1)') 'bmad_com[electric_dipole_moment] = ', bmad_com%electric_dipole_moment

ele => lat%ele(0) 

if (ele%floor%r(1) /= 0)   write (iu, '(2a)') 'beginning[x_position]     = ', trim(re_str(ele%floor%r(1)))
if (ele%floor%r(2) /= 0)   write (iu, '(2a)') 'beginning[y_position]     = ', trim(re_str(ele%floor%r(2)))
if (ele%floor%r(3) /= 0)   write (iu, '(2a)') 'beginning[z_position]     = ', trim(re_str(ele%floor%r(3)))
if (ele%floor%theta /= 0)  write (iu, '(2a)') 'beginning[theta_position] = ', trim(re_str(ele%floor%theta))
if (ele%floor%phi /= 0)    write (iu, '(2a)') 'beginning[phi_position]   = ', trim(re_str(ele%floor%phi))
if (ele%floor%psi /= 0)    write (iu, '(2a)') 'beginning[psi_position]   = ', trim(re_str(ele%floor%psi))

if (ele%s /= 0)            write (iu, '(2a)') 'beginning[s]        = ', trim(re_str(ele%s))
if (ele%ref_time /= 0)     write (iu, '(2a)') 'beginning[ref_time] = ', trim(re_str(ele%ref_time))

if (lat%param%geometry /= closed$) then
  write (iu, '(2a)')
  if (ele%a%beta /= 0)     write (iu, '(2a)') 'beginning[beta_a]   = ', trim(re_str(ele%a%beta))
  if (ele%a%alpha /= 0)    write (iu, '(2a)') 'beginning[alpha_a]  = ', trim(re_str(ele%a%alpha))
  if (ele%a%phi /= 0)      write (iu, '(2a)') 'beginning[phi_a]    = ', trim(re_str(ele%a%phi))
  if (ele%x%eta /= 0)      write (iu, '(2a)') 'beginning[eta_x]    = ', trim(re_str(ele%x%eta))
  if (ele%x%etap /= 0)     write (iu, '(2a)') 'beginning[etap_x]   = ', trim(re_str(ele%x%etap))
  if (ele%b%beta /= 0)     write (iu, '(2a)') 'beginning[beta_b]   = ', trim(re_str(ele%b%beta))
  if (ele%b%alpha /= 0)    write (iu, '(2a)') 'beginning[alpha_b]  = ', trim(re_str(ele%b%alpha))
  if (ele%b%phi /= 0)      write (iu, '(2a)') 'beginning[phi_b]    = ', trim(re_str(ele%b%phi))
  if (ele%y%eta /= 0)      write (iu, '(2a)') 'beginning[eta_y]    = ', trim(re_str(ele%y%eta))
  if (ele%y%etap /= 0)     write (iu, '(2a)') 'beginning[etap_y]   = ', trim(re_str(ele%y%etap))
  if (ele%c_mat(1,1) /= 0) write (iu, '(2a)') 'beginning[cmat_11]  = ', trim(re_str(ele%c_mat(1,1)))
  if (ele%c_mat(1,2) /= 0) write (iu, '(2a)') 'beginning[cmat_12]  = ', trim(re_str(ele%c_mat(1,2)))
  if (ele%c_mat(2,1) /= 0) write (iu, '(2a)') 'beginning[cmat_21]  = ', trim(re_str(ele%c_mat(2,1)))
  if (ele%c_mat(2,2) /= 0) write (iu, '(2a)') 'beginning[cmat_22]  = ', trim(re_str(ele%c_mat(2,2)))
endif

! beam_start

if (lat%beam_start%vec(1) /= 0) write (iu, '(2a)') 'beam_start[x]  = ', trim(re_str(lat%beam_start%vec(1)))
if (lat%beam_start%vec(2) /= 0) write (iu, '(2a)') 'beam_start[px] = ', trim(re_str(lat%beam_start%vec(2)))
if (lat%beam_start%vec(3) /= 0) write (iu, '(2a)') 'beam_start[y]  = ', trim(re_str(lat%beam_start%vec(3)))
if (lat%beam_start%vec(4) /= 0) write (iu, '(2a)') 'beam_start[py] = ', trim(re_str(lat%beam_start%vec(4)))
if (lat%beam_start%vec(5) /= 0) write (iu, '(2a)') 'beam_start[z]  = ', trim(re_str(lat%beam_start%vec(5)))
if (lat%beam_start%vec(6) /= 0) write (iu, '(2a)') 'beam_start[pz] = ', trim(re_str(lat%beam_start%vec(6)))

if (lat%beam_start%spin(1) /= 0) write (iu, '(2a)') 'beam_start[spin_x] = ', trim(re_str(lat%beam_start%spin(1)))
if (lat%beam_start%spin(2) /= 0) write (iu, '(2a)') 'beam_start[spin_y] = ', trim(re_str(lat%beam_start%spin(2)))
if (lat%beam_start%spin(3) /= 1) write (iu, '(2a)') 'beam_start[spin_z] = ', trim(re_str(lat%beam_start%spin(3)))

! Named constants

write (iu, '(a)')

do i = 1, lat%n_control_max
  if (.not. allocated(lat%control(i)%stack)) cycle
  stack => lat%control(i)%stack
  do j = 1, size(stack)
    if (stack(j)%type == end_stack$) exit
    if (stack(j)%type /= variable$) cycle
    if (stack(j)%name == '') cycle
    if (any(stack(j)%name == physical_const_list%name)) cycle
    call find_indexx(stack(j)%name, str_index, ix, add_to_list = .true., has_been_added = has_been_added)
    if (.not. (has_been_added)) cycle  ! Avoid duuplicates
    write (iu, '(3a)') trim(stack(j)%name), ' = ', trim(re_str(stack(j)%value))
  enddo
enddo

! Element stuff

write (iu, '(a)')
write (iu, '(a)') '!-------------------------------------------------------'
write (iu, '(a)')

ixs = 0
n_names = 0
n = lat%n_ele_max
allocate (names(n), an_indexx(n), named_eles(n))

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  if (ib > 0) then
    write (iu, '(a)')
    write (iu, '(a)')  '!-------------------------------------------------------'
    write (iu, '(2a)') '! Branch: ', trim(branch%name)
    write (iu, '(a)')
  endif

  branch%ele%iyy = 0

  ele_loop: do ie = 1, branch%n_ele_max

    if (ie == branch%n_ele_track+1) then
      write (iu, '(a)')
      write (iu, '(a)') '!-------------------------------------------------------'
      write (iu, '(a)') '! Overlays, groups, etc.'
      write (iu, '(a)')
    endif

    ele => branch%ele(ie)
    if (ie == ele%branch%n_ele_track .and. ele%name == 'END' .and. ele%key == marker$) cycle

    ele_dflt => ele_default(ele%key)

    multi_lord => pointer_to_multipass_lord (ele, ix_pass) 

    ! Superposition stragegy: There is a problem where add_superposition will be confused by having a 
    ! superposition in a region next to an element (eg patch) with a negative length. 
    ! To get around this, try to minimize the number of superpositions needed.
    ! So don't superimpose any super_lord if the first slave of the lord does not have another super_lord.
    ! Otherwise create a dummy drift to be used for superposition. 

    ! ele%iyy = 1  -> ele is the first slave to a superlord that is *not* superimposed in the new lattice.
    ! ele%iyy = 2  -> ele is slave (but not the first) to a superlord that is *not* superimposed in the new lattice.
    ! ele%iyy = 3  -> ele has zero length and is not superimposed but it will be superimposed so ignore it in the
    !                    element line list.
    ! ele%iyy = 10 -> ele is a super_lord that is *not* superimposed in the new lattice.
    ! ele%iyy = 20 -> ele is a multipass_lord that is *not* superimposed but needs to be superimposed.

    ! A zero length element that has superposition on both sides gets superimposed to avoid problems
    ! if the element is within the region of a super_lord that is *not* being superimposed.

    if (ele%slave_status == multipass_slave$) then
      if (ele%value(l$) == 0 .and. ele%ix_ele /= branch%n_ele_max) then
        if (branch%ele(ie-1)%slave_status == super_slave$ .and. branch%ele(ie+1)%slave_status == super_slave$) then
          ele%iyy = 3   ! Mark so that element will not be in the line list of elements.
          lord => pointer_to_lord(ele, 1)
          lord%iyy = 20          
        endif
      endif

      cycle 
    endif
 
    !

    if (ele%key == null_ele$) cycle
    if (ele%lord_status == super_lord$ .and. ix_pass > 0) cycle
    if (ele%slave_status == super_slave$ .and. ix_pass > 1) cycle
    if (ele%slave_status == super_slave$ .and. ele%iyy == 2) cycle

    ! super_slave bookkeeping

    if (ele%slave_status == super_slave$) then
      lord => pointer_to_lord(ele, 1)
      slave => pointer_to_slave(lord, 1)
      if (num_lords(slave, super_lord$) == 1) then
        lord%iyy = 10
        slave%iyy = 1
        do i = 2, lord%n_slave
          slave => pointer_to_slave(lord, i)
          slave%iyy = 2
        enddo
      elseif (ele%iyy == 0) then
        ixs = ixs + 1
        ele%ixx = ixs
        write (iu, '(a, i0, 2a)') 'slave_drift_', ixs, ': drift, l = ', trim(re_str(ele%value(l$)))
      endif
      cycle
    endif

    ! Do not write anything for elements that have a duplicate name.

    call find_indexx (ele%name, names, an_indexx, n_names, ix_match)
    if (ix_match > 0) then
      ! Only check match elements for now.
      if (ele%key == match$) then
        if (any(named_eles(ix_match)%ele%value /= ele%value)) then
          call out_io (s_fatal$, r_name, 'TWO MATCH ELEMENTS WITH THE SAME NAME HAVE DIFFERENT PARAMETERS! ' // ele%name)
          return
        endif
      endif
      cycle
    endif

    if (size(names) < n_names + 1) then
      call re_allocate(names, 2*size(names))
      call re_allocate(an_indexx, 2*size(names))
      call re_allocate_eles(named_eles, 2*size(names))
    endif
    call find_indexx (ele%name, names, an_indexx, n_names, ix_match, add_to_list = .true.)
    named_eles(n_names)%ele => ele

    ! Overlays and groups

    if (ele%key == overlay$ .or. ele%key == group$) then
      if (ele%key == overlay$) then
        write (line, '(2a)') trim(ele%name), ': overlay = {'
      else
        write (line, '(2a)') trim(ele%name), ': group = {'
      endif
      j_loop: do j = 1, ele%n_slave
        slave => pointer_to_slave(ele, j, ctl)
        ! do not use elements w/ duplicate names & attributes
        do k = 1, j-1 
          slave2 => pointer_to_slave(ele, k, ctl2)
          if (slave2%name == slave%name .and. ctl2%ix_attrib == ctl%ix_attrib) cycle j_loop
        enddo
        ! Now write the slave info
        if (j == 1) then
          write (line, '(3a)') trim(line), trim(slave%name)
        else
          write (line, '(3a)') trim(line), ', ', trim(slave%name)
        endif
        name = attribute_name(slave, ctl%ix_attrib)  
        if (name /= ele%control_var(1)%name) line = trim(line) // '[' // trim(name) // ']'
        string = expression_stack_to_string(ctl%stack)
        if (string /= ele%control_var(1)%name) write (line, '(3a)') trim(line), ':', trim(string)

        if (len_trim(line) > len(line)/2) call write_lat_line(line, iu, .false.)
      enddo j_loop
      line = trim(line) // '}, var = {' // ele%control_var(1)%name

      do j = 2, size(ele%control_var)
        line = trim(line) // ', ' // ele%control_var(j)%name
      enddo

      line = trim(line) // '}'

      do j = 1, size(ele%control_var)
        if (ele%control_var(j)%value /= 0) then
          line = trim(line) // ', ' // trim(ele%control_var(j)%name) // ' = ' // trim(re_str(ele%control_var(j)%value))
        endif
        if (ele%control_var(j)%old_value /= 0) then
          line = trim(line) // ', old_' // trim(ele%control_var(j)%name) // ' = ' // trim(re_str(ele%control_var(j)%value))
        endif
      enddo

      if (is_false(ele%value(gang$))) line = trim(line) // ', gang = False'
      if (ele%type /= ' ') line = trim(line) // ', type = "' // trim(ele%type) // '"'
      if (ele%alias /= ' ') line = trim(line) // ', alias = "' // trim(ele%alias) // '"'
      if (associated(ele%descrip)) line = trim(line) // ', descrip = "' // trim(ele%descrip) // '"'
      call write_lat_line (line, iu, .true.)
      cycle
    endif

    ! Girder

    if (ele%key == girder$) then
      write (line, '(2a)') trim(ele%name), ': girder = {'
      do j = 1, ele%n_slave
        slave => pointer_to_slave(ele, j)
        if (j == ele%n_slave) then
          write (line, '(3a)') trim(line), trim(slave%name), '}'
        else
          write (line, '(3a)') trim(line), trim(slave%name), ', '
        endif
      enddo

    ! Everything but a girder

    else
      line = trim(ele%name) // ': ' // key_name(ele%key)
    endif

    ! Branch

    if (ele%key == fork$ .or. ele%key == photon_fork$) then
      n = nint(ele%value(ix_to_branch$))
      line = trim(line) // ', to_line = ' // trim(lat%branch(n)%name)
      if (ele%value(ix_to_element$) > 0) then
        i = nint(ele%value(ix_to_element$))
        line = trim(line) // ', to_element = ' // trim(lat%branch(n)%ele(i)%name)
      endif
    endif

    ! Other elements

    if (ele%type /= ' ') line = trim(line) // ', type = "' // trim(ele%type) // '"'
    if (ele%alias /= ' ') line = trim(line) // ', alias = "' // trim(ele%alias) // '"'
    if (associated(ele%descrip)) line = trim(line) // ', descrip = "' // trim(ele%descrip) // '"'

    ! Fill in the superposition information.

    if (ele%lord_status == multipass_lord$) then
      slave => pointer_to_slave(ele, 1)
      if (slave%lord_status == super_lord$) then
        slave2 => pointer_to_slave(slave, 1)
        if (num_lords(slave2, super_lord$) /= 1) then
          do i = 1, slave2%n_lord
            super => pointer_to_lord(slave2, i)
            if (super%iyy == 10) exit
          enddo
          lord2 => pointer_to_lord(super, 1) ! Multipass_lord
          line = trim(line) // ', superimpose, ele_origin = beginning, ref_origin = beginning, ref = ' // &
                trim(lord2%name) // ', offset = ' // trim(re_str(slave%s_start-super%s_start))
        endif
      endif
    endif

    if (ele%lord_status == super_lord$) then
      slave => pointer_to_slave(ele, 1)
      if (num_lords(slave, super_lord$) /= 1) then
        do i = 1, slave%n_lord
          super => pointer_to_lord(slave, i)
          if (super%iyy == 10) exit
        enddo
        line = trim(line) // ', superimpose, ele_origin = beginning, ref_origin = beginning, ref = ' // &
                trim(super%name) // ', offset = ' // trim(re_str(ele%s_start-super%s_start))
      endif
    endif

    ! A zero length element that has superposition on both sides gets superimposed to avoid problems
    ! if the element is within the region of a super_lord that is *not* being superimposed.

    if (ele%value(l$) == 0 .and. ele%ix_ele /= branch%n_ele_max .and. ele%slave_status /= super_slave$) then
      if (branch%ele(ie-1)%slave_status == super_slave$ .and. branch%ele(ie+1)%slave_status == super_slave$) then
        ele%iyy = 3   ! Mark so that element will not be in the line list of elements.
        ele2 => branch%ele(ie-1)
        do i = 1, ele2%n_lord
          lord => pointer_to_lord(ele2, i)
          if (lord%iyy /= 10) cycle
          line = trim(line) // ', superimpose, ref_origin = beginning, ref = ' // &
                trim(lord%name) // ', offset = ' // trim(re_str(ele%s_start-lord%s_start))
          exit
        enddo
      endif
    endif

    if (ele%iyy == 20) then
      slave => pointer_to_slave(ele, 1)
      ele2 => pointer_to_next_ele(slave, -1)
      do i = 1, ele2%n_lord
        lord => pointer_to_lord(ele2, i)
        if (lord%iyy /= 10) cycle
        lord2 => lord
        if (lord2%slave_status == multipass_slave$) lord2 => pointer_to_lord(lord2, 1)
        line = trim(line) // ', superimpose, ref_origin = beginning, ref = ' // &
              trim(lord2%name) // ', offset = ' // trim(re_str(slave%s_start-lord%s_start))
        exit
      enddo
    endif

    ! AC_Kicker

    if (associated(ele%ac_kick)) then
      ac => ele%ac_kick
      if (allocated(ac%amp_vs_time)) then
        line = trim(line) // ', amp_vs_time = {(' // trim(re_str(ac%amp_vs_time(1)%time)) // &
                             ', ' // trim(re_str(ac%amp_vs_time(1)%amp)) // ')'
        do i = 2, size(ac%amp_vs_time)
          line = trim(line) // ', (' // trim(re_str(ac%amp_vs_time(i)%time)) // &
                             ', ' // trim(re_str(ac%amp_vs_time(i)%amp)) // ')'
        enddo
        line = trim(line) // '}'

      else
        line = trim(line) // ', frequencies = {(' // trim(re_str(ac%frequencies(1)%f)) // &
                  ', ' // trim(re_str(ac%frequencies(1)%amp)) // ', ' // trim(re_str(ac%frequencies(1)%phi))  // ')'
        do i = 2, size(ac%frequencies)
          line = trim(line) // ', (' // trim(re_str(ac%frequencies(i)%f)) // &
                  ', ' // trim(re_str(ac%frequencies(i)%amp)) // ', ' // trim(re_str(ac%frequencies(i)%phi)) // ')'
        enddo
        line = trim(line) // '}'
      endif
    endif

    ! Wall3d

    if (associated(ele%wall3d)) then

      ! First find out out if a wall file has been written
      found = .false.
      do ibr = 0, ubound(lat%branch, 1)
        branch2 => lat%branch(ibr)
        imax = branch2%n_ele_max
        if (ibr == branch%ix_branch) imax = ie-1
        do ie2 = 1, imax
          ele2 => branch2%ele(ie2)
          if (ele2%slave_status == multipass_slave$) cycle
          if (.not. associated(ele2%wall3d)) cycle
          if (.not. associated(ele2%wall3d, ele%wall3d)) cycle
          found = .true.
          exit
        enddo
        if (found) exit
      enddo

      if (found) then
        call str_downcase(name, ele2%name)
        line = trim(line) // ', wall = call::wall_' // trim(name)
      else
        call str_downcase(name, ele%name)
        line = trim(line) // ', wall = call::wall_' // trim(name)
        iu2 = lunget()
 
        open (iu2, file = trim(path) // '/wall_' // trim(name))
        write (iu2, '(a)') '{'
        write (iu2, '(2x, 3a)') 'ele_anchor_pt = ', trim(anchor_pt_name(ele%wall3d(1)%ele_anchor_pt)), ','
        do i = 1, size(ele%wall3d(1)%section)
          section => ele%wall3d(1)%section(i)
          write (iu2, '(2x, a)')   'section = {'

          if (ele%key == diffraction_plate$ .or. ele%key == mask$) then
            write (iu2, '(4x, 3a)') 'type = ', trim(wall3d_section_type_name(section%type)), ','
          else
            write (iu2, '(4x, 3a)')  's     = ', trim(re_str(section%s)), ','
            if (section%dr_ds /= real_garbage$) write (iu2, '(4x, 3a)')  'dr_ds = ', trim(re_str(section%s)), ','
          endif

          end_str = ','
          do j = 1, size(section%v)
            if (j == size(section%v)) then
              end_str = '},'
              if (i == size(ele%wall3d(1)%section)) end_str = '}}'
            endif
            v => section%v(j)
            if (v%tilt /= 0) then
              write (iu2, '(4x, a, i0, 3a)') 'v(', j, ') = ', &
                    trim(array_re_str([v%x, v%y, v%radius_x, v%radius_y, v%tilt], '{}')), end_str
            elseif (v%radius_y /= 0) then
              write (iu2, '(4x, a, i0, 3a)') 'v(', j, ') = ', &
                    trim(array_re_str([v%x, v%y, v%radius_x, v%radius_y], '{}')), end_str
            elseif (v%radius_x /= 0) then
              write (iu2, '(4x, a, i0, 3a)') 'v(', j, ') = ', &
                    trim(array_re_str([v%x, v%y, v%radius_x], '{}')), end_str
            else
              write (iu2, '(4x, a, i0, 3a)') 'v(', j, ') = ', &
                    trim(array_re_str([v%x, v%y], '{}')), end_str
            endif
          enddo
        enddo
        close (iu2)

      endif
    endif

    ! field overlap

    if (ele%n_slave_field /= 0) then
      slave => pointer_to_slave (ele, 1, field_overlap_ptr = .true.)
      line = trim(line) // ', field_overlaps = {' // slave%name
      do n = 2, ele%n_slave_field
        slave => pointer_to_slave (ele, n, field_overlap_ptr = .true.)
        line = trim(line) // ', ' // slave%name
      enddo
      line = trim(line) // '}'
    endif

    ! cartesian_map. Periodic style wigglers do not have their map written out.

    if (associated(ele%cartesian_map) .and. .not. &
                      ((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%sub_key == periodic_type$)) then
      do im = 1, size(ele%cartesian_map)
        ct_map => ele%cartesian_map(im)

        ! First find out out if an file has been written
        call find_matching_fieldmap (ct_map%ptr%file, ele, cartesian_map$, ele2, ix_ptr) 

        if (ix_ptr > 0) then
          call form_this_field_map_name(string, '.cartesian_map', ele2, ix_ptr, output_form)
          write (line, '(3a)')  trim(line), ', cartesian_map = call::', trim(string)

        else
          call form_this_field_map_name(string, '.cartesian_map', ele, im, output_form)
          line = trim(line) // ', cartesian_map = call::' // trim(string)
          string = trim(path) // '/' // trim(string)

          if (integer_option(binary$, output_form) == binary$) then
            call write_binary_cartesian_map(string, ct_map, err_flag)
          else
            iu2 = lunget()
            open (iu2, file = string)

            write (iu2, '(a)') '{'
            if (ct_map%master_parameter > 0) write (iu2, '(2x, 3a)') &
                                    'master_parameter  = ', trim(attribute_name(ele, ct_map%master_parameter)), ','
            write (iu2, '(2x, 3a)') 'field_scale       = ', trim(re_str(ct_map%field_scale)), ','
            write (iu2, '(2x, 4a)') 'r0                = ', trim(array_re_str(ct_map%r0)), ','
            write (iu2, '(2x, 3a)') 'ele_anchor_pt     = ', trim(anchor_pt_name(ct_map%ele_anchor_pt)), ','
            write (iu2, '(2x, 3a)') 'field_type        = ', trim(em_field_type_name(ct_map%field_type)), ','

            do j = 1, size(ct_map%ptr%term)
              ct_term => ct_map%ptr%term(j)
              last = '}, &'
              if (j == size(ct_map%ptr%term)) last = '} &'
              select case (ct_term%type)
              case (hyper_y_family_x$, hyper_xy_family_x$, hyper_x_family_x$)
                name = 'X'
              case (hyper_y_family_y$, hyper_xy_family_y$, hyper_x_family_y$)
                name = 'Y'
              case (hyper_y_family_qu$, hyper_xy_family_qu$, hyper_x_family_qu$)
                name = 'QU'
              case (hyper_y_family_sq$, hyper_xy_family_sq$, hyper_x_family_sq$)
                name = 'SQ'
              end select
              write (iu2, '(17a)') '  term = {', trim(re_str(ct_term%coef)), ', ', &
                trim(re_str(ct_term%kx)), ', ', trim(re_str(ct_term%ky)), ', ', trim(re_str(ct_term%kz)), &
                ', ', trim(re_str(ct_term%x0)), ', ', trim(re_str(ct_term%y0)), ', ', trim(re_str(ct_term%phi_z)), ', ', trim(name), trim(last)
            enddo

            write (iu2, '(a)') '}'
            close (iu2)
          endif
        endif
      enddo
    endif

    ! cylindrical_map

    if (associated(ele%cylindrical_map)) then
      do im = 1, size(ele%cylindrical_map)
        cl_map => ele%cylindrical_map(im)

        call find_matching_fieldmap (cl_map%ptr%file, ele, cylindrical_map$, ele2, ix_ptr) 

        if (ix_ptr > 0) then
          call form_this_field_map_name(string, '.cylindrical_map', ele2, ix_ptr, output_form)
          write (line, '(3a)')  trim(line), ', cylindrical_map = call::', trim(string)

        else
          call form_this_field_map_name(string, '.cylindrical_map', ele, im, output_form)
          line = trim(line) // ', cylindrical_map = call::' // trim(string)
          string = trim(path) // '/' // trim(string)

          if (integer_option(binary$, output_form) == binary$) then
            call write_binary_cylindrical_map(string, cl_map, err_flag)
          else
            iu2 = lunget()
            open (iu2, file = string)

            write (iu2, '(a)') '{'
            if (cl_map%master_parameter > 0) write (iu2, '(2x, 3a)') &
                                          'master_parameter  = ', trim(attribute_name(ele, cl_map%master_parameter)), ','
            write (iu2, '(2x, 3a)')       'field_scale       = ', trim(re_str(cl_map%field_scale)), ','
            write (iu2, '(2x, 3a)')       'ele_anchor_pt     = ', trim(anchor_pt_name(cl_map%ele_anchor_pt)), ','
            write (iu2, '(2x, a, i0, a)') 'm                 = ', cl_map%m, ','
            write (iu2, '(2x, a, i0, a)') 'harmonic          = ', cl_map%harmonic, ','
            write (iu2, '(2x, 3a)')       'dz                = ', trim(re_str(cl_map%dz)), ','
            write (iu2, '(2x, 4a)')       'r0                = ', trim(array_re_str(cl_map%r0)), ','
            write (iu2, '(2x, 3a)')       'phi0_fieldmap     = ', trim(re_str(cl_map%phi0_fieldmap)), ','
            write (iu2, '(2x, 3a)', advance = 'NO') 'theta0_azimuth      = ', trim(re_str(cl_map%theta0_azimuth))

            if (any(real(cl_map%ptr%term%e_coef) /= 0)) call write_map_coef ('E_coef_re', real(cl_map%ptr%term%e_coef))
            if (any(aimag(cl_map%ptr%term%e_coef) /= 0)) call write_map_coef ('E_coef_im', aimag(cl_map%ptr%term%e_coef))
            if (any(real(cl_map%ptr%term%b_coef) /= 0)) call write_map_coef ('B_coef_re', real(cl_map%ptr%term%b_coef))
            if (any(aimag(cl_map%ptr%term%b_coef) /= 0)) call write_map_coef ('B_coef_im', aimag(cl_map%ptr%term%b_coef))

            write (iu2, '(2x, a)') '}'
            close (iu2)
          endif
        endif
      enddo
    endif

    ! grid_field

    if (associated(ele%grid_field)) then
      do im = 1, size(ele%grid_field)
        g_field => ele%grid_field(im)

        ! First find out out if an file has been written
        call find_matching_fieldmap (g_field%ptr%file, ele, grid_field$, ele2, ix_ptr) 

        if (ix_ptr > 0) then
          call form_this_field_map_name(string, '.grid_field', ele2, ix_ptr, output_form)
          write (line, '(3a)')  trim(line), ', grid_field = call::', trim(string)

        else
          call form_this_field_map_name(string, '.grid_field', ele, im, output_form)
          line = trim(line) // ', grid_field = call::' // trim(string)
          string = trim(path) // '/' // trim(string)

          if (integer_option(binary$, output_form) == binary$) then
            call write_binary_grid_field (string, g_field, err_flag)
          else
            iu2 = lunget()
            open (iu2, file = string)

            write (iu2, '(a)') '{'
            n = grid_field_dimension(g_field%geometry)
            write (iu2, '(2x, 3a)')       'geometry          = ', trim(grid_field_geometry_name(g_field%geometry)), ','
            if (g_field%master_parameter > 0) write (iu2, '(2x, 3a)') &
                                          'master_parameter  = ', trim(attribute_name(ele, g_field%master_parameter)), ','
            write (iu2, '(2x, 3a)')       'field_scale       = ', trim(re_str(g_field%field_scale)), ','
            write (iu2, '(2x, 3a)')       'ele_anchor_pt     = ', trim(anchor_pt_name(g_field%ele_anchor_pt)), ','
            write (iu2, '(2x, 3a)')       'field_type        = ', trim(em_field_type_name(g_field%field_type)), ','
            write (iu2, '(2x, a, i0, a)') 'harmonic          = ', g_field%harmonic, ','
            write (iu2, '(2x, 3a)')       'phi0_fieldmap     = ', trim(re_str(g_field%phi0_fieldmap)), ','
            write (iu2, '(2x, 4a)')       'dr                = ', trim(array_re_str(g_field%dr(1:n))), ','
            write (iu2, '(2x, 4a)')       'r0                = ', trim(array_re_str(g_field%r0(1:n))), ','
            write (iu2, '(2x, a, l1, a)') 'curved_ref_frame  = ', g_field%curved_ref_frame, ','

            end_str = '),'

            do id1 = lbound(g_field%ptr%pt, 1), ubound(g_field%ptr%pt, 1)          
            do id2 = lbound(g_field%ptr%pt, 2), ubound(g_field%ptr%pt, 2)
            do id3 = lbound(g_field%ptr%pt, 3), ubound(g_field%ptr%pt, 3)

              if (all([id1, id2, id3] == ubound(g_field%ptr%pt))) end_str = ') &'

              select case (grid_field_dimension(g_field%geometry))
              case (1)
                write (string, '(2x, a, i0, 13a)') 'pt(', id1, ') = ('
              case (2)
                write (string, '(2x, 2(a, i0), 13a)') 'pt(', id1, ',', id2, ') = ('
              case (3)
                write (string, '(2x, 3(a, i0), 13a)') 'pt(', id1, ',', id2, ',', id3, ') = ('
              end select

              select case (g_field%field_type)
              case (mixed$)
                write (iu2, '(2x, a, 13a)') trim(string), &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(1))), ',', &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(2))), ',', &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(3))), ',', &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(1))), ',', &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(2))), ',', &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(3))), end_str
              case (electric$)
                write (iu2, '(2x, a, 13a)') trim(string), &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(1))), ',', &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(2))), ',', &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(3))), end_str
              case (magnetic$)
                write (iu2, '(2x, a, 13a)') trim(string), &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(1))), ',', &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(2))), ',', &
                                                   trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(3))), end_str
              end select

            enddo
            enddo
            enddo

            write (iu2, '(4x, a)') '}'
            close (iu2)
          endif
        endif
      enddo
    endif

    ! Taylor_field

    if (associated(ele%taylor_field)) then
      do im = 1, size(ele%taylor_field)
        t_field => ele%taylor_field(im)

        ! First find out out if an file has been written
        call find_matching_fieldmap (t_field%ptr%file, ele, taylor_field$, ele2, ix_ptr) 

        if (ix_ptr > 0) then
          call form_this_field_map_name(string, '.taylor_field', ele2, ix_ptr, output_form)
          write (line, '(3a)')  trim(line), ', taylor_field = call::', trim(string)

        else
          call form_this_field_map_name(string, '.taylor_field', ele, im, output_form)
          line = trim(line) // ', taylor_field = call::' // trim(string)
          string = trim(path) // '/' // trim(string)

          if (integer_option(binary$, output_form) == binary$) then
            call write_binary_taylor_field (string, t_field, err_flag)
          else
            iu2 = lunget()
            open (iu2, file = string)

            write (iu2, '(a)') '{'
            if (t_field%master_parameter > 0) write (iu2, '(2x, 3a)') &
                                          'master_parameter   = ', trim(attribute_name(ele, t_field%master_parameter)), ','
            write (iu2, '(2x, 3a)')       'field_scale        = ', trim(re_str(t_field%field_scale)), ','
            write (iu2, '(2x, 3a)')       'ele_anchor_pt      = ', trim(anchor_pt_name(t_field%ele_anchor_pt)), ','
            write (iu2, '(2x, 3a)')       'field_type         = ', trim(em_field_type_name(t_field%field_type)), ','
            write (iu2, '(2x, 3a)')       'dz                 = ', trim(re_str(t_field%dz)), ','
            write (iu2, '(2x, 4a)')       'r0                 = ', trim(array_re_str(t_field%r0)), ','
            write (iu2, '(2x, a, l1, a)') 'curved_ref_frame   = ', t_field%curved_ref_frame, ','
            write (iu2, '(2x, a, l1, a)') 'canonical_tracking = ', t_field%canonical_tracking, ','

            do ip = lbound(t_field%ptr%plane, 1), ubound(t_field%ptr%plane, 1)
              write (iu2, '(2x, a, i0, a)') 'plane(', ip, ') = {'
              line2 = ''
              do k = 1, 3
                do it = 1, size(t_field%ptr%plane(ip)%field(k)%term)
                  t_term => t_field%ptr%plane(ip)%field(k)%term(it)
                  if (line2 == '') then
                    write (line2, '(4x, 5a, 2i2, a)') '{', field_plane_name(k), ': ', trim(re_str(t_term%coef)), ',', t_term%expn, '}'
                  else
                    write (line2, '(6a, 2i2, a)') trim(line2), ', {', field_plane_name(k), ': ', trim(re_str(t_term%coef)), ',', t_term%expn, '}'
                  endif
                enddo
              enddo
              line2 = trim(line2) // ' }'
              if (ip == ubound(t_field%ptr%plane, 1)) then
                line2 = trim(line2) // ' &'
              else
                line2 = trim(line2) // ','
              endif
              call write_lat_line (line2, iu2, .true.)
            enddo

            write (iu2, '(a)') '}'
            close (iu2)
          endif
        endif
      enddo
    endif

    ! If the wake file is not BMAD Format (Eg: XSIF format) then create a new wake file.
    ! If first three characters of the file name are '...' then it is a foreign file.

    if (associated(ele%wake)) then

      ! Short-range

      if (ele%wake%sr_file /= ' ') then

        wake_name = ele%wake%sr_file

        if (wake_name(1:3) == '...') then
          found = .false.
          do n = 1, n_sr
            if (wake_name == sr_wake_name(n)) then
              found = .true.
              exit
            endif
          enddo
          if (.not. found) then
            n = n_sr + 1
            n_sr = n
            sr_wake_name(n_sr) = wake_name
          endif
          write (wake_name, '(a, i0, a)') 'sr_wake_file_', n, '.bmad'
          if (.not. found) then
            call out_io (s_info$, r_name, 'Creating SR Wake file: ' // trim(wake_name))
            iuw = lunget()
            open (iuw, file = wake_name)
            write (iuw, '(a)') '! Pseudo Wake modes:'
            write (iuw, '(a)') '!                 Amp	  Damp    K   Phase    Polarization  Transverse_dependence'
            write (iuw, '(a)') '! Longitudinal: [V/C/m] [1/m] [1/m] [rad]'
            write (iuw, '(a)') '! Transverse: [V/C/m^2] [1/m] [1/m] [rad]'
            write (iuw, '(a)') ''
            write (iuw, '(a)') '&short_range_modes'
            do n = 1, size(ele%wake%sr_long%mode)
              sr => ele%wake%sr_long%mode(n)
              dependence = '' ! Use for default
              if (sr%transverse_dependence /= none$) dependence = sr_transverse_dependence_name(sr%transverse_dependence)
              polar = ''      ! Use for default
              if (sr%polarization /= sr0%polarization .or. dependence /= '') polar = sr_polarization_name(sr%polarization)
              write (iuw, '(a, i0, a, 4es16.8, 2x, 2a)') 'logitudinal(', n, ') =', sr%amp, sr%damp, sr%k, sr%phi, polar, dependence
            enddo

            write (iuw, '(a)') ''
            do n = 1, size(ele%wake%sr_trans%mode)
              sr => ele%wake%sr_trans%mode(n)
              dependence = '' ! Use for default
              if (sr%transverse_dependence /= linear_leading$) dependence = sr_transverse_dependence_name(sr%transverse_dependence)
              polar = ''      ! Use for default
              if (sr%polarization /= sr0%polarization .or. dependence /= '') polar = sr_polarization_name(sr%polarization)
              write (iuw, '(a, i0, a, 4es16.8, 2x, 2a)') 'transverse(', n, ') =', sr%amp, sr%damp, sr%k, sr%phi, polar, dependence
            enddo
            write (iuw, '(a)') '/'
            close(iuw)
          endif
        endif

        line = trim(line) // ',  sr_wake_file = "' // trim(wake_name) // '"'

      endif

      ! Long-range

      if (ele%wake%lr_file /= ' ') then

        wake_name = ele%wake%lr_file

        if (wake_name(1:3) == '...') then
          found = .false.
          do n = 1, n_lr
            if (wake_name == lr_wake_name(n)) then
              found = .true.
              exit
            endif
          enddo
          if (.not. found) then
            n = n_lr + 1
            n_lr = n
            lr_wake_name(n_lr) = wake_name
          endif
          write (wake_name, '(a, i0, a)') 'lr_wake_file_', n, '.bmad'
          if (.not. found) then
            call out_io (s_info$, r_name, 'Creating LR Wake file: ' // trim(wake_name))
            iuw = lunget()
            open (iuw, file = wake_name)
            write (iuw, '(14x, a)') &
   'Freq         R/Q        Q       m    Polarization     b_sin         b_cos         a_sin         a_cos         t_ref'
            write (iuw, '(14x, a)') &
              '[Hz]  [Ohm/m^(2m)]             [Rad/2pi]'
            do n = lbound(ele%wake%lr_mode, 1), ubound(ele%wake%lr_mode, 1)
              lr => ele%wake%lr_mode(n)
              if (lr%polarized) then
                angle = re_str(lr%angle)
              else
                angle = '     unpol'
              endif
              if (any ( [lr%b_sin, lr%b_cos, lr%a_sin, lr%a_cos, lr%t_ref ]/= 0)) then
                write (iuw, '(a, i0, a, 3es16.7, i6, a, 5es14.6)') 'lr(', n, ') =', lr%freq_in, &
                      lr%R_over_Q, lr%Q, lr%m, angle, lr%b_sin, lr%b_cos, lr%a_sin, lr%a_cos, lr%t_ref
              else
                write (iuw, '(a, i0, a, 3es16.8, i6, 2x, a)') 'lr(', n, ') =', &
                      lr%freq_in, lr%R_over_Q, lr%Q, lr%m, trim(angle)
              endif
            enddo
            close(iuw)
          endif
        endif

        line = trim(line) // ',  lr_wake_file = "' // trim(wake_name) // '"'

      endif

    endif

    ! Decide if x1_limit, etc. are to be output directly or combined. 

    x_lim = ele%value(x1_limit$) 
    x_lim_good = .false.
    if (x_lim /=0 .and. ele%value(x2_limit$) == x_lim) x_lim_good = .true.

    y_lim = ele%value(y1_limit$) 
    y_lim_good = .false.
    if (y_lim /=0 .and. ele%value(y2_limit$) == y_lim) y_lim_good = .true.

    !----------------------------------------------------------------------------
    ! Print the element attributes.

    attribute_loop: do j = 1, num_ele_attrib$
      attrib = attribute_info(ele, j)
      val = ele%value(j)
      if (val == ele_dflt%value(j)) cycle
      if (ele%key == sbend$) then
        if (j == fintx$ .and. ele%value(fintx$) == ele%value(fint$)) cycle
        if (j == hgapx$ .and. ele%value(hgapx$) == ele%value(hgap$)) cycle
      endif
      if (j == check_sum$) cycle
      if (x_lim_good .and. (j == x1_limit$ .or. j == x2_limit$)) cycle
      if (y_lim_good .and. (j == y1_limit$ .or. j == y2_limit$)) cycle
      if (.not. attribute_free (ele, attrib%name, .false., .true.)) cycle
      if ((attrib%name == 'P0C' .or. attrib%name == 'P0C_START') .and. &
                          (ele%lord_status /= multipass_lord$ .or. ele%value(n_ref_pass$) /= 0)) cycle
      if (attrib%name == 'DS_STEP' .and. val == bmad_com%default_ds_step) cycle
      if (attrib%name == 'E_TOT') cycle        ! Will use p0c instead.
      if (attrib%name == 'E_TOT_START') cycle  ! Will use p0c_start instead.
      if (attrib%name == null_name$) then
        call out_io (s_error$, r_name, 'ELEMENT: ' // ele%name, 'HAS AN UNKNOWN ATTRIBUTE INDEX: \i0\ ', i_array = [j])
        if (global_com%exit_on_error) call err_exit
        return
      endif

      if (attrib%name == 'COUPLER_AT') then
        if (nint(val) /= downstream_end$) then
          line = trim(line) // ', coupler_at = ' // end_at_name(nint(val))
        endif
        cycle
      endif

      select case (attribute_type(attrib%name))
      case (is_logical$)
        write (line, '(4a, l1)') trim(line), ', ', trim(attrib%name), ' = ', (val /= 0)
      case (is_integer$)
        write (line, '(4a, i0)') trim(line), ', ', trim(attrib%name), ' = ', int(val)
      case (is_real$)
        line = trim(line) // ', ' // trim(attrib%name) // ' = ' // re_str(val)
      case (is_switch$)
        name = switch_attrib_value_name (attrib%name, val, ele, is_default)
          if (.not. is_default) then
            line = trim(line) // ', ' // trim(attrib%name) // ' = ' // name
          endif
      end select

    enddo attribute_loop ! attribute loop

    !----------------------------------------------------------------------------
    ! Print the combined limits if needed.

    if (x_lim_good .and. y_lim_good .and. x_lim == y_lim) then
      line = trim(line) // ', aperture = ' // re_str(x_lim)
    else
      if (x_lim_good) line = trim(line) // ', x_limit = ' // re_str(x_lim)
      if (y_lim_good) line = trim(line) // ', y_limit = ' // re_str(y_lim)
    endif

    ! Encode methods, etc.

    if (has_attribute (ele, 'MAT6_CALC_METHOD') .and. (ele%mat6_calc_method /= ele_dflt%mat6_calc_method)) &
                                      line = trim(line) // ', mat6_calc_method = ' // mat6_calc_method_name(ele%mat6_calc_method)
    if (has_attribute (ele, 'TRACKING_METHOD') .and. (ele%tracking_method /= ele_dflt%tracking_method)) &
                                      line = trim(line) // ', tracking_method = ' // tracking_method_name(ele%tracking_method)
    if (has_attribute (ele, 'SPIN_TRACKING_METHOD') .and. (ele%spin_tracking_method /= ele_dflt%spin_tracking_method)) &
                                      line = trim(line) // ', spin_tracking_method = ' // spin_tracking_method_name(ele%spin_tracking_method)
    if (has_attribute (ele, 'PTC_INTEGRATION_TYPE') .and. (ele%ptc_integration_type /= ele_dflt%ptc_integration_type)) &
                                      line = trim(line) // ', ptc_integration_type = ' // ptc_integration_type_name(ele%ptc_integration_type)
    if (has_attribute (ele, 'FIELD_CALC') .and. (ele%field_calc /= ele_dflt%field_calc)) &
                                      line = trim(line) // ', field_calc = ' // field_calc_name(ele%field_calc)

    if (has_attribute (ele, 'APERTURE_AT') .and. (ele%aperture_at /= ele_dflt%aperture_at)) &
                                      line = trim(line) // ', aperture_at = ' // aperture_at_name(ele%aperture_at)
    if (has_attribute (ele, 'APERTURE_TYPE') .and. (ele%aperture_type /= ele_dflt%aperture_type)) &
                                      line = trim(line) // ', aperture_type = ' // aperture_type_name(ele%aperture_type)

    if (has_attribute (ele, 'SYMPLECTIFY') .and. ele%symplectify) line = trim(line) // ', symplectify'

    if (has_attribute (ele, 'FIELD_MASTER') .and. (ele%field_master .neqv. ele_dflt%field_master)) &
                                      write (line, '(2a, l1)') trim(line), ', field_master = ', ele%field_master
    if (has_attribute (ele, 'IS_ON') .and. (ele%is_on .neqv. ele_dflt%is_on)) &
                                      write (line, '(2a, l1)') trim(line), ', is_on = ', ele%is_on
    if (has_attribute (ele, 'SCALE_MULTIPOLES') .and. (ele%scale_multipoles .neqv. ele_dflt%scale_multipoles)) &
                                      write (line, '(2a, l1)') trim(line), ', scale_multipoles = ', ele%scale_multipoles
    if (has_attribute (ele, 'MULTIPOLES_ON') .and. (ele%multipoles_on .neqv. ele_dflt%multipoles_on)) &
                                      write (line, '(2a, l1)') trim(line), ', multipoles_on = ', ele%multipoles_on
    if (has_attribute (ele, 'TAYLOR_MAP_INCLUDES_OFFSETS') .and. (ele%taylor_map_includes_offsets .neqv. ele_dflt%taylor_map_includes_offsets)) &
                                      write (line, '(2a, l1)') trim(line), ', taylor_map_includes_offsets = ', ele%taylor_map_includes_offsets
    if (has_attribute (ele, 'CSR_CALC_ON') .and. (ele%csr_calc_on .neqv. ele_dflt%csr_calc_on)) &
                                      write (line, '(2a, l1)') trim(line), ', csr_calc_on = ', ele%csr_calc_on
    if (has_attribute (ele, 'OFFSET_MOVES_APERTURE') .and. (ele%offset_moves_aperture .neqv. ele_dflt%offset_moves_aperture)) &
                                      write (line, '(2a, l1)') trim(line), ', offset_moves_aperture = ', ele%offset_moves_aperture

    if (has_attribute (ele, 'ORIGIN_ELE') .and. ele%component_name /= '')      line = trim(line) // ', origin_ele = ' // ele%component_name 
    if (has_attribute (ele, 'CRYSTAL_TYPE') .and. ele%component_name /= '')    line = trim(line) // ', crystal_type = ' // ele%component_name 
    if (has_attribute (ele, 'MATERIAL_TYPE') .and. ele%component_name /= '')   line = trim(line) // ', material_type = ' // ele%component_name 
    if (has_attribute (ele, 'PHYSICAL_SOURCE') .and. ele%component_name /= '') line = trim(line) // ', physical_source = ' // ele%component_name 


    call write_lat_line (line, iu, .false.)  

    ! Encode taylor map. Hybrid elements do not have default terms.

    if (ele%key == taylor$ .or. (ele%key == hybrid$ .and. ele%tracking_method == taylor$)) then
      do j = 1, 6
        unit_found = .false.
        unit = 0
        unit(j:j) = 1
        do k = 1, size(ele%taylor(j)%term)
          tm = ele%taylor(j)%term(k)
          write_term = .false.
          if (all(tm%expn == unit)) then
            unit_found = .true.
            if (tm%coef /= 1) write_term = .true.
          else
            write_term = .true.
          endif

          if (write_term .or. ele%key == hybrid$) then
            if (sum(tm%expn) < 6) then
              name = ''
              do ix = 1, 6
                do n = 1, tm%expn(ix)
                  write (name, '(a, i1)') trim(name), ix
                enddo
              enddo
              write (line, '(2a, i0, 5a)') trim(line), ', {', j, ': ', trim(re_str(tm%coef)), ' | ', trim(name), '}'
            else
              write (line, '(2a, i0, 3a, 6(1x, i0), a)') trim(line), ', {', j, ': ', trim(re_str(tm%coef)), ',', tm%expn, '}'
            endif
          endif
        enddo

        if (ele%key == taylor$ .and. .not. unit_found) write (line, '(2a, i0, a, 6i2, a)') trim(line), ', {', j, ': 0,', tm%expn, '}'
      enddo

      do j1 = 1, 3;  do j2 = 1, 3
        do k = 1, size(ele%spin_taylor(j1,j2)%term)
          tm = ele%spin_taylor(j1,j2)%term(k)
          write (line, '(7a, 6i2, a)') trim(line), ', {', xyz(j1), xyz(j2), ': ', trim(re_str(tm%coef)), ',', tm%expn, '}'
        enddo
      enddo;  enddo

      if (any(ele%taylor%ref /= 0)) then
        write (line, '(16a)') trim(line), ', ref_orbit = (', &
                trim(re_str(ele%taylor(1)%ref)), ', ', trim(re_str(ele%taylor(2)%ref)), ', ', &
                trim(re_str(ele%taylor(3)%ref)), ', ', trim(re_str(ele%taylor(4)%ref)), ', ', &
                trim(re_str(ele%taylor(5)%ref)), ', ', trim(re_str(ele%taylor(6)%ref)), ')'
      endif
    endif

    ! Encode linear hybrid

    if (ele%key == hybrid$ .and. ele%tracking_method == linear$) then
      do i = 1, 6
        if (ele%vec0(i) == 0) cycle
        write (line, '(2a, i0, 2a)') trim(line), ', kick', i, ' = ', trim(re_str(ele%vec0(i)))
      enddo 
      do i = 1, 6
        do j = 1, 6
          if (ele%mat6(i,j) == 0) cycle
          write (line, '(2a, 2i0, 2a)') trim(line), ', mat', i, j, ' = ', trim(re_str(ele%mat6(i,j)))
        enddo
      enddo
    endif

    ! Encode multipoles

    if (associated(ele%a_pole)) then
      do j = 0, ubound(ele%a_pole, 1)
        if (ele%a_pole(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+a0$)) // ' = ' // re_str(ele%a_pole(j))
        if (ele%b_pole(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+b0$)) // ' = ' // re_str(ele%b_pole(j))
      enddo
    endif
    
    if (associated(ele%a_pole_elec)) then
      do j = 0, ubound(ele%a_pole_elec, 1)
        if (ele%a_pole_elec(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+a0_elec$)) // ' = ' // re_str(ele%a_pole_elec(j))
        if (ele%b_pole_elec(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+b0_elec$)) // ' = ' // re_str(ele%b_pole_elec(j))
      enddo
    endif
    
    call write_lat_line (line, iu, .true.)  

  enddo ele_loop
enddo  ! branch loop

!----------------------------------------------------------
! Lattice Layout...

! Multipass stuff...

allocate (m_region%branch(0:ubound(lat%branch, 1)))
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  allocate (m_region%branch(ib)%ele(branch%n_ele_max))
  m_region%branch(ib)%ele(:)%ix_region = 0
  m_region%branch(ib)%ele(:)%region_start_pt = .false.
  m_region%branch(ib)%ele(:)%region_stop_pt   = .false.
enddo

call multipass_all_info (lat, m_info)

if (size(m_info%lord) /= 0) then

  ! Go through and mark all 1st pass regions
  ! In theory the original lattice file could have something like:
  !   lat: line = (..., m1, m2, ..., m1, -m2, ...)
  ! where m1 and m2 are multipass lines. The first pass region (m1, m2) looks 
  ! like this is one big region but the later (m1, -m2) signals that this 
  ! is not so.
  ! We thus go through all the first pass regions and compare them to the
  ! corresponding higher pass regions. If we find two elements that are contiguous
  ! in the first pass region but not contiguous in some higher pass region, 
  ! we need to break the first pass region into two.

  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)
    e_region => m_region%branch(ib)%ele

    ix_r = 0
    in_multi_region = .false.

    do ie = 1, branch%n_ele_track
      ele => branch%ele(ie)
      e_info => m_info%branch(ib)%ele(ie)
      ix_pass = e_info%ix_pass
      if (ix_pass /= 1) then  ! Not a first pass region
        if (in_multi_region) e_region(ie-1)%region_stop_pt = .true.
        in_multi_region = .false.
        cycle
      endif
      ! If start of a new region...
      if (.not. in_multi_region) then  
        ix_r = ix_r + 1
        e_region(ie)%ix_region = ix_r
        e_region(ie)%region_start_pt = .true.
        in_multi_region = .true.
        ix_lord = e_info%ix_lord(1)
        ix_super = e_info%ix_super(1)
        ss1 => m_info%lord(ix_lord)%slave(:,ix_super)
        cycle
      endif
      ix_lord = e_info%ix_lord(1)
      ix_super = e_info%ix_super(1)
      ss2 => m_info%lord(ix_lord)%slave(:, ix_super)

      need_new_region = .false.
      if (size(ss1) /= size(ss2)) then
        need_new_region = .true.
      else
        do ix_pass = 2, size(ss1)
          if (abs(ss1(ix_pass)%ele%ix_ele - ss2(ix_pass)%ele%ix_ele) == 1) cycle
          ! not contiguous then need a new region
          need_new_region = .true.
          exit
        enddo
      endif

      if (need_new_region) then
        ix_r = ix_r + 1
        e_region(ie-1)%region_stop_pt = .true.
        e_region(ie)%region_start_pt = .true.
      endif

      ss1 => ss2
      e_region(ie)%ix_region = ix_r
    enddo

  enddo

  if (in_multi_region) e_region(branch%n_ele_track)%region_stop_pt = .true.

  ! Each 1st pass region is now a valid multipass line.
  ! Write out this info.

  write (iu, '(a)')
  write (iu, '(a)') '!-------------------------------------------------------'

  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)
    e_region => m_region%branch(ib)%ele

    ix_r = 0
    in_multi_region = .false.

    do ie = 1, branch%n_ele_track

      ele => branch%ele(ie)
      if (ie == ele%branch%n_ele_track .and. ele%name == 'END' .and. ele%key == marker$) cycle

      ix_pass = m_info%branch(ib)%ele(ie)%ix_pass
      if (ix_pass /= 1) cycle 

      if (m_region%branch(ib)%ele(ie)%region_start_pt) then
        if (ix_r > 0) then
          line = line(:len_trim(line)-1) // ')'
          call write_lat_line (line, iu, .true.)
        endif
        ix_r = ix_r + 1
        write (iu, '(a)')
        write (line, '(a, i2.2, a)') 'multi_line_', ix_r, ': line[multipass] = ('
      endif

      call write_line_element (line, iu, ele, lat)

    enddo
  enddo

  line = line(:len_trim(line)-1) // ')'
  call write_lat_line (line, iu, .true.)

end if

! Lines for all the branches.
! If we get into a multipass region then name in the main_line list is "multi_line_nn".
! But only write this once.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  write (iu, '(a)')
  name = branch%name
  if (name == '') name = 'lat_line'
  line = trim(name) // ': line = ('

  in_multi_region = .false.
  do ie = 1, branch%n_ele_track
    e_info => m_info%branch(ib)%ele(ie)
    ele => branch%ele(ie)
    if (ie == ele%branch%n_ele_track .and. ele%name == 'END' .and. ele%key == marker$) cycle

    if (.not. e_info%multipass) then
      call write_line_element (line, iu, ele, lat)
      cycle
    endif

    ix_lord = e_info%ix_lord(1)
    ix_super = e_info%ix_super(1)
    do j = 1, ubound(m_info%lord(ix_lord)%slave, 1)
      if (m_info%lord(ix_lord)%slave(j,ix_super)%ele%ix_branch /= ib) cycle
      ix1 = m_info%lord(ix_lord)%slave(j,ix_super)%ele%ix_ele
      exit
    enddo
    e_region => m_region%branch(ib)%ele
    ix_r = e_region(ix1)%ix_region

    ! If entering new multipass region
    if (.not. in_multi_region) then
      in_multi_region = .true.
      if (e_region(ix1)%region_start_pt) then
        write (line, '(2a, i2.2, a)') trim(line), ' multi_line_', ix_r, ','
        look_for = 'stop'
      else
        write (line, '(2a, i2.2, a)') trim(line), ' -multi_line_', ix_r, ','
        look_for = 'start'
      endif
    endif

    if (look_for == 'start' .and. e_region(ix1)%region_start_pt .or. &
        look_for == 'stop' .and. e_region(ix1)%region_stop_pt) then 
      in_multi_region = .false.
    endif

  enddo

  line = line(:len_trim(line)-1) // ')'
  call write_lat_line (line, iu, .true.)

  ! Branch line info

  if (ib == 0) cycle

  write (iu, '(a)')
  write (iu, '(3a)') trim(branch%name), '[geometry] = ', trim(geometry_name(branch%param%geometry))
  if (branch%param%default_tracking_species /= ref_particle$) write (iu, '(3a)') trim(branch%name), &
                        '[default_tracking_species] = ', trim(species_name(branch%param%default_tracking_species))
 
  if (branch%ix_from_branch > -1) then
    branch2 => lat%branch(branch%ix_from_branch)
    if (branch2%param%particle == branch%param%particle) cycle
  endif

  ele0 => branch%ele(0)

  write (iu, '(3a)') trim(branch%name), '[particle] = ', trim(species_name(branch%param%particle))
  write (iu, '(3a)') trim(branch%name), '[p0c]      = ', trim(re_str(ele0%value(p0c$)))

  if (is_false (ele0%value(floor_set$))) then
    if (ele0%floor%r(1) /= 0)   write (iu, '(3a)') trim(branch%name), '[x_position]     = ', trim(re_str(ele0%floor%r(1)))
    if (ele0%floor%r(2) /= 0)   write (iu, '(3a)') trim(branch%name), '[y_position]     = ', trim(re_str(ele0%floor%r(2)))
    if (ele0%floor%r(3) /= 0)   write (iu, '(3a)') trim(branch%name), '[z_position]     = ', trim(re_str(ele0%floor%r(3)))
    if (ele0%floor%theta /= 0)  write (iu, '(3a)') trim(branch%name), '[theta_position] = ', trim(re_str(ele0%floor%theta))
    if (ele0%floor%phi /= 0)    write (iu, '(3a)') trim(branch%name), '[phi_position]   = ', trim(re_str(ele0%floor%phi))
    if (ele0%floor%psi /= 0)    write (iu, '(3a)') trim(branch%name), '[psi_position]   = ', trim(re_str(ele0%floor%psi))
  endif

  if (ele0%s /= 0)              write (iu, '(3a)') trim(branch%name), '[s]              = ', trim(re_str(ele0%s))
  if (ele0%ref_time /= 0)       write (iu, '(3a)') trim(branch%name), '[ref_time]       = ', trim(re_str(ele0%ref_time))
  if (branch%param%n_part /= 0) write (iu, '(3a)') trim(branch%name), '[n_part]         = ', trim(re_str(lat%param%n_part))

  if (branch%param%geometry == open$) then
    write (iu, '(3a)')
    if (ele0%a%beta /= 0)     write (iu, '(3a)') trim(branch%name), '[beta_a]   = ', trim(re_str(ele0%a%beta))
    if (ele0%a%alpha /= 0)    write (iu, '(3a)') trim(branch%name), '[alpha_a]  = ', trim(re_str(ele0%a%alpha))
    if (ele0%a%phi /= 0)      write (iu, '(3a)') trim(branch%name), '[phi_a]    = ', trim(re_str(ele0%a%phi))
    if (ele0%x%eta /= 0)      write (iu, '(3a)') trim(branch%name), '[eta_x]    = ', trim(re_str(ele0%x%eta))
    if (ele0%x%etap /= 0)     write (iu, '(3a)') trim(branch%name), '[etap_x]   = ', trim(re_str(ele0%x%etap))
    if (ele0%b%beta /= 0)     write (iu, '(3a)') trim(branch%name), '[beta_b]   = ', trim(re_str(ele0%b%beta))
    if (ele0%b%alpha /= 0)    write (iu, '(3a)') trim(branch%name), '[alpha_b]  = ', trim(re_str(ele0%b%alpha))
    if (ele0%b%phi /= 0)      write (iu, '(3a)') trim(branch%name), '[phi_b]    = ', trim(re_str(ele0%b%phi))
    if (ele0%y%eta /= 0)      write (iu, '(3a)') trim(branch%name), '[eta_y]    = ', trim(re_str(ele0%y%eta))
    if (ele0%y%etap /= 0)     write (iu, '(3a)') trim(branch%name), '[etap_y]   = ', trim(re_str(ele0%y%etap))
    if (ele0%c_mat(1,1) /= 0) write (iu, '(3a)') trim(branch%name), '[cmat_11]  = ', trim(re_str(ele0%c_mat(1,1)))
    if (ele0%c_mat(1,2) /= 0) write (iu, '(3a)') trim(branch%name), '[cmat_12]  = ', trim(re_str(ele0%c_mat(1,2)))
    if (ele0%c_mat(2,1) /= 0) write (iu, '(3a)') trim(branch%name), '[cmat_21]  = ', trim(re_str(ele0%c_mat(2,1)))
    if (ele0%c_mat(2,2) /= 0) write (iu, '(3a)') trim(branch%name), '[cmat_22]  = ', trim(re_str(ele0%c_mat(2,2)))
  endif

enddo

! Use line

line = 'use'
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%ix_from_branch > -1) cycle
  name = branch%name
  if (name == '') name = 'lat_line'
  line = trim(line) // ', ' // name
enddo

write (iu, '(a)')
write (iu, '(a)') trim(line)

! If there are multipass lines then expand the lattice and write out
! the post-expand info as needed.

expand_branch_out = .false.
do ie = 1, lat%n_ele_max
  ele => lat%ele(ie)
  if (ele%slave_status == super_slave$) cycle

  if (ele%key == lcavity$ .or. ele%key == rfcavity$) then
    if (ele%value(phi0_multipass$) == 0) cycle
    if (.not. expand_branch_out) call write_expand_lat_header
    write (iu, '(3a)') trim(ele%name), '[phi0_multipass] = ', trim(re_str(ele%value(phi0_multipass$)))
  endif

enddo

! cleanup

close(iu)
deallocate (names, an_indexx)
deallocate (m_region%branch)
call deallocate_multipass_all_info_struct (m_info)

if (present(err)) err = .false.

!--------------------------------------------------------------------------------
contains

subroutine write_expand_lat_header ()

write (iu, '(a)')
write (iu, '(a)') '!-------------------------------------------------------'
write (iu, '(a)')
write (iu, '(a)') 'expand_lattice'
write (iu, '(a)')
expand_branch_out = .true.

end subroutine write_expand_lat_header

!--------------------------------------------------------------------------------
! contains

subroutine form_this_field_map_name(string, field_type, ele, ix_map, output_form)

type (ele_struct) ele
integer ix_map
integer, optional :: output_form
character(*) string, field_type
character(40) name

!

call str_downcase(name, ele%name)
if (ix_map == 1) then
  write (string, '(2a)') trim(name), field_type
else
  write (string, '(2a, i0, a)') trim(name), '_', ix_map, field_type
endif

if (integer_option(binary$, output_form) == binary$) string = trim(string) // '.bin'

end subroutine form_this_field_map_name

!--------------------------------------------------------------------------------
! contains

subroutine write_map_coef (tag, array)

real(rp) :: array(:)
character(*) tag
integer j, k0, k, n, n5

!

write (iu2, '(a)') ','
write (iu2, '(2x, 2a)', advance = 'NO') tag, ' = ('

n = size(array)
n5 = (n - 1) / 5

do j = 1, n5
  k0 = 5 * (j - 1)
  write (iu2, '(5(es14.6, a))') (array(k0+k), ',', k = 1, 5) 
  write (iu2, '(15x)', advance = 'NO')
enddo

if (n5 == 0) write (iu2, '(15x)', advance = 'NO')
write (iu2, '(5(es14.6, a))', advance = 'NO') (array(k), ',', k = 5*n5+1, n-1), array(n), ')' 

end subroutine write_map_coef

end subroutine write_bmad_lattice_file

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

subroutine write_line_element (line, iu, ele, lat)

implicit none

type (lat_struct), target :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, m_lord, slave

character(*) line
character(40) lord_name

integer iu
integer j, ix

!

if (ele%iyy == 3) return

if (ele%slave_status == super_slave$) then

  if (ele%iyy == 1) then
    lord => pointer_to_lord(ele, 1)
    if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(lord, 1)

    if (lord%orientation == 1) then
      write (line, '(4a)') trim(line), ' ', trim(lord%name), ','
    else
      write (line, '(4a)') trim(line), ' --', trim(lord%name), ','
    endif

    return
  endif  

  if (ele%iyy == 2) return

  write (line, '(2a, i0, a)') trim(line), ' slave_drift_', ele%ixx, ','

elseif (ele%slave_status == multipass_slave$) then
  lord => pointer_to_lord(ele, 1)
  write (line, '(4a)') trim(line), ' ', trim(lord%name), ','

else
  if (ele%orientation == 1) then
    write (line, '(4a)') trim(line), ' ', trim(ele%name), ','
  else
    write (line, '(4a)') trim(line), ' --', trim(ele%name), ','
  endif
endif

if (len_trim(line) > 80) call write_lat_line(line, iu, .false.)

end subroutine write_line_element

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function re_str(rel) result (str_out)

implicit none

real(rp) rel
integer pl
character(24) str_out
character(16) fmt

!

if (rel == 0) then
  str_out = '0'
  return
endif

pl = floor(log10(abs(rel)))

if (pl > 5) then
  fmt = '(2a, i0)'
  write (str_out, fmt) trim(rchomp(rel/10.0_rp**pl, 0)), 'E', pl

elseif (pl > -3) then
  str_out = rchomp(rel, pl)

else
  fmt = '(2a, i0)'
  write (str_out, fmt) trim(rchomp(rel*10.0_rp**(-pl), 0)), 'E', pl

endif

end function re_str

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function array_re_str(arr, parens_in) result (str_out)

real(rp) arr(:)
integer i
character(100) str_out
character(*), optional :: parens_in
character(2) parens

!

parens = '()'
if (present(parens_in)) parens = parens_in

str_out = parens(1:1) // re_str(arr(1))
do i = 2, size(arr)
  str_out = trim(str_out) // ', ' // re_str(arr(i))
enddo
str_out = trim(str_out) // parens(2:2)

end function array_re_str

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function cmplx_re_str(cmp) result (str_out)

complex(rp) cmp
character(40) str_out

!

if (imag(cmp) == 0) then
  str_out = trim(re_str(real(cmp)))
else
  str_out = '(' // trim(re_str(real(cmp))) // ', ' // trim(re_str(imag(cmp))) // ')'
endif

end function cmplx_re_str

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function rchomp (rel, plc) result (out)

implicit none

real(rp) rel
character(24) out
character(8) :: fmt = '(f24.xx)'
integer it, plc, ix

!

write (fmt(6:7), '(i2.2)') 13-plc  ! 14 digits of accuracy
write (out, fmt) rel
do it = len(out), 1, -1
  if (out(it:it) == ' ') cycle
  if (out(it:it) == '0') then
    out(it:it) = ' '
    cycle
  endif
  if (out(it:it) == '.') out(it:it) = ' '
  call string_trim(out, out, ix)
  return
enddo

end function rchomp

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!+
! Subroutine write_lat_line (line, iu, end_is_neigh, continue_char)
!
! Routine to write strings to a lattice file.
! This routine will break the string up into multiple lines
! if the string is too long and add a continuation character if needed.
!
! If the "line" arg does not represent a full "sentence" (end_is_neigh = False), 
! then only part of the line may be written and the part not writen will be returned.
!
! Module needed:
!   use write_lat_file_mod
!
! Input:
!   line          -- character(*): String of text.
!   iu            -- Integer: Unit number to write to.
!   end_is_neigh  -- Logical: If true then write out everything.
!                      Otherwise wait for a full line of max_char characters or so.
!   continue_char -- character(1), optional. Default is not to use any continue character 
!                      and put line breaks after commas.
!
! Output:
!   line          -- Character(*): part of the string not writen. 
!                       If end_is_neigh = T then line will be blank.
!-

subroutine write_lat_line (line, iu, end_is_neigh, continue_char)

implicit none

character(*) line
integer i, iu
logical end_is_neigh
logical, save :: init = .true.
integer, save :: max_char = 90
character(1), optional :: continue_char

!

outer_loop: do 

  if (len_trim(line) <= max_char) then
    if (end_is_neigh) then
      call write_this (line)
      init = .true.
    endif
    return
  endif
      
  do i = max_char, 1, -1
    if (line(i:i) == ',') then
      call write_this (line(:i), continue_char)
      line = line(i+1:)
      cycle outer_loop
    endif
  enddo

  do i = max_char+1, len_trim(line)
    if (line(i:i) == ',') then
      call write_this (line(:i), continue_char)
      line = line(i+1:)
      cycle outer_loop
    endif
  enddo

  if (end_is_neigh) then
    call write_this (line)
    init = .true.
    return
  endif

enddo outer_loop

contains

subroutine write_this (line2, c_char)

character(*) line2
character(1), optional :: c_char
character(20) fmt

!

if (init) then
  fmt = '(a, 1x, a)'
  init = .false.
else
  fmt = '(2x, a, 1x, a)'
endif

if (present(c_char)) then
  write (iu, fmt) trim(line2), c_char
else
  write (iu, fmt) trim(line2)
endif

end subroutine write_this

end subroutine write_lat_line

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+ 
! Subroutine write_lattice_in_foreign_format (out_type, out_file_name, lat, ref_orbit, &
!        use_matrix_model, include_apertures, dr12_drift_max, ix_start, ix_end, ix_branch, converted_lat, err)
!
! Subroutine to write a MAD-8, MAD-X, OPAL, SAD, or XSIF lattice file using the 
! information in a lat_struct. Optionally, only part of the lattice can be generated.
!
! Also see: write_bmad_lattice_file
!
! NOTE: When translating to XSIF or MAD: sad_mult and patch element are translated
!  to a XSIF/MAD matrix element (which is a 2nd order map). In this case, the ref_orbit orbit is
!  used as the reference orbit for construction of the 2nd order map.
!
! If a sad_mult or patch element is translated to a matrix element, and the referece orbit
! is non-zero, the calculation must use 2nd order maps thourghout in order to avoid "feed down".
! If the PTC map order is different from 2, PTC will be temperarily switched to 2. 
!
! The MAD drift model is approximate and this can be a problem if the reference orbit is large.
! For a drift, the value of transfer matrix element R12 is equal to L/(1+pz) for small
! deviations of the ref_orbit from zero. dr12_drift_max sets the maximum deviation of R12 beyound 
! which an extra matrix element is inserted to make the MAD model better agree with Bmad.
!
! Note: sol_quad elements are replaced by a drift-matrix-drift or solenoid-quad model.
! Note: wiggler elements are replaced by a drift-matrix-drift or drift-bend model.
!
! Modules needed:
!   use write_lat_file_mod
!
! Input:
!   out_type          -- character(*): Either 'XSIF', 'MAD-8', 'MAD-X', 'SAD', or 'OPAL-T'.
!   out_file_name     -- character(*): Name of the mad output lattice file.
!   lat               -- lat_struct: Holds the lattice information.
!   ref_orbit(0:)     -- coord_struct, allocatable, optional: Referece orbit for sad_mult and patch elements.
!                          This argument must be present if the lattice has sad_mult or patch elements.
!   use_matrix_model  -- logical, optional: Use a drift-matrix_drift model for wigglers/undulators?
!                           [A MAD "matrix" is a 2nd order Taylor map.] This switch is ignored for SAD conversion.
!                           Default is False -> Use a bend-drift-bend model. 
!                           Note: sol_quad elements always use a drift-matrix-drift model.
!   include_apertures -- logical, optional: Include aperture info? Ignored for SAD convertion. Default is True.
!   dr12_drift_max    -- real(rp), optional: Max deviation for drifts allowed before a correction matrix element
!                           is added. Default value is 1d-5.
!   ix_start          -- integer, optional: Starting index of lat%ele(i)
!                           used for output.
!   ix_end            -- integer, optional: Ending index of lat%ele(i)
!                           used for output.
!   ix_branch         -- Integer, optional: Index of lattice branch to use. Default = 0.
!
! Output:
!   converted_lat     -- lat_struct, optional: Equivalent Bmad lattice with wiggler and 
!                           sol_quad elements replaced by their respective models.
!                           This is only valid for MAD-8, MAD-X, and XSIF conversions.
!   err               -- logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_lattice_in_foreign_format (out_type, out_file_name, lat, ref_orbit, &
      use_matrix_model, include_apertures, dr12_drift_max, ix_start, ix_end, ix_branch, converted_lat, err)

use mad_mod

implicit none

type (lat_struct), target :: lat, lat_model, lat_out
type (lat_struct), optional, target :: converted_lat
type (ele_struct), pointer :: ele, ele1, ele2, lord, sol_ele, first_sol_edge
type (ele_struct), save :: drift_ele, ab_ele, taylor_ele, col_ele, kicker_ele, null_ele, bend_ele, quad_ele
type (coord_struct) orb_start, orb_end, orb_center
type (coord_struct), allocatable, optional :: ref_orbit(:)
type (coord_struct), allocatable :: orbit_out(:)
type (taylor_term_struct) :: term
type (branch_struct), pointer :: branch, branch_out
type (mad_energy_struct) energy
type (mad_map_struct) mad_map
type (ptc_parameter_struct) ptc_param
type (taylor_struct) taylor_a(6), taylor_b(6)

real(rp), optional :: dr12_drift_max
real(rp) field, hk, vk, tilt, limit(2), length, a, b, f, e2
real(rp), pointer :: val(:)
real(rp) knl(0:n_pole_maxx), tilts(0:n_pole_maxx), a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer, optional :: ix_start, ix_end, ix_branch
integer, allocatable :: n_repeat(:), an_indexx(:)
integer i, j, ib, j2, k, n, ix, i_unique, i_line, iout, iu, n_names, j_count, f_count, ix_ele
integer ie, ie1, ie2, ios, t_count, s_count, a_count, ix_lord, ix_match, iv, ifa, ix_pole_max
integer ix1, ix2, n_lord, aperture_at, n_name_change_warn, n_elsep_warn, n_taylor_order_saved
integer :: ix_line_min, ix_line_max, n_warn_max = 10

character(*), parameter :: r_name = "write_lattice_in_foreign_format"
character(*) out_type, out_file_name
character(300) line, knl_str, ksl_str
character(40) orig_name, str
character(40), allocatable :: names(:)
character(4000) line_out   ! Can be this large for taylor maps.
character(2) continue_char, eol_char, comment_char, separator_char

logical, optional :: use_matrix_model, include_apertures, err
logical init_needed, mad_out
logical parsing, warn_printed, converted

! SAD translation

if (out_type == 'SAD') then
  call write_lat_in_sad_format (out_file_name, lat, include_apertures, ix_start, ix_end, ix_branch, converted_lat, err)
  return
endif

! open file

if (present(err)) err = .true.
n_taylor_order_saved = ptc_com%taylor_order_ptc

iu = lunget()
call fullfilename (out_file_name, line)
open (iu, file = line, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(out_file_name))
  return
endif

! Use ptc exact_model = True since this is needed to get the drift nonlinear terms

call get_ptc_params(ptc_param)
call set_ptc (exact_modeling = .true.)

! Init

ix = integer_option(0, ix_branch)
if (ix < 0 .or. ix > ubound(lat%branch, 1)) then
  call out_io (s_error$, r_name, 'BRANCH INDEX OUT OF RANGE: /i0/ ', i_array = [ix])
  return
endif

branch => lat%branch(ix)

if (out_type == 'MAD-X' .or. out_type == 'OPAL-T') then
  comment_char = '//'
  continue_char = ''
  eol_char = ';'
  separator_char = ','
  ix_line_max = 100

elseif (out_type == 'MAD-8' .or. out_type == 'XSIF') then
  comment_char = '!'
  continue_char = ' &'
  eol_char = ''
  separator_char = ','
  ix_line_max = 80

else
  call out_io (s_error$, r_name, 'BAD OUT_TYPE: ' // out_type)
  return
endif

mad_out = .false.
if (out_type == 'MAD-X' .or. out_type == 'MAD-8') mad_out = .true.

ix_line_min = ix_line_max - 20

call init_ele (col_ele)
call init_ele (drift_ele, drift$)
call init_ele (taylor_ele, taylor$)
call init_ele (ab_ele, ab_multipole$)
call init_ele (kicker_ele, kicker$) 
call init_ele (quad_ele, quadrupole$)
call init_ele (bend_ele, sbend$)
call multipole_init (ab_ele, magnetic$)
null_ele%key = null_ele$

ie1 = integer_option(1, ix_start)
ie2 = integer_option(branch%n_ele_track, ix_end)

allocate (names(branch%n_ele_max), an_indexx(branch%n_ele_max)) ! list of element names

call out_io (s_info$, r_name, &
      'Note: Bmad lattice elements have attributes that cannot be translated. ', &
      '      For example, higher order terms in a Taylor element.', &
      '      Please use caution when using a translated lattice.')

!-----------------------------------------------------------------------------
! Translation is a two step process:
!   1) Create a new lattice called lat_out making substitutions for sol_quad and wiggler elements, etc..
!   2) Use lat_out to create the lattice file.

lat_out = lat
call allocate_lat_ele_array(lat_out, 2*branch%n_ele_max, branch%ix_branch)
branch_out => lat_out%branch(branch%ix_branch)

if (present(ref_orbit)) then
  call reallocate_coord(orbit_out, size(ref_orbit))
  orbit_out = ref_orbit
else
  call reallocate_coord(orbit_out, branch%n_ele_max)
endif

f_count = 0    ! fringe around bends and quads. Also drift nonlinearities.
j_count = 0    ! drift around solenoid or sol_quad index. Also z shift count.
t_count = 0    ! taylor element count.
a_count = 0    ! Aperture count
i_unique = 1000

! Loop over all input elements

nullify(first_sol_edge)
n_name_change_warn = 0
n_elsep_warn = 0
ix_ele = ie1 - 1

do

  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => branch_out%ele(ix_ele)
  val => ele%value

  ! If the name has more than 16 characters then replace the name by something shorter and unique.

  orig_name = ele%name

  if (len_trim(ele%name) > 16) then
    i_unique = i_unique + 1
    write (ele%name, '(a, i0)') ele%name(1:11), i_unique
  endif

  ! Replace element name containing "/" or "#" with "_"

  do
    j = index (ele%name, '\')         ! '
    j = index (ele%name, '#')   
    if (j == 0) exit
    ele%name(j:j) = '_'
  enddo

  if (ele%name /= orig_name .and. n_name_change_warn <= n_warn_max) then
    call out_io (s_info$, r_name, 'Element name changed from: ' // trim(orig_name) // ' to: ' // ele%name)
    if (n_name_change_warn == n_warn_max) call out_io (s_info$, r_name, &
                           'Enough name change warnings. Will stop issuing them now.')
    n_name_change_warn = n_name_change_warn + 1
  endif

  ! If there is an aperture with an element that is not an ecoll or rcoll then need to make a separate
  ! element with the aperture info. Exception: MAD-X can handle apertures on non collimator elements.

  if ((val(x1_limit$) /= 0 .or. val(x2_limit$) /= 0 .or. val(y1_limit$) /= 0 .or. val(y2_limit$) /= 0) .and. &
      ele%key /= ecollimator$ .and. ele%key /= rcollimator$ .and. logic_option(.true., include_apertures) .and. &
      (ele%key == drift$ .or. out_type /= 'MAD-X')) then

    if (val(x1_limit$) /= val(x2_limit$)) then
      call out_io (s_warn$, r_name, 'Asymmetric x_limits cannot be converted for: ' // ele%name, &
                                    'Will use largest limit here.')
      val(x1_limit$) = max(val(x1_limit$), val(x2_limit$))
    endif

    if (val(y1_limit$) /= val(y2_limit$)) then
      call out_io (s_warn$, r_name, 'Asymmetric y_limits cannot be converted for: ' // ele%name, &
                                    'Will use largest limit here.')
      val(y1_limit$) = max(val(y1_limit$), val(y2_limit$))
    endif

    ! create ecoll and rcoll elements.

    if (ele%aperture_type == rectangular$) then
      col_ele%key = rcollimator$
    else
      col_ele%key = ecollimator$
    endif
    a_count = a_count + 1
    write (col_ele%name, '(a, i0)')  'COLLIMATOR_N', a_count
    col_ele%value = val
    col_ele%value(l$) = 0
    val(x1_limit$) = 0; val(x2_limit$) = 0; val(y1_limit$) = 0; val(y2_limit$) = 0; 
    aperture_at = ele%aperture_at  ! Save since ele pointer will be invalid after the insert
    if (aperture_at == both_ends$ .or. aperture_at == downstream_end$ .or. aperture_at == continuous$) then
      call insert_element (lat_out, col_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
      ie2 = ie2 + 1
    endif
    if (aperture_at == both_ends$ .or. aperture_at == upstream_end$ .or. aperture_at == continuous$) then
      call insert_element (lat_out, col_ele, ix_ele, branch_out%ix_branch, orbit_out)
      ie2 = ie2 + 1
    endif
    ix_ele = ix_ele - 1 ! Want to process the element again on the next loop.

    cycle ! cycle since ele pointer is invalid

  endif

  ! If the bend has a roll then put kicker elements just before and just after

  if (ele%key == sbend$ .and. val(roll$) /= 0) then
    j_count = j_count + 1
    write (kicker_ele%name,   '(a, i0)') 'ROLL_Z', j_count
    kicker_ele%value(hkick$) =  val(angle$) * (1 - cos(val(roll$))) / 2
    kicker_ele%value(vkick$) = -val(angle$) * sin(val(roll$)) / 2
    val(roll$) = 0   ! So on next iteration will not create extra kickers.
    call insert_element (lat_out, kicker_ele, ix_ele, branch_out%ix_branch, orbit_out)
    call insert_element (lat_out, kicker_ele, ix_ele+2, branch_out%ix_branch, orbit_out)
    ie2 = ie2 + 2
    cycle
  endif

  ! If there is a multipole component then put multipole elements at half strength 
  ! just before and just after the element.

  if (ele%key /= multipole$ .and. ele%key /= ab_multipole$ .and. ele%key /= null_ele$ .and. ele%key /= sad_mult$) then
    call multipole_ele_to_ab (ele, .true., ix_pole_max, ab_ele%a_pole, ab_ele%b_pole)
    if (ix_pole_max > -1) then
      ab_ele%a_pole = ab_ele%a_pole / 2
      ab_ele%b_pole = ab_ele%b_pole / 2
      if (associated(ele%a_pole)) deallocate (ele%a_pole, ele%b_pole)
      j_count = j_count + 1
      write (ab_ele%name, '(a1, a, i0)') key_name(ele%key), 'MULTIPOLE_', j_count
      call insert_element (lat_out, ab_ele, ix_ele, branch_out%ix_branch, orbit_out)
      call insert_element (lat_out, ab_ele, ix_ele+2, branch_out%ix_branch, orbit_out)
      ie2 = ie2 + 2
      cycle
    endif
  endif

  ! If there are nonzero kick values and this is not a kick type element then put
  ! kicker elements at half strength just before and just after the element.
  ! Also add a matrix element to get the change in z correct.
  ! A sad_mult gets translated to a matrix element which has kick components so no extra kickers needed here.

  if (has_hkick_attributes(ele%key)) then
    if (ele%key /= kicker$ .and. ele%key /= hkicker$ .and. ele%key /= vkicker$ .and. ele%key /= sad_mult$) then
      if (val(hkick$) /= 0 .or. val(vkick$) /= 0) then
        j_count = j_count + 1
        write (kicker_ele%name, '(a1, a, i0)') key_name(ele%key), '_KICKER_', j_count
        kicker_ele%value(hkick$) = val(hkick$) / 2
        kicker_ele%value(vkick$) = val(vkick$) / 2
        val(hkick$) = 0; val(vkick$) = 0
        if (ele%key == sbend$) then
          f = val(g_err$) * val(l$) / 2
          kicker_ele%value(hkick$) = kicker_ele%value(hkick$) - cos(ele%value(ref_tilt_tot$)) * f
          kicker_ele%value(vkick$) = kicker_ele%value(vkick$) - sin(ele%value(ref_tilt_tot$)) * f
          val(g_err$) = 0
        endif
        write (taylor_ele%name, '(a, i0)') 'Z_SHIFTER', j_count 
        call taylor_make_unit(taylor_ele%taylor)
        orb_start = orbit_out(ix_ele-1)
        orb_start%vec(2) = orb_start%vec(2) - kicker_ele%value(hkick$)
        orb_start%vec(4) = orb_start%vec(4) - kicker_ele%value(vkick$)
        call track1 (orb_start, ele, branch_out%param, orb_end) 
        f = (ele%map_ref_orb_out%vec(5) - ele%map_ref_orb_in%vec(5)) - (orb_end%vec(5) - orb_start%vec(5))
        call add_taylor_term (taylor_ele%taylor(5), f, [0, 0, 0, 0, 0, 0])
        call insert_element (lat_out, kicker_ele, ix_ele, branch_out%ix_branch, orbit_out)
        call insert_element (lat_out, kicker_ele, ix_ele+2, branch_out%ix_branch, orbit_out)
        call insert_element (lat_out, taylor_ele, ix_ele+3, branch_out%ix_branch, orbit_out)
        ie2 = ie2 + 3
        cycle
      endif
    endif
  endif

  ! A quadrupole with fringe = full or soft_edge_only has its fringe kicks modeled as a 2nd order map.

  iv = nint(ele%value(fringe_type$))
  if (mad_out .and. ele%key == quadrupole$ .and. (iv == full$ .or. iv == soft_edge_only$)) then
    quad_ele = ele
    ele%value(fringe_type$) = none$

    if (ptc_com%taylor_order_ptc /= 2) call set_ptc (taylor_order = 2) 

    f_count = f_count + 1
    ie = ix_ele

    ifa = nint(ele%value(fringe_at$))
    if (ifa == entrance_end$ .or. ifa == both_ends$) then
      quad_ele%value(fringe_at$) = entrance_end$
      quad_ele%value(l$) = 1d-30
      call ele_to_taylor (quad_ele, branch_out%param, taylor_ele%taylor, orbit_out(ie-1))
      write (taylor_ele%name, '(a, i0)') 'Q_FRINGE_IN', f_count
      call insert_element (lat_out, taylor_ele, ie, branch_out%ix_branch, orbit_out)
      ie = ie + 1
      ie2 = ie2 + 1
    endif

    if (ifa == exit_end$ .or. ifa == both_ends$) then
      quad_ele%value(fringe_at$) = exit_end$
      quad_ele%value(l$) = 1d-30
      call ele_to_taylor (quad_ele, branch_out%param, taylor_ele%taylor, orbit_out(ie))
      write (taylor_ele%name, '(a, i0)') 'Q_FRINGE_OUT', f_count
      call insert_element (lat_out, taylor_ele, ie+1, branch_out%ix_branch, orbit_out)
      ie2 = ie2 + 1
    endif

    cycle
  endif

  ! A bend with fringe = sad_full or has non-zero g_err has its fringe kicks modeled as a 1st order map.

  iv = nint(ele%value(fringe_type$))
  if (mad_out .and. ele%key == sbend$ .and. (iv == sad_full$ .or. ele%value(g_err$) /= 0)) then

    if (ptc_com%taylor_order_ptc /= 1) call set_ptc (taylor_order = 1)

    f_count = f_count + 1
    ie = ix_ele

    bend_ele = ele
    bend_ele%value(l$) = ele%value(l$)/2
    bend_ele%value(angle$) = ele%value(angle$)/2
    bend_ele%value(e2$) = 0
    call set_fringe_on_off (bend_ele%value(fringe_at$), exit_end$, off$)
    call track1 (orbit_out(ie-1), bend_ele, branch_out%param, orb_center)

    if (at_this_ele_end(entrance_end$, nint(ele%value(fringe_at$))) .or. ele%value(g_err$) /= 0) then
      call ele_to_taylor (bend_ele, branch_out%param, taylor_a, orbit_out(ie-1))

      bend_ele%value(fringe_type$) = basic_bend$
      bend_ele%value(g_err$) = 0
      orb_start = orb_center
      orb_start%direction = -1
      orb_start%species = antiparticle(orb_center%species)
      call track1 (orb_start, bend_ele, branch_out%param, orb_start)  ! bactrack to entrance end
      call ele_to_taylor (bend_ele, branch_out%param, taylor_b, orb_start)

      call taylor_inverse (taylor_b, taylor_b)
      call concat_taylor (taylor_a, taylor_b, taylor_ele%taylor)
      write (taylor_ele%name, '(a, i0)') 'B_FRINGE_IN', f_count
      call insert_element (lat_out, taylor_ele, ie, branch_out%ix_branch, orbit_out)
      ele => branch_out%ele(ix_ele+1)
      call kill_taylor (taylor_a)
      call kill_taylor (taylor_b)
      ie = ie + 1
      ie2 = ie2 + 1
    endif

    if (at_this_ele_end(exit_end$, nint(ele%value(fringe_at$))) .or. ele%value(g_err$) /= 0) then
      bend_ele = ele
      bend_ele%value(l$) = ele%value(l$)/2
      bend_ele%value(angle$) = ele%value(angle$)/2
      bend_ele%value(e1$) = 0
      call set_fringe_on_off (bend_ele%value(fringe_at$), entrance_end$, off$)

      call ele_to_taylor (bend_ele, branch_out%param, taylor_a, orb_center)

      bend_ele%value(fringe_type$) = basic_bend$
      bend_ele%value(g_err$) = 0
      call ele_to_taylor (bend_ele, branch_out%param, taylor_b, orb_center)
      call taylor_inverse (taylor_b, taylor_b)

      call concat_taylor (taylor_b, taylor_a, taylor_ele%taylor)
      write (taylor_ele%name, '(a, i0)') 'B_FRINGE_OUT', f_count
      call insert_element (lat_out, taylor_ele, ie+1, branch_out%ix_branch, orbit_out)
      call kill_taylor (taylor_a)
      call kill_taylor (taylor_b)
      ie2 = ie2 + 1
    endif

    ele%value(fringe_type$) = basic_bend$
    ele%value(g_err$) = 0
    cycle
  endif

  ! A drift where the ref orbit is too large needs an added 1st order matrix element 

  f = ele%value(l$) / (1 + orbit_out(ele%ix_ele)%vec(6))
  if (mad_out .and. ele%key == drift$ .and. abs(ele%mat6(1,2) - f) > real_option(1d-5, dr12_drift_max)) then

    if (ptc_com%taylor_order_ptc /= 1) call set_ptc (taylor_order = 1) 

    drift_ele = ele
    drift_ele%value(l$) = -ele%value(l$)
    call make_mat6_mad (drift_ele, branch_out%param, orbit_out(ix_ele), orb_end)
    call mat6_to_taylor (drift_ele%vec0, drift_ele%mat6, taylor_a)

    drift_ele%value(l$) = ele%value(l$)
    call ele_to_taylor (drift_ele, branch_out%param, taylor_b, orbit_out(ix_ele-1))
    call concat_taylor (taylor_a, taylor_b, taylor_ele%taylor)
    call kill_taylor (taylor_a)
    call kill_taylor (taylor_b)

    f_count = f_count + 1
    write (taylor_ele%name, '(a, i0)') 'D_NONLIN', f_count
    call insert_element (lat_out, taylor_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
    ie2 = ie2 + 1
    ix_ele = ix_ele + 1
    cycle
  endif

  ! Convert sol_quad_and wiggler elements to an "equivalent" set of elements.
  ! NOTE: FOR NOW, SOL_QUAD USES DRIFT-MATRIX-DRIFT MODEL!

  if (ele%key == wiggler$ .or. ele%key == undulator$ .or. ele%key == sol_quad$) then
    if (logic_option(.false., use_matrix_model) .or. ele%key == sol_quad$) then
      call out_io (s_warn$, r_name, 'Converting element to drift-matrix-drift model: ' // ele%name)
      drift_ele%value(l$) = -val(l$) / 2
      call make_mat6 (drift_ele, branch_out%param)
      taylor_ele%mat6 = matmul(matmul(drift_ele%mat6, ele%mat6), drift_ele%mat6)
      call mat6_to_taylor (taylor_ele%vec0, taylor_ele%mat6, taylor_ele%taylor)

      ! Add drifts before and after wigglers and sol_quads so total length is invariant
      j_count = j_count + 1
      t_count = t_count + 1
      write (drift_ele%name, '(a, i0)') 'DRIFT_Z', j_count
      write (taylor_ele%name, '(a, i0)') 'SOL_QUAD', j_count
      drift_ele%value(l$) = val(l$) / 2
      ele%key = -1 ! Mark for deletion
      call remove_eles_from_lat (lat_out)
      call insert_element (lat_out, drift_ele, ix_ele, branch_out%ix_branch, orbit_out)
      call insert_element (lat_out, taylor_ele, ix_ele+1, branch_out%ix_branch, orbit_out)
      call insert_element (lat_out, drift_ele, ix_ele+2, branch_out%ix_branch, orbit_out)
      ie2 = ie2 + 2
      cycle

    ! Non matrix model...
    ! If the wiggler has been sliced due to superposition, throw 
    ! out the markers that caused the slicing.

    else
      if (ele%key == wiggler$ .or. ele%key == undulator$) then
        call out_io (s_warn$, r_name, 'Converting element to drift-bend-drift model: ' // ele%name)
        if (ele%slave_status == super_slave$) then
          ! Create the wiggler model using the super_lord
          lord => pointer_to_lord(ele, 1)
          call create_wiggler_model (lord, lat_model)
          ! Remove all the slave elements and markers in between.
          call out_io (s_warn$, r_name, &
              'Note: Not translating to MAD/XSIF the markers within wiggler: ' // lord%name)
          lord%key = -1 ! mark for deletion
          call find_element_ends (lord, ele1, ele2)
          ix1 = ele1%ix_ele; ix2 = ele2%ix_ele
          ! If the wiggler wraps around the origin we are in trouble.
          if (ix2 < ix1) then 
            call out_io (s_fatal$, r_name, 'Wiggler wraps around origin. Cannot translate this!')
            if (global_com%exit_on_error) call err_exit
          endif
          do i = ix1+1, ix2
            branch_out%ele(i)%key = -1  ! mark for deletion
          enddo
          ie2 = ie2 - (ix2 - ix1 - 1)
        else
          call create_wiggler_model (ele, lat_model)
        endif
      else
        call create_sol_quad_model (ele, lat_model)  ! NOT YET IMPLEMENTED!
      endif
      ele%key = -1 ! Mark for deletion
      call remove_eles_from_lat (lat_out)
      do j = 1, lat_model%n_ele_track
        call insert_element (lat_out, lat_model%ele(j), ix_ele+j-1, branch_out%ix_branch, orbit_out)
      enddo
      ie2 = ie2 + lat_model%n_ele_track - 1
      cycle
    endif
  endif

enddo

! For a patch that is *not* associated with the edge of a solenoid: A z_offset must be split into a drift + patch

ix_ele = ie1 - 1

do
  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => branch_out%ele(ix_ele)
  val => ele%value

  if (ele%key == patch$ .and. ele%value(z_offset$) /= 0) then
    drift_ele%name = 'DRIFT_' // ele%name
    drift_ele%value(l$) = val(z_offset$)
    call insert_element (lat_out, drift_ele, ix_ele, branch_out%ix_branch, orbit_out)
    ix_ele = ix_ele + 1
    ele => branch_out%ele(ix_ele)
    val => ele%value
    val(z_offset$) = 0
  endif
enddo

!-------------------------------------------------------------------------------------------------
! Now write info to the output file...
! lat lattice name

write (iu, '(3a)') comment_char, ' File generated by: write_lattice_in_foreign_format', trim(eol_char)
write (iu, '(4a)') comment_char, ' Bmad Lattice File: ', trim(lat%input_file_name), trim(eol_char)
if (lat%lattice /= '') write (iu, '(4a)') comment_char, ' Bmad Lattice: ', trim(lat%lattice), trim(eol_char)
write (iu, '(a)')

! beam definition

select case (out_type)
case ('MAD-8', 'MAD-X', 'XSIF')
  ele => branch_out%ele(ie1-1)

  write (line_out, '(7a)') 'beam_def: Beam, Particle = ', trim(species_name(branch_out%param%particle)),  &
        ', Energy = ', trim(re_str(1d-9*ele%value(E_TOT$))), ', Npart = ', trim(re_str(branch_out%param%n_part)), trim(eol_char)
  call write_line (line_out)
  write (iu, '(a)')

end select

! write element parameters

n_names = 0                          ! number of names stored in the list

do ix_ele = ie1, ie2

  ele => branch_out%ele(ix_ele)
  val => ele%value

  if (out_type == 'XSIF') then
    if (ele%key == elseparator$) then 
      n_elsep_warn = n_elsep_warn + 1
      ele%key = drift$  ! XSIF does not have elsep elements.
      call out_io (s_info$, r_name, 'Elseparator being converted into a drift for XSIF conversion: ' // ele%name)  
    endif
  endif

  ! do not make duplicate specs

  call find_indexx (ele%name, names, an_indexx, n_names, ix_match)
  if (ix_match > 0) cycle

  ! Add to the list of elements

  if (size(names) < n_names + 1) then
    call re_allocate(names, 2*size(names))
    call re_allocate(an_indexx, 2*size(names))
  endif

  call find_indexx (ele%name, names, an_indexx, n_names, ix_match, add_to_list = .true.)

  !----------
  ! OPAL case
  
  if (out_type == 'OPAL-T') then

    select case (ele%key)

    ! OPAL-T
    case (marker$)
      write (line_out, '(a)') trim(ele%name) // ': marker'
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)


    ! OPAL-T
    case (drift$, instrument$, pipe$, detector$, monitor$)
      write (line_out, '(2a)') trim(ele%name) // ': drift, l = ', re_str(val(l$))
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)

    ! OPAL-T
    case (sbend$)
      write (line_out, '(2a)') trim(ele%name) // ': sbend, l = ', re_str(val(l$))
      call value_to_line (line_out, val(b_field$), 'k0', 'R')
      call value_to_line (line_out, val(e_tot$), 'designenergy', 'R')
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)

    ! OPAL-T
    case (quadrupole$)
      write (line_out, '(2a)') trim(ele%name) // ': quadrupole, l = ', re_str(val(l$))
      !Note that OPAL-T has k1 = dBy/dx, and that bmad needs a -1 sign for electrons
      call value_to_line (line_out, -1*val(b1_gradient$), 'k1', 'R')
      !elemedge The edge of the field is specifieda bsolute (floor space co-ordinates) in m.
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)

    ! OPAL-T
    case default
      call out_io (s_error$, r_name, 'UNKNOWN ELEMENT TYPE: ' // key_name(ele%key), &
             'CONVERTING TO DRIFT')
      write (line_out, '(2a)') trim(ele%name) // ': drift, l = ', re_str(val(l$))
      call value_to_line (line_out, ele%s - val(L$), 'elemedge', 'R', .false.)

    end select

    call write_line(line_out)
    cycle
  endif

  !-----------------------------------
  ! For anything else but OPAL

  select case (ele%key)

  ! drift MAD

  case (drift$, instrument$, pipe$, detector$, monitor$)

    write (line_out, '(2a)') trim(ele%name) // ': drift, l = ', re_str(val(l$))
  
  ! beambeam MAD

  case (beambeam$)

    line_out = trim(ele%name) // ': beambeam'
    call value_to_line (line_out, val(sig_x$), 'sigx', 'R')
    call value_to_line (line_out, val(sig_y$), 'sigy', 'R')
    call value_to_line (line_out, val(x_offset$), 'xma', 'R')
    call value_to_line (line_out, val(y_offset$), 'yma', 'R')
    call value_to_line (line_out, val(charge$), 'charge', 'R')


  ! r/ecollimator MAD

  case (ecollimator$, rcollimator$)

    if (out_type == 'MAD-X') then
      write (line_out, '(2a)') trim(ele%name) // ': collimator, l = ', re_str(val(l$))
    else
      write (line_out, '(2a)') trim(ele%name) // ': ' // trim(key_name(ele%key)) // ', l = ', re_str(val(l$))
      call value_to_line (line_out, val(x1_limit$), 'xsize', 'R')
      call value_to_line (line_out, val(y1_limit$), 'ysize', 'R')
    endif

  ! elseparator MAD

  case (elseparator$)

    write (line_out, '(2a)') trim(ele%name) // ': elseparator, l = ', re_str(val(l$))
    hk = val(hkick$)
    vk = val(vkick$)

    if (hk /= 0 .or. vk /= 0) then

      ix = len_trim(line_out) + 1
      field = 1.0d3 * sqrt(hk**2 + vk**2) * val(E_TOT$) / val(l$)
      if (out_type == 'MAD-X') then
        write (line_out(ix:), '(2a)') ', ey = ', re_str(field)
      else
        write (line_out(ix:), '(2a)') ', e = ',re_str(field)
      endif

      if (branch_out%param%particle == positron$) then
        tilt = -atan2(hk, vk) + val(tilt$)
      else
        tilt = -atan2(hk, vk) + val(tilt$) + pi
      endif
      ix = len_trim(line_out) + 1
      write (line_out(ix:), '(2a)') ', tilt = ', re_str(tilt)

    endif

  ! hkicker MAD

  case (hkicker$)

    write (line_out, '(2a)') trim(ele%name) // ': hkicker, l = ', re_str(val(l$))

    call value_to_line (line_out, val(hkick$), 'kick', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! kicker MAD

  case (kicker$)

    write (line_out, '(2a)') trim(ele%name) // ': kicker, l = ', re_str(val(l$))

    call value_to_line (line_out, val(hkick$), 'hkick', 'R')
    call value_to_line (line_out, val(vkick$), 'vkick', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! vkicker MAD

  case (vkicker$)

    write (line_out, '(2a)') trim(ele%name) // ': vkicker, l = ', re_str(val(l$))

    call value_to_line (line_out, val(vkick$), 'kick', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! marker MAD

  case (marker$, fork$, photon_fork$)

    line_out = trim(ele%name) // ': marker'

  ! octupole MAD

  case (octupole$)

    write (line_out, '(2a)') trim(ele%name) // ': octupole, l = ', re_str(val(l$))

    call value_to_line (line_out, val(k3$), 'k3', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! quadrupole MAD

  case (quadrupole$)

    write (line_out, '(2a)') trim(ele%name) // ': quadrupole, l = ', re_str(val(l$))
    call value_to_line (line_out, val(k1$), 'k1', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! sbend MAD

  case (sbend$)

    write (line_out, '(2a)') trim(ele%name) // ': sbend, l = ', re_str(val(l$))

    call value_to_line (line_out, val(angle$), 'angle', 'R')
    call value_to_line (line_out, val(e1$), 'e1', 'R')
    call value_to_line (line_out, val(e2$), 'e2', 'R')
    call value_to_line (line_out, val(k1$), 'k1', 'R')
    call value_to_line (line_out, val(ref_tilt$), 'tilt', 'R')
    if (out_type == 'MAD-X') then
      call value_to_line (line_out, val(fint$), 'fint', 'R')
      call value_to_line (line_out, val(fintx$), 'fintx', 'R')
      call value_to_line (line_out, val(hgap$), 'hgap', 'R')
    else
      if (val(fintx$) /= val(fint$)) then
        call out_io (s_info$, r_name, 'FINTX != FINT FOR BEND' // ele%name, 'CANNOT TRANSLATE FINTX')
      endif
      call value_to_line (line_out, val(fint$), 'fint', 'R')
      call value_to_line (line_out, val(hgap$), 'hgap', 'R')
    endif

  ! sextupole MAD

  case (sextupole$)

    write (line_out, '(2a)') trim(ele%name) // ': sextupole, l = ', re_str(val(l$))
    call value_to_line (line_out, val(k2$), 'k2', 'R')
    call value_to_line (line_out, val(tilt$), 'tilt', 'R')

  ! taylor MAD

  case (taylor$, sad_mult$, patch$)

    if (.not. associated (ele%taylor(1)%term)) then
      if (.not. present(ref_orbit)) then
        call out_io (s_error$, r_name, &
                      'ORBIT ARGUMENT NEEDS TO BE PRESENT WHEN TRANSLATING', &
                      'A LATTICE WITH A SAD_MULT OR PATCH ELEMENT')           
        cycle
      endif
      if (ptc_com%taylor_order_ptc /= 2) call set_ptc (taylor_order = 2) 
      call ele_to_taylor (ele, branch%param, ele%taylor, orbit_out(ix_ele-1), .true.)
    endif

    line_out = trim(ele%name) // ': matrix'
    warn_printed = .false.
    call value_to_line (line_out, ele%value(l$), 'l', 'R')

    do i = 1, 6
      do k = 1, size(ele%taylor(i)%term)
        term = ele%taylor(i)%term(k)

        select case (sum(term%expn))
        case (0)
          select case (out_type)
          case ('MAD-8') 
            write (str, '(a, i0, a)') 'kick(', i, ')'
          case ('MAD-X') 
            write (str, '(a, i0)') 'kick', i
          case ('XSIF') 
            call out_io (s_error$, r_name, 'XSIF DOES NOT HAVE A CONSTRUCT FOR ZEROTH ORDER TAYLOR TERMS NEEDED FOR: ' // ele%name)
            cycle
          end select
          call value_to_line (line_out, term%coef, str, 'R')

        case (1)
          j = maxloc(term%expn, 1)
          select case (out_type)
          case ('MAD-8')
            write (str, '(a, i0, a, i0, a)') 'rm(', i, ',', j, ')'
          case ('MAD-X')
            write (str, '(a, 2i0)') 'rm', i, j
          case ('XSIF')
            write (str, '(a, 2i0)') 'r', i, j
          end select

          if (j == i) then
            call value_to_line (line_out, term%coef, str, 'R', .false.)
          else
            call value_to_line (line_out, term%coef, str, 'R')
          endif

        case (2)
          j = maxloc(term%expn, 1)
          term%expn(j) = term%expn(j) - 1
          j2 = maxloc(term%expn, 1)
          select case (out_type)
          case ('MAD-8')
            write (str, '(a, 3(i0, a))') 'tm(', i, ',', j, ',', j2, ')'
          case ('MAD-X')
            write (str, '(a, 3i0)') 'tm', i, j, j2
          case ('XSIF')
            write (str, '(a, 3i0)') 't', i, j, j2
          end select
          call value_to_line (line_out, term%coef, str, 'R')

        case default
          if (.not. warn_printed .and. ele%key == taylor$) then
            call out_io (s_warn$, r_name, &
                  'Higher order taylor term(s) in: ' // trim(ele%name) // &
                  ' cannot be converted to mad matrix term')
            warn_printed = .true.
          endif  
        end select
      enddo

    enddo

  ! rfcavity MAD

  case (rfcavity$)

    write (line_out, '(2a)') trim(ele%name) // ': rfcavity, l = ', re_str(val(l$))
    call value_to_line (line_out, val(voltage$)/1E6, 'volt', 'R')
    call value_to_line (line_out, val(phi0$)+val(phi0_multipass$)+0.5, 'lag', 'R')
    call value_to_line (line_out, val(harmon$), 'harmon', 'I')

  ! lcavity MAD

  case (lcavity$)

    write (line_out, '(2a)') trim(ele%name) // ': lcavity, l = ', re_str(val(l$))
    call value_to_line (line_out, val(gradient$)*val(l$)/1d6, 'deltae', 'R')
    call value_to_line (line_out, val(rf_frequency$)/1d6, 'freq', 'R')
    call value_to_line (line_out, val(phi0$)+val(phi0_multipass$), 'phi0', 'R')

  ! solenoid MAD

  case (solenoid$)

    write (line_out, '(2a)') trim(ele%name) // ': solenoid, l = ', re_str(val(l$))
    call value_to_line (line_out, val(ks$), 'ks', 'R')

  ! multipole MAD

  case (multipole$, ab_multipole$)

    knl = 0; tilts = 0
    call multipole_ele_to_kt (ele, .true., ix_pole_max, knl, tilts)
    write (line_out, '(2a)') trim(ele%name) // ': multipole'  

    if (out_type == 'MAD-X') then
      knl_str = ''; ksl_str = ''
      call multipole_ele_to_ab (ele, .true., ix_pole_max, a_pole, b_pole)
      do i = 0, 9
        if (all(knl(i:) == 0)) exit
        if (abs(a_pole(i)) < 1d-12 * abs(b_pole(i))) a_pole(i) = 0  ! Round to zero insignificant value
        if (abs(b_pole(i)) < 1d-12 * abs(a_pole(i))) b_pole(i) = 0  ! Round to zero insignificant value
        call value_to_line (knl_str,  b_pole(i) * factorial(i), '', 'R', .false.)
        call value_to_line (ksl_str, -a_pole(i) * factorial(i), '', 'R', .false.)
      enddo
      if (any(b_pole /= 0)) line_out = trim(line_out) // ', knl = {' // trim(knl_str(3:)) // '}'
      if (any(a_pole /= 0)) line_out = trim(line_out) // ', ksl = {' // trim(ksl_str(3:)) // '}'

    else
      do i = 0, 9
        write (str, '(a, i0, a)') 'K', i, 'L'
        call value_to_line (line_out, knl(i), str, 'R')
        write (str, '(a, i0)') 'T', i
        call value_to_line (line_out, tilts(i), str, 'R')
      enddo
    endif

  ! unknown MAD

  case default

    call out_io (s_error$, r_name, 'UNKNOWN ELEMENT TYPE: ' // key_name(ele%key), &
                                  'CONVERTING TO MARKER')
    line_out = trim(ele%name) // ': marker'

  end select

  ! Add apertures for mad-x. Use 1 meter for unset apertures

  if (out_type == 'MAD-X' .and. logic_option(.true., include_apertures)) then
    if (val(x1_limit$) /= 0 .or. val(y1_limit$) /= 0) then
      limit = [val(x1_limit$), val(y1_limit$)]
      where (limit == 0) limit = 1
      if (ele%aperture_type == rectangular$) then
        line_out = trim(line_out) // ', apertype = rectangle'
      else
        line_out = trim(line_out) // ', apertype = ellipse'
      endif
      write (line_out, '(6a)') trim(line_out), ', aperture = (', trim(re_str(limit(1))), ', ', trim(re_str(limit(2))), ')'
    endif
  endif

  ! write element spec to file

  call write_line(line_out)

enddo

!---------------------------------------------------------------------------------------
! Write the lattice line
! MAD has a limit of 4000 characters so we may need to break the lat into pieces.

i_unique = 1000
i_line = 0
init_needed = .true.
line = ' '

do n = ie1, ie2

  ele => branch_out%ele(n)

  if (init_needed) then
    write (iu, '(a)')
    write (iu, '(3a)') comment_char, '---------------------------------', trim(eol_char)
    write (iu, '(a)')
    i_line = i_line + 1
    write (line_out, '(a, i0, 2a)') 'line_', i_line, ': line = (', ele%name
    iout = 0
    init_needed = .false.

  else

    ix = len_trim(line_out) + len_trim(ele%name)

    if (ix > 75) then
      write (iu, '(3a)') trim(line_out), trim(separator_char), trim(continue_char)
      iout = iout + 1
      line_out = '   ' // ele%name
    else
      line_out = trim(line_out) // trim(separator_char) // ' ' // ele%name
    endif
  endif

  ! Output line if long enough or at end

  if (n == ie2 .or. iout > 48) then
    line_out = trim(line_out) // ')'
    write (iu, '(2a)') trim(line_out), trim(eol_char)
    line_out = ' '
    init_needed = .true.
  endif

enddo

! Write twiss parameters for a non-closed lattice.

ele => branch_out%ele(ie1-1)
if (branch_out%param%geometry == open$ .and. &
                  (out_type == 'MAD-8' .or. out_type == 'MAD-X' .or. out_type == 'XSIF')) then
  write (iu, '(a)')
  write (iu, '(3a)') comment_char, '---------------------------------', trim(eol_char)
  write (iu, '(a)')
  write (iu, '(6a)') 'TWISS, betx = ', trim(re_str(ele%a%beta)), ', bety = ', trim(re_str(ele%b%beta)), ', ', trim(continue_char)
  write (iu, '(5x, 6a)') 'alfx = ', trim(re_str(ele%a%alpha)), ', alfy = ', trim(re_str(ele%b%alpha)), ', ', trim(continue_char)
  write (iu, '(5x, 6a)') 'dx = ', trim(re_str(ele%a%eta)), ', dpx = ', trim(re_str(ele%a%etap)), ', ', trim(continue_char)
  write (iu, '(5x, 6a)') 'dy = ', trim(re_str(ele%b%eta)), ', dpy = ', trim(re_str(ele%b%etap)), trim(eol_char)
endif

!------------------------------------------
! Use statement

write (iu, '(a)')
write (iu, '(3a)') comment_char, '---------------------------------', trim(eol_char)
write (iu, '(a)')

line_out = 'lat: line = (line_1'

do i = 2, i_line
  write (line_out, '(3a, i0)') trim(line_out), trim(separator_char), ' line_', i
enddo

line_out = trim(line_out) // ')'
call write_line (line_out)

if (out_type == 'MAD-X') then
  write (iu, '(a)') 'use, period = lat;'
elseif (out_type /= 'OPAL-T') then
  write (iu, '(a)') 'use, lat'
endif

!---------------------------------------------------
! Element offsets for MAD.
! This must come after use statement.

if (out_type(1:3) == 'MAD') then

  write (iu, '(a)')
  write (iu, '(3a)') comment_char, '---------------------------------', trim(eol_char)
  write (iu, '(a)')

  allocate (n_repeat(n_names))
  n_repeat = 0

  do ix_ele = ie1, ie2

    ele => branch_out%ele(ix_ele)
    val => ele%value

    ! sad_mult and patch elements are translated to a matrix which does not have offsets.
    ! And marker like elements also do not have offsets

    if (ele%key == sad_mult$ .or. ele%key == patch$) cycle
    if (ele%key == marker$ .or. ele%key == fork$ .or. ele%key == photon_fork$) cycle

    !

    call find_indexx (ele%name, names, an_indexx, n_names, ix_match)
    n_repeat(ix_match) = n_repeat(ix_match) + 1
    
    if (val(x_pitch$) == 0 .and. val(y_pitch$) == 0 .and. &
        val(x_offset_tot$) == 0 .and. val(y_offset_tot$) == 0 .and. val(z_offset_tot$) == 0) cycle

    write (iu, '(3a)') 'select, flag = error, clear', trim(eol_char)
    write (iu, '(3a, i0, 2a)') 'select, flag = error, range = ', trim(ele%name), &
                                    '[', n_repeat(ix_match), ']', trim(eol_char)

    line_out = 'ealign'
    call value_to_line (line_out,  val(x_pitch$), 'dtheta', 'R')
    call value_to_line (line_out, -val(y_pitch$), 'dphi', 'R')
    call value_to_line (line_out, val(x_offset$) - val(x_pitch$) * val(l$) / 2, 'dx', 'R')
    call value_to_line (line_out, val(y_offset$) - val(y_pitch$) * val(l$) / 2, 'dy', 'R')
    call value_to_line (line_out, val(z_offset$), 'ds', 'R')
    call write_line (line_out)

  enddo

  deallocate (n_repeat)

endif

! End stuff

call out_io (s_info$, r_name, 'Written ' // trim(out_type) // ' lattice file: ' // trim(out_file_name))

deallocate (names)
if (present(err)) err = .false.

if (present(converted_lat)) then
  converted_lat = lat
  converted_lat%branch(branch%ix_branch) = branch_out
  converted_lat%n_ele_max = converted_lat%n_ele_track
  do ib = 0, ubound(converted_lat%branch, 1)
    branch => converted_lat%branch(ib)
    do i = 1, branch%n_ele_track
      branch%ele(i)%slave_status = free$
      branch%ele(i)%n_lord = 0
    enddo
  enddo
  converted_lat%n_control_max = 0
  converted_lat%n_ic_max = 0
endif

call deallocate_lat_pointers (lat_out)
call deallocate_lat_pointers (lat_model)

! Restore ptc settings

if (n_taylor_order_saved /= ptc_com%taylor_order_ptc) call set_ptc (taylor_order = n_taylor_order_saved) 
call set_ptc (exact_modeling = ptc_param%exact_model)

!------------------------------------------------------------------------
contains

subroutine write_line (line_out)

implicit none

character(*) line_out
integer ix, ix1, ix2, ix3

! Prefer to breakup a line after a comma

do
  if (len_trim(line_out) < ix_line_max) exit
  ix1 = index(line_out(ix_line_min+1:), ',')
  ix2 = index(line_out(ix_line_min+1:), '=')
  ix3 = index(line_out(ix_line_min+1:), ' ')

  if (ix1 /= 0 .and. ix1+ix_line_min < ix_line_max) then
    ix = ix1 + ix_line_min
  elseif (ix2 /= 0 .and. ix2+ix_line_min < ix_line_max) then
    ix = ix2 + ix_line_min
  elseif (ix3 /= 0 .and. ix3+ix_line_min < ix_line_max) then
    ix = ix3 + ix_line_min
  elseif (ix1 /= 0) then
    ix = ix1 + ix_line_min
  elseif (ix2 /= 0) then
    ix = ix2 + ix_line_min
  else
    ix = ix3 + ix_line_min
  endif

  write (iu, '(2a)') line_out(:ix), trim(continue_char)
  line_out = '    ' // line_out(ix+1:)
enddo

write (iu, '(2a)') trim(line_out), trim(eol_char)

end subroutine write_line

end subroutine write_lattice_in_foreign_format

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

subroutine value_to_line (line, value, str, typ, ignore_if_zero, use_comma)

use precision_def

implicit none

character(*) line, str
character(40) fmt, val_str
character(*) typ

real(rp) value

integer ix

logical, optional :: ignore_if_zero, use_comma

!

if (value == 0 .and. logic_option(.true., ignore_if_zero)) return

if (logic_option(.true., use_comma)) then
  if (str == '') then
    line = trim(line) // ','
  else
    line = trim(line) // ', ' // trim(str) // ' ='
  endif
else
  if (str /= '') then
    line = trim(line) // ' ' // trim(str) // ' ='
  endif
endif

if (value == 0) then
  line = trim(line) // ' 0'
  return
endif

if (typ == 'R') then
  val_str = re_str(value)
elseif (typ == 'I') then
  write (val_str, '(i0)') nint(value)
else
  print *, 'ERROR IN VALUE_TO_LINE. BAD "TYP": ', typ 
  if (global_com%exit_on_error) call err_exit
endif

call string_trim(val_str, val_str, ix)
line = trim(line) // ' ' // trim(val_str)

end subroutine value_to_line

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+ 
! Subroutine write_lat_in_sad_format (out_file_name, lat, include_apertures, ix_start, ix_end, ix_branch, converted_lat, err)
!
! Private routine used by write_lat_in_foreign_format and not for general use. 
! See write_lat_in_foreign_format for details about the argument list.
!-

subroutine write_lat_in_sad_format (out_file_name, lat, include_apertures, ix_start, ix_end, ix_branch, converted_lat, err)

implicit none

type (lat_struct), target :: lat, lat_model, lat_out
type (lat_struct), optional, target :: converted_lat
type (ele_struct), pointer :: ele, ele1, ele2, lord, sol_ele, first_sol_edge
type (branch_struct), pointer :: branch, branch_out
type (ele_struct), save :: drift_ele, ab_ele, taylor_ele, col_ele, kicker_ele, null_ele, bend_ele, quad_ele
type (ptc_parameter_struct) ptc_param

real(rp), pointer :: val(:)
real(rp) knl(0:n_pole_maxx), tilts(0:n_pole_maxx), a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) bs_field, old_bs_field, a, b
real(rp) xp, yp, xo, yo, zo, sad_fshift
real(rp) hk, vk, tilt


integer, parameter :: sad_f1$  = custom_attribute1$
integer, parameter :: sad_geo_bound$ = custom_attribute2$
integer, parameter :: sad_bz$  = custom_attribute3$
integer, parameter :: sad_fshift$ = custom_attribute4$
integer, parameter :: sad_mark_offset$  = custom_attribute5$

integer, optional :: ix_start, ix_end, ix_branch
integer, allocatable :: n_repeat(:), an_indexx(:)
integer i, j, n, ib, iout, iu, ix, ix1, ix2, ios, ie1, ie2, n_taylor_order_saved, ix_ele, n_elsep_warn, n_name_change_warn
integer s_count, t_count, j_count
integer ix_manch, n_names, aperture_at, ix_pole_max, ix_match
integer :: ix_line_min, ix_line_max, n_warn_max = 10

character(*), parameter :: r_name = "write_lat_in_sad_format"
character(*) out_file_name
character(40), allocatable :: names(:)
character(40) str, orig_name
character(300) line, knl_str, ksl_str
character(2000) line_out

logical, optional :: include_apertures, err
logical converted, init_needed, in_solenoid

! open file

if (present(err)) err = .true.
n_taylor_order_saved = ptc_com%taylor_order_ptc

iu = lunget()
call fullfilename (out_file_name, line)
open (iu, file = line, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(out_file_name))
  return
endif

! Use ptc exact_model = True since this is needed to get the drift nonlinear terms

call get_ptc_params(ptc_param)
call set_ptc (exact_modeling = .true.)

! Init

ix = integer_option(0, ix_branch)
if (ix < 0 .or. ix > ubound(lat%branch, 1)) then
  call out_io (s_error$, r_name, 'BRANCH INDEX OUT OF RANGE: /i0/ ', i_array = [ix])
  return
endif

branch => lat%branch(ix)

ix_line_max = 120
ix_line_min = ix_line_max - 20

call init_ele (col_ele)
call init_ele (drift_ele, drift$)
call init_ele (ab_ele, ab_multipole$)
call init_ele (kicker_ele, kicker$) 
call init_ele (quad_ele, quadrupole$)
call init_ele (bend_ele, sbend$)
call multipole_init (ab_ele, magnetic$)
null_ele%key = null_ele$

ie1 = integer_option(1, ix_start)
ie2 = integer_option(branch%n_ele_track, ix_end)

allocate (names(branch%n_ele_max), an_indexx(branch%n_ele_max)) ! list of element names

call out_io (s_info$, r_name, &
      'Note: Bmad lattice elements have attributes that cannot be translated. ', &
      '      For example, higher order terms in a Taylor element.', &
      '      Please use caution when using a translated lattice.')

!-----------------------------------------------------------------------------
! Translation is a two step process:
!   1) Create a new lattice called lat_out making substitutions for sol_quad and wiggler elements, etc..
!   2) Use lat_out to create the lattice file.

lat_out = lat
call allocate_lat_ele_array(lat_out, 2*branch%n_ele_max, branch%ix_branch)
branch_out => lat_out%branch(branch%ix_branch)

j_count = 0    ! drift around solenoid or sol_quad index. Also z shift count.
t_count = 0    ! taylor element count.
s_count = 0    ! SAD solenoid count
sad_fshift = 0
in_solenoid = .false.

! Loop over all input elements

branch_out%ele%value(sad_mark_offset$) = 0  ! SAD mark offset
nullify(first_sol_edge)
old_bs_field = 0
n_name_change_warn = 0
n_elsep_warn = 0
ix_ele = ie1 - 1

do

  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => branch_out%ele(ix_ele)
  val => ele%value

  if (ele%sub_key == not_set$) cycle ! See below

  ! Patch element that was inserted in SAD to Bmad translation for fshift is to be removed 
  ! when translating back.

  if (ele%key == patch$ .and. ele%value(sad_fshift$) /= 0) then
    sad_fshift = ele%value(sad_fshift$)
    ele%sub_key = not_set$
    cycle
  endif

  ! Must create SAD sol elements as needed. 
  ! If there is a marker or patch element that can be used, convert it.
  ! Otherwise, insert a new element. 
  ! Bmad does not have a sol element so use null_ele to designate the sol element in the Bmad lattice.
  ! This works since there cannot be any actual null_eles in the lattice.

  if ((ele%key == patch$ .or. ele%key == marker$) .and. ele%value(sad_geo_bound$) /= 0) then 
    ele%key = null_ele$  ! Convert to SOL
    if (ele%value(sad_geo_bound$) > 0) then
      if (in_solenoid) then
        ele%iyy = exit_end$
      else
        ele%iyy = entrance_end$
      endif
      in_solenoid = .not. in_solenoid
    endif

  elseif (ele%key == patch$ .or. ele%value(l$) /= 0 .or. ix_ele == branch_out%n_ele_max) then
    bs_field = 0
    if (has_attribute (ele, 'BS_FIELD')) bs_field = ele%value(bs_field$)

    if (bs_field /= old_bs_field) then
      if (ele%key == marker$ .or. ele%key == patch$) then
        sol_ele => ele
      else
        sol_ele => pointer_to_next_ele(ele, -1)
        ! Look to see if there is a marker or patch element that can be converted to a SAD SOL.
        do
          if (sol_ele%ix_ele == 1) exit
          if (sol_ele%key == marker$ .or. sol_ele%key == patch$) exit
          if (sol_ele%value(l$) /= 0) exit  ! No suitable marker or patch found
          sol_ele => pointer_to_next_ele(sol_ele, -1)
        enddo
        sol_ele%value(sad_geo_bound$) = 0
      endif

      if (sol_ele%key /= marker$ .and. sol_ele%key /= patch$) then
        s_count = s_count + 1
        write (sol_ele%name, '(a, i0)') 'SOL_', s_count  
        call insert_element (lat_out, sol_ele, ix_ele, branch_out%ix_branch)
        sol_ele => branch_out%ele(ix_ele)
        ie2 = ie2 + 1
        ix_ele = ix_ele + 1
        ele => branch_out%ele(ix_ele)
        val => ele%value
      endif

      if (old_bs_field == 0) then
        sol_ele%value(sad_geo_bound$) = 2
        sol_ele%iyy = entrance_end$  ! Entering solenoid
        in_solenoid = .true.
        first_sol_edge => sol_ele
      elseif (bs_field == 0) then
        if (nint(ele%value(sad_geo_bound$)) == 2) then
          sol_ele%value(sad_geo_bound$) = 2
          first_sol_edge%value(sad_geo_bound$) = 1
        else
          sol_ele%value(sad_geo_bound$) = 1
        endif
        sol_ele%iyy = exit_end$  ! Entering solenoid
        in_solenoid = .false.
      else
        sol_ele%value(sad_geo_bound$) = 0
      endif

      sol_ele%value(bs_field$) = bs_field
      if (bs_field == 0) sol_ele%value(bs_field$) = sol_ele%value(sad_bz$)
      sol_ele%key = null_ele$

      old_bs_field = bs_field
    endif
  endif

  ! With an element superimposed with a marker create a whole element
  ! plus a marker with an offset

  if (ele%slave_status == super_slave$) then
    lord => pointer_to_lord(ele, 1, ix_slave = ix)
    ele1 => pointer_to_next_ele(ele)
    ele2 => pointer_to_slave(lord, lord%n_slave)
    if (lord%n_slave == 2 .and. ele1%key == marker$ .and. &
          num_lords(ele, super_lord$) == 1 .and. num_lords(ele2, super_lord$) == 1) then
      ele1%value(sad_mark_offset$) = -ele2%value(l$)/lord%value(l$) ! marker offset
      ele = lord                              ! Super_slave piece becomes the entire element.
      ele2%sub_key = not_set$                 ! Ignore this super_slave piece.
    endif
  endif

  ! Replace element name containing "/" or "#" with "_"

  orig_name = ele%name

  do
    j = index (ele%name, '\')         ! '
    j = index (ele%name, '#')   
    if (j == 0) exit
    ele%name(j:j) = '_'
  enddo

  if (ele%name /= orig_name .and. n_name_change_warn <= n_warn_max) then
    call out_io (s_info$, r_name, 'Element name changed from: ' // trim(orig_name) // ' to: ' // ele%name)
    if (n_name_change_warn == n_warn_max) call out_io (s_info$, r_name, &
                           'Enough name change warnings. Will stop issuing them now.')
    n_name_change_warn = n_name_change_warn + 1
  endif

  ! SAD: If there is an aperture with an element that is not an ecoll or rcoll then need to make a separate
  ! element with the aperture info. 

  if ((val(x1_limit$) /= 0 .or. val(x2_limit$) /= 0 .or. val(y1_limit$) /= 0 .or. val(y2_limit$) /= 0) .and. &
      .not. (((ele%key == ecollimator$ .or. ele%key == rcollimator$) .and. ele%value(l$) == 0) .or. ele%key == sad_mult$) .and. &
      logic_option(.true., include_apertures)) then

    if (val(x1_limit$) /= val(x2_limit$)) then
      call out_io (s_warn$, r_name, 'Asymmetric x_limits cannot be converted for: ' // ele%name, &
                                    'Will use largest limit here.')
      val(x1_limit$) = max(val(x1_limit$), val(x2_limit$))
    endif

    if (val(y1_limit$) /= val(y2_limit$)) then
      call out_io (s_warn$, r_name, 'Asymmetric y_limits cannot be converted for: ' // ele%name, &
                                    'Will use largest limit here.')
      val(y1_limit$) = max(val(y1_limit$), val(y2_limit$))
    endif

    ! Create ecoll and rcoll elements.
    ! If original element is itself a collimator, turn it into a drift.

    if (ele%aperture_type == rectangular$) then
      col_ele%key = rcollimator$
    else
      col_ele%key = ecollimator$
    endif

    if (ele%key == ecollimator$ .or. ele%key == rcollimator$) then
      col_ele%name = ele%name
      ele%key = drift$
      ele%name = 'DRIFT_' // trim(ele%name)
    else
      col_ele%name = 'COLLIMATOR_' // trim(ele%name)
    endif

    col_ele%value = val
    col_ele%value(l$) = 0
    val(x1_limit$) = 0; val(x2_limit$) = 0; val(y1_limit$) = 0; val(y2_limit$) = 0; 
    aperture_at = ele%aperture_at  ! Save since ele pointer will be invalid after the insert
    if (aperture_at == both_ends$ .or. aperture_at == downstream_end$ .or. aperture_at == continuous$) then
      call insert_element (lat_out, col_ele, ix_ele+1, branch_out%ix_branch)
      ie2 = ie2 + 1
    endif
    if (aperture_at == both_ends$ .or. aperture_at == upstream_end$ .or. aperture_at == continuous$) then
      call insert_element (lat_out, col_ele, ix_ele, branch_out%ix_branch)
      ie2 = ie2 + 1
    endif
    ix_ele = ix_ele - 1 ! Want to process the element again on the next loop.

    cycle ! cycle since ele pointer is invalid
  endif

  ! If the bend has a roll then put kicker elements just before and just after

  if (ele%key == sbend$ .and. val(roll$) /= 0) then
    j_count = j_count + 1
    write (kicker_ele%name,   '(a, i0)') 'ROLL_Z', j_count
    kicker_ele%value(hkick$) =  val(angle$) * (1 - cos(val(roll$))) / 2
    kicker_ele%value(vkick$) = -val(angle$) * sin(val(roll$)) / 2
    val(roll$) = 0   ! So on next iteration will not create extra kickers.
    call insert_element (lat_out, kicker_ele, ix_ele, branch_out%ix_branch)
    call insert_element (lat_out, kicker_ele, ix_ele+2, branch_out%ix_branch)
    ie2 = ie2 + 2
    cycle
  endif

  ! If there is a multipole component then put multipole elements at half strength 
  ! just before and just after the element.

  if (ele%key /= multipole$ .and. ele%key /= ab_multipole$ .and. ele%key /= null_ele$ .and. ele%key /= sad_mult$) then
    call multipole_ele_to_ab (ele, .true., ix_pole_max, ab_ele%a_pole, ab_ele%b_pole)
    if (ix_pole_max > -1) then
      ab_ele%a_pole = ab_ele%a_pole / 2
      ab_ele%b_pole = ab_ele%b_pole / 2
      if (associated(ele%a_pole)) deallocate (ele%a_pole, ele%b_pole)
      j_count = j_count + 1
      write (ab_ele%name, '(a1, a, i0)') key_name(ele%key), 'MULTIPOLE_', j_count
      call insert_element (lat_out, ab_ele, ix_ele, branch_out%ix_branch)
      call insert_element (lat_out, ab_ele, ix_ele+2, branch_out%ix_branch)
      ie2 = ie2 + 2
      cycle
    endif
  endif

  ! Convert wiggler elements to an "equivalent" set of elements.
  ! If the wiggler has been sliced due to superposition, throw 
  ! out the markers that caused the slicing.


  if (ele%key == wiggler$ .or. ele%key == undulator$) then

    call out_io (s_warn$, r_name, 'Converting element to drift-bend-drift model: ' // ele%name)
    if (ele%slave_status == super_slave$) then
      ! Create the wiggler model using the super_lord
      lord => pointer_to_lord(ele, 1)
      call create_wiggler_model (lord, lat_model)
      ! Remove all the slave elements and markers in between.
      call out_io (s_warn$, r_name, 'Note: Not translating the markers within wiggler: ' // lord%name)
      lord%key = -1 ! mark for deletion
      call find_element_ends (lord, ele1, ele2)
      ix1 = ele1%ix_ele; ix2 = ele2%ix_ele
      ! If the wiggler wraps around the origin we are in trouble.
      if (ix2 < ix1) then 
        call out_io (s_fatal$, r_name, 'Wiggler wraps around origin. Cannot translate this!')
        if (global_com%exit_on_error) call err_exit
      endif
      do i = ix1+1, ix2
        branch_out%ele(i)%key = -1  ! mark for deletion
      enddo
      ie2 = ie2 - (ix2 - ix1 - 1)
    else
      call create_wiggler_model (ele, lat_model)
    endif

    ele%key = -1 ! Mark for deletion
    call remove_eles_from_lat (lat_out)
    do j = 1, lat_model%n_ele_track
      call insert_element (lat_out, lat_model%ele(j), ix_ele+j-1, branch_out%ix_branch)
    enddo
    ie2 = ie2 + lat_model%n_ele_track - 1
    cycle
  endif

enddo

! If there is a finite bs_field then create a final null_ele element

if (bs_field /= 0) then
  ele => branch_out%ele(ie2)
  if (ele%key == marker$) then
    ele%key = null_ele$
    ele%value(bs_field$) = bs_field
  else
    s_count = s_count + 1
    write (null_ele%name, '(a, i0)') 'SOL_', s_count  
    null_ele%value(bs_field$) = bs_field
    ie2 = ie2 + 1
    call insert_element (lat_out, null_ele, ie2, branch_out%ix_branch)
    ele => branch_out%ele(ie2)
  endif

  ele%value(sad_geo_bound$) = 1
endif

! For a patch that is *not* associated with the edge of a solenoid: A z_offset must be split into a drift + patch

ix_ele = ie1 - 1

do
  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => branch_out%ele(ix_ele)
  val => ele%value

  if (ele%key == patch$ .and. ele%value(z_offset$) /= 0) then
    drift_ele%name = 'DRIFT_' // ele%name
    drift_ele%value(l$) = val(z_offset$)
    call insert_element (lat_out, drift_ele, ix_ele, branch_out%ix_branch)
    ix_ele = ix_ele + 1
    ele => branch_out%ele(ix_ele)
    val => ele%value
    val(z_offset$) = 0
  endif
enddo

!-------------------------------------------------------------------------------------------------
! SAD: Now write info to the output file...
! lat lattice name

write (iu, '(3a)') '! File generated by: write_lattice_in_foreign_format;'
write (iu, '(4a)') '! Bmad Lattice File: ', trim(lat%input_file_name), ';'
if (lat%lattice /= '') write (iu, '(4a)') '! Bmad Lattice: ', trim(lat%lattice), ';'
write (iu, '(a)')

write (iu, '(3a)') 'MOMENTUM = ',  trim(re_str(ele%value(p0c$))), ';'
if (sad_fshift /= 0) write (iu, '(3a)') 'FSHIFT = ', trim(re_str(sad_fshift)), ';'

! write element parameters

n_names = 0                          ! number of names stored in the list
old_bs_field = 0

do ix_ele = ie1, ie2

  ele => branch_out%ele(ix_ele)
  val => ele%value

  if (ele%sub_key == not_set$) cycle 

  ! do not make duplicate specs

  call find_indexx (ele%name, names, an_indexx, n_names, ix_match)
  if (ix_match > 0) cycle

  ! Add to the list of elements

  if (size(names) < n_names + 1) then
    call re_allocate(names, 2*size(names))
    call re_allocate(an_indexx, 2*size(names))
  endif

  call find_indexx (ele%name, names, an_indexx, n_names, ix_match, add_to_list = .true.)

  converted = .false.
  if (.not. associated (ele%a_pole) .and. ele%value(hkick$) == 0 .and. ele%value(vkick$) == 0) then
    select case (ele%key)

    ! SAD
    case (octupole$)
      write (line_out, '(4a)') 'OCT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(k3$)*val(l$), 'K3', 'R', .true., .false.)
      converted = .true.

    ! SAD
    case (quadrupole$)
      write (line_out, '(4a)') 'QUAD ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(k1$)*val(l$), 'K1', 'R', .true., .false.)
      call value_to_line (line_out, -sign_of(val(fq1$)) * sqrt(24*abs(val(fq1$))), 'F1', 'R', .true., .false.)
      call value_to_line (line_out, val(fq2$), 'F2', 'R', .true., .false.)
      converted = .true.

    ! SAD
    case (sextupole$)
      write (line_out, '(4a)') 'SEXT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(k2$)*val(l$), 'K2', 'R', .true., .false.)
      converted = .true.
    end select
  endif

  ! If not yet converted due to the presence of multipoles or a kick

  if (.not. converted) then

    a_pole = 0; b_pole = 0
    if (ele%key /= null_ele$) call multipole_ele_to_ab (ele, .false., ix_pole_max, a_pole, b_pole)

    select case (ele%key)

    ! SAD
    case (drift$, pipe$)
      write (line_out, '(4a)') 'DRIFT ', trim(ele%name), ' = (L = ', re_str(val(l$))

    ! SAD
    case (instrument$, detector$, monitor$)
      if (ele%value(l$) == 0) then
        write (line_out, '(4a)') 'MONI ', trim(ele%name), ' = ('
      else
        write (line_out, '(4a)') 'DRIFT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      endif

    ! SAD
    case (ab_multipole$, multipole$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = ('
      do i = 0, ubound(a_pole, 1)
        write (str, '(a, i0)') 'K', i
        call value_to_line (line_out, a_pole(i) * factorial(i), 'S' // str, 'R', .true., .false.)
        call value_to_line (line_out, b_pole(i) * factorial(i), str, 'R', .true., .false.)
      enddo

    ! SAD
    case (bend_sol_quad$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call multipole1_kt_to_ab (val(angle$), val(bend_tilt$), 0, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b
      call multipole1_kt_to_ab (val(k1$)*val(l$), val(quad_tilt$), 1, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b
      if (val(dks_ds$) /= 0) call out_io (s_error$, r_name, &
                        'DKS_DS OF BEND_SOL_QUAD CANNOT BE CONVERTED FOR: ' // ele%name)
      if (val(x_quad$) /= 0 .or. val(y_quad$) /= 0) call out_io (s_error$, r_name, &
                        'X/Y_QUAD OF BEND_SOL_QUAD CANNOT BE CONVERTED FOR: ' // ele%name)

    ! SAD
    case (ecollimator$)
      write (line_out, '(4a)') 'APERT ', trim(ele%name), ' = ('
      call value_to_line (line_out, val(x_offset$), 'DX', 'R', .true., .false.)
      call value_to_line (line_out, val(y_offset$), 'DY', 'R', .true., .false.)
      call value_to_line (line_out, val(x1_limit$), 'AX', 'R', .true., .false.)
      call value_to_line (line_out, val(y1_limit$), 'AY', 'R', .true., .false.)
      
    ! SAD
    case (rcollimator$)
      write (line_out, '(4a)') 'APERT ', trim(ele%name), ' = ('
      call value_to_line (line_out, -val(x1_limit$), 'DX1', 'R', .true., .false.)
      call value_to_line (line_out, -val(y1_limit$), 'DY1', 'R', .true., .false.)
      call value_to_line (line_out,  val(x2_limit$), 'DX2', 'R', .true., .false.)
      call value_to_line (line_out,  val(y2_limit$), 'DY2', 'R', .true., .false.)

    ! SAD
    case (elseparator$)
      call out_io (s_warn$, r_name, 'Elseparator will be converted into a mult: ' // ele%name)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call multipole1_kt_to_ab (-val(hkick$), 0.0_rp, 0, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b
      call multipole1_kt_to_ab (-val(vkick$), pi/2, 0, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    ! SAD
    case (hkicker$)
      write (line_out, '(4a)') 'BEND ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, -val(kick$), 'K0', 'R', .true., .false.)
!      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
!      call multipole1_kt_to_ab (-val(kick$), 0.0_rp, 0, a, b)
!      a_pole = a_pole + a;  b_pole = b_pole + b

    ! SAD
    case (vkicker$)
      tilt = -val(tilt$) - pi/2
      write (line_out, '(4a)') 'BEND ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, -val(kick$), 'K0', 'R', .true., .false.)
      call value_to_line (line_out, tilt, 'ROTATE', 'R', .true., .false.)
!      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
!      call multipole1_kt_to_ab (-val(kick$), pi/2, 0, a, b)
!      a_pole = a_pole + a;  b_pole = b_pole + b

    ! SAD
    case (kicker$)
      write (line_out, '(4a)') 'BEND ', trim(ele%name), ' = (L = ', re_str(val(l$))
      hk = -val(hkick$)
      vk = -val(vkick$)
      tilt = atan2(hk, vk) + pi/2 - val(tilt$)
      if (hk /= 0 .or. vk /= 0) then
        call value_to_line (line_out, -sqrt(hk**2 + vk**2), 'K0', 'R', .true., .false.)
        call value_to_line (line_out, tilt, 'ROTATE', 'R', .true., .false.)
      endif
!      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
!      call multipole1_kt_to_ab (-val(hkick$), 0.0_rp, 0, a, b)
!      a_pole = a_pole + a;  b_pole = b_pole + b
!      call multipole1_kt_to_ab (-val(vkick$), pi/2, 0, a, b)
!      a_pole = a_pole + a;  b_pole = b_pole + b

    ! SAD
    case (lcavity$)
      write (line_out, '(4a)') 'CAVI ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(rf_frequency$), 'FREQ', 'R', .true., .false.)
      call value_to_line (line_out, val(voltage$), 'VOLT', 'R', .true., .false.)
      call value_to_line (line_out, 0.25 - val(phi0$), 'PHI', 'R', .true., .false.)
      call value_to_line (line_out, -val(phi0_err$), 'DPHI', 'R', .true., .false.)

    ! SAD
    case (marker$)
      write (line_out, '(4a)') 'MARK ', trim(ele%name), ' = ('
      call value_to_line (line_out, val(sad_mark_offset$), 'OFFSET', 'R', .true., .false.)
      call value_to_line (line_out, val(sad_geo_bound$), 'GEO', 'R', .true., .false.)
      if (branch_out%param%geometry == open$ .and. ix_ele == 1) then
        call value_to_line (line_out, ele%a%beta, 'BX', 'R', .true., .false.)
        call value_to_line (line_out, ele%b%beta, 'BY', 'R', .true., .false.)
        call value_to_line (line_out, ele%a%alpha, 'AX', 'R', .true., .false.)
        call value_to_line (line_out, ele%b%alpha, 'AY', 'R', .true., .false.)
        call value_to_line (line_out, ele%x%eta, 'PEX', 'R', .true., .false.)
        call value_to_line (line_out, ele%y%eta, 'PEY', 'R', .true., .false.)
        call value_to_line (line_out, ele%x%etap, 'PEPX', 'R', .true., .false.)
        call value_to_line (line_out, ele%y%etap, 'PEPY', 'R', .true., .false.)
        call value_to_line (line_out, lat%a%emit, 'EMITX', 'R', .true., .false.)
        call value_to_line (line_out, lat%b%emit, 'EMITy', 'R', .true., .false.)
      endif

    ! SAD
    case (null_ele$)
      write (line_out, '(4a)') 'SOL ', trim(ele%name), ' = ('
      call value_to_line (line_out, val(bs_field$), 'BZ', 'R', .true., .false.)
      call value_to_line (line_out, ele%value(sad_f1$), 'F1', 'R', .true., .false.)
      select case (nint(ele%value(sad_geo_bound$)))
      case (1)
        line_out = trim(line_out) // ' BOUND = 1'
      case (2)
        line_out = trim(line_out) // ' BOUND = 1'
        line_out = trim(line_out) // ' GEO = 1'
      end select

    ! SAD with nonzero multipoles or kick
    case (octupole$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call multipole1_kt_to_ab (val(k3$), 0.0_rp, 3, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    ! SAD
    case (patch$)
      write (line_out, '(4a)') 'COORD ', trim(ele%name), ' = ('

    ! SAD with nonzero multipoles or kick
    case (quadrupole$) 
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, -sign_of(val(fq1$)) * sqrt(24*abs(val(fq1$))), 'F1', 'R', .true., .false.)
      call value_to_line (line_out, val(fq2$), 'F2', 'R', .true., .false.)

      call multipole1_kt_to_ab (val(k1$), 0.0_rp, 1, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    ! SAD
    case (rfcavity$)
      write (line_out, '(4a)') 'CAVI ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(rf_frequency$), 'FREQ', 'R', .true., .false.)
      call value_to_line (line_out, val(voltage$), 'VOLT', 'R', .true., .false.)
      call value_to_line (line_out, twopi * val(phi0$), 'DPHI', 'R', .true., .false.)

      select case (nint(val(fringe_at$)))
      case (entrance_end$)
        line_out = trim(line_out) // ' fringe = 1'
      case (exit_end$)
        line_out = trim(line_out) // ' fringe = 2'
      end select

    ! SAD
    case (sad_mult$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, -sign_of(val(fq1$)) * sqrt(24*abs(val(fq1$))), 'F1', 'R', .true., .false.)
      call value_to_line (line_out, val(fq2$), 'F2', 'R', .true., .false.)
      if (val(eps_step_scale$) /= 1) call value_to_line (line_out, val(eps_step_scale$), 'EPS', 'R', .true., .false.)
      call value_to_line (line_out, val(x_offset_mult$), 'DX', 'R', .true., .false.)
      call value_to_line (line_out, val(y_offset_mult$), 'DY', 'R', .true., .false.)
      call value_to_line (line_out, val(x_pitch_mult$), 'CHI1', 'R', .true., .false.)
      call value_to_line (line_out, val(y_pitch_mult$), 'CHI2', 'R', .true., .false.)
      if (val(x1_limit$) == val(y1_limit$)) then
        call value_to_line (line_out, val(x1_limit$), 'RADIUS', 'R', .true., .false.)
      else
        call out_io (s_warn$, r_name, 'Asymmetric x_limit vs y_limit cannot be converted for: ' // ele%name, &
                                  'Will use largest limit here.')
        if (val(x1_limit$) /= 0 .and. val(y1_limit$) /= 0) then
          call value_to_line (line_out, max(val(x1_limit$), val(y1_limit$)), 'RADIUS', 'R', .true., .false.)
        endif
      endif

    ! SAD
    case (sbend$)
      write (line_out, '(4a)') 'BEND ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(angle$), 'ANGLE', 'R', .true., .false.)
      call value_to_line (line_out, val(g_err$)*val(l$), 'K0', 'R', .true., .false.)
      call value_to_line (line_out, val(k1$)*val(l$), 'K1', 'R', .true., .false.)
      call value_to_line (line_out, -val(ref_tilt$), 'ROTATE', 'R', .true., .false.)
      if (val(fintx$)*val(hgapx$) == val(fint$)*val(hgap$)) then
        call value_to_line (line_out, 12*val(fint$)*val(hgap$), 'F1', 'R', .true., .false.)
      else
        call value_to_line (line_out, 12*val(fint$)*val(hgap$), 'FB1', 'R', .true., .false.)
        call value_to_line (line_out, 12*val(fintx$)*val(hgapx$), 'FB2', 'R', .true., .false.)
      endif
      if (val(angle$) == 0) then
        call value_to_line (line_out, val(e1$), 'AE1', 'R', .true., .false.)
        call value_to_line (line_out, val(e2$), 'AE2', 'R', .true., .false.)
      else
        call value_to_line (line_out, val(e1$)/val(angle$), 'E1', 'R', .true., .false.)
        call value_to_line (line_out, val(e2$)/val(angle$), 'E2', 'R', .true., .false.)
      endif

      select case (nint(val(fringe_type$)))
      case (hard_edge_only$)
        ! Nothing to be done
      case (soft_edge_only$)
        line_out = trim(line_out) // ' FRINGE = 1 DISFRIN = 1'
      case (sad_full$)
        line_out = trim(line_out) // ' FRINGE = 1'
      end select

    ! SAD with nonzero multipoles or kick
    case (sextupole$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call multipole1_kt_to_ab (val(k3$), 0.0_rp, 1, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    ! SAD
    case (solenoid$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))

    ! SAD
    case (sol_quad$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call multipole1_kt_to_ab (val(k1$), 0.0_rp, 1, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    ! SAD
    case default
      call out_io (s_error$, r_name, 'UNKNOWN ELEMENT TYPE: ' // key_name(ele%key), &
             'CONVERTING TO DRIFT')
      write (line_out, '(4a)') 'DRIFT ', trim(ele%name), ' = (L = ', re_str(val(l$))
    end select

    if (line_out(1:4) == 'MULT') then
      if (has_attribute (ele, 'HKICK') .and. ele%key /= kicker$) then
        call multipole1_kt_to_ab (-val(hkick$), -val(tilt_tot$), 0, a, b)
        a_pole = a_pole + a;  b_pole = b_pole + b
        call multipole1_kt_to_ab (-val(vkick$), pi/2-val(tilt_tot$), 0, a, b)
        a_pole = a_pole + a;  b_pole = b_pole + b
      endif

      do i = 0, 21
        write (str, '(i0)') i
        call value_to_line (line_out, b_pole(i)*factorial(i), 'K'//trim(str), 'R', .true., .false.)
        call value_to_line (line_out, a_pole(i)*factorial(i), 'SK'//trim(str), 'R', .true., .false.)
      enddo
    endif
  endif     ! Not converted

  ! fringe

  if (ele%key == quadrupole$ .or. ele%key == sad_mult$) then
    select case (nint(val(fringe_type$)))
    case (soft_edge_only$, full$)
      select case (nint(val(fringe_at$)))
      case (entrance_end$);         line_out = trim(line_out) // ' FRINGE = 1'
      case (exit_end$);             line_out = trim(line_out) // ' FRINGE = 2'
      case (both_ends$);            line_out = trim(line_out) // ' FRINGE = 3'
      end select
    end select

    select case (nint(val(fringe_type$)))
    case (none$, soft_edge_only$)
      line_out = trim(line_out) // ' DISFRIN = 1'
    case (hard_edge_only$)
      select case (nint(val(fringe_at$)))
      case (no_end$)
        line_out = trim(line_out) // ' DISFRIN = 1'
      case (entrance_end$, exit_end$)
        line_out = trim(line_out) // ' CANNOT TRANSLATE BMAD FRINGE_TYPE/FRINGE_AT!'
      end select
    case (full$)
      if (nint(val(fringe_at$)) == no_end$) line_out = trim(line_out) // ' DISFRIN = 1'
    end select
  endif

  ! misalignments
  ! Note: SAD applies pitches and offsets in reverse order to Bmad.

  xp = val(x_pitch$);  yp = val(y_pitch$)

  if (xp /= 0) then
    zo =  val(z_offset$) * cos(xp) + val(x_offset$) * sin(xp)
    xo = -val(z_offset$) * sin(xp) + val(x_offset$) * cos(xp)
    val(z_offset$) = zo
    val(x_offset$) = xo
  endif

  if (yp /= 0) then
    zo =  val(z_offset$) * cos(yp) + val(y_offset$) * sin(yp)
    yo = -val(z_offset$) * sin(yp) + val(y_offset$) * cos(yp)
    val(z_offset$) = zo
    val(y_offset$) = yo
  endif

  if (ele%key == null_ele$) then ! patch -> SOL
    if (ele%iyy == entrance_end$) then
      val(x_offset$) = -val(x_offset$)
      val(y_offset$) = -val(y_offset$)
      val(z_offset$) = -val(z_offset$)
    else
      val(z_offset$) = -val(z_offset$)
    endif
    val(z_offset$) = val(z_offset$) + ele%value(t_offset$) * c_light
  endif

  call value_to_line (line_out, val(x_offset$), 'DX', 'R', .true., .false.)
  call value_to_line (line_out, val(y_offset$), 'DY', 'R', .true., .false.)
  call value_to_line (line_out, val(z_offset$), 'DZ', 'R', .true., .false.)

  if (ele%key /= marker$) then
    call value_to_line (line_out, -val(x_pitch$), 'CHI1', 'R', .true., .false.)
    call value_to_line (line_out, -val(y_pitch$), 'CHI2', 'R', .true., .false.)
  endif

  if (ele%key == patch$ .or. ele%key == null_ele$) then   ! null_ele -> SOL
    call value_to_line (line_out, -val(tilt$),    'CHI3', 'R', .true., .false.)
  elseif (ele%key /= sbend$ .and. ele%key /= kicker$ .and. ele%key /= vkicker$) then
    call value_to_line (line_out, -val(tilt$),    'ROTATE', 'R', .true., .false.)
  endif

  if (ele%key == null_ele$ .and. nint(ele%value(sad_geo_bound$)) == 2) then ! At geo = 1 end
    if (ele%iyy == entrance_end$) then
      call value_to_line (line_out, val(x_pitch$), 'DPX', 'R', .true., .false.)
      call value_to_line (line_out, val(y_pitch$), 'DPY', 'R', .true., .false.)
    else
      call value_to_line (line_out, -val(x_pitch$), 'DPX', 'R', .true., .false.)
      call value_to_line (line_out, -val(y_pitch$), 'DPY', 'R', .true., .false.)
    endif
  endif      

  !

  line_out = trim(line_out) // ')'
  call write_line(line_out)
  cycle

  ! write element spec to file

  call write_line(line_out)

enddo

!---------------------------------------------------------------------------------------
! Write the lattice line

write (iu, '(a)')
write (iu, '(3a)') '!---------------------------------;'
write (iu, '(a)')
write (line_out, '(2a)') 'LINE ASC = ('

do n = ie1, ie2

  ele => branch_out%ele(n)
  if (ele%sub_key == not_set$) cycle

  ix = len_trim(line_out) + len_trim(ele%name)

  if (ix > 75) then
    write (iu, '(3a)') trim(line_out)
    line_out = '   ' // ele%name
  else
    line_out = trim(line_out) // ' ' // ele%name
  endif

  ! Output line if long enough or at end

  if (n == ie2) then
    line_out = trim(line_out) // ')'
    write (iu, '(2a)') trim(line_out), ';'
    line_out = ' '
    init_needed = .true.
  endif

enddo

!------------------------------------------
! Use statement

write (iu, '(a)')
write (iu, '(3a)') '!---------------------------------;'
write (iu, '(a)')

write (iu, '(a)') 'FFS USE ASC;'

! End stuff

call out_io (s_info$, r_name, 'Written SAD lattice file: ' // trim(out_file_name))

deallocate (names)
if (present(err)) err = .false.

if (present(converted_lat)) then
  converted_lat = lat
  converted_lat%branch(branch%ix_branch) = branch_out
  converted_lat%n_ele_max = converted_lat%n_ele_track
  do ib = 0, ubound(converted_lat%branch, 1)
    branch => converted_lat%branch(ib)
    do i = 1, branch%n_ele_track
      branch%ele(i)%slave_status = free$
      branch%ele(i)%n_lord = 0
    enddo
  enddo
  converted_lat%n_control_max = 0
  converted_lat%n_ic_max = 0
endif

call deallocate_lat_pointers (lat_out)
call deallocate_lat_pointers (lat_model)

! Restore ptc settings

if (n_taylor_order_saved /= ptc_com%taylor_order_ptc) call set_ptc (taylor_order = n_taylor_order_saved) 
call set_ptc (exact_modeling = ptc_param%exact_model)

!------------------------------------------------------------------------
contains

subroutine write_line (line_out)

implicit none

character(*) line_out
integer ix, ix1, ix2, ix3

! Prefer to breakup a line after a comma

do
  if (len_trim(line_out) < ix_line_max) exit
  ix1 = index(line_out(ix_line_min+1:), ',')
  ix2 = index(line_out(ix_line_min+1:), '=')
  ix3 = index(line_out(ix_line_min+1:), ' ')

  if (ix1 /= 0 .and. ix1+ix_line_min < ix_line_max) then
    ix = ix1 + ix_line_min
  elseif (ix2 /= 0 .and. ix2+ix_line_min < ix_line_max) then
    ix = ix2 + ix_line_min
  elseif (ix3 /= 0 .and. ix3+ix_line_min < ix_line_max) then
    ix = ix3 + ix_line_min
  elseif (ix1 /= 0) then
    ix = ix1 + ix_line_min
  elseif (ix2 /= 0) then
    ix = ix2 + ix_line_min
  else
    ix = ix3 + ix_line_min
  endif

  write (iu, '(2a)') line_out(:ix)
  line_out = '    ' // line_out(ix+1:)
enddo

write (iu, '(2a)') trim(line_out), ';'

end subroutine write_line

end subroutine write_lat_in_sad_format


end module

