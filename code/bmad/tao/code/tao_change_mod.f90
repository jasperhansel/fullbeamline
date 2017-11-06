module tao_change_mod

use tao_mod
use tao_data_and_eval_mod
use lat_ele_loc_mod

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_change_var (name, num_str, silent)
!
! Routine to change a variable in the model lattice.
!
! Input:
!   name     -- Character(*): Name of variable or element.
!   num_str  -- Character(*): Change in value. 
!                                A '@' signifies a absolute set.
!                                A 'd' signifies a set relative design.        
!   silent   -- Logical: If True then do not print any info.
!
! Output:
!    %u(s%com%default_universe)%model -- model lattice where the variable lives.
!-

subroutine tao_change_var (name, num_str, silent)

implicit none

type (tao_universe_struct), pointer :: u
type (tao_var_array_struct), allocatable, save :: v_array(:)
type (tao_var_struct), pointer :: var

real(rp), allocatable :: change_number(:)
real(rp) model_value, old_merit, new_merit, max_val
real(rp) old_value, new_value, design_value, delta

integer nl, i, ixa, ix, err_num, n

character(*) name, num_str
character(20) :: r_name = 'tao_change_var'
character(20) abs_or_rel, component, str
character(100) l1, num, fmt
character(200), allocatable, save :: lines(:)

logical err, exists, silent

!-------------------------------------------------

call tao_find_var (err, name, v_array = v_array, component = component)
if (err) return
if (.not. allocated(v_array)) then
  call out_io (s_error$, r_name, 'BAD VARIABLE NAME: ' // name)
  return
endif
if (component /= "" .and. component /= "model") then
  call out_io (s_error$, r_name, &
            '"change var" ONLY MODIFIES THE "model" COMPONENT.', &
            'USE "set var" INSTEAD.')
  return
endif

! find change value(s)

call tao_to_change_number (num_str, size(v_array), change_number, abs_or_rel, err);  if (err) return
old_merit = tao_merit()

! We need at least one variable to exist.

exists = .false.
do i = 1, size(v_array)
  if (v_array(i)%v%exists) exists = .true.
enddo

if (.not. exists) then
  call out_io (s_error$, r_name, 'VARIABLE DOES NOT EXIST')
  return
endif

! now change all desired variables

max_val = 0
do i = 1, size(v_array)
  var => v_array(i)%v
  if (.not. var%exists) cycle
  var%old_value = var%model_value
  if (abs_or_rel == '@') then
    call tao_set_var_model_value (var, change_number(i))
  elseif (abs_or_rel == 'd') then
    call tao_set_var_model_value (var, var%design_value + change_number(i))
  elseif (abs_or_rel == '%') then
    call tao_set_var_model_value (var, var%model_value * (1 + 0.01 * change_number(i)))
  else
    call tao_set_var_model_value (var, var%model_value + change_number(i))
  endif
  max_val = max(max_val, abs(var%old_value))
  max_val = max(max_val, abs(var%design_value)) 
  max_val = max(max_val, abs(var%model_value))
enddo

! print results

new_merit = tao_merit()

str = real_num_fortran_format(max_val, 14, 2)
fmt = '(5x, i5, 2x, ' // trim(str) // ', a, 4' // trim(str) // ')'

call re_allocate (lines, 200)
nl = 0
l1 = '     Index       Old_Model         New_Model         Delta    Old-Design    New-Design'
nl=nl+1; lines(nl) = l1

n = size(v_array)
call re_allocate (lines, n+100)

do i = 1, size(v_array)
  var => v_array(i)%v
  if (.not. var%exists) cycle
  delta = var%model_value - var%old_value
  nl=nl+1; write (lines(nl), fmt) i, var%old_value, '  ->', &
       var%model_value, delta, var%old_value - var%design_value, &
	     var%model_value - var%design_value
enddo
nl=nl+1; lines(nl) = l1

str = real_num_fortran_format(max(abs(old_merit), abs(new_merit)), 14, 1)
fmt = '(5x, 2(a, ' // trim(str) // '), ' // trim(str) // ')'

nl=nl+1; lines(nl) = ' '
nl=nl+1; write (lines(nl), fmt) 'Merit:      ', old_merit, '  ->', new_merit
nl=nl+1; write (lines(nl), fmt) 'dMerit:     ', new_merit - old_merit
if (delta /= 0) then
  nl=nl+1; write (lines(nl), '(5x, a, es12.3)') 'dMerit/dVar:', (new_merit-old_merit) / delta
endif

if (.not. silent) call out_io (s_blank$, r_name, lines(1:nl))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_change_ele (ele_name, attrib_name, num_str)
!
! Routine to change a variable in the model lattice.
!
! Input:
!   ele_name    -- Character(*): Name of variable or element.
!   attrib_name -- Character(*): Attribute name of element.
!   num_str     -- Character(*): Change in value. 
!                                A '@' signifies a absolute set.
!                                A 'd' signifies a set relative design.        
!
! Output:
!    %u(s%com%default_universe)%model -- model lattice where the variable lives.
!-

subroutine tao_change_ele (ele_name, attrib_name, num_str)

implicit none

type (tao_universe_struct), pointer :: u
type (all_pointer_struct), allocatable, save :: d_ptr(:), m_ptr(:)
type (ele_pointer_struct), allocatable, save :: eles(:)
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d2_data_struct), pointer :: d2_dat

real(rp), allocatable, save :: change_number(:), old_value(:)
real(rp) new_merit, old_merit, new_value, delta, max_val

integer i, ix, iu, nl, len_name, nd
integer, parameter :: len_lines = 200

character(*) ele_name
character(*) attrib_name
character(*) num_str
character(40) e_name, a_name, fmt, name
character(20) :: r_name = 'tao_change_ele'
character(len_lines), allocatable, save :: lines(:)
character(20) abs_or_rel, str

logical err, etc_added
logical, allocatable :: this_u(:), free(:)

!-------------------------------------------------

if (s%com%common_lattice) then
  call re_allocate2 (this_u, lbound(s%u, 1), ubound(s%u, 1))
  this_u = .false.
  this_u(ix_common_uni$) = .true.
  e_name = ele_name
else
  call tao_pick_universe (ele_name, e_name, this_u, err)
  if (err) return
endif

! 

call string_trim (e_name, e_name, ix)
call str_upcase (e_name, e_name)

call string_trim (attrib_name, a_name, ix)
call str_upcase (a_name, a_name)

etc_added = .false.
nl = 0
call re_allocate (lines, 100)
nl=nl+1;write (lines(nl), '(11x, a)') &
                  'Old           New    Old-Design    New-Design         Delta'

do iu = lbound(s%u, 1), ubound(s%u, 1)

  if (.not. this_u(iu)) cycle
  u => s%u(iu)

  call pointers_to_attribute (u%design%lat, e_name, a_name, .true., d_ptr, err, .true., eles)
  if (err) return

  call pointers_to_attribute (u%model%lat, e_name, a_name, .true., m_ptr, err, .true., eles)
  if (err) return

  ! Make sure attributes are free to vary. 
  ! With something like "change ele quad::* x_offset ..." then need to ignore any super_slave elements.
  ! So things are OK if at least one free attribute exists.

  nd = size(d_ptr)
  call re_allocate (old_value, nd)
  call re_allocate (free, nd)

  free = .true.

  ! Can vary beam energy with RF off even in rings.
  ! Or when doing multi_turn_orbit data taking.

  if (e_name == 'BEAM_START') then
    if (u%model%lat%param%geometry == closed$) then
      free = .false.
      write (name, '(i0, a)') iu, '@multi_turn_orbit'
      call tao_find_data (err, name, d2_array, print_err = .false.)
      if (a_name == 'PZ' .and. .not. s%global%rf_on) free = .true.
      if (size(d2_array) > 0) free = .true.
    endif
  else
    do i = 1, nd
      free(i) = attribute_free (eles(i)%ele, a_name, .false.)
    end do
  endif

  if (all(.not. free)) then
    call out_io (s_error$, r_name, 'ATTRIBUTE NOT FREE TO VARY. NOTHING DONE')
    return
  endif

  ! Find change value(s)

  call tao_to_change_number (num_str, nd, change_number, abs_or_rel, err);  if (err) return
  old_merit = tao_merit()

  ! put in change

  do i = 1, nd

    if (.not. free(i)) cycle

    old_value(i) = m_ptr(i)%r

    if (abs_or_rel == '@') then
      m_ptr(i)%r = change_number(i)
    elseif (abs_or_rel == 'd') then
      m_ptr(i)%r = d_ptr(i)%r + change_number(i)
    elseif (abs_or_rel == '%') then
      m_ptr(i)%r = m_ptr(i)%r * (1 + 0.01 * change_number(i))
    else
      m_ptr(i)%r = m_ptr(i)%r + change_number(i)
    endif

    delta = m_ptr(i)%r - old_value(i)

    call tao_set_flags_for_changed_attribute(u, e_name, eles(i)%ele, m_ptr(i)%r)

    max_val = max(abs(old_value(i)), abs(m_ptr(i)%r), abs(d_ptr(1)%r)) 
    str = real_num_fortran_format(max_val, 14, 2)
    fmt = '(5' // trim(str) // ', 4x, a)'

    ! Record change but only for the first 10 variables.

    if (nl < 11) then
      name = 'BEAM_START'
      if (associated(eles(i)%ele)) name = eles(i)%ele%name
      nl=nl+1; write (lines(nl), fmt) old_value(i), m_ptr(i)%r, &
                              old_value(i)-d_ptr(i)%r, m_ptr(i)%r-d_ptr(i)%r, &
                              m_ptr(i)%r-old_value(i), trim(name)
    else
      if (.not. etc_added) then
        nl=nl+1; lines(nl) = '   ... etc ...'
        etc_added = .true.
      endif
    endif

  enddo

enddo

!----------------------------------
! print results

new_merit = tao_merit()

str = real_num_fortran_format(max(abs(old_merit), abs(new_merit)), 14, 1)
fmt = '(2(a, ' // trim(str) // '), a, ' // trim(str) // ')'

nl=nl+1; lines(nl) = ' '
nl=nl+1;write (lines(nl), fmt) 'Merit:      ', old_merit, '  ->', new_merit
nl=nl+1;write (lines(nl), fmt) 'dMerit:     ', new_merit - old_merit
if (delta /= 0) then
  nl=nl+1; write (lines(nl), '(a, es12.3)') 'dMerit/dValue:  ', &
                                        (new_merit-old_merit) / delta
endif

call out_io (s_blank$, r_name, lines(1:nl))

end subroutine tao_change_ele

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_to_change_number (num_str, n_size, change_number, abs_or_rel, err)

implicit none

type (tao_expression_info_struct), allocatable, save :: info(:)

real(rp), allocatable :: change_number(:)
integer ix, ios, n_size

character(*) num_str
character(*) abs_or_rel
character(len(num_str)) number_str
character(*), parameter :: r_name = 'tao_to_change_number'

logical err

!

call string_trim(num_str, number_str, ix)
abs_or_rel = ' ' 

select case (number_str(1:1))
case ('@', 'd', '%')
  abs_or_rel = number_str(1:1)
  number_str(1:1) = ' '
end select

call tao_evaluate_expression (number_str, n_size, .false., change_number, info, err)

end subroutine tao_to_change_number

end module
