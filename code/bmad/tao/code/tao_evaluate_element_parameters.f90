!+
! Subroutine tao_evaluate_element_parameters (err, param_name, values, print_err, dflt_source, dflt_component, dflt_uni)
!
! Routine to evaluate a lattice element parameter of the form 
!     <universe>@ele::{<class>}::<ele_name_or_num>[<parameter>]{|<component>}
! or to evaluate at the middle of the element
!     <universe>@ele_mid::{<class>}::<ele_name_or_num>[<parameter>]{|<component>}
! Note: size(values) can be zero without an error
! 
! Input:
!   param_name      -- character(*): parameter name.
!   print_err       -- logical: Print error message? 
!   dflt_source     -- character(*): Default source
!   dflt_component  -- character(*), optional: Default component
!   dflt_uni        -- integer, optional :: Default universe to use.
!
! Output:
!   err       -- Logical: True if there is an error in syntax. False otherwise
!   values(:) -- Real(rp), allocatable: Array of datum valuse.
!-

subroutine tao_evaluate_element_parameters (err, param_name, values, print_err, dflt_source, dflt_component, dflt_uni)

use tao_utils, except_dummy => tao_evaluate_element_parameters

implicit none

type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (ele_struct) ele3
type (coord_struct), pointer :: this_orb(:)
type (coord_struct) orb
type (branch_struct), pointer :: branch
type (all_pointer_struct) a_ptr

character(*) param_name
character(*) dflt_source
character(*), optional :: dflt_component
character(60) name, class_ele, parameter, component
character(*), parameter :: r_name = 'tao_evaluate_element_parameters'

real(rp), allocatable :: values(:)
real(rp) :: real_val

integer, optional :: dflt_uni
integer i, j, ix, num, ixe, ix1, ios, n_tot

logical err, valid, middle
logical :: print_err

!

call tao_pick_universe (param_name, name, scratch%this_u, err, dflt_uni = dflt_uni)
if (err) return

err = .true.

if (name(1:5) == 'ele::') then
  name = name(6:)  ! Strip off 'ele::'
  middle = .false.
elseif (name(1:9) == 'ele_mid::') then   
  name = name(10:)  ! Strip off 'ele_mid::'
  middle = .true.
elseif (dflt_source /= 'element') then
  return
endif

! Get component

ix = index(name, '|')
if (ix == 0) then
  component = 'model'
  if (present(dflt_component)) then
    if (dflt_component /= '') component = dflt_component
  endif
else
  component = name(ix+1:)
  name = name(1:ix-1)
endif

! Get class:name

ix1 = index(name, '[');  
if (ix1 == 0) return

ix1 = index(name, '[');  if (ix1 == 0) return
class_ele = name(1:ix1-1)
name = name(ix1+1:)
if (class_ele(1:2) == '::') class_ele = class_ele(3:)
ix1 = index(name, ']');  if (ix1 == 0) return
parameter = name(1:ix1-1)

select case (parameter)
case ('l', 'angle');    middle = .false.
end select

! Evaluate

n_tot = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. scratch%this_u(i)) cycle
  u => s%u(i)
  call tao_locate_elements (class_ele, u%ix_uni, scratch%eles, err)
  if (err) return
  call re_allocate (values, n_tot + size(scratch%eles))

  do j = 1, size(scratch%eles)
    ixe = scratch%eles(j)%ele%ix_ele

    if (parameter == 'index') then
      values(n_tot+j) = ixe
      cycle
    endif

    select case (component)
    case ('model')   
      lat => u%model%lat
      this_orb => u%model%tao_branch(scratch%eles(j)%ele%ix_branch)%orbit
      branch => u%model%lat%branch(scratch%eles(j)%ele%ix_branch)
    case ('base')  
      lat => u%base%lat
      this_orb => u%base%tao_branch(scratch%eles(j)%ele%ix_branch)%orbit
      branch => u%base%lat%branch(scratch%eles(j)%ele%ix_branch)
    case ('design')
      lat => u%design%lat
      this_orb => u%design%tao_branch(scratch%eles(j)%ele%ix_branch)%orbit
      branch => u%design%lat%branch(scratch%eles(j)%ele%ix_branch)
    case default
      call out_io (s_error$, r_name, 'BAD DATUM COMPONENT FOR: ' // param_name)
      return
    end select

    if (middle .and. ixe /= 0) then
      call twiss_and_track_intra_ele (branch%ele(ixe), lat%param, 0.0_rp, branch%ele(ixe)%value(l$)/2, &
                .true., .false., this_orb(ixe-1), orb, branch%ele(ixe-1), ele3, err)
      call tao_orbit_value (parameter, orb, values(n_tot+j), err)
      if (err) then
        call pointer_to_attribute (ele3, parameter, .true., a_ptr, err, print_err)
        if (err) return
        values(n_tot+j) = value_of_all_ptr(a_ptr)
      endif

    else
      call tao_orbit_value (parameter, this_orb(ixe), values(n_tot+j), err)
      if (err) then
        call pointer_to_attribute (branch%ele(ixe), parameter, .true., a_ptr, err, print_err)
        if (err) return
        values(n_tot+j) = value_of_all_ptr(a_ptr)
      endif
    endif

  enddo

  n_tot = n_tot + size(values)
enddo

err = .false.

end subroutine tao_evaluate_element_parameters
