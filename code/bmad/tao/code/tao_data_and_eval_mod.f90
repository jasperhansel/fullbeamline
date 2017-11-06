module tao_data_and_eval_mod

use tao_mod
use spin_mod
use utilities_mod
use measurement_mod
use geometry_mod

implicit none

! Used for parsing expressions

integer, parameter :: var_num$ = 101, lat_num$ = 102, data_num$ = 103, ele_num$ = 104

private tao_scratch_values_calc

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_lat_or_beam_data (err, data_name, values, print_err,
!                    default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_component, dflt_uni)
!
! Routine to evaluate data with a lat or beam source of the form:
!     <universe>@lat::<data_type>[<ix_ele_start>&<ix_ele>]|<component>
!
! Input:
!   data_name      -- character(*): data name.
!   print_err      -- logical: Print error message?
!   dflt_source    -- character(*): If not blank: Default source: 'lat' or 'beam'.
!   dflt_ele_ref   -- ele_struct, pointer, optional: Default reference element.
!   dflt_ele_start -- ele_struct, pointer, optional: Default start element.
!   dflt_ele       -- ele_struct, pointer, optional: Default element to evaluate at.
!   dflt_component -- character(*), optional: Default component: 'model', 'base', or 'design'.
!   dflt_uni       -- integer, optional: Default universe to use
!
! Output:
!   err       -- Logical: True if there is an error. False otherwise
!   values(:) -- Real(rp), allocatable: Array of datum valuse.
!-

subroutine tao_evaluate_lat_or_beam_data (err, data_name, values, print_err, &
                         default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_component, dflt_uni)

type (tao_data_struct) datum
type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable, save :: eles(:)
type (ele_struct), pointer, optional :: dflt_ele_ref, dflt_ele_start, dflt_ele

character(*) data_name
character(*) default_source
character(*), optional :: dflt_component
character(60) name, ele_name, component
character(*), parameter :: r_name = 'tao_evaluate_lat_or_beam_data'

real(rp), allocatable :: values(:)

integer, optional :: dflt_uni
integer i, j, num, ix, ix1, ios, n_tot, n_loc

logical err, valid, err_flag
logical print_err, use_dflt_ele
logical, allocatable, save :: this_u(:)

!

datum%exists = .true.
datum%ele_start_name = ''
datum%ele_ref_name = ''
datum%ix_ele_start = -1
datum%ix_ele_ref = -1

call tao_pick_universe (data_name, name, this_u, err, dflt_uni = dflt_uni)
if (err) return

err = .true.
if (name(1:5) == 'lat::') then
  datum%data_source = 'lat'
  name = name(6:)  ! Strip off 'lat:'
elseif (name(1:6) == 'beam::') then
  datum%data_source = 'beam'
  name = name(7:)  ! Strip off 'beam:'
elseif (default_source /= '') then
  datum%data_source = default_source
else
  if (print_err) call out_io (s_error$, r_name, 'DATUM NOT "LAT::" OR "BEAM::"' // data_name)
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

! Get data type

ix1 = index(name, '[');
if (ix1 == 0) then
  if (.not. present(dflt_ele) .or. .not. associated(dflt_ele)) then
    if (print_err) call out_io (s_error$, r_name, 'NO "[" FOUND IN:' // data_name)
    return
  endif
  datum%data_type = name
  name = dflt_ele%name
  use_dflt_ele = .true.

else
  datum%data_type = name(1:ix1-1)
  name = name(ix1+1:)
  ix1 = index(name, ']')
  if (ix1 == 0) then
    if (print_err) call out_io (s_error$, r_name, 'NO "]" FOUND IN:' // data_name)
    return
  endif
  name(ix1:ix1) = ''
  if (name(ix1+1:) /= '') then
    if (print_err) call out_io (s_error$, r_name, 'MANGLED CONSTRUCT:' // data_name)
    return
  endif
  use_dflt_ele = .false.

  ! Get ele_ref & ele

  ix = index(name, '&')
  if (ix /= 0) then
    if (is_integer(name)) then
      read (name(:ix-1), *, iostat = ios) datum%ix_ele_ref
      if (ios /= 0) then
        if (print_err) call out_io (s_error$, r_name, 'BAD ELE_REF: ' // data_name)
        return
      endif
    else
      datum%ele_ref_name = name(:ix-1)
    endif
    ele_name = name(ix+1:)
  else
    ele_name = name
  endif
endif

! Evaluate

n_tot = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. this_u(i)) cycle
  u => s%u(i)

  if (use_dflt_ele) then
    if (associated(dflt_ele_start)) then
      n_loc = dflt_ele%ix_ele - dflt_ele_start%ix_ele + 1
    else
      n_loc = 1
    endif
    if (associated(dflt_ele_ref)) then
      datum%ix_ele_ref = dflt_ele_ref%ix_ele
    endif

  else
    if (datum%ele_ref_name /= '') then
      call lat_ele_locator (datum%ele_ref_name, u%model%lat, eles, n_loc, err_flag)
      if (err_flag) return
      if (n_loc /= 1) then
        if (print_err) call out_io (s_error$, r_name, &
                        'MULTIPLE ELEMENTS MATCH REFERENCE NAME: ' // datum%ele_ref_name)
        return
      endif
      datum%ix_ele_ref = eles(1)%ele%ix_ele
    endif

    call lat_ele_locator (ele_name, u%model%lat, eles, n_loc, err_flag)
    if (err_flag) return
  endif

  call re_allocate (values, n_tot + n_loc)

  do j = 1, n_loc
    if (use_dflt_ele) then
      datum%ix_branch = dflt_ele%ix_branch
      if (associated(dflt_ele_start)) then
        datum%ix_ele = dflt_ele_start%ix_ele + j - 1
      else
        datum%ix_ele = dflt_ele%ix_ele
      endif
    else
      datum%ele_name = eles(j)%ele%name
      datum%ix_ele = eles(j)%ele%ix_ele
      datum%ix_branch = eles(j)%ele%ix_branch
    endif

    select case (component)
    case ('model')   
      call tao_evaluate_a_datum (datum, u, u%model, values(n_tot+j), valid)
    case ('base')  
      call tao_evaluate_a_datum (datum, u, u%base, values(n_tot+j), valid)
    case ('design')  
      call tao_evaluate_a_datum (datum, u, u%design, values(n_tot+j), valid)
    case default
      if (print_err) call out_io (s_error$, r_name, 'BAD DATUM COMPONENT FOR: ' // data_name)
      return
    end select
  enddo
  n_tot = n_tot + size(values)
enddo

if (n_tot == 0) then
  if (print_err) call out_io (s_error$, r_name, 'ELEMENT NOT FOUND: ' // ele_name)
  return
endif

err = .false.

end subroutine tao_evaluate_lat_or_beam_data

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
!
! Buffer routine for to_phase_and_coupling_reading.
!
! Input:
!   ele -- Ele_struct: The monitor.
!
! Output:
!   bpm_data     -- Bpm_phase_coupling_struct: Monitor values
!   valid_value  -- Logical: Valid data value?
!-

subroutine tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)

type (ele_struct) ele
type (bpm_phase_coupling_struct) bpm_data
type (bpm_phase_coupling_struct), save :: old_bpm_data

integer, save :: ix_ele_old = -1

logical valid_value
logical, save :: err

!

if (ix_ele_old /= ele%ix_ele) then
  call to_phase_and_coupling_reading (ele, old_bpm_data, err)
  ix_ele_old = ele%ix_ele
endif

bpm_data = old_bpm_data
valid_value = .not. err

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_get_data (data_value, data_weight, data_meas_value, dat_ix_dModel)
!
! Subroutine to get the values of the data used in optimization and put them
! in an array. The data is ordered starting with the first universe
!
! Input:
! 
! Output:
!   data_value(:)       -- Real(Rp), allocatable, optional: Data model values.
!   data_weight(:)      -- Real(Rp), allocatable, optional: Data weights in the merit function.
!   data_meas_value(:)  -- Real(Rp), allocatable, optional: Data values when the data was taken.
!   data_ix_dModel(:)   -- Integer, allocatable, optional: Data ix_dModel indices
!-

subroutine tao_get_data (data_value, data_weight, data_meas_value, data_ix_dModel)

real(rp), allocatable, optional :: data_value(:), data_meas_value(:), data_weight(:)
integer, allocatable, optional :: data_ix_dModel(:)

integer i, j, iu
integer n_data

!
  
n_data = 0
do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(iu)%is_on) cycle
  n_data  = n_data + count(s%u(iu)%data(:)%useit_opt)
enddo
if (present(data_value))      call re_allocate (data_value, n_data)
if (present(data_meas_value)) call re_allocate (data_meas_value, n_data)
if (present(data_weight))     call re_allocate (data_weight, n_data)
if (present(data_ix_dModel))  call re_allocate (data_ix_dModel, n_data)

j = 0
do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(iu)%is_on) cycle
  do i = 1, size(s%u(iu)%data)
    if (.not. s%u(iu)%data(i)%useit_opt) cycle
    j = j + 1
    if (present(data_value))        data_value(j)      = s%u(iu)%data(i)%model_value
    if (present(data_weight))       data_weight(j)     = s%u(iu)%data(i)%weight
    if (present(data_meas_value))   data_meas_value(j) = s%u(iu)%data(i)%meas_value
    if (present(data_ix_dModel))    data_ix_dModel(j)  = s%u(iu)%data(i)%ix_dModel
  enddo
enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_data_coupling_init (u)
!
! Routine to initialize the coupling structure for a lattice branch.
! This routine is called by tao_lattic_calc and is not meant for general use.
!
! Input:
!   branch -- branch_struct: New lattice branch.
!-

subroutine tao_data_coupling_init (branch)

type (branch_struct) branch
integer m

! 

m = branch%n_ele_max
if (.not. allocated(scratch%cc)) allocate (scratch%cc(0:m))
if (ubound(scratch%cc, 1) < m) then
  deallocate(scratch%cc)
  allocate(scratch%cc(0:m))
endif

scratch%cc%coupling_calc_done = .false.
scratch%cc%amp_calc_done = .false.

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, why_invalid)
!
! Subroutine to put the proper data in the specified datum
!
! Input:
!   datum         -- Tao_data_struct: What type of datum
!   u             -- Tao_universe_struct: Which universe to use.
!   tao_lat       -- Tao_lattice_struct: Lattice to use.
!     
! Output:
!   datum          -- Tao_data_struct: 
!     %ix_ele_merit   -- For max/min type constraints: Place where value is max/min. 
!   datum_value   -- Real(rp): Value of the datum.
!   valid_value   -- Logical: Set false when there is a problem. Set true otherwise.
!   why_invalid   -- Character(*), optional: Tells why datum value is invalid.
!-

recursive subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, why_invalid)

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_lattice_branch_struct), pointer :: tao_branch
type (tao_building_wall_section_struct), pointer :: section
type (tao_building_wall_point_struct), pointer :: pt
type (tao_data_struct), pointer :: dp
type (tao_data_array_struct), allocatable, save :: d_array(:)
type (tao_lattice_struct), target :: tao_lat
type (tao_expression_info_struct), allocatable :: info(:)
type (lat_struct), pointer :: lat
type (normal_modes_struct) mode
type (spin_polar_struct) polar
type (ele_struct), pointer :: ele, ele_start, ele_ref, ele2
type (ele_struct) ele_at_s
type (coord_struct), pointer :: orb0, orbit(:), orb
type (coord_struct) :: orb_at_s
type (bpm_phase_coupling_struct) bpm_data
type (taylor_struct), save :: taylor_save(6), taylor(6) ! Saved taylor map
type (floor_position_struct) floor
type (branch_struct), pointer :: branch
type (bunch_params_struct), pointer :: bunch_params(:)
type (normal_form_struct), pointer :: normal_form
type (taylor_struct), pointer :: taylor_ptr
type (complex_taylor_struct), pointer :: complex_taylor_ptr
type (all_pointer_struct) a_ptr
type (spin_polar_struct) polar_spin

real(rp) datum_value, mat6(6,6), vec0(6), angle, px, py, vec2(2)
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), a_vec(4), mc2
real(rp) gamma, one_pz, w0_mat(3,3), w_mat(3,3), vec3(3), value, s_len
real(rp) dz, dx, cos_theta, sin_theta, z_pt, x_pt, z0_pt, x0_pt
real(rp) z_center, x_center, x_wall, s_eval, s_eval_ref, phase, amp
real(rp), allocatable, save :: value_vec(:)
real(rp), allocatable, save :: expression_value_vec(:)
real(rp) theta, phi, psi

! Cf: Sands Eq 5.46 pg 124.
real(rp), parameter :: const_q_factor = 55 * h_bar_planck * c_light / (32 * sqrt_3) 

integer i, j, k, m, n, k_old, ix, ie, is, iz, ix_ele, ix_start, ix_ref
integer n_size, ix0, which, expnt(6), n_track, n_max

character(*), optional :: why_invalid
character(6) expn_str
character(16) constraint
character(20) :: r_name = 'tao_evaluate_a_datum'
character(40) head_data_type, sub_data_type, data_source, name, dflt_dat_index
character(300) str

logical found, valid_value, err, taylor_is_complex, use_real_part, compute_floor
logical, allocatable, save :: good(:)

! If does not exist

if (.not. datum%exists) then
  datum_value = real_garbage$
  valid_value = .false.
  if (present(why_invalid)) why_invalid = 'Datum does not exist.'
  return
endif

! To save time, don't evaluate if unnecessary when the running an optimizer.
! Exception: When there are datums that use expressions, things are 
!   complicated so don't try to save time in this case.

if (s%com%optimizer_running .and. .not. datum%useit_opt .and. .not. s%com%have_datums_using_expressions) then
  datum_value = 0
  valid_value = .false.
  return
endif

! See if there is a hook for this datum

call tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value, why_invalid)
if (found) return

! Set ix_ele, ix_ref, and ix_start. 
! Note: To check that a start element was set, need to look at datum%ele_start_name, not ix_start.

data_source = datum%data_source
head_data_type = datum%data_type
lat => tao_lat%lat

if (head_data_type == 'null') then
  datum_value = 0
  valid_value = .false.
  return
endif

ele => tao_pointer_to_datum_ele (lat, datum%ele_name, datum%ix_ele, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_ele = -1
if (associated(ele)) ix_ele = ele%ix_ele

ele_ref => tao_pointer_to_datum_ele (lat, datum%ele_ref_name, datum%ix_ele_ref, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_ref = -1
if (associated(ele_ref)) ix_ref = ele_ref%ix_ele

ele_start => tao_pointer_to_datum_ele (lat, datum%ele_start_name, datum%ix_ele_start, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_start = ix_ele
if (associated(ele_start)) ix_start = ele_start%ix_ele

! Some inits

valid_value = .false.

datum_value = 0           ! default
datum%ix_ele_merit = -1   ! default

branch => lat%branch(datum%ix_branch)
tao_branch => tao_lat%tao_branch(datum%ix_branch)
orbit => tao_branch%orbit
bunch_params => tao_branch%bunch_params

n_track = branch%n_ele_track
n_max   = branch%n_ele_max

call re_allocate2 (value_vec, 0, n_track, .false.) ! Make sure is large enough if used.
call re_allocate2 (good,      0, n_track, .false.) ! Make sure is large enough if used.

ix = index(head_data_type, '.')
if (head_data_type(1:11) == 'expression:') then
  head_data_type = 'expression:'
elseif (ix /= 0) then
  sub_data_type  = head_data_type(ix+1:)
  head_data_type = head_data_type(1:ix) 
endif

if (head_data_type  == 'rad_int.' .or. head_data_type == 'rad_int1.') then
  if (index(head_data_type, '_e') /= 0 .and. (ix_ref > -1 .or. ix_ele > -1)) then
    if (.not. allocated(tao_branch%rad_int%ele)) then
      call out_io (s_error$, r_name, 'tao_branch%rad_int not allocated')
      return
    endif
  endif
endif

select case (data_source)
case ('lat', 'beam')
  ! Valid data source
case default
  if ( head_data_type /= 'expression:') then
    call tao_set_invalid (datum, 'UNKNOWN DATA_SOURCE: ' // data_source, why_invalid)
    return
  endif
end select

! ele_ref must not be specified for some data types. Check this.

select case (head_data_type)
case ('wall.')
  if (datum%ele_ref_name /= '') then
    call tao_set_invalid (datum, 'SPECIFYING ELE_REF NOT VALID', why_invalid)
    return
 endif
end select


! ele_start must not be specified for some data types. Check this.

if (datum%data_type(1:11) == 'periodic.tt' .or. datum%data_type == 'sigma.pz') then
  if (datum%ele_start_name /= '') then
    call tao_set_invalid (datum, 'SPECIFYING ELE_START NOT VALID', why_invalid)
    return
  endif
endif

if (datum%data_type(1:11) == 'periodic.tt') then
  if (datum%ele_ref_name /= '') then
    call tao_set_invalid (datum, 'SPECIFYING ELE_REF NOT VALID', why_invalid)
    return
  endif
endif

if (tao_branch%track_state /= moving_forward$ .and. ix_ele >= tao_branch%track_state) then
  if ((data_source == 'beam' .and. head_data_type /= 'n_particle_loss') .or. &
                         head_data_type(1:4) == 'bpm_' .or. head_data_type == 'orbit.') then
    if (present(why_invalid)) why_invalid = 'CANNOT EVALUATE DUE TO PARTICLE LOSS.'
    return
  endif
endif

!---------------------------------------------------

if (datum%s_offset /= 0 .or. datum%eval_point == anchor_center$) then
  if (data_source /= 'lat') then
    call tao_set_invalid (datum, 'CANNOT USE A BEAM DATA_SOURCE WITH A FINITE S_OFFSET OR EVAL_POINT = CENTER.', why_invalid)
    return
  endif

  if (.not. associated(ele)) then
    call tao_set_invalid (datum, 'THERE MUST BE AN ASSOCIATED ELEMENT WHEN S_OFFSET IS NON-ZERO OR EVAL_POINT = CENTER.', why_invalid)
    return
  endif

  select case (datum%eval_point)
  case (anchor_beginning$)
    ! Here tao_pointer_to_datum_ele has pointed ele to the element before the element specified in the datum.
    s_eval = ele%s + datum%s_offset
    if (associated(ele_ref)) s_eval_ref = ele%s
  case (anchor_center$)
    s_eval = (ele%s_start + ele%s) / 2 + datum%s_offset
    if (associated(ele_ref)) s_eval_ref = (ele%s_start + ele%s) / 2
  case (anchor_end$)
    s_eval = ele%s + datum%s_offset
    if (associated(ele_ref)) s_eval_ref = ele%s
  end select

  compute_floor = (head_data_type == 'floor.')

  call twiss_and_track_at_s (branch%lat, s_eval, ele_at_s, orbit, orb_at_s, branch%ix_branch, &
                                                                      err, compute_floor_coords = compute_floor)
  if (err) then
    call tao_set_invalid (datum, 'CANNOT TRACK TO OFFSET POSITION.', why_invalid)
    return
  endif

  datum_value = tao_bmad_parameter_value (datum%data_type, ele_at_s, orb_at_s, err)
  if (err) then
    call tao_set_invalid (datum, 'CANNOT EVALUATE DATUM AT OFFSET POSITION.', why_invalid)
    return
  endif

  if (associated(ele_ref)) then
    call twiss_and_track_at_s (branch%lat, s_eval_ref, ele_at_s, orbit, orb_at_s, branch%ix_branch, &
                                                                      err, compute_floor_coords = compute_floor)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT TRACK TO REFERENCE POSITION.', why_invalid)
      return
    endif

    datum_value = datum_value - tao_bmad_parameter_value (datum%data_type, ele_at_s, orb_at_s, err)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT EVALUATE DATUM AT REFERENCE POSITION.', why_invalid)
      return
    endif
  endif
  
  valid_value = .true.
  return
endif

!---------------------------------------------------

select case (head_data_type)

!-----------

case ('alpha.')

  select case (datum%data_type)

  case ('alpha.a')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%a%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%a%norm_emit /= 0)
      if (.not. valid_value) call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif
    
  case ('alpha.b')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%b%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%b%norm_emit /= 0)
      if (.not. valid_value) call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif

  case ('alpha.z')
    if (data_source == 'lat') return
    call tao_load_this_datum (bunch_params(:)%z%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('apparent_emit.', 'norm_apparent_emit.')

  select case (datum%data_type)

  case ('apparent_emit.x', 'norm_apparent_emit.x')
    do i = ix_start, ix_ele
      if (data_source == 'lat') then
        value_vec(i) = tao_lat_emit_calc (x_plane$, apparent_emit$, branch%ele(i), tao_branch%modes)
      else
        value_vec(i) = tao_beam_emit_calc (x_plane$, apparent_emit$, branch%ele(i), bunch_params(i))
      endif
    enddo

    if (ix_ref > -1) then
      if (data_source == 'lat') then
        value_vec(ix_ref) = tao_lat_emit_calc (x_plane$, apparent_emit$, branch%ele(ix_ref), tao_branch%modes)
      else
        value_vec(ix_ref) = tao_beam_emit_calc (x_plane$, apparent_emit$, branch%ele(i), bunch_params(ix_ref))
      endif
    endif

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

    if (datum%data_type == 'norm_apparent_emit.x') then
      call convert_total_energy_to (ele%value(E_tot$), lat%param%particle, gamma)
      datum_value = datum_value * gamma
    endif

  case ('apparent_emit.y', 'norm_apparent_emit.y')
    do i = ix_start, ix_ele
      if (data_source == 'lat') then
        value_vec(i) = tao_lat_emit_calc (y_plane$, apparent_emit$, branch%ele(i), tao_branch%modes)
      else
        value_vec(i) = tao_beam_emit_calc (y_plane$, apparent_emit$, branch%ele(i), bunch_params(i))
      endif
    enddo

    if (ix_ref > -1) then
      if (data_source == 'lat') then
        value_vec(ix_ref) = tao_lat_emit_calc (y_plane$, apparent_emit$, branch%ele(ix_ref), tao_branch%modes)
      else
        value_vec(ix_ref) = tao_beam_emit_calc (y_plane$, apparent_emit$, branch%ele(i), bunch_params(ix_ref))
      endif
    endif

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

    if (datum%data_type == 'norm_apparent_emit.y') then
      call convert_total_energy_to (ele%value(E_tot$), lat%param%particle, gamma)
      datum_value = datum_value * gamma
    endif


  case default
    call tao_set_invalid (datum, 'UNKNOWN DATUM TYPE: ' // datum%data_type, why_invalid)
    return

  end select

!-----------

case ('beta.')

  select case (datum%data_type)

  case ('beta.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      if (bunch_params(ix_ele)%x%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    else
      call tao_set_invalid (datum, 'DATA_SOURCE: ' // trim(data_source) // ' INVALID WITH: beta.x DATA_TYPE', why_invalid)
    endif
    
  case ('beta.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      if (bunch_params(ix_ele)%y%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    else
      call tao_set_invalid (datum, 'DATA_SOURCE: ' // trim(data_source) // ' INVALID WITH: beta.y DATA_TYPE', why_invalid)
    endif

  case ('beta.z')
    if (data_source == 'lat') return
    call tao_load_this_datum (bunch_params(:)%z%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    if (bunch_params(ix_ele)%z%norm_emit == 0) then
      valid_value = .false.
      call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif

  case ('beta.a')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%a%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      if (bunch_params(ix_ele)%a%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    endif
    
  case ('beta.b')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%b%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      if (bunch_params(ix_ele)%b%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('bpm_orbit.')

  select case (datum%data_type)
  case ('bpm_orbit.x')
    which = x_plane$
  case ('bpm_orbit.y')
    which = y_plane$
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return
  end select

  if (data_source == 'beam') return ! bad
  call to_orbit_reading (orbit(ix_ele), ele, which, datum_value, err)
  valid_value = .not. (err .or. (tao_branch%track_state /= moving_forward$ .and. ix_ele > tao_branch%track_state))

!-----------

case ('bpm_eta.')

  select case (datum%data_type)
  case ('bpm_eta.x')
    which = x_plane$
  case ('bpm_eta.y')
    which = y_plane$
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return
  end select

  if (data_source == 'beam') return ! bad
  vec2 = [ele%x%eta, ele%y%eta]
  call to_eta_reading (vec2, ele, which, datum_value, err)
  valid_value = .not. err

!-----------

case ('bpm_phase.')

  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)

  select case (datum%data_type)
  case ('bpm_phase.a')
    datum_value = bpm_data%phi_a
  case ('bpm_phase.b')
    datum_value = bpm_data%phi_b
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    valid_value = .false.
    return
  end select

!-----------

case ('bpm_k.')

  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)

  select case (datum%data_type)
  case ('bpm_k.22a')
    datum_value = bpm_data%k_22a
  case ('bpm_k.12a')
    datum_value = bpm_data%k_12a
  case ('bpm_k.11b')
    datum_value = bpm_data%k_11b
  case ('bpm_k.12b')
    datum_value = bpm_data%k_12b
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    valid_value = .false.
    return
  end select

!-----------

case ('bpm_cbar.')

  if (data_source == 'beam') return ! bad
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)

  select case (datum%data_type)
  case ('bpm_cbar.22a')
    datum_value = bpm_data%cbar22_a
  case ('bpm_cbar.12a')
    datum_value = bpm_data%cbar12_a
  case ('bpm_cbar.11b')
    datum_value = bpm_data%cbar11_b
  case ('bpm_cbar.12b')
    datum_value = bpm_data%cbar12_b
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    valid_value = .false.
    return
  end select

!-----------
case ('bunch_max.', 'bunch_min.')
  if (data_source /= 'beam') return ! bad
  select case (datum%data_type(11:))
  case ('x'); i=1
  case ('px');i=2
  case ('y'); i=3
  case ('py');i=4
  case ('z'); i=5
  case ('pz');i=6
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return
  end select
  
  select case (datum%data_type(1:10))
  case ('bunch_max.')
    call tao_load_this_datum (bunch_params(:)%rel_max(i), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  case ('bunch_min.')
    call tao_load_this_datum (bunch_params(:)%rel_min(i), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  end select

!-----------

case ('c_mat.')

  select case (datum%data_type)

  case ('c_mat.11')
    if (data_source == 'beam') return
    call tao_load_this_datum (branch%ele(:)%c_mat(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case ('c_mat.12')
    if (data_source == 'beam') return
    call tao_load_this_datum (branch%ele(:)%c_mat(1,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case ('c_mat.21')
    if (data_source == 'beam') return
    call tao_load_this_datum (branch%ele(:)%c_mat(2,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case ('c_mat.22')
    if (data_source == 'beam') return
    call tao_load_this_datum (branch%ele(:)%c_mat(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('cbar.')

  select case (datum%data_type)

  case ('cbar.11')
    if (data_source == 'beam') return
    call tao_load_this_datum (scratch%cc%cbar(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case ('cbar.12')
    if (data_source == 'beam') return
    call tao_load_this_datum (scratch%cc%cbar(1,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case ('cbar.21')
    if (data_source == 'beam') return
    call tao_load_this_datum (scratch%cc%cbar(2,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case ('cbar.22')
    if (data_source == 'beam') return
    call tao_load_this_datum (scratch%cc%cbar(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('chrom.')
  
  if (data_source == 'beam') return

  if (.not. associated(tao_branch%low_E_lat%ele)) then
    if (branch%param%unstable_factor == 0) then
      why_invalid = 'Chrom bookkeeping problem. Please contact DCS.'
    else
      why_invalid = 'Unstable lattice.'
    endif
    return
  endif

  !----

  select case (datum%data_type)

  case ('chrom.dtune.a', 'chrom.a')
    datum_value = tao_branch%a%chrom
    valid_value = .true.

  case ('chrom.dtune.b', 'chrom.b')
    datum_value = tao_branch%b%chrom
    valid_value = .true.

   
  case ('chrom.dbeta.a')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_branch%high_E_lat%ele(i)%a%beta - tao_branch%low_E_lat%ele(i)%a%beta) / &
                        (tao_lat%lat%ele(i)%a%beta * s%global%delta_e_chrom)
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.dbeta.b')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_branch%high_E_lat%ele(i)%b%beta - tao_branch%low_E_lat%ele(i)%b%beta) / &
                        (tao_lat%lat%ele(i)%b%beta * s%global%delta_e_chrom)
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  
  case ('chrom.dphi.a')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_branch%high_E_lat%ele(i)%a%phi - tao_branch%low_E_lat%ele(i)%a%phi)/ s%global%delta_e_chrom
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.dphi.b')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_branch%high_E_lat%ele(i)%b%phi - tao_branch%low_E_lat%ele(i)%b%phi)/ s%global%delta_e_chrom
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.deta.x')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_branch%high_E_lat%ele(i)%x%eta - tao_branch%low_E_lat%ele(i)%x%eta)/ s%global%delta_e_chrom
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.deta.y')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_branch%high_E_lat%ele(i)%y%eta - tao_branch%low_E_lat%ele(i)%y%eta)/ s%global%delta_e_chrom
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.detap.x')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_branch%high_E_lat%ele(i)%x%etap - tao_branch%low_E_lat%ele(i)%x%etap)/ s%global%delta_e_chrom
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.detap.y')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_branch%high_E_lat%ele(i)%y%etap - tao_branch%low_E_lat%ele(i)%y%etap)/ s%global%delta_e_chrom
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('damp.')

  if (data_source == 'beam') return ! bad

  select case (datum%data_type)

  case ('damp.j_a')
    datum_value = tao_branch%modes%a%j_damp
    valid_value = .true.

  case ('damp.j_b')
    datum_value = tao_branch%modes%b%j_damp
    valid_value = .true.

  case ('damp.j_z')
    datum_value = tao_branch%modes%z%j_damp
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('dpx_dx') 
  if (data_source == 'lat') then
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) = tao_branch%linear(ix_ref)%sigma(1,2) / tao_branch%linear(ix_ref)%sigma(1,1)
      value_vec(ix_ele) = tao_branch%linear(ix_ele)%sigma(1,2) / tao_branch%linear(ix_ele)%sigma(1,1)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (tao_branch%linear%sigma(1,2) / tao_branch%linear%sigma(1,1), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  else
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) = bunch_params(ix_ref)%sigma(1,2) / bunch_params(ix_ref)%sigma(1,1)
      value_vec(ix_ele) = bunch_params(ix_ele)%sigma(1,2) / bunch_params(ix_ele)%sigma(1,1)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params%sigma(1,2) / bunch_params%sigma(1,1), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  endif

case ('dpy_dy') 
  if (data_source == 'lat') then
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  tao_branch%linear(ix_ref)%sigma(3,4) / tao_branch%linear(ix_ref)%sigma(3,3)
      value_vec(ix_ele) = tao_branch%linear(ix_ele)%sigma(3,4) / tao_branch%linear(ix_ele)%sigma(3,3)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (tao_branch%linear%sigma(3,4) / tao_branch%linear%sigma(3,3), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  else
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  bunch_params(ix_ref)%sigma(3,4) / bunch_params(ix_ref)%sigma(3,3)
      value_vec(ix_ele) = bunch_params(ix_ele)%sigma(3,4) / bunch_params(ix_ele)%sigma(3,3)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params%sigma(3,4) / bunch_params%sigma(3,3), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  endif

case ('dpz_dz') 
  if (data_source == 'lat') then
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  tao_branch%linear(ix_ref)%sigma(5,6) / tao_branch%linear(ix_ref)%sigma(5,5)
      value_vec(ix_ele) = tao_branch%linear(ix_ele)%sigma(5,6) / tao_branch%linear(ix_ele)%sigma(5,5)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (tao_branch%linear%sigma(5,6) / tao_branch%linear%sigma(5,5), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  else
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  bunch_params(ix_ref)%sigma(5,6) / bunch_params(ix_ref)%sigma(5,5)
      value_vec(ix_ele) = bunch_params(ix_ele)%sigma(5,6) / bunch_params(ix_ele)%sigma(5,5)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params%sigma(5,6) / bunch_params%sigma(5,5), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  endif

case ('e_tot')

  call out_io (s_warn$, r_name, '"e_tot" has been renamed to "orbit.e_tot" to avoid confusion with lattice element referece energy parameter.', &
                'Please modify your input file appropriately.')

  if (ix_ref > -1) then
    if (data_source == 'beam') then
      orb => bunch_params(ix_ref)%centroid
    else
      orb => orbit(ix_ref)
    endif
    if (orb%state == not_set$) return
    call convert_pc_to ((1 + orb%vec(6))*orb%p0c, orb%species, e_tot = value_vec(ix_ref))
  endif

  do i = ix_start, ix_ele
    if (data_source == 'beam') then
      orb => bunch_params(i)%centroid
    else
      orb => orbit(i)
    endif
    if (orb%state == not_set$) return
    call convert_pc_to ((1 + orb%vec(6))*orb%p0c, orb%species, e_tot = value_vec(i))
  enddo

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('e_tot_ref')
  if (data_source == 'beam') return
  call tao_load_this_datum (branch%ele(:)%value(e_tot$), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('element_attrib.')

  name = upcase(datum%data_type(16:))
  value_vec = 0
  good = .false.

  do i = ix_start, ix_ele
    call pointer_to_attribute (branch%ele(i), name, .false., a_ptr, err, .false.)
    if (.not. associated (a_ptr%r)) cycle
    value_vec(i) = a_ptr%r
    good(i) = .true.
  enddo

  if (ix_ref > -1) then
    call pointer_to_attribute (ele_ref, name, .false., a_ptr, err, .false.)
    if (associated (a_ptr%r)) then
      value_vec(ix_ref) = a_ptr%r
      good(ix_ref) = .true.
    endif
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, good = good)

!-----------

case ('emit.', 'norm_emit.')

  if (associated(ele)) then
    call convert_total_energy_to (ele%value(e_tot$), lat%param%particle, gamma = gamma)
  else
    gamma = 0
  endif

  select case (datum%data_type)

  case ('emit.x', 'norm_emit.x')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = tao_lat_emit_calc (x_plane$, projected_emit$, branch%ele(i), tao_branch%modes)
      enddo
      if (ix_ref > -1) then
        value_vec(ix_ref) = tao_lat_emit_calc (x_plane$, projected_emit$, branch%ele(ix_ref), tao_branch%modes)
      endif
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

    if (datum%data_type == 'norm_emit.x') datum_value = datum_value * gamma

  case ('emit.y', 'norm_emit.y')  
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = tao_lat_emit_calc (y_plane$, projected_emit$, branch%ele(i), tao_branch%modes)
      enddo
      if (ix_ref > -1) then
        value_vec(ix_ref) = tao_lat_emit_calc (y_plane$, projected_emit$, branch%ele(ix_ref), tao_branch%modes)
      endif
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

    if (datum%data_type == 'norm_emit.y') datum_value = datum_value * gamma

  case ('emit.z', 'norm_emit.z')
    if (data_source == 'lat') then
      return
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%z%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

    if (datum%data_type == 'norm_emit.y') datum_value = datum_value * gamma

  case ('emit.a', 'norm_emit.a')
    if (data_source == 'lat') then
      if (lat%param%geometry == open$) then
        if (.not. allocated(tao_branch%rad_int%ele)) then
          call out_io (s_error$, r_name, 'tao_branch%rad_int not allocated')
          return
        endif
        call tao_load_this_datum (tao_branch%rad_int%ele%lin_norm_emit_a, &
                                ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
        datum_value = datum_value / gamma
      else
        datum_value = tao_branch%modes%a%emittance
        valid_value = .true.
      endif
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

    if (datum%data_type == 'norm_emit.a') datum_value = datum_value * gamma
    
  case ('emit.b', 'norm_emit.b')  
    if (data_source == 'lat') then
      if (lat%param%geometry == open$) then
        if (.not. allocated(tao_branch%rad_int%ele)) then
          call out_io (s_error$, r_name, 'tao_branch%rad_int not allocated')
          return
        endif
        call tao_load_this_datum (tao_branch%rad_int%ele%lin_norm_emit_b, &
                                ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
        datum_value = datum_value / gamma
      else
        datum_value = tao_branch%modes%b%emittance
        valid_value = .true.
      endif
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

    if (datum%data_type == 'norm_emit.b') datum_value = datum_value * gamma

  case ('emit.c', 'norm_emit.c')  
    if (data_source == 'lat') then
      return
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%c%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

    if (datum%data_type == 'norm_emit.y') datum_value = datum_value * gamma

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('eta.')

  select case (datum%data_type)

  case ('eta.a')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%a%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.b')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%b%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%x%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%y%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.z')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%z%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%z%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('etap.')

  select case (datum%data_type)

  case ('etap.a')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%a%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('etap.b')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%b%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('etap.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%x%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('etap.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%y%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('expression:')
  ! The point here is that tao_evaluate_stack is much quicker than tao_evaluate_expression.
  ! So on the fist time through, construct datum%stack and for subsequent times, use
  ! datum%stack with tao_evaluate_stack.
  !! if (allocated (datum%stack)) then
  !!   call tao_evaluate_stack (datum%stack, 1, .false., value1, good_exp, err, .true.)
  !!   if (err) return
  !!   datum_value = value1(1)
  !!   valid_value = good_exp(1)

  !! else ! Only do this first time through...
    write (dflt_dat_index, '(i0)') datum%ix_d1
    str = datum%data_type(12:)
    do
      ix = index(str, 'ele::#[')
      if (ix == 0) exit
      if (ix_ele == -1) then
        call tao_set_invalid (datum, 'NO ASSOCIATED ELEMENT' // datum%data_type(12:))
        return
      endif
      str = str(1:ix+4) // trim(ele_loc_to_string(ele)) // str(ix+6:)
    enddo

    call tao_evaluate_expression (str, 0, .false., expression_value_vec, info, err, .true., &
               datum%stack, 'model', datum%data_source, ele_ref, ele_start, ele, dflt_dat_index, u%ix_uni)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT EVALUATE EXPRESSION: ' // datum%data_type(12:))
      return
    endif

    select case (datum%merit_type)
    case ('min')
      datum_value = minval(expression_value_vec)
    case ('max') 
      datum_value = maxval(expression_value_vec)
    case ('abs_min')
      datum_value = minval(abs(expression_value_vec))
    case ('abs_max') 
      datum_value = maxval(abs(expression_value_vec))
    case ('target')
      if (size(expression_value_vec) /= 1) then
        call out_io (s_error$, r_name, &
                  'DATUM: ' // tao_datum_name(datum), &
                  'HAS "TARGET" MERIT_TYPE BUT DOES NOT EVALUATE TO A SINGLE NUMBER!')
        return
      endif
      datum_value = expression_value_vec(1)
    case default
      call out_io (s_error$, r_name, &
                  'SINCE THIS DATUM: ' // tao_datum_name(datum), &
                  'SPECIFIES A RANGE OF ELEMENTS, THEN THIS MERIT_TYPE: ' // datum%merit_type, &
                  'IS NOT VALID. VALID MERIT_TYPES ARE MIN, MAX, ABS_MIN, AND ABS_MAX.')
      return
    end select

    ! Make sure that any datums used in the expression have already been evaluated.
    do i = 1, size(datum%stack)
      if (datum%stack(i)%name == '') cycle
      call tao_find_data (err, datum%stack(i)%name, d_array = d_array, print_err = .false.)
      if (err .or. size(d_array) == 0) cycle  ! Err -> This is not associated then not a datum.
      dp => d_array(1)%d
      if (dp%d1%d2%ix_uni < u%ix_uni) cycle ! OK
      if (dp%d1%d2%ix_uni == u%ix_uni .and. dp%ix_data < datum%ix_data) cycle
      call out_io (s_error$, r_name, 'DATUM: ' // tao_datum_name(datum), &
                                     'WHICH IS OF TYPE EXPRESSION:' // datum%data_type(12:), &
                                     'THE EXPRESSION HAS A COMPONENT: ' // datum%stack(i)%name, &
                                     'AND THIS COMPONENT IS EVALUATED AFTER THE EXPRESSION!')
      return
    enddo
    valid_value = .true.
  !! endif

!-----------

case ('floor.')

  select case (datum%data_type)

  case ('floor.x')
    call tao_load_this_datum (branch%ele(:)%floor%r(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.y')
    call tao_load_this_datum (branch%ele(:)%floor%r(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.z')
    call tao_load_this_datum (branch%ele(:)%floor%r(3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.theta')
    call tao_load_this_datum (branch%ele(:)%floor%theta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.phi')
    call tao_load_this_datum (branch%ele(:)%floor%phi, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.psi')
    call tao_load_this_datum (branch%ele(:)%floor%psi, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('gamma.')

  select case (datum%data_type)

  case ('gamma.a')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%a%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%a%norm_emit /= 0)
      call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif
    
  case ('gamma.b')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%b%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%b%norm_emit /= 0)
      call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif

  case ('gamma.z')
    if (data_source == 'lat') return
    call tao_load_this_datum (bunch_params(:)%z%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('k.')

  select case (datum%data_type)

  case ('k.11b')
    if (data_source == 'beam') return
    call tao_load_this_datum (scratch%cc%k_11a, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
  case ('k.12a')
    if (data_source == 'beam') return
    call tao_load_this_datum (scratch%cc%k_12a, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
  case ('k.12b')
    if (data_source == 'beam') return
    call tao_load_this_datum (scratch%cc%k_12b, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
  case ('k.22a')
    if (data_source == 'beam') return
    call tao_load_this_datum (scratch%cc%k_22b, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('normal.')

  ! Fetches normal_form components.
  if (datum%ix_branch /= 0) then
    call out_io (s_fatal$, r_name, 'TRANSFER MAP CALC NOT YET MODIFIED FOR BRANCHES.')
    return
  endif

  if (data_source == 'beam') return

  if (s%global%rf_on) then
    normal_form => branch%normal_form_with_rf
  else
    normal_form => branch%normal_form_no_rf
  endif
  
  ! Do nothing it the map wasn't made
  if (.not. associated(normal_form%ele_origin) ) return



  ! Expect: taylor.#.######
  ! Example: normal.dhdj.2.000001 is the b-mode chromaticity
  !          head   sub
  ! Get position of first number. 
  iz = index(sub_data_type, '.') + 1
  
  ! Component i
  i = tao_read_this_index (sub_data_type, iz)

  ! Point to taylor
  taylor_is_complex = .false.
  if (sub_data_type(1:5) == 'dhdj.') then
    taylor_ptr => normal_form%dhdj(i)
  else if (sub_data_type(1:2) == 'A.') then
    taylor_ptr => normal_form%A(i)
  else if (sub_data_type(1:2) == 'A_inv.') then
    taylor_ptr => normal_form%A_inv(i)
  else if (sub_data_type(1:2) == 'M.') then
    taylor_ptr => normal_form%M(i)
  else if (sub_data_type(1:4) == 'ReF.') then
    taylor_is_complex = .true.
    use_real_part = .true.
    complex_taylor_ptr => normal_form%f(i)   
  else if (sub_data_type(1:4) == 'ImF.') then
    taylor_is_complex = .true.
    use_real_part = .false.
    complex_taylor_ptr => normal_form%F(i)      
  else if (sub_data_type(1:4) == 'ReL.') then
    taylor_is_complex = .true.
    use_real_part = .true.
    complex_taylor_ptr => normal_form%L(i)   
  else if (sub_data_type(1:4) == 'ImL.') then
    taylor_is_complex = .true.
    use_real_part = .false.
    complex_taylor_ptr => normal_form%L(i)       
    
  endif
 
  ! Check for second dot
  if (sub_data_type(iz+1:iz+1) /= '.') then
   call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
   call out_io (s_error$, r_name, 'datum%data_type: '//trim(datum%data_type) )
   call out_io (s_error$, r_name, 'expect dot: ', sub_data_type(1:iz)//'.######' )
  endif
 
  ! Get exponent
  expn_str = sub_data_type(iz+2:iz+7)
  expnt = 0
  do j = 1, 6
    if (expn_str(j:j) == ' ') exit
    expnt(j) = index('0123456789', expn_str(j:j)) - 1
  enddo
  
  ! Coefficient
  if (taylor_is_complex) then
    if (use_real_part) then
      datum_value = real(complex_taylor_coef(complex_taylor_ptr, expnt))
    else
      datum_value = aimag(complex_taylor_coef(complex_taylor_ptr, expnt))
    endif
  else
    datum_value = taylor_coef(taylor_ptr, expnt)
  endif
  valid_value = .true.  


!-----------

case ('momentum')
  if (data_source == 'beam') return
  call tao_load_this_datum (branch%ele(0:n_track)%value(p0c$) * (1+orbit(0:n_track)%vec(6)), &
                            ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('momentum_compaction')
  if (data_source == 'beam') return

  if (ix_ref < 0) then
    ix_ref = 0
    ele_ref => lat%branch(datum%ix_branch)%ele(ix_ref)
  endif

  orb0 => orbit(ix_ref)
  call make_v_mats (ele_ref, v_mat, v_inv_mat)
  eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb0%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz

  call transfer_matrix_calc (lat, mat6, vec0, ix_ref, ix_start, datum%ix_branch)

  do i = ix_start, ix_ele
    s_len = branch%ele(i)%s - branch%ele(ix_ref)%s
    if (s_len == 0) then
      value_vec(i) = 0
    else
      value_vec(i) = -(sum(mat6(5,1:4) * eta_vec) + mat6(5,6)) / s_len
    endif
    if (i /= ix_ele) mat6 = matmul(branch%ele(i+1)%mat6, mat6)
  enddo
  call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('n_particle_loss')
  if (data_source /= 'beam') return
  if (ix_ele < 0) ix_ele = branch%n_ele_track
  datum_value = sum(bunch_params(ix_ref+1:ix_ele)%n_particle_lost_in_ele)
  valid_value = .true.

!-----------

case ('orbit.')

  if (tao_branch%track_state /= moving_forward$ .and. ix_ele > tao_branch%track_state) then
    valid_value = .false.
    call tao_set_invalid (datum, 'Particle lost.', why_invalid)
    return
  endif

  select case (datum%data_type)

  case ('orbit.e_tot')
    if (ix_ref > -1) then
      if (data_source == 'beam') then
        orb => bunch_params(ix_ref)%centroid
      else
        orb => orbit(ix_ref)
      endif
      if (orb%state == not_set$) return
      value_vec(ix_ref) = (1 + orb%vec(6)) * orb%p0c / orb%beta
    endif

    do i = ix_start, ix_ele
      if (data_source == 'beam') then
        orb => bunch_params(i)%centroid
      else
        orb => orbit(i)
      endif
      if (orb%state == not_set$) return
      call convert_pc_to ((1 + orb%vec(6))*orb%p0c, orb%species, e_tot = value_vec(i))
    enddo

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('orbit.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (orbit(:)%vec(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (orbit(:)%vec(3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.z')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (orbit(:)%vec(5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.px')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (orbit(:)%vec(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.py')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (orbit(:)%vec(4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.pz')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (orbit(:)%vec(6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.amp_a')
    if (data_source == 'beam') return ! bad
    call tao_load_this_datum (scratch%cc%amp_a, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case ('orbit.amp_b')
    if (data_source == 'beam') return ! bad
    call tao_load_this_datum (scratch%cc%amp_b, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case ('orbit.norm_amp_a')
    if (data_source == 'beam') return ! bad
    call tao_load_this_datum (scratch%cc%amp_na, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case ('orbit.norm_amp_b')
    if (data_source == 'beam') return ! bad
    call tao_load_this_datum (scratch%cc%amp_nb, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('pc')
  if (ix_ref > -1) then
    if (data_source == 'beam') then
      value_vec(ix_ref) = (1 + bunch_params(ix_ref)%centroid%vec(6)) * bunch_params(ix_ref)%centroid%p0c
    else
      value_vec(ix_ref) = (1 + orbit(ix_ref)%vec(6)) * orbit(ix_ref)%p0c
    endif
  endif

  if (data_source == 'beam') then
    do i = ix_start, ix_ele
      value_vec(i) = (1 + bunch_params(i)%centroid%vec(6)) * bunch_params(ix_ref)%centroid%p0c
    enddo
  else
    do i = ix_start, ix_ele
      value_vec(i) = (1 + orbit(i)%vec(6)) * orbit(ix_ref)%p0c
    enddo
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('periodic.')

  ix = index(datum%data_type(10:), '.') + 9
  select case (datum%data_type(1:ix))

  case ('periodic.tt.')
    if (data_source == 'beam') return
    if (lat%param%geometry /= closed$ .and. .not. associated(ele_ref)) then
      call tao_set_invalid (datum, 'LATTICE MUST BE CIRCULAR FOR A DATUM LIKE: ' // datum%data_type, why_invalid)
      call err_exit
    endif

    call transfer_map_calc (lat, taylor, err, ix_ele, ix_ele, orbit(ix_ele), branch%ix_branch, one_turn = .true.)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT INVERT MAP', why_invalid)
      return
    endif

    do i = 1, 4
      call add_taylor_term (taylor(i), -1.0_rp, taylor_expn([i]))
    enddo
    call taylor_inverse (taylor, taylor, err)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT INVERT MAP', why_invalid)
      return
    endif

    expnt = 0
    i = tao_read_this_index (datum%data_type, 13); if (i == 0) return
    do j = 14, 24
      if (datum%data_type(j:j) == ' ') exit
      k = tao_read_this_index (datum%data_type, j); if (k == 0) return
      expnt(k) = expnt(k) + 1
    enddo

    datum_value = taylor_coef (taylor(i), expnt)
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('phase.', 'phase_frac.')

  select case (datum%data_type)

  case ('phase.a', 'phase_frac.a')
    if (data_source == 'beam') return ! bad
    if (ix_ref < 0) then
      datum_value = ele%a%phi
    else
      datum_value = ele%a%phi - ele_ref%a%phi
      if (ix_ref > ix_ele) datum_value = datum_value - branch%ele(0)%a%phi + branch%ele(n_track)%a%phi 
    if (datum%data_type == 'phase_frac.a') datum_value = modulo2(datum_value, pi)
    endif
    valid_value = .true.

  case ('phase.b', 'phase_frac.b')
    if (data_source == 'beam') return ! bad
    if (ix_ref < 0) then
      datum_value = ele%b%phi
    else
      datum_value = ele%b%phi - ele_ref%b%phi
      if (ix_ref > ix_ele) datum_value = datum_value - branch%ele(0)%b%phi + branch%ele(n_track)%b%phi 
    endif
    if (datum%data_type == 'phase_frac.b') datum_value = modulo2(datum_value, pi)
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('phase_frac_diff')
  if (data_source == 'beam') return ! bad
  if (ix_ref < 0) then
    px = ele%a%phi 
    py = ele%b%phi 
  else
    px = ele%a%phi - ele_ref%a%phi
    py = ele%b%phi - ele_ref%b%phi
    if (ix_ref > ix_ele) px = px - branch%ele(0)%a%phi + branch%ele(n_track)%a%phi 
    if (ix_ref > ix_ele) py = py - branch%ele(0)%b%phi + branch%ele(n_track)%b%phi 
  endif
  datum_value = modulo2 (px - py, pi)
  valid_value = .true.

!-----------

case ('photon.')

  select case (datum%data_type)

  case ('photon.intensity_x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%centroid%field(1)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (orbit(:)%field(1)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('photon.intensity_y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%centroid%field(2)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (orbit(:)%field(2)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('photon.intensity')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%centroid%field(1)**2+bunch_params(:)%centroid%field(2)**2, &
                                                                      ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (orbit(:)%field(1)**2 + orbit(:)%field(2)**2, ele_ref, ele_start, ele, &
                                                                           datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('photon.phase_x')
    if (data_source == 'beam') return ! bad
    call tao_load_this_datum (orbit(:)%phase(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('photon.phase_y')
    if (data_source == 'beam') return ! bad
    call tao_load_this_datum (orbit(:)%phase(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  end select

!-----------

case ('ping_a.')

  if (.not. associated(ele)) return  ! Bad

  select case (datum%data_type)
  case ('ping_a.amp_x')
    datum_value = ele%gamma_c * sqrt(ele%a%beta)
    if (associated(ele_ref)) datum_value = datum_value - ele%gamma_c * sqrt(ele_ref%a%beta)
    valid_value = .true.

  case ('ping_a.phase_x')
    datum_value = ele%a%phi
    if (associated(ele_ref)) datum_value = datum_value - ele_ref%a%phi
    valid_value = .true.

  case ('ping_a.amp_y')
    call tao_scratch_values_calc (ix_ele, datum, branch, orbit)
    datum_value = sqrt(ele%b%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(2,2)**2))
    if (associated(ele_ref)) then
      call tao_scratch_values_calc (ix_ref, datum, branch, orbit)
      datum_value = datum_value - sqrt(ele_ref%b%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(2,2)**2))
    endif
    valid_value = .true.

  case ('ping_a.phase_y')
    datum_value = ele%a%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), - scratch%cc(ix_ele)%cbar(2,2))
    if (associated(ele_ref)) then
      call tao_scratch_values_calc (ix_ref, datum, branch, orbit)
      datum_value = datum_value - ele_ref%a%phi - atan2(scratch%cc(ix_ref)%cbar(1,2), - scratch%cc(ix_ref)%cbar(2,2))
    endif
    valid_value = .true.

  case ('ping_a.sin_y')
    call tao_scratch_values_calc (ix_ele, datum, branch, orbit)
    datum_value = -sqrt(ele%b%beta) * scratch%cc(ix_ele)%cbar(1,2)

    if (associated(ele_ref)) then
      call tao_scratch_values_calc (ix_ref, datum, branch, orbit)
      datum_value = datum_value + sqrt(ele_ref%b%beta) * scratch%cc(ix_ref)%cbar(1,2)
    endif
    valid_value = .true.

  case ('ping_a.cos_y')
    call tao_scratch_values_calc (ix_ele, datum, branch, orbit)
    datum_value = -sqrt(ele%b%beta) * scratch%cc(ix_ele)%cbar(2,2)

    if (associated(ele_ref)) then
      call tao_scratch_values_calc (ix_ref, datum, branch, orbit)
      datum_value = datum_value + sqrt(ele_ref%b%beta) * scratch%cc(ix_ref)%cbar(2,2)
    endif
    valid_value = .true.
  end select

!-----------

case ('ping_b.')

  if (.not. associated(ele)) return  ! Bad

  select case (datum%data_type)
  case ('ping_b.amp_y')
    datum_value = ele%gamma_c * sqrt(ele%b%beta)
    if (associated(ele_ref)) datum_value = datum_value - ele%gamma_c * sqrt(ele_ref%b%beta)
    valid_value = .true.

  case ('ping_b.phase_y')
    datum_value = ele%b%phi
    if (associated(ele_ref)) datum_value = datum_value - ele_ref%b%phi
    valid_value = .true.

  case ('ping_b.amp_x')
    call tao_scratch_values_calc (ix_ele, datum, branch, orbit)
    datum_value = sqrt(ele%a%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(1,1)**2))
    if (associated(ele_ref)) then
      call tao_scratch_values_calc (ix_ref, datum, branch, orbit)
      datum_value = datum_value - sqrt(ele_ref%a%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(1,1)**2))
    endif
    valid_value = .true.

  case ('ping_b.phase_x')
    call tao_scratch_values_calc (ix_ele, datum, branch, orbit)
    datum_value = ele%b%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), scratch%cc(ix_ele)%cbar(1,1))
    if (associated(ele_ref)) then
      call tao_scratch_values_calc (ix_ref, datum, branch, orbit)
      datum_value = datum_value - ele_ref%b%phi - atan2(scratch%cc(ix_ref)%cbar(1,2), scratch%cc(ix_ref)%cbar(1,1))
    endif
    valid_value = .true.

  case ('ping_b.sin_x')
    call tao_scratch_values_calc (ix_ele, datum, branch, orbit)
    datum_value = -sqrt(ele%a%beta) * scratch%cc(ix_ele)%cbar(1,2)

    if (associated(ele_ref)) then
      call tao_scratch_values_calc (ix_ref, datum, branch, orbit)
      datum_value = datum_value + sqrt(ele_ref%a%beta) * scratch%cc(ix_ref)%cbar(1,2)
    endif
    valid_value = .true.

  case ('ping_b.cos_x')
    call tao_scratch_values_calc (ix_ele, datum, branch, orbit)
    datum_value = sqrt(ele%a%beta) * scratch%cc(ix_ele)%cbar(1,1)

    if (associated(ele_ref)) then
      call tao_scratch_values_calc (ix_ref, datum, branch, orbit)
      datum_value = datum_value - sqrt(ele_ref%a%beta) * scratch%cc(ix_ref)%cbar(1,1)
    endif
    valid_value = .true.
  end select

!-----------

case ('r.')
  if (data_source == 'beam') return
  i = tao_read_this_index (datum%data_type, 3); if (i == 0) return
  j = tao_read_this_index (datum%data_type, 4); if (j == 0) return

  if (ix_ref < 0) ix_ref = 0
  call transfer_matrix_calc (lat, mat6, vec0, ix_ref, ix_start, datum%ix_branch)

  k = ix_start
  do 
    value_vec(k) = mat6(i, j)
    if (k == ix_ele) exit
    k = k + 1
    if (k > n_track) k = 0
    mat6 = matmul(branch%ele(k)%mat6, mat6)
  enddo

  call tao_load_this_datum (value_vec, NULL(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('r56_compaction')
  if (data_source == 'beam') return

  if (ix_ref < 0) then
    ix_ref = 0
    ele_ref => lat%branch(datum%ix_branch)%ele(ix_ref)
  endif

  orb0 => orbit(ix_ref)
  call make_v_mats (ele_ref, v_mat, v_inv_mat)
  eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb0%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz

  call transfer_matrix_calc (lat, mat6, vec0, ix_ref, ix_start, datum%ix_branch)

  do i = ix_start, ix_ele
    value_vec(i) = sum(mat6(5,1:4) * eta_vec) + mat6(5,6)
    if (i /= ix_ele) mat6 = matmul(branch%ele(i+1)%mat6, mat6)
  enddo
  call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('rad_int.')

  if (ix_ref > -1 .or. ix_ele > -1) then
    if (ix_ele < 0) ix_ele = branch%n_ele_track
    if (ix_ref < 0) ix_ref = 0
  endif

  if (data_source == 'beam') return
  if (.not. allocated(tao_branch%rad_int%ele)) return

  select case (datum%data_type(9:))
  case ('i1')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%i1)
    else
      datum_value = tao_branch%modes%synch_int(1)
    endif

  case ('i2')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%i2)
    else
      datum_value = tao_branch%modes%synch_int(2)
    endif

  case ('i2_e4')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%lin_i2_e4)
    else
      datum_value = tao_branch%modes%lin%i2_e4
    endif

  case ('i3')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%i3)
    else
      datum_value = tao_branch%modes%synch_int(3)
    endif

  case ('i3_e7')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%lin_i3_e7)
    else
      datum_value = tao_branch%modes%lin%i3_e7
    endif

  case ('i4a')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%i4a)
    else
      datum_value = tao_branch%modes%a%synch_int(4)
    endif

  case ('i4b')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%i4b)
    else
      datum_value = tao_branch%modes%b%synch_int(4)
    endif

  case ('i4z')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%i4z)
    else
      datum_value = tao_branch%modes%z%synch_int(4)
    endif

  case ('i5a')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%i5a)
    else
      datum_value = tao_branch%modes%a%synch_int(5)
    endif

  case ('i5a_e6')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%lin_i5a_e6)
    else
      datum_value = tao_branch%modes%lin%i5a_e6
    endif

  case ('i5b')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%i5b)
    else
      datum_value = tao_branch%modes%b%synch_int(5)
    endif

  case ('i5b_e6')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%lin_i5b_e6)
    else
      datum_value = tao_branch%modes%lin%i5b_e6
    endif

  case ('i6b')
    if (ix_ele > -1) then
      datum_value = sum(tao_branch%rad_int%ele(ix_ref:ix_ele)%i6b)
    else
      datum_value = tao_branch%modes%b%synch_int(6)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

  valid_value = .true.

!-----------

case ('rad_int1.')

  if (data_source == 'beam') return
  if (ix_ele < 0) return
  if (.not. allocated(tao_branch%rad_int%ele)) return

  select case (datum%data_type(10:))
  case ('i1')
    datum_value = tao_branch%rad_int%ele(ix_ele)%i1
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%i1

  case ('i2')
    datum_value = tao_branch%rad_int%ele(ix_ele)%i2
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%i2

  case ('i2_e4')
    datum_value = tao_branch%rad_int%ele(ix_ele)%lin_i2_e4
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%lin_i2_e4

  case ('i3')
    datum_value = tao_branch%rad_int%ele(ix_ele)%i3
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%i3

  case ('i3_e7')
    datum_value = tao_branch%rad_int%ele(ix_ele)%lin_i3_e7
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%lin_i3_e7

  case ('i4a')
    datum_value = tao_branch%rad_int%ele(ix_ele)%i4a
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%i4a

  case ('i5a')
    datum_value = tao_branch%rad_int%ele(ix_ele)%i5a
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%i5a

  case ('i5a_e6')
    datum_value = tao_branch%rad_int%ele(ix_ele)%lin_i5a_e6
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%lin_i5a_e6

  case ('i4b')
    datum_value = tao_branch%rad_int%ele(ix_ele)%i4b
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%i4b

  case ('i5b')
    datum_value = tao_branch%rad_int%ele(ix_ele)%i5b
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%i5b

  case ('i5b_e6')
    datum_value = tao_branch%rad_int%ele(ix_ele)%lin_i5b_e6
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%lin_i5b_e6

  case ('i6b')
    datum_value = tao_branch%rad_int%ele(ix_ele)%i6b
    if (ix_ref > -1) datum_value = datum_value - tao_branch%rad_int%ele(ix_ref)%i6b

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

  valid_value = .true.

!-----------

case ('ref_time')
    call tao_load_this_datum (branch%ele(:)%ref_time, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    
!-----------

case ('rel_floor.')

  select case (datum%data_type)

  case ('rel_floor.x', 'rel_floor.y', 'rel_floor.z')

    if (ix_ref < 0) ele_ref => lat%branch(datum%ix_branch)%ele(0)

    call floor_angles_to_w_mat (-ele_ref%floor%theta, -ele_ref%floor%phi, -ele_ref%floor%psi, w0_mat)

    i = ix_start
    do 
      ele2 => branch%ele(i)
      vec3 = ele2%floor%r - ele_ref%floor%r
      vec3 = matmul (w0_mat, vec3)
      select case (datum%data_type)
      case ('rel_floor.x')
        value_vec(i) = vec3(1)
      case ('rel_floor.y')
        value_vec(i) = vec3(2)
      case ('rel_floor.z')
        value_vec(i) = vec3(3)
      end select
      if (i == ix_ele) exit
      i = i + 1
      if (i > n_track) i = 0
    enddo

    call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('rel_floor.theta', 'rel_floor.phi', 'rel_floor.psi')

    if (ix_ref < 0) ele_ref => lat%branch(datum%ix_branch)%ele(0)

    call floor_angles_to_w_mat (-ele_ref%floor%theta, -ele_ref%floor%phi, -ele_ref%floor%psi, w0_mat)

    i = ix_start
    do 
      ele2 => branch%ele(i)
      call floor_angles_to_w_mat (ele2%floor%theta, ele2%floor%phi, ele2%floor%psi, w_mat)
      w_mat = matmul (w0_mat, w_mat)
      call floor_w_mat_to_angles (w_mat, theta, phi, psi)

      select case (datum%data_type)
      case ('rel_floor.theta')
        value_vec(i) = theta
      case ('rel_floor.phi')
        value_vec(i) = phi
      case ('rel_floor.psi')
        value_vec(i) = psi
      end select
      if (i == ix_ele) exit
      i = i + 1
      if (i > n_track) i = 0
    enddo

    call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('sigma.')

  ! Looks for numbers: e.g. sigma.13
  i = index('123456', datum%data_type(7:7))
  j = index('123456', datum%data_type(8:8))
  if (i > 0 .and. j > 0) then
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%linear%sigma(i,j), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(i,j), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
    return
  endif

  select case (datum%data_type)

  case ('sigma.x')
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%linear%sigma(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    endif

  case ('sigma.px')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%linear%sigma(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    endif

  case ('sigma.y')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%linear%sigma(3,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(3,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    endif
    
  case ('sigma.py')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%linear%sigma(4,4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(4,4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    endif
    
  case ('sigma.z')
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%linear%sigma(5,5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(5,5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    endif

  case ('sigma.pz')  
    if (data_source == 'lat') then
      if (lat%param%geometry == closed$) then
        call tao_load_this_datum (tao_branch%linear%sigma(6,6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit)
        datum_value = sqrt(datum_value)
      else
        if (ix_ele == -1) ix_ele = branch%n_ele_track
        datum_value = sqrt(4 * const_q_factor * classical_radius_factor * &
                                 sum(tao_branch%rad_int%ele(ix_ref+1:ix_ele)%lin_i3_e7) / 3) / mass_of(lat%param%particle)
        valid_value = .true.
      endif
    else
      call tao_load_this_datum (bunch_params(:)%sigma(6,6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    endif
    
  case ('sigma.xy')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%linear%sigma(1,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(1,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('sigma.Lxy')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%linear%sigma(1,4) - tao_branch%linear%sigma(2,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(1,4) - bunch_params(:)%sigma(2,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('spin.')

  if (.not. bmad_com%spin_tracking_on) then
    call tao_set_invalid (datum, 'NO SPIN TRACKING WHEN BMAD_COM%SPIN_TRACKING_ON = FALSE!', why_invalid)
    return
  endif

  select case (datum%data_type)

  case ('spin.x', 'spin.y', 'spin.z')
    do i = ix_start, ix_ele
      if (data_source == 'beam') then
        vec3 = polar_to_vec(bunch_params(i)%spin)
      else
        vec3 = orbit(i)%spin
      endif

      select case (datum%data_type)
      case ('spin.x');  value_vec(i) = vec3(1)
      case ('spin.y');  value_vec(i) = vec3(2)
      case ('spin.z');  value_vec(i) = vec3(3)
      end select
    enddo

    if (ix_ref > -1) then
      if (data_source == 'beam') then
        vec3 = polar_to_vec(bunch_params(ix_ref)%spin)
      else
        vec3 = orbit(ix_ref)%spin
      endif

      select case (datum%data_type)
      case ('spin.x');  value_vec(ix_ref) = vec3(1)
      case ('spin.y');  value_vec(ix_ref) = vec3(2)
      case ('spin.z');  value_vec(ix_ref) = vec3(3)
      end select
    endif

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('spin.theta', 'spin.phi', 'spin.amplitude')
    do i = ix_start, ix_ele
      if (data_source == 'beam') then
        polar_spin = bunch_params(i)%spin
      else
        polar_spin = vec_to_polar(orbit(i)%spin)
      endif

      select case (datum%data_type)
      case ('spin.theta');   value_vec(i) = polar_spin%theta
      case ('spin.phi');     value_vec(i) = polar_spin%phi
      case ('spin.amp');     value_vec(i) = polar_spin%polarization
      end select
    enddo

    if (ix_ref > -1) then
      if (data_source == 'beam') then
        polar_spin = bunch_params(ix_ref)%spin
      else
        polar_spin = vec_to_polar(orbit(ix_ref)%spin)
      endif

      select case (datum%data_type)
      case ('spin.theta');   value_vec(ix_ref) = polar_spin%theta
      case ('spin.phi');     value_vec(ix_ref) = polar_spin%phi
      case ('spin.amp');     value_vec(ix_ref) = polar_spin%polarization
      end select
    endif

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('s_position') 
  if (data_source == 'beam') return
  if (ix_ref >= 0) then
    datum_value = ele%s - ele_ref%s
  else
    datum_value = ele%s 
  endif
  valid_value = .true.

!-----------

case ('time')
  if (data_source == 'beam') then
    call tao_load_this_datum (bunch_params%centroid%t, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  else
    call tao_load_this_datum (orbit(:)%t, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  endif

!-----------

case ('tune.')

  select case (datum%data_type)

  case ('tune.a')
    if (data_source == 'beam') return ! bad
    datum_value = branch%ele(branch%n_ele_track)%a%phi
    valid_value = .true.

  case ('tune.b')
    if (data_source == 'beam') return ! bad
    datum_value = branch%ele(branch%n_ele_track)%b%phi
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('t.', 'tt.')
  if (datum%ix_branch /= 0) then
    call out_io (s_fatal$, r_name, 'TRANSFER MAP CALC NOT YET MODIFIED FOR BRANCHES.')
    return
  endif
  if (data_source == 'beam') return

  expnt = 0
  if (head_data_type == 't.') then
    i = tao_read_this_index (datum%data_type, 3)
    do j = 4, 5
      k = tao_read_this_index (datum%data_type, j); if (k == 0) exit
      expnt(k) = expnt(k) + 1
    enddo
  else
    i = tao_read_this_index (datum%data_type, 4)
    do j = 5, 15
      if (datum%data_type(j:j) == ' ') exit
      k = tao_read_this_index (datum%data_type, j); if (k == 0) exit
      expnt(k) = expnt(k) + 1
    enddo
  endif

  if (i == 0 .or. k == 0) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return
  endif

  if (ix_ref < 0) ix_ref = 0

  ! Computation if there is no range

  if (ix_start == ix_ele) then
    if (s%com%ix_ref_taylor /= ix_ref .or. s%com%ix_ele_taylor /= ix_ele) then
      ix0 = s%com%ix_ele_taylor
      if (s%com%ix_ref_taylor == ix_ref .and. ix_ele > ix0) then
        call transfer_map_calc (lat, taylor_save, err, ix0, ix_ele, orbit(ix0), unit_start = .false.)
      else
        call transfer_map_calc (lat, taylor_save, err, ix_ref, ix_ele, orbit(ix_ref))
      endif

      if (err) then
        call tao_set_invalid (datum, 'MAP TERM OVERFLOW', why_invalid)
        return
      endif

      s%com%ix_ref_taylor = ix_ref
      s%com%ix_ele_taylor = ix_ele
    endif
    datum_value = taylor_coef (taylor_save(i), expnt)
    valid_value = .true.

  ! Here if there is a range.
  else
    k = ix_start
    call transfer_map_calc (lat, taylor, err, ix_ref, k, orbit(ix_ref))
    if (err) then
      call tao_set_invalid (datum, 'MAP TERM OVERFLOW', why_invalid)
      return
    endif

    do
      value_vec(k) = taylor_coef (taylor(i), expnt)
      if (k == ix_ele) exit
      k_old = k
      k = k + 1
      if (k > branch%n_ele_track) k = 0
      call transfer_map_calc (lat, taylor, err, k_old, k, unit_start = .false.)
      if (err) then
       call tao_set_invalid (datum, 'MAP TERM OVERFLOW', why_invalid)
        return
      endif
    enddo
    call tao_load_this_datum (value_vec, NULL(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  endif

!-----------

case ('unstable.')

  select case (datum%data_type)

  case ('unstable.orbit')

    if (lat%param%geometry /= open$) return
    if (datum%ele_name == '') ix_ele = branch%n_ele_track

    if (data_source == 'beam') then
      datum_value = 0
      do i = 1, ix_ele
        datum_value = datum_value + (1 + ix_ele - i) * bunch_params(i)%n_particle_lost_in_ele
      enddo
      datum_value = datum_value / tao_branch%bunch_params(ix_ele)%n_particle_tot
      datum%ix_ele_merit = -1

    else
      iz = tao_branch%track_state
      if (iz /= moving_forward$ .and. iz <= ix_ele) then
        datum_value = 1 + ix_ele - iz
        if (orbit(iz)%s < branch%ele(iz)%s) then
          orb => orbit(iz-1)
          datum%ix_ele_merit = iz - 1
          if (branch%ele(iz)%value(L$) /= 0) then
            ! Add s_rel/L 
            datum_value = datum_value + (branch%ele(iz)%s - orbit(iz)%s)/branch%ele(iz)%value(L$)
          endif
        else
          datum_value = datum_value - 0.5
          orb => orbit(iz)
          datum%ix_ele_merit = iz
        endif

        datum_value = datum_value + 0.5 * tanh(lat%param%unstable_factor)

      endif
    endif
    valid_value = .true.

  case ('unstable.ring')
    if (data_source == 'beam') return
    datum_value = lat%param%unstable_factor
    ! unstable_penalty is needed since at the meta stable borderline the growth rate is zero.
    if (.not. lat%param%stable) datum_value = datum_value + s%global%unstable_penalty
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
    return

  end select

!-----------

case ('wall.')
  if (data_source == 'beam') return

  constraint = datum%data_type(6:)
  z0_pt = 0
  datum_value = 1e10   ! Something large

  do i = 1, size(s%building_wall%section)
    section => s%building_wall%section(i)
    if (section%constraint /= constraint) cycle
    do ie = ix_start, ix_ele
      ele => branch%ele(ie)
      
      do is = 1, size(section%point)
        pt => section%point(is)
        dz = pt%z - ele%floor%r(3); dx = pt%x - ele%floor%r(1)
        cos_theta = cos(ele%floor%theta); sin_theta = sin(ele%floor%theta)
        z_pt =  dz * cos_theta + dx * sin_theta
        x_pt = -dz * sin_theta + dx * cos_theta

        ! The perpendicular to the machine line intersects the segment if
        ! z_pt and z0_pt have a different sign.

        if (is > 1 .and. z_pt * z0_pt <= 0) then
          if (pt%radius == 0 .or. z_pt == 0 .or. z0_pt == 0) then
            x_wall = (x0_pt * z_pt - x_pt * z0_pt) / (z_pt - z0_pt)
          else
            dz = pt%z_center - ele%floor%r(3); dx = pt%x_center - ele%floor%r(1)
            z_center =  dz * cos_theta + dx * sin_theta
            x_center = -dz * sin_theta + dx * cos_theta
            if (pt%radius * z_pt > 0) then
              x_wall = x_center - sqrt(pt%radius**2 - z_center**2)
            else
              x_wall = x_center + sqrt(pt%radius**2 - z_center**2)
            endif
          endif
          if (datum%data_type =='wall.right_side')  x_wall = -x_wall
          datum_value = min(datum_value, x_wall)
        endif

        z0_pt = z_pt
        x0_pt = x_pt
      enddo
    enddo

  enddo

  valid_value = .true.

!-----------

case ('wire.')  
  if (data_source == 'lat') return
  read (datum%data_type(6:), '(a)') angle
  datum_value = tao_do_wire_scan (ele, angle, u%uni_branch(datum%ix_branch)%ele(ix_ele)%beam)
  valid_value = .true.
  
case default
  call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" NOT VALID', why_invalid)
  return

end select

! Set s position

if (datum%ix_ele_merit > -1) then
  datum%s = lat%branch(datum%ix_branch)%ele(datum%ix_ele_merit)%s
elseif (associated(ele)) then
  datum%s = ele%s
endif

end subroutine tao_evaluate_a_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine tao_load_this_datum (vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, orbit, good)

type (tao_data_struct) datum
type (branch_struct) branch
type (ele_struct), pointer :: ele_ref, ele_start, ele
type (coord_struct), optional :: orbit(0:)

real(rp), target :: vec(0:)
real(rp) datum_value, ref_value
real(rp), pointer :: vec_ptr(:)

character(20) :: r_name = 'tao_evaluate_a_datum'
character(*), optional :: why_invalid

integer ix_m, i, n_track, ix_m2, ix_ref, ix_start, ix_ele

logical valid_value
logical, optional :: good(0:)

!

valid_value = .true.

n_track = branch%n_ele_track
ix_start = -1; ix_ref = -1; ix_ele = -1
if (associated(ele)) ix_ele = ele%ix_ele
if (associated(ele_ref)) ix_ref = ele_ref%ix_ele
if (associated(ele_start)) ix_start = ele_start%ix_ele

if (ix_ele < 0) then
  datum%exists = .false.
  valid_value = .false.
  call tao_set_invalid (datum, 'ELEMENT INDEX FOR DATUM IS NEGATIVE!', why_invalid)
  return
endif

if (ix_ele < ix_start .and. branch%param%geometry == open$) then
  if (datum%useit_opt) call out_io (s_error$, r_name, &
                'ERROR: ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum), &
                'STARTING ELEMENT: ' // ele_start%name, &
                'IS AFTER ENDING ELEMENT: ' // ele%name)
  call tao_set_invalid (datum, 'ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum), why_invalid)
  valid_value = .false.
  return
endif

! Set up refrence value

if (ix_ref > -1) then
  if (present(orbit)) call tao_scratch_values_calc (ix_ref, datum, branch, orbit)
  ref_value = vec(ix_ref)
  if (present(good)) then
    if (.not. good(ix_ref)) then
      call tao_set_invalid (datum, 'DATA AT REFERENCE ELEMENT NOT VALID', why_invalid)
      valid_value = .false.
      return
    endif
  endif
else
  ref_value = 0
endif

 
! If ele_start does not exist

if (datum%ele_start_name == '' .or. ix_start == ix_ele) then
  if (present(orbit)) call tao_scratch_values_calc (ix_ele, datum, branch, orbit)
  datum_value = vec(ix_ele) - ref_value
  if (datum%merit_type(1:4) == 'abs_') datum_value = abs(vec(ele%ix_ele))
  if (present(good)) valid_value = good(ix_ele)
  if (.not. valid_value) call tao_set_invalid (datum, 'DATA AT START ELEMENT NOT VALID.', why_invalid)

  return
endif

! Set up the vector of values with the reference subtracted off

if (ref_value == 0) then
  vec_ptr => vec
else
  allocate(vec_ptr(0:ubound(vec,1)))
  vec_ptr = vec - ref_value
endif

!------------------------
! If there is a range

if (ix_ele < ix_start) then   ! wrap around

  if (present(orbit)) then
    do i = ix_start, n_track
      call tao_scratch_values_calc (i, datum, branch, orbit)
    enddo
    do i = 0, ix_ele
      call tao_scratch_values_calc (i, datum, branch, orbit)
    enddo
  endif

  select case (datum%merit_type)
  case ('min')
    ix_m = minloc (vec_ptr(0:ix_ele), 1) - 1
    ix_m2 = minloc (vec_ptr(ix_start:n_track), 1) + ix_start - 1
    if (vec_ptr(ix_m2) < vec_ptr(ix_m2)) ix_m = ix_m2
    datum_value = vec_ptr(ix_m)
    if (present(good)) valid_value = good(ix_m)

  case ('max')
    ix_m = maxloc (vec_ptr(0:ix_ele), 1) - 1
    ix_m2 = maxloc (vec_ptr(ix_start:n_track), 1) + ix_start - 1
    if (vec_ptr(ix_m2) > vec_ptr(ix_m2)) ix_m = ix_m2
    datum_value = vec_ptr(ix_m)
    if (present(good)) valid_value = good(ix_m)

  case ('abs_min')
    ix_m = minloc (abs(vec_ptr(0:ix_ele)), 1) - 1
    ix_m2 = minloc (abs(vec_ptr(ix_start:n_track)), 1) + ix_start - 1
    if (abs(vec_ptr(ix_m2)) < abs(vec_ptr(ix_m2))) ix_m = ix_m2
    datum_value = abs(vec_ptr(ix_m))
    if (present(good)) valid_value = good(ix_m)

  case ('abs_max')
    ix_m = maxloc (abs(vec_ptr(0:ix_ele)), 1) - 1
    ix_m2 = maxloc (abs(vec_ptr(ix_start:n_track)), 1) + ix_start - 1
    if (abs(vec_ptr(ix_m2)) > abs(vec_ptr(ix_m2))) ix_m = ix_m2
    datum_value = abs(vec_ptr(ix_m))
    if (present(good)) valid_value = good(ix_m)

  case ('int_min')
    datum_value = 0; ix_m = -1
    call integrate_min (ix_ele, n_track, datum_value, ix_m, branch, vec_ptr, datum)
    call integrate_min (0, ix_start, datum_value, ix_m, branch, vec_ptr, datum)
    if (present(good)) valid_value = all(good(ix_ele:n_track)) .and. all(good(0:ix_start))

  case ('int_max')
    datum_value = 0; ix_m = -1
    call integrate_max (ix_ele, n_track, datum_value, ix_m, branch, vec_ptr, datum)
    call integrate_max (0, ix_start, datum_value, ix_m, branch, vec_ptr, datum)
    if (present(good)) valid_value = all(good(ix_ele:n_track)) .and. all(good(0:ix_start))

  case default
    call tao_set_invalid (datum, 'BAD MERIT_TYPE WHEN THERE IS A RANGE OF ELEMENTS: ' // datum%merit_type, why_invalid)
    valid_value = .false.
    return
  end select

  if (.not. valid_value) call tao_set_invalid (datum, 'INVALID DATA IN RANGE FROM ELE_START TO ELE_REF', why_invalid)

! no wrap case
else
  if (present(orbit)) then
    do i = ix_start, ix_ele
      call tao_scratch_values_calc (i, datum, branch, orbit)
    enddo
  endif

  select case (datum%merit_type)
  case ('min')
    ix_m = minloc (vec_ptr(ix_start:ix_ele), 1) + ix_start - 1
    datum_value = vec_ptr(ix_m)
    if (present(good)) valid_value = good(ix_m)

  case ('max')
    ix_m = maxloc (vec_ptr(ix_start:ix_ele), 1) + ix_start - 1
    datum_value = vec_ptr(ix_m)
    if (present(good)) valid_value = good(ix_m)

  case ('abs_min')
    ix_m = minloc (abs(vec_ptr(ix_start:ix_ele)), 1) + ix_start - 1
    datum_value = abs(vec_ptr(ix_m))
    if (present(good)) valid_value = good(ix_m)

  case ('abs_max')
    ix_m = maxloc (abs(vec_ptr(ix_start:ix_ele)), 1) + ix_start - 1
    datum_value = abs(vec_ptr(ix_m))
    if (present(good)) valid_value = good(ix_m)

  case ('int_min')
    datum_value = 0; ix_m = -1
    call integrate_min (ix_start, ix_ele, datum_value, ix_m, branch, vec_ptr, datum)
    if (present(good)) valid_value = all(good(ix_start:ix_ele))

  case ('int_max')
    datum_value = 0; ix_m = -1
    call integrate_max (ix_start, ix_ele, datum_value, ix_m, branch, vec_ptr, datum)
    if (present(good)) valid_value = all(good(ix_start:ix_ele))

  case default
    call tao_set_invalid (datum, 'BAD MERIT_TYPE WHEN THERE IS A RANGE OF ELEMENTS: ' // datum%merit_type, why_invalid)
    valid_value = .false.
    return
  end select

  if (.not. valid_value) call tao_set_invalid (datum, 'INVALID DATA IN RANGE FROM ELE_START TO ELE_REF', why_invalid)

endif

datum%ix_ele_merit = ix_m
if (ref_value /= 0) deallocate (vec_ptr)

!

end subroutine tao_load_this_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine integrate_min (ix_start, ix_ele, datum_value, ix_m, branch, vec, datum)

type (branch_struct) branch
type (tao_data_struct) datum

real(rp) vec(0:), val0, dv0, dv1, ds
real(rp) datum_value

integer ix_m, i, ix_start, ix_ele

!

val0 = datum%meas_value

do i = ix_start, ix_ele
  if (ix_m < 0) ix_m = i
  if (vec(i) < vec(ix_m)) ix_m = i
  if (i == ix_start) cycle
  dv0 = val0 - vec(i-1) 
  dv1 = val0 - vec(i)
  ds = branch%ele(i)%s - branch%ele(i)%s_start
  if (dv0 < 0 .and. dv1 < 0) cycle
  if (dv0 > 0 .and. dv1 > 0) then
    datum_value = datum_value + 0.5 * ds * (dv0 + dv1)
  elseif (dv0 > 0) then
    datum_value = datum_value + 0.5 * ds * dv0**2 / (dv0 - dv1)
  else
    datum_value = datum_value + 0.5 * ds * dv1**2 / (dv1 - dv0)
  endif
enddo

end subroutine integrate_min

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine integrate_max (ix_start, ix_ele, datum_value, ix_m, branch, vec, datum)

type (branch_struct) branch
type (tao_data_struct) datum

real(rp) vec(0:), val0, dv0, dv1, ds
real(rp) datum_value

integer ix_m, i, ix_start, ix_ele

!

val0 = datum%meas_value

do i = ix_start, ix_ele
  if (ix_m < 0) ix_m = i
  if (vec(i) > vec(ix_m)) ix_m = i
  if (i == ix_start) cycle
  dv0 = vec(i-1) - val0
  dv1 = vec(i) - val0
  ds = branch%ele(i)%s - branch%ele(i)%s_start
  if (dv0 < 0 .and. dv1 < 0) cycle
  if (dv0 > 0 .and. dv1 > 0) then
    datum_value = datum_value + 0.5 * ds * (dv0 + dv1)
  elseif (dv0 > 0) then
    datum_value = datum_value + 0.5 * ds * dv0**2 / (dv0 - dv1)
  else
    datum_value = datum_value + 0.5 * ds * dv1**2 / (dv1 - dv0)
  endif
enddo

end subroutine integrate_max

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine tao_scratch_values_calc (ix_ele, datum, branch, orbit)

type (ele_struct), pointer :: ele
type (this_array_struct), pointer :: cc_p
type (tao_data_struct) datum
type (branch_struct) branch
type (coord_struct) orbit(0:)

integer ix_ele, ip
real(rp) f, f1, f2, gam_c, E_ratio, a_emit, b_emit, z_emit, px_b1, px_b2, py_a1, py_a2

!

cc_p => scratch%cc(ix_ele)
ele => branch%ele(ix_ele)
ip = index(datum%data_type, '.')
if (ip == 0) return

select case (datum%data_type(1:ip))

case ('k.', 'cbar.', 'ping_a.', 'ping_b.')
  if (cc_p%coupling_calc_done) return
  cc_p%coupling_calc_done = .true.

  call c_to_cbar (ele, cc_p%cbar)
  f = sqrt(ele%a%beta/ele%b%beta) 
  f1 = f / ele%gamma_c
  f2 = 1 / (f * ele%gamma_c)

  cc_p%k_11a = cc_p%cbar(1,1) * f1
  cc_p%k_12a = cc_p%cbar(1,2) * f2
  cc_p%k_12b = cc_p%cbar(1,2) * f1
  cc_p%k_22b = cc_p%cbar(2,2) * f2

! Amplitude calc
! 'orbit.amp_a', 'orbit.amp_b', 'orbit.norm_amp_a', 'orbit.norm_amp_b'

case ('orbit.')
  if (index(datum%data_type, 'amp_') == 0) return
  if (cc_p%amp_calc_done) return
  cc_p%amp_calc_done = .true.

  call orbit_amplitude_calc (ele, orbit(ix_ele), cc_p%amp_a, cc_p%amp_b, cc_p%amp_na, cc_p%amp_nb)

end select

end subroutine tao_scratch_values_calc

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine tao_do_wire_scan (ele, wire_params, theta, beam) result (moment)
!
! Returns the beam's second moment using the wire along the specified angle.
! Keep in mind that the actual correlation axis is 90 degrees off of the 
! wire angle
!
! This simulates a fast wire scanner that performs the scan over only one
! bunch. Obviously, this isn't realistic. Any dynamic effects will not be
! accounted for!
!
! Input:
!  ele         -- Element_struct: 
!    %value(noise$) -- relative wire resolution RMS 
!    %value(tilt$)  -- wire angle error in radians rms.
!  theta       -- Real(rp): wire angle wrt x axis (in degrees)
!  beam        -- Beam_struct: contains the beam distribution
!
! Output:
!   moment  -- Real(rp): second moment along axis specified by angle.
!-

function tao_do_wire_scan (ele, theta, beam) result (moment)

use random_mod

type (ele_struct) ele
type (beam_struct) beam

real(rp), allocatable, save :: dist(:)
real(rp) theta, theta_rad, moment, ran_num(2)
real(rp) avg

!

call ran_gauss (ran_num)

! angle in radians and correlation angle is 90 off from wire angle
theta_rad = ele%value(tilt_tot$)+ (theta - 90) * (2.0*pi / 360.0)

if (.not. allocated (dist)) then
  allocate (dist(size(beam%bunch(1)%particle)))
elseif (size(dist) /= size(beam%bunch(1)%particle)) then
  deallocate (dist)
  allocate (dist(size(beam%bunch(1)%particle)))
endif

! Rotating the wire scanner is equivalent to rotating the beam by -theta

dist =  beam%bunch(1)%particle%vec(1) * cos(-theta_rad ) &
        + beam%bunch(1)%particle%vec(3) * sin(-theta_rad)
  
avg = sum (dist, mask = (beam%bunch(1)%particle%state == alive$)) &
          / count (beam%bunch(1)%particle%state == alive$)
        
moment = (1 + ele%value(noise$)*ran_num(2)) * sum ((dist-avg)*(dist-avg), &
                 mask = (beam%bunch(1)%particle%state == alive$)) &
          / count (beam%bunch(1)%particle%state == alive$)

end function tao_do_wire_scan

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Function tao_pointer_to_datum_ele (lat, ix_ele, datum, valid) result (ele)
! 
! Routine to see if an element index corresponds to an element with a definite 
! location such as an overlay or multipass element.
!
! If the element is a super_lord then the super_slave element at the exit end
! of the lord will be returned. Otherwise ix_loc will be set to ix_ele.
!
! Input:
!   lat    -- Lat_struct: Lattice
!   ix_ele -- Integer: Index of element.
!   datum  -- Tao_data_struct: Used for error messages and gives branch index.
!   valid  -- Logical: Set False if element does not have a definite location.
!               Set True otherwise
!
! Output:
!   ele   -- Ele_struct, pointer :: Pointer to the element. Set to NULL if not valid
!-

function tao_pointer_to_datum_ele (lat, ele_name, ix_ele, datum, valid, why_invalid) result (ele)

type (lat_struct) lat
type (tao_data_struct) datum
type (ele_struct), pointer :: ele

integer ix_ele, ixc, n_track, n_max

logical valid

character(*) ele_name
character(*), optional :: why_invalid
character(*), parameter :: r_name = 'tao_pointer_to_datum_ele'

! If ele_name is blank but ix_ele is not -1 then this datum was constructed internally
! by Tao (as opposed being constructed from tao init file info).
! If this is the case, assume that everything is OK and just return a pointer to the element.

valid = .true.
nullify (ele)

if (ele_name == '') then
  if (ix_ele /= -1) ele => pointer_to_ele (lat, ix_ele, datum%ix_branch)
  return
endif

! Here if ele_name is not blank.
! Do some checking...

n_track = lat%branch(datum%ix_branch)%n_ele_track
n_max   = lat%branch(datum%ix_branch)%n_ele_max

if (ix_ele < 0 .or. ix_ele > n_max) then
  call out_io (s_error$, r_name, 'ELEMENT INDEX OUT OF RANGE! \i5\ ', ix_ele)
  call tao_set_invalid (datum, 'ELEMENT INDEX OUT OF RANGE FOR: ' // tao_datum_name(datum), why_invalid)
  valid = .false.
  return
endif

ele => pointer_to_ele (lat, ix_ele, datum%ix_branch)

if (ix_ele <= n_track) then
  if (datum%eval_point == anchor_beginning$) ele => pointer_to_next_ele(ele, -1)
  return
endif

if (ele%lord_status == super_lord$ .or. ele%lord_status == overlay_lord$ .or. ele%lord_status == group_lord$) then
  if (datum%data_type(1:8) == 'rad_int1') return  ! Since rad_int1 is integraged radiation integral over the element.
  if (datum%eval_point == anchor_beginning$) then
    ele => pointer_to_slave(ele, 1)
    ele => pointer_to_next_ele(ele, -1)
  else
    ele => pointer_to_slave(ele, ele%n_slave)
  endif
  return
endif

valid = .false.
call out_io (s_error$, r_name, &
            'ELEMENT: ' // trim(lat%ele(ix_ele)%name) // &
            '    WHICH IS A: ' // control_name(lat%ele(ix_ele)%lord_status), &
            'CANNOT BE USED IN DEFINING A DATUM SINCE IT DOES NOT HAVE ', &
            '   A DEFINITE LOCATION IN THE LATTICE.', &
            'FOR DATUM: ' // tao_datum_name(datum) )
call tao_set_invalid (datum, 'NO DEFINITE LOCATION IN LATTICE FOR: ' // tao_datum_name(datum), why_invalid)

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_to_real (expression, value, err_flag, use_good_user, good, stack, &
!                             print_err, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele)
!
! Mathematically evaluates an expression.
!
! Input:
!   expression    -- character(*): arithmetic expression
!
! Output:
!   value        -- real(rp): Value of arithmetic expression.
!   err_flag     -- Logical: TRUE on error.
!-

subroutine tao_to_real (expression, value, err_flag)

character(*) :: expression

type (tao_expression_info_struct), allocatable, save :: info(:)

real(rp) value
real(rp), allocatable, save :: vec(:)

logical err_flag

!

call tao_evaluate_expression (expression, 1, .false., vec, info, err_flag)
if (err_flag) return
value = vec(1)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_expression (expression, n_size, use_good_user, &
!      value, info, err_flag, print_err, stack, dflt_component, dflt_source, &
!      dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, dflt_uni)
!
! Mathematically evaluates a character expression.
!
! Input:
!   expression     -- Character(*): Arithmetic expression.
!   n_size         -- Integer: Size of the value array. If the expression evaluates to a
!                      a scalar, each value in the value array will get this value.
!                      If n_size = 0, the natural size is determined by the expression itself.
!   use_good_user  -- Logical: Use the good_user logical in evaluating good(:)
!   print_err      -- Logical, optional: If False then supress evaluation error messages.
!                      This does not affect syntax error messages. Default is True.
!   dflt_component -- Character(*), optional: Component to use if not specified in the expression. 'model', 'base', or 'design'.
!   dflt_source    -- Character(*), optional: Default source ('lat', 'data', etc.). Default is ''.
!   dflt_ele_ref   -- Ele_struct, pointer, optional: Default reference element.
!   dflt_ele_start -- Ele_struct, pointer, optional: Default start element for ranges.
!   dflt_ele       -- Ele_struct, pointer, optional: Default element to evaluate at.
!   dflt_dat_or_var_index -- Character(*), optional: Default datum or variable index to use.
!   dflt_uni       -- Integer, optional: Default universe to use. If 0 or not present, use viewed universe.
!
! Output:
!   value(:)  -- Real(rp), allocatable: Value of arithmetic expression.
!   info(:)    -- tao_expression_info_struct, allocatable: Is the value valid?, etc.
!                  Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                  orbit.x[23]|good_user is False.
!   err_flag  -- Logical: True on an error. EG: Invalid expression.
!                  A divide by zero is not an error but good(:) will be set to False.
!   stack(:)  -- Tao_eval_stack1_struct, allocatable, optional: Evaluation stack for the
!                  expression. This is useful to save if the same expression is
!                  to be evaluated repeatedly. 
!                  With this, tao_evaluate_stack can be called directly.
!-

subroutine tao_evaluate_expression (expression, n_size, use_good_user, value, &
          info, err_flag, print_err, stack, dflt_component, dflt_source, &
          dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, dflt_uni)

use random_mod

type (tao_eval_stack1_struct), save :: stk(100)
type (tao_eval_stack1_struct), allocatable, optional :: stack(:)
type (ele_struct), optional, pointer :: dflt_ele_ref, dflt_ele_start, dflt_ele
type (tao_expression_info_struct), allocatable :: info(:)

integer, optional :: dflt_uni
integer i_lev, i_op, i, ios, n, n_size, n__size
integer op(200), ix_word, i_delim, i2, ix, ix_word2, ixb

real(rp), allocatable :: value(:)

character(*) :: expression
character(*), optional :: dflt_component, dflt_source
character(*), optional :: dflt_dat_or_var_index

character(:), allocatable :: phrase
character(1) delim
character(80) word, word2, default_source
character(*), parameter :: r_name = "tao_evaluate_expression"
character(40) saved_prefix

logical delim_found, split, ran_function_pending, use_good_user
logical err_flag, err, wild, printit
logical, optional :: print_err

! Don't destroy the input expression

err_flag = .true.
saved_prefix = ''
printit = logic_option(.true., print_err)
default_source = ''
if (present(dflt_source)) default_source = dflt_source

phrase = expression

! if phrase is blank then return 0.0

call string_trim (phrase, phrase, ios)
if (ios == 0) then
  call out_io (s_warn$, r_name, "Expression is blank")
  call re_allocate (value, max(1, n_size))
  value = 0.0
  return
endif
 
! General idea: Create a reverse polish stack that represents the expression.
! Reverse polish means that the operand goes last so that 2 * 3 is writen 
! on the stack as: [2, 3, *]

! The stack is called: stk
! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stk.

! init

i_lev = 0
i_op = 0
ran_function_pending = .false.

do i = 1, size(stk)
  stk(i)%name = ''
  if (allocated(stk(i)%info)) deallocate (stk(i)%info)
  if (allocated(stk(i)%value_ptr)) deallocate (stk(i)%value_ptr)
enddo

! parsing loop to build up the stack.

parsing_loop: do

  ! get a word

  call word_read (phrase, '+-*/()^,}[ ', word, ix_word, delim, delim_found, phrase)

  if (ran_function_pending .and. (ix_word /= 0 .or. delim /= ')')) then
    call out_io (s_warn$, r_name, 'RAN AND RAN_GAUSS DO NOT TAKE AN ARGUMENT')
    return
  endif

  !--------------------------
  ! Preliminary: If we have split up something that should have not been split
  ! then put it back together again...

  ! just make sure we are not chopping a number in two, e.g. "3.5e-7" should not
  ! get split at the "-" even though "-" is a delimiter

  split = .true.         ! assume initially that we have a split number
  if (ix_word == 0) then
    split = .false.
  elseif (word(ix_word:ix_word) /= 'E' .and. word(ix_word:ix_word) /= 'e' .and. &
          word(ix_word:ix_word) /= 'D' .and. word(ix_word:ix_word) /= 'd'  ) then
    split = .false.
  endif
  if (delim /= '-' .and. delim /= '+') split = .false.
  do i = 1, ix_word-1
    if (index('.0123456789', word(i:i)) == 0) split = .false.
  enddo

  ! If still SPLIT = .TRUE. then we need to unsplit

  if (split) then
    word = word(:ix_word) // delim
    call word_read (phrase, '+-*/()^,}', word2, ix_word2, delim, delim_found, phrase)
    word = word(:ix_word+1) // word2
    ix_word = ix_word + ix_word2
  endif

  ! Something like "lcav[lr(2).freq]" or "[2,4]@orbit.x[1,4] will get split on the "["

  do
    if (delim /= '[') exit

    call word_read (phrase, ']', word2, ix_word2, delim, delim_found, phrase)
    if (.not. delim_found) then
      call out_io (s_warn$, r_name, "NO MATCHING ']' FOR OPENING '[':" // expression)
      return
    endif
    word = word(:ix_word) // '[' // trim(word2) // ']'
    ix_word = ix_word + ix_word2 + 2
    if (phrase == ' ') then  
      delim_found = .false.
      delim = ' '
    elseif (phrase(1:1) == ' ') then  
      call string_trim (phrase, phrase, ix)
      if (index('+-*/()^,}[', phrase(1:1)) == 0) then
        delim = ' '
      else
        delim = phrase(1:1)
        phrase = phrase(2:)
      endif
    else          ! even more...
      call word_read (phrase, '[+-*/()^,}', word2, ix_word2, delim, delim_found, phrase)
      word = word(:ix_word) // trim(word2)       
      ix_word = ix_word + ix_word2 
    endif

  enddo

  ! If delim = "*" then see if this is being used as a wildcard
  ! Examples: "[*]|", "*.*|", "*.x|", "*@orbit.x|", "*@*|", "orbit.*[3]|"
  ! If so, we have split in the wrong place and we need to correct this.

  if (delim == '*') then

    wild = .false.

    select case (phrase(1:1))
    case ( ']', '[', '|', '@')
      wild = .true.
    case ('.')
      if (index('0123456789', phrase(2:2)) == 0) wild = .true.
    end select
    ixb = index(phrase, '|')
    if (ixb == 0) wild = .false.

    if (wild) then
      word = word(:ix_word) // '*' // phrase(1:ixb) 
      phrase = phrase(ixb+1:)
      call word_read (phrase, '+-*/()^,}', word2, ix_word2, delim, delim_found, phrase)
      word = trim(word) // trim(word2)       
      ix_word = len_trim(word)
    endif

  endif

  !---------------------------
  ! Now see what we got...

  ! For a word ending in '|' then must be a construct like 'orbit.x|-model'.
  ! So store the 'orbit.x|' prefix

  if (ix_word /= 0) then
    if (word(ix_word:ix_word) == '|') then
      saved_prefix = word
      word = ''
      ix_word = 0
    endif
  endif

  ! if the word is a datum without an "|", and dflt_component is present, 
  ! Then use the dflt_component.

  if (present(dflt_component) .and. index(word, '|') == 0) then
    call tao_find_data (err, word, print_err = .false.)
    if (.not. err) then
      phrase = trim(word) // '|' // trim(dflt_component) // delim // trim(phrase)
      cycle   ! Try again
    endif
  endif

  ! For a "(" delim we must have a function

  if (delim == '(') then

    ran_function_pending = .false.
    if (ix_word /= 0) then
      word2 = word
      call downcase_string (word2)
      select case (word2)
      case ('sin')
        call pushit (op, i_op, sin$)
      case ('sinc')
        call pushit (op, i_op, sinc$)
      case ('cos')
        call pushit (op, i_op, cos$)
      case ('tan') 
        call pushit (op, i_op, tan$)
      case ('asin') 
        call pushit (op, i_op, asin$)
      case ('acos') 
        call pushit (op, i_op, acos$)
      case ('atan') 
        call pushit (op, i_op, atan$)
      case ('atan2') 
        call pushit (op, i_op, atan2$)
      case ('abs') 
        call pushit (op, i_op, abs$)
      case ('sqrt') 
        call pushit (op, i_op, sqrt$)
      case ('log') 
        call pushit (op, i_op, log$)
      case ('exp') 
        call pushit (op, i_op, exp$)
      case ('factorial') 
        call pushit (op, i_op, factorial$)
      case ('ran') 
        call pushit (op, i_op, ran$)
        ran_function_pending = .true.
      case ('ran_gauss')
        call pushit (op, i_op, ran_gauss$)
        ran_function_pending = .true.
      case ('int')
        call pushit (op, i_op, int$)
      case ('nint')
        call pushit (op, i_op, nint$)
      case ('floor')
        call pushit (op, i_op, floor$)
      case ('ceiling')
        call pushit (op, i_op, ceiling$)
      case default
        call out_io (s_warn$, r_name, 'UNEXPECTED CHARACTERS BEFORE "(": ', 'IN EXPRESSION: ' // expression)
        return
      end select
    endif

    call pushit (op, i_op, l_parens$)
    cycle parsing_loop

  ! for a unary "-"

  elseif (delim == '-' .and. ix_word == 0) then
    call pushit (op, i_op, unary_minus$)
    cycle parsing_loop

  ! for a unary "+"

  elseif (delim == '+' .and. ix_word == 0) then
    call pushit (op, i_op, unary_plus$)
    cycle parsing_loop

  ! for a ")" delim

  elseif (delim == ')') then
    if (ix_word == 0) then
      if (.not. ran_function_pending) then
        call out_io (s_warn$, r_name, 'CONSTANT OR VARIABLE MISSING BEFORE ")"', 'IN EXPRESSION: ' // expression)
        return
      endif
    else
      call pushit (stk%type, i_lev, numeric$)
      call tao_param_value_routine (word, saved_prefix, stk(i_lev), err, printit, &
             dflt_component, default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, dflt_uni)
      if (err) then
        if (printit) call out_io (s_error$, r_name, &
                        'ERROR IN EVALUATING EXPRESSION: ' // expression, &
                        'CANNOT EVALUATE: ' // word)
        return
      endif
    endif

    ran_function_pending = .false.

    do
      do i = i_op, 1, -1       ! release pending ops
        if (op(i) == l_parens$) exit            ! break do loop
        call pushit (stk%type, i_lev, op(i))
      enddo

      if (i == 0) then
        call out_io (s_warn$, r_name, 'UNMATCHED ")" IN EXPRESSION: ' // expression)
        return
      endif

      i_op = i - 1

      call word_read (phrase, '+-*/()^,}', word, ix_word, delim, delim_found, phrase)
      if (ix_word /= 0) then
        call out_io (s_warn$, r_name, 'UNEXPECTED CHARACTERS AFTER ")" IN EXPRESSION: ' // expression)
        return
      endif

      if (delim /= ')') exit    ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      call out_io (s_warn$, r_name, '")(" CONSTRUCT DOES NOT MAKE SENSE IN EXPRESSION: ' // expression)
      return
    endif

  ! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      call out_io (s_warn$, r_name, 'CONSTANT OR VARIABLE MISSING IN EXPRESSION: ' // expression)
      return
    endif
    call pushit (stk%type, i_lev, numeric$)
    call tao_param_value_routine (word, saved_prefix, stk(i_lev), err, printit, &
            dflt_component, default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, dflt_uni)
    if (err) then
      if (printit) call out_io (s_error$, r_name, &
                        'ERROR IN EXPRESSION: ' // expression, &
                        'CANNOT EVALUATE: ' // word)
      return
    endif
  endif

  ! If we are here then we have an operation that is waiting to be identified

  select case (delim)
  case ('+')
    i_delim = plus$
  case ('-')
    i_delim = minus$
  case ('*')
    i_delim = times$
  case ('/')
    i_delim = divide$
  case (')')
    i_delim = r_parens$
  case ('^')
    i_delim = power$
  case (',', '}')
    i_delim = no_delim$
    call out_io (s_error$, r_name, &
                      'DELIMITOR FOUND OUT OF PLACE: ' // delim, &
                      'IN EXPRESSION: ' // expression)
      return
  case default
    if (delim_found) then
      if (delim == ' ') then
        call out_io (s_error$, r_name, 'MALFORMED EXPRESSION: ' // expression)
        return
      endif

      call out_io (s_error$, r_name, 'INTERNAL ERROR')
      call err_exit
    endif
    i_delim = no_delim$
  end select

  ! now see if there are operations on the OP stack that need to be transferred
  ! to the STK stack

  do i = i_op, 1, -1
    if (expression_eval_level(op(i)) >= expression_eval_level(i_delim)) then
      if (op(i) == l_parens$) then
        if (i > 1 .and. op(max(1,i-1)) == atan2$ .and. delim == ',') cycle parsing_loop
        call out_io (s_warn$, r_name, 'UNMATCHED "(" IN EXPRESSION: ' // expression)
        return
      endif
      call pushit (stk%type, i_lev, op(i))
    else
      exit
    endif
  enddo

  ! put the pending operation on the OP stack

  i_op = i
  if (i_delim == no_delim$) then
    exit parsing_loop
  else
    call pushit (op, i_op, i_delim)
  endif

enddo parsing_loop

!------------------------------------------------------------------
! Some error checks

if (i_op /= 0) then
  call out_io (s_warn$, r_name, 'UNMATCHED "(" IN EXPRESSION: ' // expression)
  return
endif

if (i_lev == 0) then
  call out_io (s_warn$, r_name, 'NO VALUE FOUND IN EXPRESSION: ' // expression)
  return
endif

n__size = 1
do i = 1, i_lev
  if (stk(i)%type < numeric$) cycle
  n = size(stk(i)%value)
  if (n == 1) cycle
  if (n__size == 1) n__size = n
  if (n /= n__size) then
    call out_io (s_warn$, r_name, 'ARRAY SIZE MISMATCH IN EXPRESSION: ' // expression)
    return
  endif
enddo

if (n_size /= 0) then
  if (n__size /= 1 .and. n_size /= n__size) then
    call out_io (s_warn$, r_name, 'ARRAY SIZE MISMATCH IN EXPRESSION: ' // expression)
    return
  endif
  n__size = n_size
endif

call tao_evaluate_stack (stk(1:i_lev), n__size, use_good_user, value, info, err_flag, printit)

! If the stack argument is present then copy stk to stack

if (present(stack)) then
  if (allocated(stack)) deallocate(stack)
  allocate (stack(i_lev))
  do i = 1, i_lev
    if (allocated (stk(i)%value)) then
      n = size(stk(i)%value)
      allocate (stack(i)%value(n), stack(i)%info(n))
      if (allocated (stack(i)%value_ptr)) allocate (stack(i)%value_ptr(n))
    endif
    stack(i) = stk(i)
  enddo
endif

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
contains

subroutine pushit (stack, i_lev, value)

integer stack(:), i_lev, value

character(6) :: r_name = "pushit"

!

i_lev = i_lev + 1

if (i_lev > size(stack)) then
  call out_io (s_warn$, r_name, 'STACK OVERFLOW.')
  call err_exit
endif

stack(i_lev) = value

end subroutine pushit
                       
end subroutine tao_evaluate_expression

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

recursive &
subroutine tao_param_value_routine (str, saved_prefix, stack, err_flag, print_err, &
      dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, dflt_uni)

type (tao_eval_stack1_struct) stack, stack2
type (tao_real_pointer_struct), allocatable :: re_array(:)
type (tao_data_array_struct), allocatable :: d_array(:)
type (tao_integer_array_struct), allocatable :: int_array(:)
type (tao_var_array_struct), allocatable :: v_array(:)
type (tao_data_struct) datum
type (tao_data_struct), pointer :: d
type (ele_struct), pointer, optional :: dflt_ele_ref, dflt_ele_start, dflt_ele

integer, optional :: dflt_uni
integer ios, i, m, n, ix, ix2, ix_word, ix_uni

character(*) str, saved_prefix
character(*), optional :: dflt_source, dflt_component
character(*), optional :: dflt_dat_or_var_index

character(1) delim
character(16) s_str, source
character(60) name, word2
character(200) str2
character(*), parameter :: r_name = 'tao_param_value_routine'

logical err_flag, print_err, print_error, delim_found

! pi, etc

err_flag = .false.

call match_word(str, physical_const_list%name, ix, .false., .false.)
if (ix > 0) then
  call re_allocate(stack%value, 1)
  stack%value(1) = physical_const_list(ix)%value
  return
endif

! An array "[...]"

if (str(1:1) == '[') then
  n = len_trim(str)
  if (str(n:n) /= ']') then
    if (print_err) call out_io (s_warn$, r_name, "Malformed array: " // str)
    err_flag = .true.
    return
  endif

  str2 = str(2:n-1)
  call re_allocate (stack%value, 0)

  do
    call word_read (str2, ',', word2, ix_word, delim, delim_found, str2, ignore_interior = .true.)
    call tao_param_value_routine (word2, saved_prefix, stack2, err_flag, print_err, &
      dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, dflt_uni)
    if (err_flag) return
    m = size(stack%value)
    n = size(stack2%value)
    call re_allocate(stack%value, m+n)
    stack%value(m+1:m+n) = stack2%value
    if (.not. delim_found) return
  enddo

endif

! Case where str represents a number.

if (is_real(str)) then
  call re_allocate(stack%value, 1)
  read (str, *, iostat = ios) stack%value(1)
  if (ios /= 0) then
    if (print_err) call out_io (s_warn$, r_name, "This doesn't seem to be a number: " // str)
    err_flag = .true.
  endif
  return
endif

! Case where str is a variable name.
! Remember the last string so 'orbit.x|meas-ref' translates to 'orbit.x|meas - orbit.x|ref.'

ix = index(str, '|')
stack%name = str

if (ix == 0) then
  select case (str)
  case ('model', 'design', 'base', 'meas', 'ref', 'old', 'fit')
    stack%name = trim(saved_prefix) // str
  end select
else
  saved_prefix = str(1:ix)
endif

name = stack%name

allocate (re_array(0))
allocate (d_array(0))
allocate (int_array(0))
allocate (v_array(0))

! Decide data source

source = dflt_source
ix = index(name, '::')
if (ix /= 0) then
  s_str = name(max(1,ix-4):ix-1)
  if (s_str(1:1) == '@') s_str = s_str(2:)
  if (s_str == 'lat' .or. s_str == 'ele' .or. s_str == 'data' .or. s_str == 'var') then
    source = s_str
  else if (name(max(1,ix-7):ix-1) == 'ele_mid') then
    source = 'ele'
  elseif (name(max(1, ix-4):ix-1) == 'beam') then
    source = 'beam'
  endif
endif


! Look for a lat datum.

if (source == 'lat' .or. source == 'beam') then
  call tao_evaluate_lat_or_beam_data (err_flag, name, stack%value, print_err, &
                                dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_component, dflt_uni)
  call tao_re_allocate_expression_info (stack%info, size(stack%value))
  stack%info%good = .not. err_flag
  stack%type = lat_num$
  return

! Look for a lattice element parameter 

elseif (source == 'ele') then
  call tao_evaluate_element_parameters (err_flag, name, stack%value, print_err, dflt_source, dflt_component)
  call tao_re_allocate_expression_info (stack%info, size(stack%value))
  stack%info%good = .not. err_flag
  stack%type = ele_num$
  return

! Look for variable or data values

else

  err_flag = .true.
  print_error = print_err
  if (source == '') print_error = .false. ! Don't generate unnecessary messages

  if (index(name, '|') == 0 .and. present(dflt_component)) name = trim(name) // '|' // trim(dflt_component)
  
  if (source == 'var' .or. source == '') then
    call tao_find_var (err_flag, name, v_array = v_array,  re_array = re_array, print_err = print_error, dflt_var_index = dflt_dat_or_var_index)
    stack%type = var_num$
  endif

  if (source == 'data' .or. (err_flag .and. source == '')) then
    call tao_find_data (err_flag, name, d_array = d_array, re_array = re_array, int_array = int_array, &
                        dflt_index = dflt_dat_or_var_index, print_err = print_error, ix_uni = dflt_uni)
    stack%type = data_num$
  endif

  if (err_flag) then
    if (source == '') call out_io (s_error$, r_name, 'CANNOT EVALUATE AS DATUM OR VARIABLE VALUE: ' // name)
    return
  endif

endif

! Now transfer the information to the stack

if (size(re_array) /= 0) then
  n = size(re_array)
  if (allocated(stack%value_ptr)) then
    if (size(stack%value_ptr) /= n) deallocate (stack%value_ptr)
  endif

  if (.not. allocated(stack%value_ptr)) allocate (stack%value_ptr(n))
  call re_allocate (stack%value, n)
  call tao_re_allocate_expression_info (stack%info, n)

  stack%scale =  1

  if (index(name, 'ping_a.') /= 0 .and. index(name, 'ping_a.phase') == 0) then
    ix_uni = d_array(1)%d%d1%d2%ix_uni
    if (index(name, '|meas') /= 0) then
      stack%scale = s%u(ix_uni)%ping_scale%a_mode_meas
    elseif (index(name, '|ref') /= 0) then
      stack%scale = s%u(ix_uni)%ping_scale%a_mode_ref
    endif
  endif

  if (index(name, 'ping_b.') /= 0 .and. index(name, 'ping_b.phase') == 0) then
    ix_uni = d_array(1)%d%d1%d2%ix_uni
    if (index(name, '|meas') /= 0) then
      stack%scale = s%u(ix_uni)%ping_scale%b_mode_meas
    elseif (index(name, '|ref') /= 0) then
      stack%scale = s%u(ix_uni)%ping_scale%b_mode_ref
    endif
  endif

  do i = 1, n
    stack%value(i) = re_array(i)%r
    stack%value_ptr(i)%r => re_array(i)%r
    ! good is only used with data and not variables
    if (associated(re_array(i)%good_value)) then
      stack%value_ptr(i)%good_user => re_array(i)%good_user
      stack%value_ptr(i)%good_value => re_array(i)%good_value
    else
      stack%info(i)%good = .true.
    endif

    select case (stack%type)
    case (var_num$)
      if (v_array(i)%v%exists) then
        stack%info(i)%s      = v_array(i)%v%s
        stack%info(i)%ix_ele = v_array(i)%v%slave(1)%ix_ele
        stack%info(i)%good   = v_array(i)%v%exists
        stack%value_ptr(i)%good_user => v_array(i)%v%good_user
      endif
    case (data_num$)
      d => d_array(i)%d
      if (d%exists) then
        stack%info(i)%s      = d%s
        stack%info(i)%ix_ele = d%ix_ele
        stack%info(i)%good   = d%exists
        stack%value_ptr(i)%good_user => d%good_user
      endif
    end select
  enddo

elseif (size(int_array) /= 0) then
  n = size(int_array)
  call re_allocate (stack%value, n)
  call tao_re_allocate_expression_info (stack%info, n)
  do i = 1, n
    stack%value(i) = int_array(i)%i
    stack%info(i)%good  = .true.
  enddo

else
  if (print_err) call out_io (s_warn$, r_name, &
               'THIS IS NOT A DATUM OR A VARIABLE VALUE: ' // name, &
               '[PERHAPS MISSING "|<component>" SUFFIX.]')
  err_flag = .true.
  return
endif

end subroutine tao_param_value_routine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_stack (stack, n_size, use_good_user, value, info, err_flag, print_err)
!
! Routine to evaluate an expression stack.
!
! Input:
!   stack(:)      -- Tao_eval_stack1_struct: Expression stack
!   n_size        -- Integer: Result array size.
!   use_good_user -- Logical: Use the good_user logical in evaluating good(:)
!   print_err     -- Logical: If False then supress evaluation error messages.
!                      This does not affect syntax error messages. Default is True.
!
! Output:
!   value(:)     -- Real(rp), allocatable: Value of arithmetic expression.
!   info(:)      -- tao_expression_info_struct, allocatable: Is the value valid? 
!                     Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                     orbit.x[23]|good_user is False.
!   err_flag     -- Logical: True on error. False otherwise
!-

subroutine tao_evaluate_stack (stack, n_size, use_good_user, value, info, err_flag, print_err)

type (tao_eval_stack1_struct), target :: stack(:)
type (tao_eval_stack1_struct), pointer :: s(:)
type (tao_eval_stack1_struct) stk2(20)
type (tao_expression_info_struct), allocatable :: info(:)

real(rp), allocatable :: value(:)

integer i, i2, j, n
integer n_size

logical err_flag, use_good_user, print_err

character(20) :: r_name = 'tao_evaluate_stack'

! Calculate good

s => stack   ! For debugging purposes
err_flag = .true.

call re_allocate (value, n_size)
call tao_re_allocate_expression_info (info, n_size)
info = tao_expression_info_struct()

do i = 1, size(stack)

  if (allocated(stack(i)%value_ptr)) then
    if (associated(stack(i)%value_ptr(1)%good_value)) then    
      do j = 1, size(stack(i)%value_ptr)
        if (use_good_user) then
          stack(i)%info(j)%good = stack(i)%value_ptr(j)%good_value .and. stack(i)%value_ptr(j)%good_user
        else
          stack(i)%info(j)%good = stack(i)%value_ptr(j)%good_value
        endif
      enddo
    endif
  endif

  if (.not. allocated(stack(i)%info)) cycle

  if (size(stack(i)%info) == 1) then
    info%good = info%good .and. stack(i)%info(1:1)%good
  else
   info%good = info%good .and. stack(i)%info%good
  endif

  if (size(stack(i)%info) == size(info)) then
    do j = 1, size(info)
      if (stack(i)%info(j)%s /= real_garbage$) info(j)%s = stack(i)%info(j)%s
      if (stack(i)%info(j)%ix_ele > -1) info(j)%ix_ele = stack(i)%info(j)%ix_ele
    enddo 
  endif

enddo

! Go through the stack and perform the operations...

i2 = 0  ! stack pointer
do i = 1, size(stack)

  select case (stack(i)%type)
  case (numeric$) 
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (lat_num$, ele_num$)
    !!! This needs to be fixed to include default stuff
    !!! call tao_param_value_routine (stack(i)%name, '', stack(i), err_flag, print_err)
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (var_num$, data_num$)
    do j = 1, size(stack(i)%value)
      stack(i)%value(j) = stack(i)%value_ptr(j)%r
    enddo
    stack(i)%value = stack(i)%value * stack(i)%scale
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (unary_minus$) 
    stk2(i2)%value = -stk2(i2)%value

  case (unary_plus$) 
    ! Nothing to do

  case (plus$) 
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value + stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) + stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value + stk2(i2)%value
    endif
    i2 = i2 - 1

  case (minus$) 
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value - stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) - stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value - stk2(i2)%value
    endif
    i2 = i2 - 1

  case (times$) 
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value * stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) * stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value * stk2(i2)%value
    endif
    i2 = i2 - 1

  case (divide$) 
    n = size(stk2(i2)%value)
    do j = 1, n
      if (stk2(i2)%value(j) == 0) then  ! Divide by 0 error!
        stk2(i2)%value(j) = 1
        if (n == 1) then
          info%good = .false.  ! All are false
        else
          info(j)%good = .false.
        endif
      endif
    enddo

    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value / stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) / stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value / stk2(i2)%value
    endif
    i2 = i2 - 1

  case (power$) 
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value ** stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) ** stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value ** stk2(i2)%value
    endif
    i2 = i2 - 1

  case (sin$) 
    stk2(i2)%value = sin(stk2(i2)%value)

  case (sinc$) 
    stk2(i2)%value = sinc(stk2(i2)%value)

  case (cos$) 
    stk2(i2)%value = cos(stk2(i2)%value)

  case (tan$) 
    stk2(i2)%value = tan(stk2(i2)%value)

  case (asin$) 
    stk2(i2)%value = asin(stk2(i2)%value)

  case (acos$) 
    stk2(i2)%value = acos(stk2(i2)%value)

  case (atan$) 
    stk2(i2)%value = atan(stk2(i2)%value)

  case (atan2$) 
    stk2(i2-1)%value = atan2(stk2(i2-1)%value, stk2(i2)%value)
    i2 = i2 - 1

  case (abs$) 
    stk2(i2)%value = abs(stk2(i2)%value)

  case (sqrt$) 
    stk2(i2)%value = sqrt(stk2(i2)%value)

  case (log$) 
    stk2(i2)%value = log(stk2(i2)%value)

  case (exp$) 
    stk2(i2)%value = exp(stk2(i2)%value)

  case (factorial$) 
    do n = 1, size(stk2(i2)%value)
      stk2(i2)%value(n) = factorial(nint(stk2(i2)%value(n)))
    enddo

  case (ran$) 
    i2 = i2 + 1
    call re_allocate(stk2(i2)%value, n_size)
    call ran_uniform(stk2(i2)%value)

  case (ran_gauss$) 
    i2 = i2 + 1
    call re_allocate(stk2(i2)%value, n_size)
    call ran_gauss(stk2(i2)%value)

  case (int$)
    stk2(i2)%value = int(stk2(i2)%value)

  case (nint$)
    stk2(i2)%value = nint(stk2(i2)%value)

  case (floor$)
    stk2(i2)%value = floor(stk2(i2)%value)

  case (ceiling$)
    stk2(i2)%value = ceiling(stk2(i2)%value)

  case default
    call out_io (s_warn$, r_name, 'INTERNAL ERROR')
    call err_exit
  end select
enddo

if (i2 /= 1) then
  call out_io (s_warn$, r_name, 'INTERNAL ERROR')
  call err_exit
endif

if (size(stk2(1)%value) == 1 .and. n_size > 1) then
  value = stk2(1)%value(1)
else
  call value_transfer (value, stk2(1)%value)
endif

where (.not. info%good) value = 0

err_flag = .false.

!-------------------------------------------------------------------------
contains

subroutine value_transfer (to_array, from_array)

real(rp), allocatable :: to_array(:)
real(rp) from_array(:)

!

call re_allocate (to_array, size(from_array))
to_array = from_array

end subroutine

end subroutine tao_evaluate_stack 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_int (str, i_int, err)
! 
! Converts a string to an integer
!
! If the string str is blank then i_int = 0
!-

subroutine tao_to_int (str, i_int, err)

character(*) str
integer ios, i_int
logical err
character(12) :: r_name = "tao_to_int"

!

call string_trim (str, str, ios)
if (ios .eq. 0) then
  i_int = 0
  return
endif
 
err = .false.
read (str, *, iostat = ios) i_int

if (ios /= 0) then
  call out_io (s_error$, r_name, 'EXPECTING INTEGER: ' // str)
  err = .true.
  return
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_set_invalid (datum, message, why_invalid)
!
! Routine to either print an error message to the terminal (if why_invalid
! is not present) or set the why_invalid string to the error message.
!
! Input:
!   datum       -- tao_data_struct: Bad datum.
!   message     -- character(*): Error message
!
! Output:
!   why_invalid -- character(*), optional: Set to message if present.
!-

Subroutine tao_set_invalid (datum, message, why_invalid)

type (tao_data_struct) datum
character(*) message
character(*), optional :: why_invalid
character(*), parameter :: r_name = 'tao_set_invalid'

!

if (present(why_invalid)) then
  why_invalid = message
else
  call out_io (s_error$, r_name, message, 'FOR DATUM: ' // tao_datum_name(datum))
endif

datum%exists = .false. ! So error message does not get generated again.

end subroutine tao_set_invalid

end module tao_data_and_eval_mod
