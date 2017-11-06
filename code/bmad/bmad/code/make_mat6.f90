!+
! Subroutine make_mat6 (ele, param, start_orb, end_orb, err_flag)
!
! Subroutine to make the 6x6 1st order transfer matrix for an element 
! along with the 0th order transfer vector.
!
! Note: Radiation fluctuations (but not damping) is turned off for the calculation.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element holding the transfer matrix.
!   param     -- lat_param_struct: Lattice global parameters.
!   start_orb -- Coord_struct, optional: Reference coordinates at the beginning of element. 
!                  If not present, default is to use the zero orbit.
!
! Output:
!   ele      -- Ele_struct: Element
!     %mat6    -- Real(rp): 1st order 6x6 transfer matrix.
!     %vec0    -- Real(rp): 0th order transfer vector.
!   end_orb  -- Coord_struct, optional: Reference coordinates at the end of element.
!   err_flag -- Logical, optional: Set True if there is an error. False otherwise.
!-

recursive subroutine make_mat6 (ele, param, start_orb, end_orb, err_flag)

use symp_lie_mod, only: symp_lie_bmad
use bookkeeper_mod, only: attribute_bookkeeper
use mad_mod, only: make_mat6_mad
use space_charge_mod, except_dummy => make_mat6
use equality_mod, only: operator(==)

implicit none

type (ele_struct), target :: ele
type (coord_struct), optional :: start_orb, end_orb
type (lat_param_struct)  param
type (coord_struct) a_start_orb, a_end_orb

real(rp), parameter :: zero_vec(6) = 0
integer mat6_calc_method, species

logical, optional :: err_flag
logical rad_fluct_save, err, finished

character(*), parameter :: r_name = 'make_mat6'

!--------------------------------------------------------
! Some init.
! If start_orb is in its not_set state (can happen if a particle is lost in 
! tracking and ele is downstream from the loss point), init the orbit to zero.

if (present(err_flag)) err_flag = .true.

if (.not. present(start_orb)) then
  call init_coord (a_start_orb, zero_vec, ele, upstream_end$, default_tracking_species(param))
else if (start_orb%state == not_set$ .or. start_orb%p0c /= ele%value(p0c_start$)) then
  call init_coord(a_start_orb, start_orb, ele, upstream_end$, default_tracking_species(param))
else
  a_start_orb = start_orb
endif

if (a_start_orb%direction == -1) then  ! Can only happen if start_orb is present
  call out_io (s_fatal$, r_name, 'TRANSFER MATRIX CALCULATION NOT ALLOWED WITH BACKWARD TRACKING.')
  if (global_com%exit_on_error) call err_exit
  return
endif

! init

if (bmad_com%auto_bookkeeper) call attribute_bookkeeper (ele, param)

mat6_calc_method = ele%mat6_calc_method
if (.not. ele%is_on) mat6_calc_method = bmad_standard$

if (any(ele%map_ref_orb_in%vec /= a_start_orb%vec)) then
  if (associated(ele%rad_int_cache)) ele%rad_int_cache%stale = .true.
endif

ele%map_ref_orb_in = a_start_orb

rad_fluct_save = bmad_com%radiation_fluctuations_on
bmad_com%radiation_fluctuations_on = .false.

! Compute matrix

err = .false.

select case (mat6_calc_method)

case (custom$)
  call make_mat6_custom (ele, param, a_start_orb, a_end_orb, err)

case (taylor$)
  call make_mat6_taylor (ele, param, a_start_orb, a_end_orb, err)

case (bmad_standard$)
  if (a_start_orb%species == photon$) then
    call make_mat6_bmad_photon (ele, param, a_start_orb, a_end_orb, err)
  else
    call make_mat6_bmad (ele, param, a_start_orb, a_end_orb, err)
  endif

case (symp_lie_ptc$)
  call make_mat6_symp_lie_ptc (ele, param, a_start_orb, a_end_orb)

case (symp_lie_bmad$)
  call symp_lie_bmad (ele, param, a_start_orb, a_end_orb, .true.)

case (tracking$)
  call make_mat6_tracking (ele, param, a_start_orb, a_end_orb)

case (mad$)
  call make_mat6_mad (ele, param, a_start_orb, a_end_orb)

! Static is used with hybrid elements since, in this case, the transfer map is not recomputable.

case (static$)
  if (present(err_flag)) err_flag = .false.
  if (present(end_orb)) call track1 (a_end_orb, ele, param, end_orb)
  if (ele%bookkeeping_state%mat6 == stale$) ele%bookkeeping_state%mat6 = ok$
  return

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN MAT6_CALC_METHOD: ' // mat6_calc_method_name(ele%mat6_calc_method))
  if (global_com%exit_on_error) call err_exit
  return
end select

if (err) then
  if (present(end_orb)) end_orb = a_end_orb
  return
endif

! Add space charge effects

if (bmad_com%space_charge_on) call make_mat6_ultra_rel_space_charge (ele, param)

! symplectify if wanted

if (ele%symplectify) call mat_symplectify (ele%mat6, ele%mat6, ele%value(p0c$)/ele%value(p0c_start$))

! Finish up

ele%map_ref_orb_out = a_end_orb

if (present(end_orb)) then
  end_orb = a_end_orb

  if (end_orb%state /= alive$) then
    end_orb%location = inside$
  elseif (end_orb%direction == 1) then
    end_orb%location = downstream_end$
  else
    end_orb%location = upstream_end$
  endif
endif

bmad_com%radiation_fluctuations_on = rad_fluct_save

if (ele%bookkeeping_state%mat6 == stale$) ele%bookkeeping_state%mat6 = ok$
if (present(err_flag)) err_flag = .false.

end subroutine

