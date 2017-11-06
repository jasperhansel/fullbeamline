module symp_lie_mod

use bmad_struct
use bmad_interface
use make_mat6_mod
use em_field_mod
use random_mod
use track1_mod

type wiggler_computations_struct
  real(rp) :: c_x = 0, s_x = 0, c_y = 0, s_y = 0, c_z = 0, s_z = 0
  real(rp) :: coef_Ax = 0, coef_Ay = 0, coef_Az = 0
  real(rp) :: sx_over_kx = 0, sy_over_ky = 0, sz_over_kz = 0
  real(rp) :: integral_sx = 0, integral_sy = 0, integral_cx = 0, integral_cy = 0
  integer :: trig_x = 0, trig_y = 0, family= 0
end type

private wiggler_computations_struct

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine symp_lie_bmad (ele, param, start_orb, end_orb, make_matrix, track, offset_ele)
!
! Subroutine to track through an element (which gives the 0th order map) 
! and optionally make the 6x6 transfer matrix (1st order map) as well.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele          -- Ele_struct: Element with transfer matrix
!   param        -- lat_param_struct: Parameters are needed for some elements.
!   start_orb    -- Coord_struct: Coordinates at the beginning of element. 
!   make_matrix  -- Logical: If True then make the 6x6 transfer matrix.
!   offset_ele   -- Logical, optional: Offset the element using ele%value(x_offset$), etc.
!                     Default is True.
!
! Output:
!   ele        -- Ele_struct: Element with transfer matrix.
!     %mat6(6,6)  -- 6x6 transfer matrix.
!     %vec0(6)    -- 0th order part of the transfer matrix.
!   end_orb    -- Coord_struct: Coordinates at the end of element.
!   track      -- Track_struct, optional: Structure holding the track information.
!                   When tracking through multiple elements, the trajectory in an element
!                   is appended to the existing trajectory. To reset: Set track%n_pt = -1.
!-

subroutine symp_lie_bmad (ele, param, start_orb, end_orb, make_matrix, track, offset_ele)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: field_ele
type (ele_pointer_struct), allocatable :: field_eles(:)
type (coord_struct) :: start_orb, end_orb, start_orb_saved
type (lat_param_struct)  param
type (track_struct), optional :: track
type (cartesian_map_struct), pointer :: ct_map
type (cartesian_map_term1_struct), pointer :: wig_term(:)
type (cartesian_map_term1_struct), pointer :: wt
type (wiggler_computations_struct), allocatable, target :: tm(:)

real(rp) rel_E, rel_E2, rel_E3, ds, ds2, s, m6(6,6), kmat6(6,6), Ax_saved, Ay_saved
real(rp) g_x, g_y, k1, k1_norm, k1_skew, x_q, y_q, ks_tot_2, ks, dks_ds, z_patch
real(rp), parameter :: z0 = 0, z1 = 1
real(rp) gamma_0, fact_d, fact_f, this_ran, g2, g3, ddAz__dx_dy
real(rp) dE_p, dpx, dpy, mc2, z_offset, orient_dir, rel_tracking_charge, charge_dir
real(rp), parameter :: rad_fluct_const = 55 * classical_radius_factor * h_bar_planck * c_light / (24 * sqrt_3)
real(rp), allocatable :: dz_offset(:)
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer i, n_step, n_field, key
integer ix_pole_max, ix_elec_max

integer num_wig_terms  ! number of wiggler terms

logical make_matrix, calculate_mat6, err, do_offset
logical, optional :: offset_ele

character(16) :: r_name = 'symp_lie_bmad'

! init

calculate_mat6 = (make_matrix .or. synch_rad_com%i_calc_on)

start_orb_saved = start_orb
end_orb = start_orb
if (end_orb%direction == 1) then
  end_orb%s = ele%s_start
  s = 0   ! longitudianal position
else
  end_orb%s = ele%s
  s = ele%value(l$)
endif

do_offset = logic_option (.true., offset_ele)
rel_E = (1 + end_orb%vec(6))
rel_E2 = rel_E**2
rel_E3 = rel_E**3

if (calculate_mat6) then
  call mat_make_unit(ele%mat6)
endif

err = .false.

orient_dir = ele%orientation * end_orb%direction
rel_tracking_charge = rel_tracking_charge_to_mass(end_orb, param)
charge_dir = rel_tracking_charge * orient_dir

! element offset 

call multipole_ele_to_ab (ele, .false., ix_pole_max, an,      bn,      magnetic$, include_kicks = .true.)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

if (do_offset) call offset_particle (ele, param, set$, end_orb, mat6 = ele%mat6, make_matrix = calculate_mat6)

if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, end_orb, magnetic$, 1.0_rp/2,   ele%mat6, calculate_mat6)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, end_orb, electric$, ele%value(l$)/2, ele%mat6, calculate_mat6)

! init

call compute_even_steps (ele%value(ds_step$), ele%value(l$), bmad_com%default_ds_step, ds, n_step)
ds2 = ds / 2

! radiation damping and fluctuations...
! The same kick is applied over the entire wiggler to save time.

if ((bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on)) then

  mc2 = mass_of(param%particle)
  gamma_0 = ele%value(e_tot$) / mass_of(param%particle)

  fact_d = 0
  if (bmad_com%radiation_damping_on) fact_d = 2 * classical_radius_factor * gamma_0**3 * ds / (3 * mc2)

  fact_f = 0
  if (bmad_com%radiation_fluctuations_on) then
    fact_f = sqrt(rad_fluct_const * ds * gamma_0**5) / mc2 
  endif

endif

!------------------------------------------------------------------
! select the element

key = ele%key
if (associated(ele%cartesian_map)) key = wiggler$

select case (ele%key)

!------------------------------------------------------------------
! Wiggler

Case (wiggler$, undulator$)

  call get_field_ele_list (ele, field_eles, dz_offset, n_field)
  do i = 1, n_field
    field_ele => field_eles(i)%ele
    if (field_ele%key == ele%key) exit
  enddo

  if (size(field_ele%cartesian_map) > 1) then
    call out_io (s_fatal$, r_name, 'ELEMENT: ' // ele%name, 'HAS MULTIPLE CARTESIAN MAPS. SYMP_LIE_PTC CAN ONLY HANDLE ONE.')
    if (global_com%exit_on_error) call err_exit
    return
  endif

  ct_map => field_ele%cartesian_map(1)
  if (ct_map%field_type /= magnetic$) then
    call out_io (s_fatal$, r_name, 'ELEMENT: ' // ele%name, 'HAS A CARTESIAN MAP THAT IS NOT MAGNETIC. SYMP_LIE_PTC CAN NOT HANDLE THIS.')
    if (global_com%exit_on_error) call err_exit
    return
  endif

  wig_term => ct_map%ptr%term
  z_offset = dz_offset(i)

  num_wig_terms = size(wig_term)

  allocate (tm(num_wig_terms))

  call calc_wig_coefs (calculate_mat6)
  call update_wig_y_terms (err); if (err) return
  call update_wig_x_terms (err); if (err) return
  call update_wig_s_terms

  if (present(track)) then
    call save_this_track_pt ()
  endif

  ! Symplectic integration gives a kick to the physical momentum (but not the canonical momentum) at the ends of an 
  ! element if there is a finite vector potential. This is necessary for symplecticity but can complicate comparisons
  ! with runge_kutta. That is, only set bmad_com%convert_to_kinetic_momentum = True for testing purposes.

  if (bmad_com%convert_to_kinetic_momentum) then
    end_orb%vec(2) = end_orb%vec(2) + Ax()
    end_orb%vec(4) = end_orb%vec(4) + Ay()
  endif

  ! loop over all steps

  do i = 1, n_step

    ! s half step

    s = s + ds2 * end_orb%direction
    call update_wig_s_terms

    ! Drift_1 = (P_x - Ax)^2 / (2 * (1 + dE))

    call apply_wig_exp_int_ax (-1, calculate_mat6)

    call apply_p_x (calculate_mat6)
    call update_wig_x_terms (err); if (err) return

    call apply_wig_exp_int_ax (+1, calculate_mat6)

    ! Drift_2 = (P_y - Ay)^2 / (2 * (1 + dE))

    call apply_wig_exp_int_ay (-1, calculate_mat6)

    call apply_p_y (calculate_mat6)
    call update_wig_y_terms (err); if (err) return

    call apply_wig_exp_int_ay (+1, calculate_mat6)

    ! Kick = Az

    dpx = dAz_dx()
    dpy = dAz_dy()
    end_orb%vec(2) = end_orb%vec(2) + ds * dpx
    end_orb%vec(4) = end_orb%vec(4) + ds * dpy

    call radiation_kick()

    if (calculate_mat6) then
      ddAz__dx_dy = ddAz_dx_dy()
      ele%mat6(2,1:6) = ele%mat6(2,1:6) + ds * (ddAz_dx_dx() * ele%mat6(1,1:6) + ddAz__dx_dy  * ele%mat6(3,1:6))
      ele%mat6(4,1:6) = ele%mat6(4,1:6) + ds * (ddAz__dx_dy  * ele%mat6(1,1:6) + ddAz_dy_dy() * ele%mat6(3,1:6))
    endif 

    ! Drift_2

    call apply_wig_exp_int_ay (-1, calculate_mat6)

    call apply_p_y (calculate_mat6)
    call update_wig_y_terms (err); if (err) return

    call apply_wig_exp_int_ay (+1, calculate_mat6)

    ! Drift_1

    call apply_wig_exp_int_ax (-1, calculate_mat6)

    call apply_p_x (calculate_mat6)
    call update_wig_x_terms (err); if (err) return

    call apply_wig_exp_int_ax (+1, calculate_mat6)

    ! s half step

    s = s + ds2 * end_orb%direction
    call update_wig_s_terms

    if (present(track)) call save_this_track_pt ()

  enddo

  ! Correction for finite vector potential at exit end.

  if (bmad_com%convert_to_kinetic_momentum) then
    end_orb%vec(2) = end_orb%vec(2) - Ax()
    end_orb%vec(4) = end_orb%vec(4) - Ay()
  endif

!----------------------------------------------------------------------------
! rf cavity

case (lcavity$, rfcavity$)

  ! loop over all steps

  do i = 1, n_step

    ! s half step

    s = s + ds2 * end_orb%direction
!    call rf_drift1 (calculate_mat6)
!    call rf_drift2 (calculate_mat6)
!    call rf_kick (calculate_mat6)
    call radiation_kick()
!    call rf_drift2 (calculate_mat6)
!    call rf_drift1 (calculate_mat6)

    s = s + ds2 * end_orb%direction

    if (present(track)) call save_this_track_pt ()

  enddo

!----------------------------------------------------------------------------
! solenoid, quadrupole, sol_quad, or bend_sol_quad
! Notice that we don't have to worry about a solenoid hard edge since we are
! using canonical coords so the hard edge is automatically included.

case (bend_sol_quad$, solenoid$, quadrupole$, sol_quad$)

  g_x = 0
  g_y = 0
  x_q = 0
  y_q = 0
  dks_ds = 0
  k1_norm = 0
  k1_skew = 0
  ks = 0

  select case (ele%key)
  case (bend_sol_quad$)
    g_x = ele%value(g$) * cos (ele%value(bend_tilt$)) * charge_dir
    g_y = ele%value(g$) * sin (ele%value(bend_tilt$)) * charge_dir
    k1_norm = ele%value(k1$) * cos (2 * ele%value(quad_tilt$)) * charge_dir
    k1_skew = ele%value(k1$) * sin (2 * ele%value(quad_tilt$)) * charge_dir
    x_q = ele%value(x_quad$)
    y_q = ele%value(y_quad$)
    ks = ele%value(ks$) * rel_tracking_charge
    dks_ds = ele%value(dks_ds$) * charge_dir
  case (solenoid$)
    ks = ele%value(ks$) * rel_tracking_charge
  case (quadrupole$)
    k1_norm = ele%value(k1$) * charge_dir 
  case (sol_quad$)
    k1_norm = ele%value(k1$) * charge_dir
    ks = ele%value(ks$) * rel_tracking_charge
  end select

  call hard_multipole_edge_kick (ele, param, first_track_edge$, end_orb, ele%mat6, calculate_mat6)
  call soft_quadrupole_edge_kick (ele, param, first_track_edge$, end_orb, ele%mat6, calculate_mat6)

  ! loop over all steps

  do i = 1, n_step

    s = s + ds2 * end_orb%direction
    ks_tot_2 = (ks + dks_ds * s) / 2

    call bsq_drift1 (calculate_mat6)
    call bsq_drift2 (calculate_mat6)
    call bsq_kick (calculate_mat6)
    call radiation_kick()
    call bsq_drift2 (calculate_mat6)
    call bsq_drift1 (calculate_mat6)

    s = s + ds2 * end_orb%direction
    ks_tot_2 = (ks + dks_ds * s) / 2

    if (present(track)) call save_this_track_pt ()

  enddo

  call soft_quadrupole_edge_kick (ele, param, second_track_edge$, end_orb, ele%mat6, calculate_mat6)
  call hard_multipole_edge_kick (ele, param, second_track_edge$, end_orb, ele%mat6, calculate_mat6)

!----------------------------------------------------------------------------
! unknown element

case default

  call out_io (s_fatal$, r_name, 'TRACKING NOT YET IMPLEMENTED FOR: ' // key_name(ele%key), &
                                 'FOR ELEMENT: ', ele%name)

end select

!----------------------------------------------------------------------------
! element offset

if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, end_orb, magnetic$, 1.0_rp/2,        ele%mat6, calculate_mat6)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, end_orb, electric$, ele%value(l$)/2, ele%mat6, calculate_mat6)

if (do_offset) call offset_particle (ele, param, unset$, end_orb, mat6 = ele%mat6, make_matrix = calculate_mat6)

! Correct z-position for wigglers, etc. 

z_patch = ele%value(delta_ref_time$) * c_light * end_orb%beta - ele%value(l$)
end_orb%vec(5) = end_orb%vec(5) + z_patch
end_orb%t = start_orb_saved%t + ele%value(delta_ref_time$) + (start_orb_saved%vec(5) - end_orb%vec(5)) / (end_orb%beta * c_light)

if (calculate_mat6) then
  ele%mat6(5,6) = ele%mat6(5,6) + ele%value(delta_ref_time$) * c_light * &
                (mass_of(param%particle)/ele%value(p0c$))**2 * (end_orb%beta / (1 + end_orb%vec(6)))**3
endif

! calc vec0

if (calculate_mat6) then
  ele%vec0(1:5) = end_orb%vec(1:5) - matmul (ele%mat6(1:5,1:6), start_orb_saved%vec)
  ele%vec0(6) = 0
endif

!

if (end_orb%direction == 1) then
  end_orb%s = ele%s
else
  end_orb%s = ele%s_start
endif

end_orb%p0c = ele%value(p0c$)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
contains

subroutine err_set (err, plane)

logical err
integer plane

!

call out_io (s_warn$, r_name, &
                'FLOATING OVERFLOW IN TRACKING:' // ele%name, &
                'PARTICLE WILL BE TAGGED AS LOST.')

if (plane == x_plane$) then
  end_orb%vec(1) = sign(2 * bmad_com%max_aperture_limit, end_orb%vec(1))
  if (end_orb%vec(1) > 0) then ; end_orb%state = lost_pos_x_aperture$
  else;                          end_orb%state = lost_neg_x_aperture$
  endif
else
  end_orb%vec(3) = sign(2 * bmad_com%max_aperture_limit, end_orb%vec(3))
  if (end_orb%vec(3) > 0) then ; end_orb%state = lost_pos_y_aperture$
  else;                          end_orb%state = lost_neg_y_aperture$
  endif
endif

err = .true.

end subroutine err_set

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine save_this_track_pt ()

integer ix

!

call save_a_step (track, ele, param, .true., end_orb, s)

if (calculate_mat6) then
  ix = track%n_pt
  track%map(ix)%mat6 = ele%mat6
  if (ele%value(tilt_tot$) /= 0) call tilt_mat6 (track%map(ix)%mat6, ele%value(tilt_tot$))
  call mat6_add_pitch (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%orientation, track%map(ix)%mat6)
  track%map(ix)%vec0 = track%orb(ix)%vec - matmul (track%map(ix)%mat6, start_orb_saved%vec)
endif

if (bmad_com%convert_to_kinetic_momentum) then
  ix = track%n_pt
  track%orb(ix)%vec(2) = track%orb(ix)%vec(2) - Ax()
  track%orb(ix)%vec(4) = track%orb(ix)%vec(4) - Ay()
endif

end subroutine save_this_track_pt

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_p_x (do_mat6)

logical do_mat6

end_orb%vec(1) = end_orb%vec(1) + ds2 * end_orb%vec(2) / rel_E
end_orb%vec(5) = end_orb%vec(5) - ds2 * end_orb%vec(2)**2 / (2*rel_E2)
end_orb%s = end_orb%s + ds2 * end_orb%direction

if (do_mat6) then
  ele%mat6(1,1:6) = ele%mat6(1,1:6) + (ds2 / rel_E)               * ele%mat6(2,1:6) - (ds2*end_orb%vec(2)/rel_E2)    * ele%mat6(6,1:6) 
  ele%mat6(5,1:6) = ele%mat6(5,1:6) - (ds2*end_orb%vec(2)/rel_E2) * ele%mat6(2,1:6) + (ds2*end_orb%vec(2)**2/rel_E3) * ele%mat6(6,1:6)
endif

end subroutine apply_p_x 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_p_y (do_mat6)

logical do_mat6

end_orb%vec(3) = end_orb%vec(3) + ds2 * end_orb%vec(4) / rel_E
end_orb%vec(5) = end_orb%vec(5) - ds2 * end_orb%vec(4)**2 / (2*rel_E2)

if (do_mat6) then
  ele%mat6(3,1:6) = ele%mat6(3,1:6) + (ds2 / rel_E)               * ele%mat6(4,1:6) - (ds2*end_orb%vec(4)/rel_E2)    * ele%mat6(6,1:6) 
  ele%mat6(5,1:6) = ele%mat6(5,1:6) - (ds2*end_orb%vec(4)/rel_E2) * ele%mat6(4,1:6) + (ds2*end_orb%vec(4)**2/rel_E3) * ele%mat6(6,1:6)
endif      

end subroutine apply_p_y

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine rf_drift1 (do_mat6)

logical do_mat6

! Drift_1 = (P_x - Ax)**2 / (2 * (1 + dE))

end_orb%vec(2) = end_orb%vec(2) 
end_orb%vec(4) = end_orb%vec(4) 

if (do_mat6) then
  ele%mat6(2,1:6) = ele%mat6(2,1:6) 
  ele%mat6(4,1:6) = ele%mat6(4,1:6) 
endif      

!

call apply_p_x (do_mat6)

!

end_orb%vec(2) = end_orb%vec(2) 
end_orb%vec(4) = end_orb%vec(4) 

if (do_mat6) then
  ele%mat6(2,1:6) = ele%mat6(2,1:6)
  ele%mat6(4,1:6) = ele%mat6(4,1:6)
endif  

end subroutine rf_drift1

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine bsq_drift1 (do_mat6)

logical do_mat6

! Drift_1 = (P_x - Ax)**2 / (2 * (1 + dE))

end_orb%vec(2) = end_orb%vec(2) + end_orb%vec(3) * ks_tot_2   !  vec(2) - Ax
end_orb%vec(4) = end_orb%vec(4) + end_orb%vec(1) * ks_tot_2   !  vec(4) - dint_Ax_dy

if (do_mat6) then
  ele%mat6(2,1:6) = ele%mat6(2,1:6) + ks_tot_2 * ele%mat6(3,1:6)
  ele%mat6(4,1:6) = ele%mat6(4,1:6) + ks_tot_2 * ele%mat6(1,1:6)
endif      

!

call apply_p_x (do_mat6)

!

end_orb%vec(2) = end_orb%vec(2) - end_orb%vec(3) * ks_tot_2   !  vec(2) + Ax
end_orb%vec(4) = end_orb%vec(4) - end_orb%vec(1) * ks_tot_2   !  vec(4) + dint_Ax_dy

if (do_mat6) then
  ele%mat6(2,1:6) = ele%mat6(2,1:6) - ks_tot_2 * ele%mat6(3,1:6)
  ele%mat6(4,1:6) = ele%mat6(4,1:6) - ks_tot_2 * ele%mat6(1,1:6)
endif  

end subroutine bsq_drift1

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine bsq_drift2 (do_mat6)

logical do_mat6

! Drift_2 = (P_y - Ay)**2 / (2 * (1 + dE))

end_orb%vec(2) = end_orb%vec(2) - end_orb%vec(3) * ks_tot_2   !  vec(2) - dint_Ay_dx
end_orb%vec(4) = end_orb%vec(4) - end_orb%vec(1) * ks_tot_2   !  vec(4) - Ay

if (do_mat6) then
  ele%mat6(2,1:6) = ele%mat6(2,1:6) - ks_tot_2 * ele%mat6(3,1:6)
  ele%mat6(4,1:6) = ele%mat6(4,1:6) - ks_tot_2 * ele%mat6(1,1:6)
endif      

!

call apply_p_y (do_mat6)

!

end_orb%vec(2) = end_orb%vec(2) + end_orb%vec(3) * ks_tot_2   !  vec(2) + dint_Ay_dx
end_orb%vec(4) = end_orb%vec(4) + end_orb%vec(1) * ks_tot_2   !  vec(4) + Ay

if (do_mat6) then
  ele%mat6(2,1:6) = ele%mat6(2,1:6) + ks_tot_2 * ele%mat6(3,1:6)
  ele%mat6(4,1:6) = ele%mat6(4,1:6) + ks_tot_2 * ele%mat6(1,1:6)
endif  

end subroutine bsq_drift2

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine bsq_kick (do_mat6)

logical do_mat6

dpx = k1_norm * (x_q - end_orb%vec(1)) - k1_skew * end_orb%vec(3) - g_x
end_orb%vec(2) = end_orb%vec(2) + ds * dpx  ! dAz_dx
              
dpy = k1_norm * (end_orb%vec(3) - y_q) - k1_skew * end_orb%vec(1) - g_y
end_orb%vec(4) = end_orb%vec(4) + ds * dpy  ! dAz_dy

if (do_mat6) then
  ele%mat6(2,1:6) = ele%mat6(2,1:6) - ds * k1_norm * ele%mat6(1,1:6) - ds * k1_skew * ele%mat6(3,1:6)
  ele%mat6(4,1:6) = ele%mat6(4,1:6) - ds * k1_skew * ele%mat6(1,1:6) + ds * k1_norm * ele%mat6(3,1:6)
endif 

end subroutine bsq_kick

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_wig_exp_int_ax (sgn, do_mat6)

integer sgn
logical do_mat6
real(rp) dAx__dy

Ax_saved = Ax()
end_orb%vec(2) = end_orb%vec(2) + sgn * Ax_saved
end_orb%vec(4) = end_orb%vec(4) + sgn * dint_Ax_dy()

if (do_mat6) then
  dAx__dy = dAx_dy()
  ele%mat6(2,1:6) = ele%mat6(2,1:6) + sgn * (dAx_dx() * ele%mat6(1,1:6) + dAx__dy          * ele%mat6(3,1:6))
  ele%mat6(4,1:6) = ele%mat6(4,1:6) + sgn * (dAx__dy  * ele%mat6(1,1:6) + ddint_Ax_dy_dy() * ele%mat6(3,1:6))
endif      

end subroutine apply_wig_exp_int_ax

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_wig_exp_int_ay (sgn, do_mat6)

integer sgn
logical do_mat6
real(rp) dAy__dx

Ay_saved = Ay()
end_orb%vec(2) = end_orb%vec(2) + sgn * dint_Ay_dx()
end_orb%vec(4) = end_orb%vec(4) + sgn * Ay_saved

if (do_mat6) then
  dAy__dx = dAy_dx()
  ele%mat6(2,1:6) = ele%mat6(2,1:6) + sgn * (ddint_Ay_dx_dx() * ele%mat6(1,1:6) + dAy__dx  * ele%mat6(3,1:6))
  ele%mat6(4,1:6) = ele%mat6(4,1:6) + sgn * (dAy__dx          * ele%mat6(1,1:6) + dAy_dy() * ele%mat6(3,1:6))
endif      

end subroutine apply_wig_exp_int_ay

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine calc_wig_coefs (do_mat6)

type (wiggler_computations_struct), pointer :: tmj
real(rp) factor, coef
integer j
logical do_mat6

factor = charge_of(start_orb_saved%species) * c_light * field_ele%value(polarity$) / ele%value(p0c$)

do j = 1, num_wig_terms
  wt => wig_term(j)
  tmj => tm(j)
  tmj = wiggler_computations_struct()
  coef = factor * wt%coef 

  select case (wt%type)
  case (hyper_y_family_x$)
    tmj%trig_x = -1
    tmj%trig_y =  1
    tmj%family = x_family$
    tmj%coef_Ax =  coef * wt%kz / wt%ky                ! Missing factor of: 1/k_y
    tmj%coef_Az =  coef * wt%kx / wt%ky * orient_dir   ! Missing factor of: 1/k_y

  case (hyper_y_family_y$)
    tmj%trig_x = -1
    tmj%trig_y =  1
    tmj%family = y_family$
    tmj%coef_Ay = -coef * wt%kz / wt%ky                ! Missing factor of: 1/k_x
    tmj%coef_Az = -coef * orient_dir                   ! Missing factor of: 1/k_x

  case (hyper_y_family_qu$)
    tmj%trig_x = -1
    tmj%trig_y =  1
    tmj%family = qu_family$
    tmj%coef_Ax =  coef                                ! Missing factor of: 1/k_z
    tmj%coef_Ay = -coef * wt%kx / wt%ky                ! Missing factor of: 1/k_z

  case (hyper_y_family_sq$)
    tmj%trig_x = -1
    tmj%trig_y =  1
    tmj%family = sq_family$
    tmj%coef_Ax =  coef                                ! Missing factor of: 1/k_z
    tmj%coef_Ay =  coef * wt%kx / wt%ky                ! Missing factor of: 1/k_z


  case (hyper_xy_family_x$)
    tmj%trig_x =  1
    tmj%trig_y =  1
    tmj%family = x_family$
    tmj%coef_Ax =  coef                                ! Missing factor of: 1/k_y
    tmj%coef_Az =  coef * wt%kx / wt%kz * orient_dir   ! Missing factor of: 1/k_y

  case (hyper_xy_family_y$)
    tmj%trig_x =  1
    tmj%trig_y =  1
    tmj%family = y_family$
    tmj%coef_Ay = -coef                                ! Missing factor of: 1/k_x
    tmj%coef_Az = -coef * wt%ky / wt%kz * orient_dir   ! Missing factor of: 1/k_x

  case (hyper_xy_family_qu$)
    tmj%trig_x =  1
    tmj%trig_y =  1
    tmj%family = qu_family$
    tmj%coef_Ax =  coef * wt%ky / wt%kz                ! Missing factor of: 1/k_z
    tmj%coef_Ay = -coef * wt%kx / wt%kz                ! Missing factor of: 1/k_z

  case (hyper_xy_family_sq$)
    tmj%trig_x =  1
    tmj%trig_y =  1
    tmj%family = sq_family$
    tmj%coef_Ax =  coef * wt%ky / wt%kz                ! Missing factor of: 1/k_z
    tmj%coef_Ay = -coef * wt%kx / wt%kz                ! Missing factor of: 1/k_z


  case (hyper_x_family_x$)
    tmj%trig_x =  1
    tmj%trig_y = -1
    tmj%family = x_family$
    tmj%coef_Ax =  coef * wt%kz / wt%kx                ! Missing factor of: 1/k_y
    tmj%coef_Az =  coef * orient_dir                   ! Missing factor of: 1/k_y

  case (hyper_x_family_y$)
    tmj%trig_x =  1
    tmj%trig_y = -1
    tmj%family = y_family$
    tmj%coef_Ay = -coef * wt%kz / wt%kx                ! Missing factor of: 1/k_x
    tmj%coef_Az = -coef * wt%ky / wt%kx * orient_dir   ! Missing factor of: 1/k_x

  case (hyper_x_family_qu$)
    tmj%trig_x =  1
    tmj%trig_y = -1
    tmj%family = qu_family$
    tmj%coef_Ax =  coef * wt%ky / wt%kx                ! Missing factor of: 1/k_z
    tmj%coef_Ay = -coef                                ! Missing factor of: 1/k_z

  case (hyper_x_family_sq$)
    tmj%trig_x =  1
    tmj%trig_y = -1
    tmj%family = sq_family$
    tmj%coef_Ax =  coef * wt%ky / wt%kx                ! Missing factor of: 1/k_z
    tmj%coef_Ay =  coef                                ! Missing factor of: 1/k_z

  case default
    call err_exit
  end select
enddo

end subroutine calc_wig_coefs

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_x_terms (err)

type (wiggler_computations_struct), pointer :: tmj

real(rp) arg0, darg, arg, x
integer j
logical err

!

x = end_orb%vec(1) - ct_map%r0(1)

do j = 1, num_wig_terms
  wt => wig_term(j)
  tmj => tm(j)

  arg0 = wt%kx * wt%x0
  darg = wt%kx * x
  arg = arg0 + darg

  select case (tmj%trig_x)
  case (-1)
    tmj%c_x = cos(arg)
    tmj%s_x = sin(arg)

    if (abs(darg) < 1d-4) then
      tmj%integral_sx = x * ((darg/2 - darg**3/24) * cos(arg0) + (1 - darg**2/6) * sin(arg0))
    else
      tmj%integral_sx = (cos(arg0) - tmj%c_x) / wt%kx
    endif

    if (tmj%family == sq_family$) then
      if (abs(darg) < 1d-4) then
        tmj%integral_cx = x * ((1 - darg**2/6) * cos(arg0) + (-darg/2 + darg**3/24) * sin(arg0))
      else
        tmj%integral_cx = (tmj%s_x - sin(arg0)) / wt%kx
      endif
    endif

  case (1)
    if (abs(arg) > 30 .or. abs(arg0) > 30) then
      call err_set (err, x_plane$)
      return
    endif

    tmj%c_x = cosh(arg)
    tmj%s_x = sinh(arg)

    if (abs(darg) < 1d-4) then
      tmj%integral_sx = x * ((darg/2 + darg**3/24) * cosh(arg0) + (1 + darg**2/6) * sinh(arg0))
    else
      tmj%integral_sx = (tmj%c_x - cosh(arg0)) / wt%kx
    endif

    if (tmj%family == sq_family$) then
      if (abs(darg) < 1d-4) then
        tmj%integral_cx = x * ((1 + darg**2/6) * cosh(arg0) + (darg/2 + darg**3/24) * sinh(arg0))
      else
        tmj%integral_cx = (tmj%s_x - sinh(arg0)) / wt%kx
      endif
    endif

  end select

  if (tmj%family == y_family$) then
    if (abs(arg) < 1d-10) then
      tmj%sx_over_kx = x + wt%x0
    else
      tmj%sx_over_kx = tmj%s_x / wt%kx
    endif
  endif

enddo

end subroutine update_wig_x_terms

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_y_terms (err)

type (wiggler_computations_struct), pointer :: tmj

real(rp) arg0, darg, arg, y
integer j
logical err

!

y = end_orb%vec(3) - ct_map%r0(2)

do j = 1, num_wig_terms
  wt => wig_term(j)
  tmj => tm(j)

  arg0 = wt%ky * wt%y0
  darg = wt%ky * y
  arg = arg0 + darg

  select case (tmj%trig_y)
  case (-1)
    tmj%c_y = cos(arg)
    tmj%s_y = sin(arg)

    if (abs(darg) < 1d-4) then
      tmj%integral_sy = y * ((darg/2 - darg**3/24) * cos(arg0) + (1 - darg**2/6) * sin(arg0))
    else
      tmj%integral_sy = (cos(arg0) - tmj%c_y) / wt%ky
    endif

    if (tmj%family == sq_family$) then
      if (abs(darg) < 1d-4) then
        tmj%integral_cy = y * ((1 - darg**2/6) * cos(arg0) + (-darg/2 + darg**3/24) * sin(arg0))
      else
        tmj%integral_cy = (tmj%s_y - sin(arg0)) / wt%ky
      endif
    endif

  case (1)
    if (abs(arg) > 30 .or. abs(arg0) > 30) then
      call err_set (err, y_plane$)
      return
    endif

    tmj%c_y = cosh(arg)
    tmj%s_y = sinh(arg)

    if (abs(darg) < 1d-4) then
      tmj%integral_sy = y * ((darg/2 + darg**3/24) * cosh(arg0) + (1 + darg**2/6) * sinh(arg0))
    else
      tmj%integral_sy = (tmj%c_y - cosh(arg0)) / wt%ky
    endif

    if (tmj%family == sq_family$) then
      if (abs(darg) < 1d-4) then
        tmj%integral_cy = y * ((1 + darg**2/6) * cosh(arg0) + (darg/2 + darg**3/24) * sinh(arg0))
      else
        tmj%integral_cy = (tmj%s_y - sinh(arg0)) / wt%ky
      endif
    endif

  end select

  if (tmj%family == x_family$) then
    if (abs(arg) < 1d-10) then
      tmj%sy_over_ky = y + wt%y0
    else
      tmj%sy_over_ky = tmj%s_y / wt%ky
    endif
  endif


enddo

end subroutine update_wig_y_terms

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_s_terms

real(rp) spz_offset, s0
real(rp) kzz(1:num_wig_terms)

integer j

!

select case (ct_map%ele_anchor_pt)
case (anchor_beginning$); s0 = 0
case (anchor_center$);    s0 = ele%value(l$) / 2
case (anchor_end$);       s0 = ele%value(l$)
end select

spz_offset = s + z_offset - s0 - ct_map%r0(3)

if (ele%orientation == -1) spz_offset = ele%value(l$) - spz_offset
kzz(1:num_wig_terms) = wig_term(1:num_wig_terms)%kz * spz_offset + wig_term(1:num_wig_terms)%phi_z

tm(1:num_wig_terms)%c_z = cos(kzz(1:num_wig_terms))
tm(1:num_wig_terms)%s_z = sin(kzz(1:num_wig_terms))

do j = 1, num_wig_terms
  select case (tm(j)%family)
  case (qu_family$, sq_family$)
    wt => wig_term(j)
    if (abs(kzz(j)) < 1d-5) then
      if (wt%kz == 0) then
        tm(j)%sz_over_kz = spz_offset 
      else
        tm(j)%sz_over_kz = spz_offset + wt%phi_z / wt%kz   ! phi_z must be small
      endif
    else
      tm(j)%sz_over_kz = tm(j)%s_z / wt%kz
    endif
  end select
end do
end subroutine update_wig_s_terms

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function Ax() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Ax * tm(j)%s_x * tm(j)%sy_over_ky * tm(j)%s_z
  case (qu_family$)
    value = value + tm(j)%coef_Ax * tm(j)%s_x * tm(j)%c_y * tm(j)%sz_over_kz
  case (sq_family$)
    value = value + tm(j)%coef_Ax * tm(j)%c_x * tm(j)%s_y * tm(j)%sz_over_kz
  end select
enddo

end function Ax

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAx_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  wt => wig_term(j)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Ax * tm(j)%c_x * tm(j)%sy_over_ky * tm(j)%s_z * wt%kx
  case (qu_family$)
    value = value + tm(j)%coef_Ax * tm(j)%c_x * tm(j)%c_y * tm(j)%sz_over_kz * wt%kx
  case (sq_family$)
    value = value + tm(j)%coef_Ax * tm(j)%s_x * tm(j)%s_y * tm(j)%sz_over_kz * wt%kx * tm(j)%trig_x
  end select
enddo

end function dAx_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAx_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  wt => wig_term(j)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Ax * tm(j)%s_x * tm(j)%c_y * tm(j)%s_z
  case (qu_family$)
    value = value + tm(j)%coef_Ax * tm(j)%s_x * tm(j)%s_y * tm(j)%sz_over_kz * wt%ky * tm(j)%trig_y
  case (sq_family$)
    value = value + tm(j)%coef_Ax * tm(j)%c_x * tm(j)%c_y * tm(j)%sz_over_kz * wt%ky
  end select
enddo

end function dAx_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_Ax_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  wt => wig_term(j)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Ax * tm(j)%integral_sx * tm(j)%c_y * tm(j)%s_z
  case (qu_family$)
    value = value + tm(j)%coef_Ax * tm(j)%integral_sx * tm(j)%s_y * tm(j)%sz_over_kz * wt%ky * tm(j)%trig_y
  case (sq_family$)
    value = value + tm(j)%coef_Ax * tm(j)%integral_cx * tm(j)%c_y * tm(j)%sz_over_kz * wt%ky
  end select
enddo

end function dint_Ax_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function ddint_Ax_dy_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  wt => wig_term(j)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Ax * tm(j)%integral_sx * tm(j)%s_y * tm(j)%s_z * wt%ky * tm(j)%trig_y
  case (qu_family$)
    value = value + tm(j)%coef_Ax * tm(j)%integral_sx * tm(j)%c_y * tm(j)%sz_over_kz * wt%ky**2 * tm(j)%trig_y
  case (sq_family$)
    value = value + tm(j)%coef_Ax * tm(j)%integral_cx * tm(j)%s_y * tm(j)%sz_over_kz * wt%ky**2 * tm(j)%trig_y
  end select
enddo

end function ddint_Ax_dy_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function Ay() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  select case (tm(j)%family)
  case (y_family$)
    value = value + tm(j)%coef_Ay * tm(j)%sx_over_kx * tm(j)%s_y * tm(j)%s_z
  case (qu_family$)
    value = value + tm(j)%coef_Ay * tm(j)%c_x * tm(j)%s_y * tm(j)%sz_over_kz
  case (sq_family$)
    value = value + tm(j)%coef_Ay * tm(j)%s_x * tm(j)%c_y * tm(j)%sz_over_kz
  end select
enddo

end function Ay

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAy_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  wt => wig_term(j)
  select case (tm(j)%family)
  case (y_family$)
    value = value + tm(j)%coef_Ay * tm(j)%c_x * tm(j)%s_y * tm(j)%s_z
  case (qu_family$)
    value = value + tm(j)%coef_Ay * tm(j)%s_x * tm(j)%s_y * tm(j)%sz_over_kz * wt%kx * tm(j)%trig_x
  case (sq_family$)
    value = value + tm(j)%coef_Ay * tm(j)%c_x * tm(j)%c_y * tm(j)%sz_over_kz * wt%kx
  end select
enddo

end function dAy_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAy_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  wt => wig_term(j)
  select case (tm(j)%family)
  case (y_family$)
    value = value + tm(j)%coef_Ay * tm(j)%sx_over_kx * tm(j)%c_y * tm(j)%s_z * wt%ky
  case (qu_family$)
    value = value + tm(j)%coef_Ay * tm(j)%c_x * tm(j)%c_y * tm(j)%sz_over_kz * wt%ky
  case (sq_family$)
    value = value + tm(j)%coef_Ay * tm(j)%s_x * tm(j)%s_y * tm(j)%sz_over_kz * wt%ky * tm(j)%trig_y
  end select
enddo

end function dAy_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_Ay_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  wt => wig_term(j)
  select case (tm(j)%family)
  case (y_family$)
    value = value + tm(j)%coef_Ay * tm(j)%c_x * tm(j)%integral_sy * tm(j)%s_z
  case (qu_family$)
    value = value + tm(j)%coef_Ay * tm(j)%s_x * tm(j)%integral_sy * tm(j)%sz_over_kz * wt%kx * tm(j)%trig_x
  case (sq_family$)
    value = value + tm(j)%coef_Ay * tm(j)%c_x * tm(j)%integral_cy * tm(j)%sz_over_kz * wt%kx
  end select
enddo

end function dint_Ay_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function ddint_Ay_dx_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  wt => wig_term(j)
  select case (tm(j)%family)
  case (y_family$)
    value = value + tm(j)%coef_Ay * tm(j)%s_x * tm(j)%integral_sy * tm(j)%s_z * wt%kx * tm(j)%trig_x
  case (qu_family$)
    value = value + tm(j)%coef_Ay * tm(j)%c_x * tm(j)%integral_sy * tm(j)%sz_over_kz * wt%kx**2 * tm(j)%trig_x
  case (sq_family$)
    value = value + tm(j)%coef_Ay * tm(j)%s_x * tm(j)%integral_cy * tm(j)%sz_over_kz * wt%kx**2 * tm(j)%trig_x
  end select
enddo

end function ddint_Ay_dx_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAz_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Az * tm(j)%s_x * tm(j)%sy_over_ky * tm(j)%c_z * wig_term(j)%kx * tm(j)%trig_x
  case (y_family$)
    value = value + tm(j)%coef_Az * tm(j)%c_x * tm(j)%c_y * tm(j)%c_z
  end select
enddo

end function dAz_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAz_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Az * tm(j)%c_x * tm(j)%c_y * tm(j)%c_z
  case (y_family$)
    value = value + tm(j)%coef_Az * tm(j)%sx_over_kx * tm(j)%s_y * tm(j)%c_z * wig_term(j)%ky * tm(j)%trig_y
  end select
enddo

end function dAz_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function ddAz_dx_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Az * tm(j)%c_x * tm(j)%sy_over_ky * tm(j)%c_z * wig_term(j)%kx**2 * tm(j)%trig_x
  case (y_family$)
    value = value + tm(j)%coef_Az * tm(j)%s_x * tm(j)%c_y * tm(j)%c_z * wig_term(j)%kx * tm(j)%trig_x
  end select
enddo

end function ddAz_dx_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function ddAz_dx_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Az * tm(j)%s_x * tm(j)%c_y * tm(j)%c_z * wig_term(j)%kx * tm(j)%trig_x
  case (y_family$)
    value = value + tm(j)%coef_Az * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z * wig_term(j)%ky * tm(j)%trig_y
  end select
enddo

end function ddAz_dx_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function ddAz_dy_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  select case (tm(j)%family)
  case (x_family$)
    value = value + tm(j)%coef_Az * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z * wig_term(j)%ky * tm(j)%trig_y
  case (y_family$)
    value = value + tm(j)%coef_Az * tm(j)%sx_over_kx * tm(j)%c_y * tm(j)%c_z * wig_term(j)%ky**2 * tm(j)%trig_y
  end select
enddo

end function ddAz_dy_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine radiation_kick()

type (ele_struct) :: temp_ele
type (ele_struct), pointer :: ele0

! Test if kick should be applied

if (.not. bmad_com%radiation_damping_on .and. .not. bmad_com%radiation_fluctuations_on) return

! g2 and g3 radiation integrals can be computed from the change in momentum.

g2 = dpx**2 + dpy**2
g3 = g2 * sqrt(g2)

! synch_rad_com%scale is normally 1 but can be set by a program for testing purposes.

call ran_gauss (this_ran)
dE_p = (1 + end_orb%vec(6)) * (fact_d * g2 + fact_f * sqrt(g3) * this_ran) * synch_rad_com%scale 

! And kick the particle.

end_orb%vec(2) = end_orb%vec(2) - (end_orb%vec(2) - Ax_saved) * dE_p
end_orb%vec(4) = end_orb%vec(4) - (end_orb%vec(4) - Ay_saved) * dE_p
end_orb%vec(6) = end_orb%vec(6) - dE_p * (1 + end_orb%vec(6))

! synch_ran_com%i_calc_on is, by default, False but a program can set this to True for testing purposes.

if (synch_rad_com%i_calc_on) then
  synch_rad_com%i2 = synch_rad_com%i2 + g2 * ds
  synch_rad_com%i3 = synch_rad_com%i3 + g3 * ds
  if (associated(ele%branch)) then
    temp_ele%mat6 = ele%mat6
    temp_ele%vec0(1:5) = end_orb%vec(1:5) - matmul (ele%mat6(1:5,1:6), start_orb_saved%vec)
    temp_ele%vec0(6) = 0
    temp_ele%map_ref_orb_in = start_orb_saved
    temp_ele%map_ref_orb_out = end_orb
    ele0 => ele%branch%ele(ele%ix_ele-1)
    call twiss_propagate1 (ele0, temp_ele)
    synch_rad_com%i5a = synch_rad_com%i5a + g3 * ds * (temp_ele%a%gamma * temp_ele%a%eta**2 + &
          2 * temp_ele%a%alpha * temp_ele%a%eta * temp_ele%a%etap + temp_ele%a%beta * temp_ele%a%etap**2)
    synch_rad_com%i5b = synch_rad_com%i5b + g3 * ds * (temp_ele%b%gamma * temp_ele%b%eta**2 + &
          2 * temp_ele%b%alpha * temp_ele%b%eta * temp_ele%b%etap + temp_ele%b%beta * temp_ele%b%etap**2)
  endif
endif

end subroutine radiation_kick

end subroutine

end module
