!+
! Module em_field_mod
!
! Module to define the electric and magnetic fields for an elemet.
!-

module em_field_mod

use bmad_struct
use bmad_interface
use spline_mod

implicit none

contains

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!+
! Function g_bend_from_em_field (B, E, orbit) result (g_bend)
!
! Routine to calculate the bending strength (1/bending_radius) for
! a given particle for a given field.
!
! Input:
!   B(3)  -- real(rp): Magnetic field.
!   E(3)  -- real(rp): Electric field
!   orbit -- coord_struct: particle orbit
!
! Output:
!   g_bend(3) -- real(rp): bending strength vector.
!-

function g_bend_from_em_field (B, E, orbit) result (g_bend)

type (coord_struct) orbit
real(rp) b(3), e(3), g_bend(3)
real(rp) vel(3), rel_pc, force(3)

! vel is normalized velocity

rel_pc = 1 + orbit%vec(6)
vel(1:2) = [orbit%vec(2), orbit%vec(4)] / rel_pc
vel(3) = sqrt(1 - vel(1)**2 - vel(2)**2) * orbit%direction

force = (E + cross_product(vel, B) * orbit%beta * c_light) * charge_of(orbit%species)
g_bend = -(force - vel * (dot_product(force, vel))) / (orbit%p0c * rel_pc)

end function g_bend_from_em_field

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine em_field_calc (ele, param, s_pos, orbit, local_ref_frame, field, calc_dfield, &
!                             err_flag, potential, use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles)
!
! Routine to calculate the E and B fields at a particular place in an element.
!
! Note: Zero field will be returned if an element is turned off.
!
! Note: The fields due to any kicks will be present. It therefore important in tracking to make sure that 
! offset_particle does not add in kicks at the beginning and end which would result in double counting the kicks.
!
! Input:
!   ele             -- Ele_struct: Lattice element.
!   param           -- lat_param_struct: Lattice parameters.
!   s_pos           -- Real(rp): Longitudinal position.
!                        If local_ref_frame = T: In Body coords relative to the entrance edge of the element.
!                        If local_ref_frame = F: In Lab coords relative to the upstream edge of the element.
!   orbit           -- Coord_struct: Transverse coordinates.
!     %vec(1), %vec(3) -- Transverse coords.
!     %t               -- Used with absolute time tracking.
!     %vec(5)          -- Used with relative time tracking (except with time Runge-Kutta).
!   local_ref_frame  -- Logical, If True then take the input coordinates and output fields 
!                         as being with respect to the frame of referene of the element (ignore misalignments). 
!   calc_dfield      -- Logical, optional: If present and True then calculate the field derivatives.
!   use_overlap      -- logical, optional: Add in overlap fields from other elements? Default is True.
!   grid_allow_s_out_of_bounds 
!                    -- logical, optional: For grids, allow s-coordinate to be grossly out of bounds 
!                          and return zero instead of an error? Default: False. Used internally for overlapping fields.
!   rf_time          -- real(rp), optional: Set the time relative to the RF clock. Normally this time is calculated using
!                          orbit%t or orbit%vec(5) but sometimes it is convenient to be able to override this.
!                          For example, time_runge_kutta uses this.
!   used_eles(:)     -- ele_pointer_struct, allocatable, optional: For internal use only when this routine is
!                          called recursively. Used to prevent double counting when there is field overlap.
!
! Output:
!   field       -- em_field_struct: E and B fields and derivatives.
!   err_flag    -- logical, optional: Set True if there is an error. False otherwise.
!   potential   -- em_potential_struct, optional: The electric and magnetic potentials.
!                   This is experimental and only implemented for wigglers at present.
!-

recursive subroutine em_field_calc (ele, param, s_pos, orbit, local_ref_frame, field, calc_dfield, &
                           err_flag, potential, use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles)

use geometry_mod

type (ele_struct), target :: ele, ele2
type (ele_pointer_struct), allocatable, optional :: used_eles(:)
type (ele_pointer_struct), allocatable :: use_list(:)
type (ele_struct), pointer :: lord
type (lat_param_struct) param
type (coord_struct) :: orbit, local_orb, lab_orb, lord_orb, this_orb
type (em_potential_struct), optional :: potential
type (em_field_struct) :: field, field1, field2, lord_field, l1_field, mode_field
type (cartesian_map_struct), pointer :: ct_map
type (cartesian_map_term1_struct), pointer :: ct_term
type (cylindrical_map_struct), pointer :: cl_map
type (cylindrical_map_term1_struct), pointer :: cl_term
type (grid_field_struct), pointer :: g_field
type (grid_field_pt1_struct) g_pt
type (taylor_field_struct), pointer :: t_field
type (taylor_field_plane1_struct), pointer :: t_plane
type (floor_position_struct) lab_position, global_position, lord_position
type (spline_struct) spline

real(rp), optional :: rf_time
real(rp) :: x, y, time, s_pos, s_body, s_lab, s_lab2, z, ff, dk(3,3), ref_charge, f_p0c
real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, coef, fd(3), s0, Ex, Ey, amp
real(rp) :: cos_ang, sin_ang, sgn_x, sgn_y, sgn_z, kx, ky, dkm(2,2), cos_ks, sin_ks
real(rp) phase, gradient, r, E_r_coef, E_s, k_wave, s_eff, a_amp
real(rp) k_t, k_zn, kappa2_n, kap_rho, s_hard_offset, beta_start
real(rp) radius, phi, t_ref, tilt, omega, freq0, freq, B_phi_coef, z_center
real(rp) sx_over_kx, sy_over_ky, sz_over_kz, rot2(2,2)
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), pot
real(rp) w_ele_mat(3,3), w_lord_mat(3,3), Er, Ep, Ez, Br, Bp, Bz
real(rp) :: fld(3), dfld(3,3), fld0(3), fld1(3), dfld0(3,3), dfld1(3,3)
real(rp) phi0_autoscale, field_autoscale, ds, beta_ref

complex(rp) exp_kz, expt, dEp, dEr, E_rho, E_phi, E_z, B_rho, B_phi, B_z
complex(rp) Im_0, Im_plus, Im_minus, Im_0_R, kappa_n, Im_plus2, cm, sm, q

integer i, j, m, n, ix, trig_x, trig_y, status, im, iz0, iz1, izp, field_calc, ix_pole_max

logical :: local_ref_frame
logical, optional :: calc_dfield, err_flag, use_overlap, grid_allow_s_out_of_bounds
logical do_df_calc, err, dfield_computed, add_kicks

character(20) :: r_name = 'em_field_calc'

! Initialize field
! If element is turned off then return zero

field = em_field_struct()
if (present(potential)) potential = em_potential_struct()

do_df_calc = logic_option (.false., calc_dfield)
dfield_computed = .false.

if (present(err_flag)) err_flag = .false.
if (.not. ele%is_on) return

! Has this element been used before?

if (present(used_eles)) then
  do j = 1, size(used_eles)
    if (.not. associated (used_eles(j)%ele)) then
      used_eles(j)%ele => ele
      exit
    endif

    if (associated(used_eles(j)%ele, ele)) return
  enddo

  if (j == size(used_eles) + 1) then
    call move_alloc(used_eles, use_list)
    allocate(used_eles(2*j))
    used_eles(1:j-1) = use_list
    used_eles(j)%ele => ele
  endif
endif

!----------------------------------------------------------------------------
! super_slave, and slice_slave, have their field info stored in the associated lord elements.
! Note: multipass_slave elements do store their own field info. This should be changed.

if (ele%field_calc == refer_to_lords$) then
  if (.not. present(used_eles)) allocate (use_list(ele%n_lord+5))

  ! The lord of an element may have independent misalignments.
  ! So use an orbit that is not in the slave's reference frame.

  lab_orb = orbit

  if (local_ref_frame) then
    call offset_particle (ele, param, unset$, lab_orb, set_hvkicks = .false., s_pos = s_pos, s_out = s_lab)
  else
    s_lab = s_pos
  endif

  z = lab_orb%vec(5)

  !

  lord_loop: do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)

    if (lord%field_calc == no_field$) cycle   ! Group, overlay and girder elements do not have fields.

    ds = ele%s_start - lord%s_start
    s_lab2 = s_lab + ds
    beta_ref = lord%value(p0c$) / lord%value(e_tot$)
    lab_orb%vec(5) = z - c_light * orbit%beta * &
        ((ele%value(ref_time_start$) - lord%value(ref_time_start$)) - ds / (beta_ref * c_light))

    if (present(used_eles)) then
      do j = 1, size(used_eles)
        if (.not. associated(used_eles(j)%ele)) exit
        if (associated(used_eles(j)%ele, lord)) cycle lord_loop
      enddo
      call em_field_calc (lord, param, s_lab2, lab_orb, .false., field2, calc_dfield, err, potential, &
                                               use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles)
    else
      call em_field_calc (lord, param, s_lab2, lab_orb, .false., field2, calc_dfield, err, potential, &
                                               use_overlap, grid_allow_s_out_of_bounds, rf_time, use_list)
    endif

    if (err) then
      if (present(err_flag)) err_flag = .true.
      return
    endif

    field%E = field%E + field2%E
    field%B = field%B + field2%B
    if (do_df_calc) then
      field%dE = field%dE + field2%dE
      field%dB = field%dB + field2%dB
    endif

  enddo lord_loop

  if (local_ref_frame) call convert_field_ele_to_lab(ele, s_lab, .false., field)
  return
endif

!----------------------------------------------------------------------------
! Custom field calc 

if (ele%field_calc == custom$) then
  call em_field_custom (ele, param, s_pos, orbit, local_ref_frame, field, calc_dfield, err_flag)
  return
end if

!-------------------------------
! Sad mult is complicated by the fact that the sad_mult defines additional multipole misalignment

if (ele%key == sad_mult$ .and. ele%value(sad_flag$) == 0) then
  call transfer_ele(ele, ele2)
  ele2%value(sad_flag$) = 1  ! To prevent infinite recursion

  ! Solenoid calc. Ignore multipoles
  nullify(ele2%a_pole)
  nullify(ele2%b_pole)
  call em_field_calc (ele2, param, s_pos, orbit, local_ref_frame, field1, calc_dfield, err_flag, potential, &
                                                     use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles)
  ! multipole calc
  ele2%value(ks$) = 0
  ele2%value(bs_field$) = 0
  ele2%value(x_pitch_tot$) = ele%value(x_pitch_tot$) + ele%value(x_pitch_mult$)
  ele2%value(y_pitch_tot$) = ele%value(y_pitch_tot$) + ele%value(y_pitch_mult$)
  ele2%value(x_offset_tot$) = ele%value(x_offset_tot$) + ele%value(x_offset_mult$)
  ele2%value(y_offset_tot$) = ele%value(y_offset_tot$) + ele%value(y_offset_mult$)
  call em_field_calc (ele2, param, s_pos, orbit, local_ref_frame, field2, calc_dfield, err_flag, potential, &
                                                     use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles)

  field%b  = field1%b + field2%b
  field%db = field1%db + field2%db 
  return
endif

!----------------------------------------------------------------------------
! convert to local coords

local_orb = orbit
if (local_ref_frame) then
  s_body = s_pos
else
  call offset_particle (ele, param, set$, local_orb, set_hvkicks = .false., s_pos = s_pos, s_out = s_body)
endif

!----------------------------------------------------------------------------
! Set up common variables for all (non-custom) methods

ref_charge = charge_of(param%particle)

x = local_orb%vec(1)
y = local_orb%vec(3)

if (ref_charge == 0) then
  f_p0c = 0
else
  f_p0c = ele%value(p0c$) / (c_light * ref_charge)
endif

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! field_calc methods

field_calc = ele%field_calc
if (ele%key == wiggler$ .and. ele%sub_key == periodic_type$) field_calc = fieldmap$

select case (field_calc)
  
!----------------------------------------------------------------------------
! Bmad_standard field calc 

case (bmad_standard$)

  ! Field outside of element is zero
  if (s_body < 0 .or. s_body > ele%value(l$)) goto 8000   ! Goto field overlap code.

  select case (ele%key)

  !------------------------------------------
  ! Drift, et. al. Note that kicks get added at the end for all elements

  case (drift$, ecollimator$, rcollimator$, instrument$, monitor$, pipe$, marker$, detector$)

  !------------------------------------------
  ! E_Gun

  case (e_gun$)
    if (ele%value(rf_frequency$) == 0) then
      field%e(3) = e_accel_field (ele, gradient$) / ref_charge
    else
      if (present(rf_time)) then
        time = rf_time
      else
        time = particle_rf_time(orbit, ele, .true., s_body)
      endif
      phase = (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + ele%value(phi0_autoscale$))
      field%e(3) = e_accel_field (ele, gradient$) * cos(twopi * (time * ele%value(rf_frequency$) + phase)) / ref_charge
    endif

  !------------------------------------------
  ! Elseparator

  case (elseparator$)
    field%e(1) = ele%value(hkick$) * ele%value(p0c$) / ele%value(l$)
    field%e(2) = ele%value(vkick$) * ele%value(p0c$) / ele%value(l$)

  !------------------------------------------
  ! HKicker

  case (hkicker$)
    field%b(2) = -ele%value(kick$) * f_p0c / ele%value(l$)

  !------------------------------------------
  ! Kicker  

  case (kicker$, ac_kicker$)
    field%b(1) =  ele%value(vkick$) * f_p0c / ele%value(l$)
    field%b(2) = -ele%value(hkick$) * f_p0c / ele%value(l$)

  !------------------------------------------
  ! RFcavity and Lcavity  bmad_standard
  !
  ! For standing wave cavity:
  ! Use N_cell half-wave pillbox formulas for TM_011 mode with infinite wall radius.
  ! See S.Y. Lee, "Accelerator Physics"
  !   E_s   = 2 * gradient * cos(k s) * cos(omega t + phase)
  !   E_r   =     gradient * k * r * sin(k s) * cos(omega t + phase)
  !   B_phi =    -gradient * k * r * cos(k s) * sin(omega t + phase) / c_light
  ! For traveling wave cavity:
  !   E_s   =  gradient * cos(omega t + phase - k s)
  !   E_r   = -gradient * k * r * sin(omega t + phase - k s) / 2
  !   B_phi = -gradient * k * r * sin(omega t + phase - k s) / c_light / 2
  ! 
  ! Note: Length of pillbox is 1/2 wavelength. Not the length of the element.  
  ! That is, the field extends to +/- c_light * freq / 2 from centerline of element.
  !
  ! Since the active (pillbox) length is different from the element length, the gradient used
  ! is different from the element-gradient = voltage / element-length stored in the element struct so that
  !   gradient-used * pillbox-length = element-gradient * element-length = voltage
  !
  ! Note: There is a discontinuity in the field at the edge. Edge focusing due to this 
  !  discontinuity can be handled in the apply_element_edge_kick routine.

  case(rfcavity$, lcavity$)

    if (ele%value(rf_frequency$) == 0) return

    phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + ele%value(phi0_autoscale$))
    if (ele%key == rfcavity$) phase = pi/2 - phase
    orbit%phase(1) = phase  ! RF phase is needed by apply_element_edge_kick when calling rf_coupler_kick.

    gradient = e_accel_field (ele, gradient$)

    if (.not. ele%is_on) gradient = 0
    gradient = (gradient + gradient_shift_sr_wake(ele, param)) / ref_charge
    gradient = gradient * ele%value(l$) / hard_edge_model_length(ele)
    omega = twopi * ele%value(rf_frequency$)
    k_wave = omega / c_light

    s_hard_offset = (ele%value(l$) - hard_edge_model_length(ele)) / 2  ! Relative to entrance end of the cavity
    s_eff = s_body - s_hard_offset
    if (s_eff < 0 .or. s_eff > hard_edge_model_length(ele)) goto 8000  ! Zero field outside

    beta_start = ele%value(p0c_start$) / ele%value(e_tot_start$)

    if (present(rf_time)) then
      time = rf_time
    else
      time = particle_rf_time(orbit, ele, .true., s_body)
    endif

    if (nint(ele%value(cavity_type$)) == traveling_wave$) then
      phi = omega * time + phase - k_wave * s_eff
      E_z        =  gradient * cos(phi)
      E_r_coef   = -gradient * k_wave * sin(phi) / 2
      B_phi_coef = -gradient * k_wave * sin(phi) / c_light / 2
    else
      E_z        = 2 * gradient * cos(k_wave * s_eff) * cos(omega * time + phase)
      E_r_coef   =     gradient * k_wave * sin(k_wave*s_eff) * cos(omega * time + phase)
      B_phi_coef =    -gradient * k_wave * cos(k_wave*s_eff) * sin(omega * time + phase) / c_light 
    endif

    field%E(1) = E_r_coef * x
    field%E(2) = E_r_coef * y
    field%E(3) = E_z
    
    field%B(1) = -B_phi_coef * y
    field%B(2) =  B_phi_coef * x

  !------------------------------------------
  ! Octupole 

  case (octupole$)

    field%b(1) = -(y**3 - 3*y*x**2) / 6 * ele%value(k3$) * f_p0c 
    field%b(2) =  (x**3 - 3*x*y**2) / 6 * ele%value(k3$) * f_p0c 

    if (do_df_calc) then
      dfield_computed = .true.
      field%dB(1,1) =  x*y * ele%value(k3$) * f_p0c
      field%dB(1,2) = (x**2 - y**2) / 2 * ele%value(k3$) * f_p0c
      field%dB(2,1) = (x**2 - y**2) / 2 * ele%value(k3$) * f_p0c
      field%dB(2,2) = -x*y * ele%value(k3$) * f_p0c
    endif

  !------------------------------------------
  ! Patch: There are no fields

  case (patch$)

  !------------------------------------------
  ! Quadrupole

  case (quadrupole$) 

    field%b(1) = y * ele%value(k1$) * f_p0c 
    field%b(2) = x * ele%value(k1$) * f_p0c 

    if (do_df_calc) then
      dfield_computed = .true.
      field%dB(1,2) =  ele%value(k1$) * f_p0c
      field%dB(2,1) =  ele%value(k1$) * f_p0c
    endif

  !------------------------------------------
  ! Sextupole 

  case (sextupole$)

    field%b(1) = x * y * ele%value(k2$) * f_p0c
    field%b(2) = (x**2 - y**2) / 2 * ele%value(k2$) * f_p0c 

    if (do_df_calc) then
      dfield_computed = .true.
      field%dB(1,1) =  y * ele%value(k2$) * f_p0c
      field%dB(1,2) =  x * ele%value(k2$) * f_p0c
      field%dB(2,1) =  x * ele%value(k2$) * f_p0c
      field%dB(2,2) = -y * ele%value(k2$) * f_p0c
    endif

  !------------------------------------------
  ! VKicker

  case (vkicker$)
    field%b(1) =  ele%value(kick$) * f_p0c / ele%value(l$)

  !------------------------------------------
  ! SBend

  case (sbend$)

    ! Finite k1 or k2 is handled with rest of multipoles
    field%b(2) = (ele%value(g$) + ele%value(g_err$)) * f_p0c 


  !------------------------------------------
  ! Sol_quad

  case (sol_quad$)

    field%b(1) = y * ele%value(k1$) * f_p0c 
    field%b(2) = x * ele%value(k1$) * f_p0c 
    field%b(3) = ele%value(ks$) * f_p0c

    if (do_df_calc) then
      dfield_computed = .true.
      field%dB(1,2) = ele%value(k1$) * f_p0c
      field%dB(2,1) = ele%value(k1$) * f_p0c
    endif

  !------------------------------------------
  ! Solenoid

  case (solenoid$, sad_mult$)

    field%b(3) = ele%value(ks$) * f_p0c

    if (do_df_calc) then
      dfield_computed = .true.
    endif

  !------------------------------------------
  ! Wiggler

  case(wiggler$, undulator$)


  !------------------------------------------
  ! Error

  case default
    call out_io (s_fatal$, r_name, 'EM FIELD NOT YET CODED FOR ELEMENT OF TYPE: ' // key_name(ele%key), &
                                   'FOR ELEMENT: ' // ele%name, 'PERHAPS "FIELD_CALC" NEEDS TO BE SET FOR THIS ELEMENT?')
    if (global_com%exit_on_error) call err_exit
    if (present(err_flag)) err_flag = .true.
    return
  end select

  !---------------------------------------------------------------------
  ! Add multipoles

  if (ele%key == sbend$ .and. nint(ele%value(exact_multipoles$)) /= off$) then
    call bend_exact_multipole_field (ele, param, orbit, local_ref_frame, field2, do_df_calc, potential)
    field%e = field%e + field2%e
    field%b = field%b + field2%b
    if (do_df_calc) then
      field%de = field2%de
      field%db = field2%db
    endif
    add_kicks = .false.  ! h/vkicks are accounted for in bend_exact_multipole_field

  ! Everything but exact bend multipoles
  else
    add_kicks = .true.
    call multipole_ele_to_ab(ele, .not. local_ref_frame, ix_pole_max, a_pole, b_pole)

    if (ele%key == sbend$) then
      if (ele%value(k1$) /= 0) then
        b_pole(1) = b_pole(1) + ele%value(k1$) * ele%value(l$)
        ix_pole_max = max(1, ix_pole_max)
      endif
      if (ele%value(k2$) /= 0) then
        b_pole(2) = b_pole(2) + ele%value(k2$) * ele%value(l$) / 2 
        ix_pole_max = max(2, ix_pole_max)
      endif
    endif

    if (ix_pole_max > -1) then

      if (ele%value(l$) == 0) then
        call out_io (s_fatal$, r_name, 'CANNOT COMPUTE FIELD OF ZERO LENGTH ELEMENT WITH MULTIPOLES. FOR: ' // ele%name)
        if (global_com%exit_on_error) call err_exit
        if (present(err_flag)) err_flag = .true.
        return
      endif

      do i = 0, ix_pole_max
        if (a_pole(i) == 0 .and. b_pole(i) == 0) cycle
        if (do_df_calc) then
          call ab_multipole_kick(a_pole(i), b_pole(i), i, local_orb%species, 0, local_orb, kx, ky, dkm)
        else
          call ab_multipole_kick(a_pole(i), b_pole(i), i, local_orb%species, 0, local_orb, kx, ky)
        endif
        field%B(1) = field%B(1) + f_p0c * ky / ele%value(l$)
        field%B(2) = field%B(2) - f_p0c * kx / ele%value(l$)
        if (do_df_calc) then
          field%dB(1,1) = field%dB(1,1) + f_p0c * dkm(2,1) / ele%value(l$)
          field%dB(1,2) = field%dB(1,2) + f_p0c * dkm(2,2) / ele%value(l$)
          field%dB(2,1) = field%dB(2,1) - f_p0c * dkm(1,1) / ele%value(l$)
          field%dB(2,2) = field%dB(2,2) - f_p0c * dkm(1,2) / ele%value(l$)
        endif
      enddo
    endif

    ! Add electric multipoles

    call multipole_ele_to_ab(ele, .not. local_ref_frame, ix_pole_max, a_pole, b_pole, electric$)
    do i = 0, ix_pole_max
      if (a_pole(i) == 0 .and. b_pole(i) == 0) cycle
      call elec_multipole_field(a_pole(i), b_pole(i), i, local_orb, Ex, Ey, dkm, do_df_calc)
      field%E(1) = field%E(1) + Ex
      field%E(2) = field%E(2) + Ey
      if (do_df_calc) field%dE(1:2,1:2) = field%dE(1:2,1:2) + dkm
    enddo

  endif

  !-------------------------------
  ! Add kicks. Since the kick field is not rotated by a tilt then we have to unrotate if in the local_ref_frame

  if (add_kicks .and. has_hkick_attributes(ele%key) .and. (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0)) then
    select case (ele%key)
    ! Kickers and elsep handled above
    case (ac_kicker$, kicker$, hkicker$, vkicker$, elseparator$)  
    ! Everything else...
    case default
      if (.not. local_ref_frame .or. ele%value(tilt_tot$) == 0) then
        field%b(1) = field%b(1) + ele%value(Vkick$) * f_p0c / ele%value(l$)
        field%b(2) = field%b(2) - ele%value(Hkick$) * f_p0c / ele%value(l$)
      else
        ! Rotate from lab to local
        tilt = ele%value(tilt_tot$)
        if (ele%key == sbend$) tilt = ele%value(ref_tilt_tot$)
        field%b(1) = field%b(1) + (ele%value(Vkick$) * cos(tilt) - ele%value(hkick$) * sin(tilt)) * f_p0c / ele%value(l$)
        field%b(2) = field%b(2) - (ele%value(Hkick$) * cos(tilt) + ele%value(vkick$) * sin(tilt)) * f_p0c / ele%value(l$)
      endif
    end select
  endif

!----------------------------------------------------------------------------
! FieldMap

case(fieldmap$)

  if (present(rf_time)) then
    time = rf_time
  else
    time = particle_rf_time(orbit, ele, .false., s_body)
  endif

  if (.not. associated(ele%cylindrical_map) .and. .not. associated(ele%cartesian_map) .and. &
      .not. associated(ele%grid_field) .and. .not. associated(ele%taylor_field)) then
    call out_io (s_fatal$, r_name, 'No associated fieldmap (cartesican_map, grid_field, etc) FOR: ' // ele%name) 
    if (global_com%exit_on_error) call err_exit
    if (present(err_flag)) err_flag = .true.
    return  
  endif

  select case (ele%key)
  case (e_gun$, em_field$, lcavity$, rfcavity$)
    phi0_autoscale = ele%value(phi0_autoscale$)
    field_autoscale = ele%value(field_autoscale$)
  case default
    phi0_autoscale = 0
    field_autoscale = 1
  end select

  !----------------------------------------------------------------------------
  ! Cartesian map field

  if (associated(ele%cartesian_map)) then
    do im = 1, size(ele%cartesian_map)
      ct_map => ele%cartesian_map(im)

      fld = 0; dfld = 0

      call to_field_map_coords (local_orb, s_body, ct_map%ele_anchor_pt, ct_map%r0, .false., x, y, z)

      n = size(ct_map%ptr%term)
      do i = 1, n
        ct_term => ct_map%ptr%term(i)
        sgn_x = 1; sgn_y = 1; sgn_z = 1

        select case (ct_term%type)
        case (hyper_y_family_x$, hyper_y_family_y$, hyper_y_family_qu$, hyper_y_family_sq$)
          coef = ct_term%coef / ct_term%ky
          c_x = cos(ct_term%kx * (x + ct_term%x0))
          s_x = sin(ct_term%kx * (x + ct_term%x0))
          c_y = cosh (ct_term%ky * (y + ct_term%y0))
          s_y = sinh (ct_term%ky * (y + ct_term%y0))
          if (ct_term%type == hyper_y_family_y$)  sgn_x = -1
          if (ct_term%type == hyper_y_family_sq$) sgn_x = -1
          if (ct_term%type == hyper_y_family_sq$) sgn_z = -1
          trig_x = -1; trig_y = 1

        case (hyper_xy_family_x$, hyper_xy_family_y$, hyper_xy_family_qu$, hyper_xy_family_sq$)
          coef = ct_term%coef / ct_term%kz
          c_x = cosh(ct_term%kx * (x + ct_term%x0))
          s_x = sinh(ct_term%kx * (x + ct_term%x0))
          c_y = cosh (ct_term%ky * (y + ct_term%y0))
          s_y = sinh (ct_term%ky * (y + ct_term%y0))
          if (ct_term%type == hyper_xy_family_sq$) sgn_z = -1
          trig_x = 1; trig_y = 1

        case (hyper_x_family_x$, hyper_x_family_y$, hyper_x_family_qu$, hyper_x_family_sq$)
          coef = ct_term%coef / ct_term%kx
          c_x = cosh(ct_term%kx * (x + ct_term%x0))
          s_x = sinh(ct_term%kx * (x + ct_term%x0))
          c_y = cos (ct_term%ky * (y + ct_term%y0))
          s_y = sin (ct_term%ky * (y + ct_term%y0))
          if (ct_term%type == hyper_x_family_x$)  sgn_y = -1
          if (ct_term%type == hyper_x_family_sq$) sgn_x = -1
          trig_x = 1; trig_y = -1
        end select

        c_z = cos (ct_term%kz * z + ct_term%phi_z)
        s_z = sin (ct_term%kz * z + ct_term%phi_z)

        select case (ct_term%type)
        case (hyper_y_family_x$, hyper_xy_family_x$, hyper_x_family_x$)
          fld(1) = fld(1) + coef  * ct_term%kx * c_x * c_y * c_z
          fld(2) = fld(2) + coef  * ct_term%ky * s_x * s_y * c_z * sgn_y
          fld(3) = fld(3) - coef  * ct_term%kz * s_x * c_y * s_z
        case (hyper_y_family_y$, hyper_xy_family_y$, hyper_x_family_y$)
          fld(1) = fld(1) + coef  * ct_term%kx * s_x * s_y * c_z * sgn_x
          fld(2) = fld(2) + coef  * ct_term%ky * c_x * c_y * c_z
          fld(3) = fld(3) - coef  * ct_term%kz * c_x * s_y * s_z
        case (hyper_y_family_qu$, hyper_xy_family_qu$, hyper_x_family_qu$)
          fld(1) = fld(1) + coef  * ct_term%kx * c_x * s_y * c_z
          fld(2) = fld(2) + coef  * ct_term%ky * s_x * c_y * c_z
          fld(3) = fld(3) - coef  * ct_term%kz * s_x * s_y * s_z
        case (hyper_y_family_sq$, hyper_xy_family_sq$, hyper_x_family_sq$)
          fld(1) = fld(1) + coef  * ct_term%kx * s_x * c_y * c_z * sgn_x
          fld(2) = fld(2) + coef  * ct_term%ky * c_x * s_y * c_z
          fld(3) = fld(3) + coef  * ct_term%kz * c_x * c_y * s_z * sgn_z
        end select

        if (do_df_calc) then
          dfield_computed = .true.
          select case (ct_term%type)
          case (hyper_y_family_x$, hyper_xy_family_x$, hyper_x_family_x$)
            ff = coef * ct_term%kx
            dfld(1,1) = dfld(1,1) + ff * ct_term%kx * s_x * c_y * c_z * trig_x
            dfld(2,1) = dfld(2,1) + ff * ct_term%ky * c_x * s_y * c_z * sgn_y
            dfld(3,1) = dfld(3,1) - ff * ct_term%kz * c_x * c_y * s_z 
            ff = coef * ct_term%ky
            dfld(1,2) = dfld(1,2) + ff * ct_term%kx * c_x * s_y * c_z * trig_y
            dfld(2,2) = dfld(2,2) + ff * ct_term%ky * s_x * c_y * c_z * sgn_y
            dfld(3,2) = dfld(3,2) - ff * ct_term%kz * s_x * s_y * s_z * trig_y
            ff = coef * ct_term%kz
            dfld(1,3) = dfld(1,3) - ff * ct_term%kx * c_x * c_y * s_z
            dfld(2,3) = dfld(2,3) - ff * ct_term%ky * s_x * s_y * s_z * sgn_y 
            dfld(3,3) = dfld(3,3) - ff * ct_term%kz * s_x * c_y * c_z
          case (hyper_y_family_y$, hyper_xy_family_y$, hyper_x_family_y$)
            ff = coef * ct_term%kx
            dfld(1,1) = dfld(1,1) + ff * ct_term%kx * c_x * s_y * c_z * sgn_x
            dfld(2,1) = dfld(2,1) + ff * ct_term%ky * s_x * c_y * c_z * trig_x
            dfld(3,1) = dfld(3,1) - ff * ct_term%kz * s_x * s_y * s_z * trig_x
            ff = coef * ct_term%ky
            dfld(1,2) = dfld(1,2) + ff * ct_term%kx * s_x * c_y * c_z * sgn_x
            dfld(2,2) = dfld(2,2) + ff * ct_term%ky * c_x * s_y * c_z * trig_y
            dfld(3,2) = dfld(3,2) - ff * ct_term%kz * c_x * c_y * s_z 
            ff = coef * ct_term%kz
            dfld(1,3) = dfld(1,3) - ff * ct_term%kx * s_x * s_y * s_z * sgn_x
            dfld(2,3) = dfld(2,3) - ff * ct_term%ky * c_x * c_y * s_z 
            dfld(3,3) = dfld(3,3) - ff * ct_term%kz * c_x * s_y * c_z 
          case (hyper_y_family_qu$, hyper_xy_family_qu$, hyper_x_family_qu$)
            ff = coef * ct_term%kx
            dfld(1,1) = dfld(1,1) + ff * ct_term%kx * s_x * s_y * c_z * trig_x
            dfld(2,1) = dfld(2,1) + ff * ct_term%ky * c_x * c_y * c_z 
            dfld(3,1) = dfld(3,1) - ff * ct_term%kz * c_x * s_y * s_z
            ff = coef * ct_term%ky
            dfld(1,2) = dfld(1,2) + ff * ct_term%kx * c_x * c_y * c_z 
            dfld(2,2) = dfld(2,2) + ff * ct_term%ky * s_x * s_y * c_z * trig_y
            dfld(3,2) = dfld(3,2) - ff * ct_term%kz * s_x * c_y * s_z 
            ff = coef * ct_term%kz
            dfld(1,3) = dfld(1,3) - ff * ct_term%kx * c_x * s_y * s_z
            dfld(2,3) = dfld(2,3) - ff * ct_term%ky * s_x * c_y * s_z 
            dfld(3,3) = dfld(3,3) - ff * ct_term%kz * s_x * s_y * c_z 
          case (hyper_y_family_sq$, hyper_xy_family_sq$, hyper_x_family_sq$)
            ff = coef * ct_term%kx
            dfld(1,1) = dfld(1,1) + ff * ct_term%kx * c_x * c_y * c_z * sgn_x
            dfld(2,1) = dfld(2,1) + ff * ct_term%ky * s_x * s_y * c_z * trig_x
            dfld(3,1) = dfld(3,1) + ff * ct_term%kz * s_x * c_y * s_z * sgn_z * trig_x
            ff = coef * ct_term%ky
            dfld(1,2) = dfld(1,2) + ff * ct_term%kx * s_x * s_y * c_z * sgn_x * trig_y
            dfld(2,2) = dfld(2,2) + ff * ct_term%ky * c_x * c_y * c_z 
            dfld(3,2) = dfld(3,2) + ff * ct_term%kz * c_x * s_y * s_z * sgn_z * trig_y
            ff = coef * ct_term%kz
            dfld(1,3) = dfld(1,3) - ff * ct_term%kx * s_x * c_y * s_z * sgn_x
            dfld(2,3) = dfld(2,3) - ff * ct_term%ky * c_x * s_y * s_z 
            dfld(3,3) = dfld(3,3) + ff * ct_term%kz * c_x * c_y * c_z * sgn_z
          end select
        endif

        if (present(potential)) then
          coef = ct_term%coef 
          select case (ct_term%type)
          case (hyper_y_family_x$, hyper_xy_family_x$, hyper_x_family_x$)
            if (abs(ct_term%ky * (y + ct_term%y0)) < 1d-10) then
              sy_over_ky = y + ct_term%y0
            else
              sy_over_ky = s_y / ct_term%ky
            endif
          case (hyper_y_family_y$, hyper_xy_family_y$, hyper_x_family_y$)
            if (abs(ct_term%kx * (x + ct_term%x0)) < 1d-10) then
              sx_over_kx = x + ct_term%x0
            else
              sx_over_kx = s_x / ct_term%kx
            endif
          case default
            if (abs(ct_term%kz * z + ct_term%phi_z) < 1d-10) then
              if (ct_term%kz == 0) then
                sz_over_kz = 0
              else
                sz_over_kz = z + ct_term%phi_z / ct_term%kz
              endif
            else
              sz_over_kz = s_z / ct_term%kz
            endif
          end select

          if (ct_map%field_type == magnetic$) then
            select case (ct_term%type)
            case (hyper_y_family_x$)
              potential%a(1) = potential%a(1) + coef * s_x * sy_over_ky * s_z * ct_term%kz / ct_term%ky
              potential%a(3) = potential%a(3) + coef * c_x * sy_over_ky * c_z * ct_term%kx / ct_term%ky
            case (hyper_xy_family_x$)
              potential%a(1) = potential%a(1) + coef * s_x * sy_over_ky * s_z
              potential%a(3) = potential%a(3) + coef * c_x * sy_over_ky * c_z * ct_term%kx / ct_term%kz
            case (hyper_x_family_x$)
              potential%a(1) = potential%a(1) + coef * s_x * sy_over_ky * s_z * ct_term%kz / ct_term%kx
              potential%a(3) = potential%a(3) + coef * c_x * sy_over_ky * c_z

            case (hyper_y_family_y$)
              potential%a(2) = potential%a(2) - coef * sx_over_kx * s_y * s_z * ct_term%kz / ct_term%ky
              potential%a(3) = potential%a(3) - coef * sx_over_kx * c_y * c_z
            case (hyper_xy_family_y$)
              potential%a(2) = potential%a(2) - coef * sx_over_kx * s_y * s_z
              potential%a(3) = potential%a(3) - coef * sx_over_kx * c_y * c_z * ct_term%ky / ct_term%kz
            case (hyper_x_family_y$)
              potential%a(2) = potential%a(2) - coef * sx_over_kx * s_y * s_z * ct_term%kz / ct_term%kx
              potential%a(3) = potential%a(3) - coef * sx_over_kx * c_y * c_z * ct_term%ky / ct_term%kx

            case (hyper_y_family_qu$)
              potential%a(1) = potential%a(1) + coef * s_x * c_y * sz_over_kz
              potential%a(2) = potential%a(2) - coef * c_x * s_y * sz_over_kz * ct_term%kx / ct_term%ky
            case (hyper_xy_family_qu$)
              potential%a(1) = potential%a(1) + coef * s_x * c_y * sz_over_kz * ct_term%ky / ct_term%kz
              potential%a(2) = potential%a(2) - coef * c_x * s_y * sz_over_kz * ct_term%kx / ct_term%kz
            case (hyper_x_family_qu$)
              potential%a(1) = potential%a(1) + coef * s_x * c_y * sz_over_kz * ct_term%ky / ct_term%kx
              potential%a(2) = potential%a(2) - coef * c_x * s_y * sz_over_kz

            case (hyper_y_family_sq$)
              potential%a(1) = potential%a(1) + coef * c_x * s_y * sz_over_kz
              potential%a(2) = potential%a(2) + coef * s_x * c_y * sz_over_kz * ct_term%kx / ct_term%ky
            case (hyper_xy_family_sq$)
              potential%a(1) = potential%a(1) + coef * c_x * s_y * sz_over_kz * ct_term%ky / ct_term%kz
              potential%a(2) = potential%a(2) - coef * s_x * c_y * sz_over_kz * ct_term%kx / ct_term%kz
            case (hyper_x_family_sq$)
              potential%a(1) = potential%a(1) + coef * c_x * s_y * sz_over_kz * ct_term%ky / ct_term%kx
              potential%a(2) = potential%a(2) + coef * s_x * c_y * sz_over_kz
            end select

          else  ! electric$
            select case (ct_term%type)
            case (hyper_y_family_x$);   potential%phi = potential%phi + coef * s_x * c_y * c_z / ct_term%ky
            case (hyper_xy_family_x$);  potential%phi = potential%phi + coef * s_y * c_y * c_z / ct_term%kz
            case (hyper_x_family_x$);   potential%phi = potential%phi + coef * s_x * c_y * c_z / ct_term%kx

            case (hyper_y_family_y$);   potential%phi = potential%phi + coef * c_x * s_y * c_z / ct_term%ky
            case (hyper_xy_family_y$);  potential%phi = potential%phi + coef * c_x * s_y * c_z / ct_term%kz
            case (hyper_x_family_y$);   potential%phi = potential%phi + coef * c_x * s_y * c_z / ct_term%kx

            case (hyper_y_family_qu$);   potential%phi = potential%phi + coef * s_x * s_y * c_z / ct_term%ky
            case (hyper_xy_family_qu$);  potential%phi = potential%phi + coef * s_x * s_y * c_z / ct_term%kz
            case (hyper_x_family_qu$);   potential%phi = potential%phi + coef * s_x * s_y * c_z / ct_term%kx

            case (hyper_y_family_sq$);   potential%phi = potential%phi + coef * c_x * c_y * c_z / ct_term%ky
            case (hyper_xy_family_sq$);  potential%phi = potential%phi + coef * c_x * c_y * c_z / ct_term%kz
            case (hyper_x_family_sq$);   potential%phi = potential%phi - coef * c_x * c_y * c_z / ct_term%kx
            end select
          endif

        endif
      enddo

      !

      fld = fld * ct_map%field_scale
      if (ct_map%master_parameter > 0) fld = fld * ele%value(ct_map%master_parameter)
      if (ele%key == sbend$) call restore_curvilinear_field(fld)

      select case (ct_map%field_type)
      case (electric$)
        field%E = field%E + fld
      case (magnetic$)
        field%B = field%B + fld
      case default
        if (global_com%exit_on_error) call err_exit
      end select

      if (do_df_calc) then
        dfld = dfld * ct_map%field_scale
        if (ct_map%master_parameter > 0) dfld = dfld * ele%value(ct_map%master_parameter)
        if (ele%key == sbend$ .and. ele%value(g$) /= 0) then
          rot2(1,:) = [ cos_ang, sin_ang]
          rot2(2,:) = [-sin_ang, cos_ang]
          dfld(1:3:2,1:3:2) = matmul(dfld(1:3:2,1:3:2), rot2)
          rot2(1,2) = -sin_ang
          rot2(2,1) =  sin_ang
          dfld(1:3:2,1:3:2) = matmul(rot2, dfld(1:3:2,1:3:2))
        endif

        select case (ct_map%field_type)
        case (electric$)
          field%dE = field%dE + dfld
        case (magnetic$)
          field%dB = field%dB + dfld
        end select
      endif

    enddo
  endif

  !----------------------------------------------------------------------------
  ! Cylindrical map field

  if (associated(ele%cylindrical_map)) then

    do i = 1, size(ele%cylindrical_map)
      cl_map => ele%cylindrical_map(i)

      if (cl_map%harmonic /= 0) then
        freq0 = ele%value(rf_frequency$)
        freq = ele%value(rf_frequency$) * cl_map%harmonic

        if (freq0 == 0) then
          call out_io (s_fatal$, r_name, 'Element frequency is zero but cylindrical_map harmonic is not in: ' // ele%name)
          if (global_com%exit_on_error) call err_exit
          if (present(err_flag)) err_flag = .true.
          return  
        endif
        t_ref = (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + &
                                             phi0_autoscale + cl_map%phi0_fieldmap) / freq0
        if (ele%key == rfcavity$) t_ref = 0.25/freq0 - t_ref
      endif

      coef = field_autoscale * cl_map%field_scale
      if (cl_map%master_parameter > 0) coef = coef * ele%value(cl_map%master_parameter)

      !

      m = cl_map%m

      if (cl_map%harmonic /= 0) k_t = twopi * freq / c_light

      call to_field_map_coords (local_orb, s_body, cl_map%ele_anchor_pt, cl_map%r0, .false., x, y, z)

      radius = sqrt(x**2 + y**2)
      phi = atan2(y, x)

      E_rho = 0; E_phi = 0; E_z = 0
      B_rho = 0; B_phi = 0; B_z = 0

      do n = 1, size(cl_map%ptr%term)

        cl_term => cl_map%ptr%term(n)
        k_zn = twopi * (n - 1) / (size(cl_map%ptr%term) * cl_map%dz)
        if (2 * n > size(cl_map%ptr%term)) k_zn = k_zn - twopi / cl_map%dz

        cos_ks = cos(k_zn * z)
        sin_ks = sin(k_zn * z)
        exp_kz = cmplx(cos_ks, sin_ks, rp)

        ! DC
        if (cl_map%harmonic == 0) then

          kap_rho = k_zn * radius
          if (m == 0) then
            Im_0    = I_bessel(0, kap_rho)
            Im_plus = I_bessel(1, kap_rho)
            E_rho = E_rho + real(cl_term%e_coef * exp_kz * Im_plus)
            E_z   = E_z   + real(cl_term%e_coef * exp_kz * Im_0 * i_imaginary)
            B_rho = B_rho + real(cl_term%b_coef * exp_kz * Im_plus)
            B_z   = B_z   + real(cl_term%b_coef * exp_kz * Im_0 * i_imaginary)
          else
            cm = exp_kz * cos(m * phi - cl_map%theta0_azimuth)
            sm = exp_kz * sin(m * phi - cl_map%theta0_azimuth)
            Im_plus  = I_bessel(m+1, kap_rho)
            Im_minus = I_bessel(m-1, kap_rho)
            Im_0     = kap_rho * (Im_minus - Im_plus) / (2 * m)

            q = cm * (Im_minus + Im_plus) / 2
            E_rho = E_rho + real(cl_term%e_coef * q)
            B_rho = B_rho + real(cl_term%b_coef * q)

            q = -sm * (Im_minus - Im_plus) / 2
            E_phi = E_phi + real(cl_term%e_coef * q)
            B_phi = B_phi + real(cl_term%b_coef * q)

            q = i_imaginary * cm * Im_0
            E_z = E_z + real(cl_term%e_coef * q)
            B_z = B_z + real(cl_term%b_coef * q)
          endif

          if (present(potential)) then
            if (k_zn == 0) then
              if (m == 0) then
                potential%phi = potential%phi - coef * real(cl_term%e_coef * z * i_imaginary)
              elseif (m == 1) then
                potential%phi = potential%phi - coef * real(cl_term%e_coef * cm * radius / 2)
              endif
            elseif (m == 0) then
              potential%phi = potential%phi - coef * real(cl_term%e_coef * exp_kz * Im_0 / k_zn)
            else
              potential%phi = potential%phi - coef * real(cl_term%e_coef * cm * Im_0 / k_zn)
            endif
          endif

        ! RF mode 
        else
          kappa2_n = k_zn**2 - k_t**2
          kappa_n = sqrt(abs(kappa2_n))
          kap_rho = kappa_n * radius
          if (kappa2_n < 0) then
            kappa_n = -i_imaginary * kappa_n
            kap_rho = -kap_rho
          endif

          if (m == 0) then
            Im_0    = I_bessel_extended(0, kap_rho)
            Im_plus = I_bessel_extended(1, kap_rho) / kappa_n

            E_rho = E_rho - cl_term%e_coef * Im_plus * exp_kz * I_imaginary * k_zn
            E_phi = E_phi + cl_term%b_coef * Im_plus * exp_kz
            E_z   = E_z   + cl_term%e_coef * Im_0    * exp_kz

            B_rho = B_rho - cl_term%b_coef * Im_plus * exp_kz * k_zn
            B_phi = B_phi - cl_term%e_coef * Im_plus * exp_kz * k_t**2 * I_imaginary
            B_z   = B_z   - cl_term%b_coef * Im_0    * exp_kz * I_imaginary

          else
            cm = exp_kz * cos(m * phi - cl_map%theta0_azimuth)
            sm = exp_kz * sin(m * phi - cl_map%theta0_azimuth)
            Im_plus  = I_bessel_extended(m+1, kap_rho) / kappa_n**(m+1)
            Im_minus = I_bessel_extended(m-1, kap_rho) / kappa_n**(m-1)

            ! Reason for computing Im_0_R like this is to avoid divide by zero when radius = 0.
            Im_0_R  = (Im_minus - Im_plus * kappa_n**2) / (2 * m) ! = Im_0 / radius
            Im_0    = radius * Im_0_R       

            E_rho = E_rho - i_imaginary * (k_zn * cl_term%e_coef * Im_plus + cl_term%b_coef * Im_0_R) * cm
            E_phi = E_phi - i_imaginary * (k_zn * cl_term%e_coef * Im_plus + cl_term%b_coef * (Im_0_R - Im_minus / m)) * sm
            E_z   = E_z +                         cl_term%e_coef * Im_0 * cm
     
            B_rho = B_rho + i_imaginary * sm * (cl_term%e_coef * (m * Im_0_R + k_zn**2 * Im_plus) + &
                                          cl_term%b_coef * k_zn * (m * Im_0_R - Im_minus / m))
            B_phi = B_phi + i_imaginary * cm * (cl_term%e_coef * (Im_minus - (k_zn**2 + k_t**2) * Im_plus) / 2 - &
                                          cl_term%b_coef * k_zn * Im_0_R)
            B_z   = B_z +                 sm * (-cl_term%e_coef * k_zn * Im_0 + cl_term%b_coef * kappa2_n * Im_0 / m)

         endif
        endif ! cl_map%harmonic /= 0
          
      enddo  ! cl_map%ptr%term

      ! Notice that phi0, phi0_multipass, and phi0_err are folded into t_ref above.

      if (cl_map%harmonic /= 0) then
        expt = exp(-I_imaginary * twopi * (freq * (time + t_ref)))
        E_rho = E_rho * expt
        E_phi = E_phi * expt
        E_z   = E_z * expt

        expt = expt / (twopi * freq)
        B_rho = B_rho * expt
        B_phi = B_phi * expt
        B_z   = B_z * expt
      endif

      Er = real(E_rho, rp); Ep = real(E_phi, rp); Ez = real(E_z, rp)
      Br = real(B_rho, rp); Bp = real(B_phi, rp); Bz = real(B_z, rp)

      mode_field%E = coef * [cos(phi) * Er - sin(phi) * Ep, sin(phi) * Er + cos(phi) * Ep, Ez]
      mode_field%B = coef * [cos(phi) * Br - sin(phi) * Bp, sin(phi) * Br + cos(phi) * Bp, Bz]

      if (ele%key == sbend$) call restore_curvilinear_field(mode_field%E, mode_field%B)

      field%E = field%E + mode_field%E
      field%B = field%B + mode_field%B

      if (present(potential)) then
      endif

    enddo

  endif

  !----------------------------------------------------------------------------
  ! Grid field calc 

  if (associated(ele%grid_field)) then
  
    ! loop over grid modes

    do i = 1, size(ele%grid_field)
      g_field => ele%grid_field(i)

      if (g_field%harmonic /= 0) then
        freq0 = ele%value(rf_frequency$)
        freq = freq0 * g_field%harmonic
        if (freq0 == 0) then
          call out_io (s_fatal$, r_name, 'ELEMENT FREQUENCY IS ZERO BUT GRID_FIELD HARMONIC IS NOT FOR: ' // ele%name)
          if (global_com%exit_on_error) call err_exit
          if (present(err_flag)) err_flag = .true.
          return  
        endif

        t_ref = (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + &
                                                  phi0_autoscale + g_field%phi0_fieldmap) / freq0
        if (ele%key == rfcavity$) t_ref = 0.25/freq0 - t_ref
      endif

      call to_field_map_coords (local_orb, s_body, g_field%ele_anchor_pt, g_field%r0, g_field%curved_ref_frame, x, y, z)

      ! DC modes should have g_field%harmonic = 0

      expt = field_autoscale * g_field%field_scale
      if (g_field%harmonic /= 0) expt = expt * exp(-I_imaginary * twopi * (freq * (time + t_ref)))
      if (g_field%master_parameter > 0) expt = expt * ele%value(g_field%master_parameter)

      ! calculate field based on grid type
      select case(g_field%geometry)

      case (xyz$)
      
        call grid_field_linear_interpolate(ele, local_orb, g_field, g_pt, err, x, y, z, &
                                allow_s_out_of_bounds = logic_option(.false., grid_allow_s_out_of_bounds))
        if (err) then
          if (present(err_flag)) err_flag = .true.
          return
        endif

        mode_field%e = real(expt * g_pt%e)
        mode_field%b = real(expt * g_pt%B)

      case(rotationally_symmetric_rz$)
        
        ! Format should be: pt (ir, iz) = ( Er, 0, Ez, 0, Bphi, 0 ) 
          
        ! Interpolate 2D (r, z) grid
        ! g_pt is a grid_field_pt_struct, which has complex E and B

        r = sqrt(x**2 + y**2)

        call grid_field_linear_interpolate(ele, local_orb, g_field, g_pt, err, r, z, &
                                allow_s_out_of_bounds = logic_option(.false., grid_allow_s_out_of_bounds))
        if (err) then
          if (global_com%exit_on_error) call err_exit
          if (present(err_flag)) err_flag = .true.
          return
        endif

        ! Transverse field is zero on axis. Otherwise:

        if (r /= 0) then
          ! Get non-rotated field
          E_rho = real(expt * g_pt%E(1))
          E_phi = real(expt * g_pt%E(2))
          B_rho = real(expt * g_pt%B(1)) 
          B_phi = real(expt * g_pt%B(2))

          ! rotate field and output Ex, Ey, Bx, By
          mode_field%e(1) = (x*E_rho - y*E_phi)/r
          mode_field%e(2) = (y*E_rho + x*E_phi)/r
          mode_field%b(1) = (x*B_rho - y*B_phi)/r
          mode_field%b(2) = (y*B_rho + x*B_phi)/r
        endif
    
        ! Ez, Bz 
        mode_field%e(3) = real(expt*g_pt%E(3))
        mode_field%b(3) = real(expt*g_pt%B(3)) 
    
      case default
        call out_io (s_fatal$, r_name, 'UNKNOWN GRID GEOMETRY: \i0\ ', &
                                       'FOR ELEMENT: ' // ele%name, i_array = [g_field%geometry])
        if (global_com%exit_on_error) call err_exit
        if (present(err_flag)) err_flag = .true.
        return
      end select
      
      if (ele%key == sbend$ .and. .not. g_field%curved_ref_frame) call restore_curvilinear_field(mode_field%E, mode_field%B)

      field%E = field%E + mode_field%E
      field%B = field%B + mode_field%B

    enddo
  endif

  !----------------------------------------------------------------------------
  ! Taylor field calc 

  if (associated(ele%taylor_field)) then
  
    ! loop over taylor modes

    do i = 1, size(ele%taylor_field)
      t_field => ele%taylor_field(i)

      fld = 0

      call to_field_map_coords (local_orb, s_body, t_field%ele_anchor_pt, t_field%r0, t_field%curved_ref_frame, x, y, z)

      iz0 = lbound(t_field%ptr%plane, 1)
      iz1 = ubound(t_field%ptr%plane, 1)
      z_center = (iz0 + iz1) * t_field%dz / 2
      if (abs(z - z_center) > (iz1 - iz0) * t_field%dz .and. &
                        .not. logic_option(.false., grid_allow_s_out_of_bounds)) then
        call out_io (s_error$, r_name, 'PARTICLE Z  \F10.3\ POSITION OUT OF BOUNDS.', &
                                       'FOR TAYLOR_FIELD IN ELEMENT: ' // ele%name, r_array = [s_body])
        return
      endif

      izp = floor(z / t_field%dz)
      if (izp < iz0 - 1 .or. izp > iz1) cycle ! Outside of defined field region field is assumed zero.

      ! Taylor upsteam of particle

      if (izp == iz0 - 1) then
        fld0 = 0
        dfld0 = 0
      else
        call evaluate_em_taylor ([x, y], t_field%ptr%plane(izp)%field, fld0, dfld0)
      endif

      ! Taylor downstream of particle

      if (izp == iz1) then
        fld1 = 0
        dfld1 = 0
      else
        call evaluate_em_taylor ([x, y], t_field%ptr%plane(izp+1)%field, fld1, dfld1)
      endif

      ! Interpolate

      do j = 1, 3
        call create_a_spline (spline, [0.0_rp, fld0(j)], [t_field%dz, fld1(j)], dfld0(j,3), dfld1(j,3))
        fld(j) = spline1 (spline, z - izp*t_field%dz)
      enddo

      !

      fld = fld * t_field%field_scale
      if (t_field%master_parameter > 0) fld = fld * ele%value(t_field%master_parameter)
      if (ele%key == sbend$ .and. .not. t_field%curved_ref_frame) call restore_curvilinear_field(fld)

      select case (t_field%field_type)
      case (electric$)
        field%E = field%E + fld
      case (magnetic$)
        field%B = field%B + fld
      case default
        if (global_com%exit_on_error) call err_exit
      end select

    enddo

  endif

! Beginning_ele, for example, has no field

case (no_field$)

  return

! Unknown field_calc

case default
  call out_io (s_fatal$, r_name, 'BAD FIELD_CALC METHOD FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  if (present(err_flag)) err_flag = .true.
  return
end select

! Scale ac_kicker element field

if (ele%key == ac_kicker$) then

  if (present(rf_time)) then
    a_amp = ac_kicker_amp (ele, rf_time)
  else
    a_amp = ac_kicker_amp (ele, particle_rf_time(orbit, ele, .true., s_body))
  endif

  field%E = a_amp * field%E
  field%B = a_amp * field%B
  field%dE = a_amp * field%dE
  field%dB = a_amp * field%dB
endif

!----------------------------------------------------------------------------
! overlapping of fields from other elements

8000 continue

if (ele%n_lord_field /= 0 .and. logic_option(.true., use_overlap)) then
  lab_orb = orbit
  if (local_ref_frame) then
    call offset_particle (ele, param, unset$, lab_orb, set_hvkicks = .false., s_pos = s_body, s_out = s_lab)
  else
    s_lab = s_body
  endif

  lab_position%r = [lab_orb%vec(1), lab_orb%vec(3), s_lab]
  global_position = coords_local_curvilinear_to_floor (lab_position, ele, w_mat = w_ele_mat, calculate_angles = .false.)

  lord_orb = lab_orb
  do i = 1, ele%n_lord_field
    lord => pointer_to_lord(ele, i, field_overlap_ptr = .true.)
    lord_position = coords_floor_to_local_curvilinear (global_position, lord, status, w_lord_mat)
    lord_orb%vec(1) = lord_position%r(1)
    lord_orb%vec(3) = lord_position%r(2)
    ! Set use_overlap = False to prevent recursion.
    call em_field_calc (lord, param, lord_position%r(3), lord_orb, .false., l1_field, calc_dfield, err, &
                        use_overlap = .false., grid_allow_s_out_of_bounds = .true., used_eles = used_eles)
    if (err) then
      if (present(err_flag)) err_flag = .true.
      return
    endif
    ! Field in lord lab coords to field in global coords
    call rotate_em_field (l1_field, transpose(w_lord_mat), w_lord_mat, calc_dfield)
    if (i == 1) then
      lord_field = l1_field
    else
      lord_field = lord_field + l1_field
    endif
  enddo

  ! Field in global coords to field in lab coords

  call rotate_em_field (lord_field, transpose(w_ele_mat), transpose(w_ele_mat))

  if (local_ref_frame) then
    call convert_field_ele_to_lab (ele, s_lab, .false., lord_field)  ! lab -> ele
    field = field + lord_field
  else
    call convert_field_ele_to_lab (ele, s_body, .true., field)
    field = field + lord_field
  endif

  return
endif

! Final

if (.not. local_ref_frame) call convert_field_ele_to_lab (ele, s_body, .true., field)

if (do_df_calc .and. .not. dfield_computed) then
  call em_field_derivatives (ele, param, s_pos, orbit, local_ref_frame, field, rf_time)
endif

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Convert fields: ele to lab coords

contains

subroutine convert_field_ele_to_lab (ele, s_here, forward_transform, field)

type (ele_struct) ele
type (em_field_struct) field

real(rp) s_here, w_mat(3,3), w_inv(3,3), w_s(3,3), w_rt(3,3), w_rt_inv(3,3)
real(rp) theta
logical forward_transform

!

if (ele%key == sbend$) then
  call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(roll$), w_mat)
  theta = ele%value(g$) * s_here - ele%value(angle$)/2
  w_s = w_mat_for_x_pitch (theta)
  if (ele%value(ref_tilt_tot$) == 0) then
    w_mat = matmul(matmul(w_s, w_mat), transpose(w_s))
  else
    w_rt = w_mat_for_tilt (ele%value(ref_tilt_tot$))
    w_rt_inv = w_mat_for_tilt (ele%value(ref_tilt_tot$), .true.)
    w_mat = matmul(matmul(matmul(matmul(matmul(w_rt, w_s), w_rt_inv), w_mat), w_rt), transpose(w_s))
  endif
  w_inv = transpose(w_mat)
else
  call floor_angles_to_w_mat (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), w_mat, w_inv)
endif

if (forward_transform) then
  call rotate_em_field (field, w_mat, w_inv, calc_dfield, potential)
else
  call rotate_em_field (field, w_inv, w_mat, calc_dfield, potential)
endif

end subroutine convert_field_ele_to_lab

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine to_field_map_coords (local_orb, s_body, ele_anchor_pt, r0, curved_ref_frame, x, y, z)

type (coord_struct) local_orb

real(rp) :: s_body, r0(3), x, y, z, x_save
integer ele_anchor_pt
logical curved_ref_frame

!

select case (ele_anchor_pt)
case (anchor_beginning$); s0 = 0
case (anchor_center$);    s0 = ele%value(l$) / 2
case (anchor_end$);       s0 = ele%value(l$)
case default
  call out_io (s_fatal$, r_name, 'BAD ELE_ANCHOR_PT FOR FIELD GRID IN ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  if (present(err_flag)) err_flag = .true.
  return
end select


!

x = local_orb%vec(1)
z = s_body - s0

!

if (ele%key == sbend$ .and. ele%value(g$) /= 0 .and. .not. curved_ref_frame) then
  cos_ang = cos(z*ele%value(g$))
  sin_ang = sin(z*ele%value(g$))

  x_save = x
  x = (x_save + ele%value(rho$) )*cos_ang - ele%value(rho$)
  z = (x_save + ele%value(rho$) )*sin_ang 
endif

!

x = x - r0(1)
y = local_orb%vec(3) - r0(2)
z = z - r0(3)

end subroutine to_field_map_coords

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! restore_curvilinear_field(field_a, field_b)
!
! For sbend with Grid calculation.

subroutine restore_curvilinear_field(field_a, field_b)

real(rp) temp, field_a(3)
real(rp), optional :: field_b(3)

! For sbend with Grid calculation Restores x and s_body, and rotates output fields.

if (ele%value(g$) == 0) return

temp       = field_a(3)*cos_ang - field_a(1)*sin_ang
field_a(1) = field_a(3)*sin_ang + field_a(1)*cos_ang
field_a(3) = temp

if (present(field_b)) then
  temp       = field_b(3)*cos_ang - field_b(1)*sin_ang
  field_b(1) = field_b(3)*sin_ang + field_b(1)*cos_ang
  field_b(3) = temp 
endif

end subroutine restore_curvilinear_field

end subroutine em_field_calc 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine rotate_em_field (field, w_mat, w_inv, calc_dfield, potential)
!
! Routine to transform the fields using the given rotation matrices.
!
! Input:
!   field       -- em_field_struct: E and B fields and derivatives.
!   w_mat(3,3)  -- real(rp): rotation matrix.
!   w_inv(3,3)  -- real(rp): rotation matrix inverse = transpose(w_mat)
!   calc_dfield -- Logical, optional: If present and True then calculate the field derivatives.
!   potential   -- em_potential_struct, optional: The electric and magnetic potentials.
!
! Output:
!   field       -- em_field_struct: E and B fields and derivatives.
!-

subroutine rotate_em_field (field, w_mat, w_inv, calc_dfield, potential)

type (em_field_struct) field
type (em_potential_struct), optional :: potential

real(rp) w_mat(3,3), w_inv(3,3)
logical, optional :: calc_dfield

!

field%B = matmul(w_mat, field%B)
field%E = matmul(w_mat, field%E)

if (present(potential)) then
  potential%A = matmul(w_mat, potential%A)
endif

if (logic_option (.false., calc_dfield)) then
  field%dB = matmul(w_mat, matmul(field%dB, w_inv))
  field%dE = matmul(w_mat, matmul(field%dE, w_inv))
endif

end subroutine rotate_em_field

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine grid_field_linear_interpolate (ele, orbit, grid, field, err_flag, x1, x2, x3, allow_s_out_of_bounds)
!
! Subroutine to interpolate the E and B fields on a rectilinear grid.
!
! Input:
!   ele      -- ele_struct: Element containing the grid.
!   orbit    -- coord_struct: Used for constructing an error message if the particle is out of bounds.
!   grid     -- grid_field_struct: Grid to interpolate.
!   err_flag -- Logical: Set to true if there is an error. False otherwise.
!   x1       -- real(rp): dimension 1 interpolation point.
!   x2       -- real(rp), optional: dimension 2 interpolation point.
!   x3       -- real(rp), optional: dimension 3 interpolation point.
!   allow_s_out_of_bounds -- logical, optional: allow s-coordinate grossly out of bounds to return
!                 zero field without an error. This is used when the field of one element overlaps
!                 the field of another.
!
! Output:
!   field    -- grid_field_pt_struct: Interpolated field (complex)
!-

subroutine grid_field_linear_interpolate (ele, orbit, grid, g_field, err_flag, x1, x2, x3, allow_s_out_of_bounds)

type (ele_struct) ele
type (coord_struct) orbit
type (grid_field_struct) :: grid
type (grid_field_pt1_struct), intent(out) :: g_field
real(rp) :: x1
real(rp), optional :: x2, x3
real(rp) rel_x1, rel_x2, rel_x3
integer i1, i2, i3, grid_dim, allow_s, lbnd, ubnd
logical err_flag
logical :: allow_s_out_of_bounds

character(32), parameter :: r_name = 'grid_field_linear_interpolate'

integer, parameter :: allow_tiny$ = 1, allow_some$ = 2, allow_all$ = 3

! Pick appropriate dimension 

err_flag = .false.

allow_s = allow_some$
if (allow_s_out_of_bounds) allow_s = allow_all$

grid_dim = grid_field_dimension(grid%geometry)

! xz grid

select case(grid_dim)
case (2)

  lbnd = lbound(grid%ptr%pt, 2); ubnd = ubound(grid%ptr%pt, 2)

  call get_this_index(x2, 2, i2, rel_x2, err_flag, allow_s); if (err_flag) return
  ! If grossly out of longitudinal bounds just return zero field. Do not test transverse position in this case.
  if (i2 < lbnd - 1 .or. i2 > ubnd) return 

  call get_this_index(x1, 1, i1, rel_x1, err_flag, allow_tiny$); if (err_flag) return

  ! Do bilinear interpolation. If just outside longitudinally, interpolate between grid edge and zero.

  if (i2 == lbnd - 1) then  ! Just outside entrance end
    g_field%E(:) = (1-rel_x1)*(rel_x2)   * grid%ptr%pt(i1,   i2+1, 1)%E(:) &
                 + (rel_x1)*(rel_x2)     * grid%ptr%pt(i1+1, i2+1, 1)%E(:) 

    g_field%B(:) = (1-rel_x1)*(rel_x2)   * grid%ptr%pt(i1,   i2+1, 1)%B(:) &
                 + (rel_x1)*(rel_x2)     * grid%ptr%pt(i1+1, i2+1, 1)%B(:)  

  elseif (i2 == ubnd) then  ! Just outside exit end
    g_field%E(:) = (1-rel_x1)*(1-rel_x2) * grid%ptr%pt(i1,   i2,   1)%E(:) &
                 + (rel_x1)*(1-rel_x2)   * grid%ptr%pt(i1+1, i2,   1)%E(:)

    g_field%B(:) = (1-rel_x1)*(1-rel_x2) * grid%ptr%pt(i1,   i2,   1)%B(:) &
                 + (rel_x1)*(1-rel_x2)   * grid%ptr%pt(i1+1, i2,   1)%B(:)

  elseif (lbnd <= i2 .and. i2 < ubnd) then   ! Inside 
    g_field%E(:) = (1-rel_x1)*(1-rel_x2) * grid%ptr%pt(i1,   i2,   1)%E(:) &
                 + (1-rel_x1)*(rel_x2)   * grid%ptr%pt(i1,   i2+1, 1)%E(:) &
                 + (rel_x1)*(1-rel_x2)   * grid%ptr%pt(i1+1, i2,   1)%E(:) &
                 + (rel_x1)*(rel_x2)     * grid%ptr%pt(i1+1, i2+1, 1)%E(:) 

    g_field%B(:) = (1-rel_x1)*(1-rel_x2) * grid%ptr%pt(i1,   i2,   1)%B(:) &
                 + (1-rel_x1)*(rel_x2)   * grid%ptr%pt(i1,   i2+1, 1)%B(:) &
                 + (rel_x1)*(1-rel_x2)   * grid%ptr%pt(i1+1, i2,   1)%B(:) &
                 + (rel_x1)*(rel_x2)     * grid%ptr%pt(i1+1, i2+1, 1)%B(:)  
  endif

! xyz grid

case (3)

  lbnd = lbound(grid%ptr%pt, 3); ubnd = ubound(grid%ptr%pt, 3)

  call get_this_index(x3, 3, i3, rel_x3, err_flag, allow_s); if (err_flag) return
  ! If grossly out of longitudinal bounds just return zero field. Do not test transverse position in this case.
  if (i3 < lbnd - 1 .or. i3 > ubnd) return 

  call get_this_index(x1, 1, i1, rel_x1, err_flag, allow_tiny$); if (err_flag) return
  call get_this_index(x2, 2, i2, rel_x2, err_flag, allow_tiny$); if (err_flag) return
    
  ! Do trilinear interpolation. If just outside longitudinally, interpolate between grid edge and zero.

  if (i3 == lbnd - 1) then  ! Just outside entrance end
    g_field%E(:) = (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1,   i2,   i3+1)%E(:) &
                 + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1,   i2+1, i3+1)%E(:) &
                 + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1+1, i2,   i3+1)%E(:) &
                 + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1+1, i2+1, i3+1)%E(:)               
               
    g_field%B(:) = (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1,   i2,   i3+1)%B(:) &
                 + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1,   i2+1, i3+1)%B(:) &
                 + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1+1, i2,   i3+1)%B(:) &
                 + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1+1, i2+1, i3+1)%B(:)

  elseif (i3 == ubnd) then  ! Just outside exit end
    g_field%E(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1,   i2,   i3  )%E(:) &
                 + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1,   i2+1, i3  )%E(:) &
                 + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1+1, i2,   i3  )%E(:) &
                 + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1+1, i2+1, i3  )%E(:)
               
    g_field%B(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1,   i2,   i3  )%B(:) &
                 + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1,   i2+1, i3  )%B(:) &
                 + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1+1, i2,   i3  )%B(:) &
                 + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1+1, i2+1, i3  )%B(:) 

  elseif (lbnd <= i3 .and. i3 < ubnd) then   ! Inside
    g_field%E(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1,   i2,   i3  )%E(:) &
                 + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1,   i2+1, i3  )%E(:) &
                 + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1+1, i2,   i3  )%E(:) &
                 + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1+1, i2+1, i3  )%E(:) &
                 + (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1,   i2,   i3+1)%E(:) &
                 + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1,   i2+1, i3+1)%E(:) &
                 + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1+1, i2,   i3+1)%E(:) &
                 + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1+1, i2+1, i3+1)%E(:)               
               
    g_field%B(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1,   i2,   i3  )%B(:) &
                 + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1,   i2+1, i3  )%B(:) &
                 + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1+1, i2,   i3  )%B(:) &
                 + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1+1, i2+1, i3  )%B(:) &
                 + (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1,   i2,   i3+1)%B(:) &
                 + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1,   i2+1, i3+1)%B(:) &
                 + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1+1, i2,   i3+1)%B(:) &
                 + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1+1, i2+1, i3+1)%B(:) 
  endif

case default
  call out_io (s_fatal$, r_name, 'BAD DIMENSION: \i0\ ', grid_field_dimension(grid%geometry))
  if (global_com%exit_on_error) call err_exit
  err_flag = .true.
  return
end select

!-------------------------------------------------------------------------------------
contains

subroutine get_this_index (x, ix_x, i0, rel_x0, err_flag, allow_out_of_bounds)

real(rp) x, rel_x0, x_norm, x_diff, x_ave
integer ix_x, i0, ig0, ig1, allow_out_of_bounds
logical err_flag

!

ig0 = lbound(grid%ptr%pt, ix_x)
ig1 = ubound(grid%ptr%pt, ix_x)

x_norm = x / grid%dr(ix_x)  ! Note that to_field_map_coords has already been called.
i0 = floor(x_norm)          ! index of lower 1 data point
rel_x0 = x_norm - i0        ! Relative distance from lower x1 grid point

! Out of bounds?

if (i0 < ig0 .or. i0 >= ig1) then
  g_field%E = 0
  g_field%B = 0

  select case (allow_out_of_bounds)
  case (allow_tiny$)
    ! Here do extrapolation is the point is within one dr/2 of the grid boundary.
    ! Otherwise it is an error.
    if (i0 == ig0 - 1 .and. rel_x0 > 0.5) then
      i0 = ig0
      rel_x0 = rel_x0 - 1
      return
    elseif (i0 == ig1 .and. rel_x0 < 0.5) then
      i0 = ig1 - 1
      rel_x0 = rel_x0 + 1
      return
    endif

  case (allow_some$)
    ! Here only generate an error message if the particle is grossly outside of the grid region.
    ! Here "gross" is defined as dOut > L_grid/2 where dOut is the distance between the
    ! particle and the grid edge and L_grid is the length of the grid.
    x_diff = (ig1 - ig0) * grid%dr(ix_x)
    x_ave = (ig1 + ig0) * grid%dr(ix_x) / 2
    if (abs(x - x_ave) < x_diff .or. i0 == ig0-1 .or. i0 == ig1) return

  case (allow_all$)
    return
  end select

  err_flag = .true.
  call out_io (s_error$, r_name, '\i0\D GRID_FIELD INTERPOLATION INDEX OUT OF BOUNDS: I\i0\ = \i0\ (POSITION = \f12.6\)', &
                                 'FOR ELEMENT: ' // ele%name // '  (' // trim(ele_loc_to_string(ele)) // ')', &
                                 'PARTICLE POSITION: \3F12.6\ ', &
                                 'SETTING FIELD TO ZERO', i_array = [grid_dim, ix_x, i0], &
                                 r_array = [x, orbit%vec(1), orbit%vec(3), orbit%s-ele%s_start])
endif

end subroutine get_this_index 

end subroutine grid_field_linear_interpolate

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Function field_interpolate_3d (position, field_mesh, deltas, position0) result (field)
!
! Function to interpolate a 3d field.
! The interpolation is such that the derivative is continuous.
!
! Note: For "interpolation" outside of the region covered by the field_mesh
! it is assumed that the field is constant, Equal to the field at the
! boundary.
!
! Module needed:
!   use em_field_mod
!
! Input:
!   position(3)       -- Real(rp): (x, y, z) position.
!   field_mesh(:,:,:) -- Real(rp): Grid of field points.
!   deltas(3)         -- Real(rp): (dx, dy, dz) distances between mesh points.
!   position0(3)      -- Real(rp), optional:  position at (ix0, iy0, iz0) where
!                            (ix0, iy0, iz0) is the lower bound of the
!                            filed_mesh(i, j, k) array. If not present then
!                            position0 is taken to be (0.0, 0.0, 0.0)
! Output:
!   field -- Real(rp): interpolated field.
!-

function field_interpolate_3d (position, field_mesh, deltas, position0) result (field)

real(rp), optional, intent(in) :: position0(3)
real(rp), intent(in) :: position(3), field_mesh(0:,0:,0:), deltas(3)
real(rp) field

real(rp) r(3), f(-1:2), g(-1:2), h(-1:2), r_frac(3)

integer i0(3), ix, iy, iz, iix, iiy, iiz

!

if (present(position0)) then
  r = (position - position0) / deltas
else
  r = position / deltas
endif

i0 = int(r)
r_frac = r - i0

do ix = -1, 2
 iix = min(max(ix + i0(1), 0), ubound(field_mesh, 1))
 do iy = -1, 2
    iiy = min(max(iy + i0(2), 0), ubound(field_mesh, 2))
    do iz = -1, 2
      iiz = min(max(iz + i0(3), 0), ubound(field_mesh, 3))
      f(iz) = field_mesh(iix, iiy, iiz)
    enddo
    g(iy) = interpolate_1d (r_frac(3), f)
  enddo
  h(ix) = interpolate_1d (r_frac(2), g)
enddo
field = interpolate_1d (r_frac(1), h)

!---------------------------------------------------------------

contains

! interpolation in 1 dimension using 4 equally spaced points: P1, P2, P3, P4.
!   x = interpolation point.
!           x = 0 -> point is at P2.
!           x = 1 -> point is at P3.
! Interpolation is done so that the derivative is continuous.
! The interpolation uses a cubic polynomial

function interpolate_1d (x, field1_in) result (field1)

real(rp) field1, x, field1_in(4), df_2, df_3
real(rp) c0, c1, c2, c3

!

df_2 = (field1_in(3) - field1_in(1)) / 2   ! derivative at P2
df_3 = (field1_in(4) - field1_in(2)) / 2   ! derivative at P3

c0 = field1_in(2)
c1 = df_2
c2 = 3 * field1_in(3) - df_3 - 3 * field1_in(2) - 2 * df_2
c3 = df_3 - 2 * field1_in(3) + 2 * field1_in(2) + df_2

field1 = c0 + c1 * x + c2 * x**2 + c3 * x**3

end function interpolate_1d

end function field_interpolate_3d 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine em_field_derivatives (ele, param, s_pos, orbit, local_ref_frame, dfield, rf_time)
!
! Routine to calculate field derivatives.
! In theory this should be handled by em_filed_calc. In practice, em_field_calc is currently incomplete.
!
! Input
!   ele             -- Ele_struct: Element
!   param           -- lat_param_struct: Lattice parameters.
!   s_pos           -- Real(rp): Longitudinal position relative to the upstream edge of the element.
!   time            -- Real(rp): Particle time.
!                       For absolute time tracking this is the absolute time.
!                       For relative time tracking this is relative to the reference particle entering the element.
!   orbit           -- Coord_struct: Transverse coordinates.
!     %vec(1), %vec(3)  -- Transverse coords. These are the only components used in the calculation.
!   local_ref_frame     -- Logical, If True then take the input coordinates and output fields 
!                                   as being with respect to the frame of referene of the element (ignore misalignments). 
!   rf_time         -- real(rp), optional: RF clock time. If not present then the time will be calculated using the standard algorithm.
!
! Output:
!   dfield       -- em_field_struct: E and B field derivatives. dfield%E and dfield%B are not touched.
!-

subroutine em_field_derivatives (ele, param, s_pos, orbit, local_ref_frame, dfield, rf_time)

type (ele_struct), target :: ele
type (lat_param_struct) param
type (em_field_struct) :: dfield, f0, f1
type (coord_struct) :: orbit, orb

real(rp), optional :: rf_time
real(rp) s_pos, s0, s1, del
logical local_ref_frame

!

orb = orbit
del = bmad_com%d_orb(1)

orb%vec(1) = orbit%vec(1) - del
call em_field_calc (ele, param, s_pos, orb, .true., f0, rf_time = rf_time)
orb%vec(1) = orbit%vec(1) + del
call em_field_calc (ele, param, s_pos, orb, .true., f1, rf_time = rf_time)

dfield%dB(:,1) = (f1%B - f0%B) / (2 * del)
dfield%dE(:,1) = (f1%E - f0%E) / (2 * del)

!

orb = orbit
del = bmad_com%d_orb(3)

orb%vec(3) = orbit%vec(3) - del
call em_field_calc (ele, param, s_pos, orb, .true., f0, rf_time = rf_time)
orb%vec(3) = orbit%vec(3) + del
call em_field_calc (ele, param, s_pos, orb, .true., f1, rf_time = rf_time)

dfield%dB(:,2) = (f1%B - f0%B) / (2 * del)
dfield%dE(:,2) = (f1%E - f0%E) / (2 * del)

!

orb = orbit
del = bmad_com%d_orb(5)

call em_field_calc (ele, param, s_pos-del, orbit, .true., f0, rf_time = rf_time)
call em_field_calc (ele, param, s_pos+del, orbit, .true., f1, rf_time = rf_time)

dfield%dB(:,3) = (f1%B - f0%B) / (2 * del)
dfield%dE(:,3) = (f1%E - f0%E) / (2 * del)

end subroutine em_field_derivatives

end module
