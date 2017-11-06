!+
! Subroutine track_a_lcavity (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through a lcavity element.
!
! Modified version of:
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
! with b_0 = b_-1 = 1. See the Bmad manual for more details.
!
! One must keep in mind that we are NOT using good canonical coordinates since
!   the energy of the reference particle is changing.
! This means that the resulting matrix will NOT be symplectic.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Thick multipole element.
!   param       -- lat_param_struct: Lattice parameters.
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before the element.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_lcavity (orbit, ele, param, mat6, make_matrix)

use track1_mod, except_dummy => track_a_lcavity

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) length, pc_start, pc_end, gradient_ref, gradient_max, dz_factor, rel_p, coef, k2
real(rp) E_start_ref, E_end_ref, pc_start_ref, pc_end_ref, alpha, sin_a, cos_a, r_mat(2,2), volt_ref
real(rp) phase, cos_phi, sin_phi, gradient_net, e_start, e_end, e_ratio, voltage_max, dp_dg, sqrt_8, f, k1
real(rp) dE_start, dE_end, dE, beta_start, beta_end, beta_start_ref, beta_end_ref
real(rp) pxy2, xp1, xp2, yp1, yp2, mc2, om, om_g, m2(2,2), kmat(6,6)
real(rp) dbeta1_dE1, dbeta2_dE2, dalpha_dt1, dalpha_dE1, dcoef_dt1, dcoef_dE1, z21, z22
real(rp) c_min, c_plu, dc_min, dc_plu, cos_term, dcos_phi, drp1_dr0, drp1_drp0, drp2_dr0, drp2_drp0
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer ix_pole_max, ix_elec_max

logical, optional :: make_matrix

!

length = ele%value(l$)
if (length == 0) return

! 

call multipole_ele_to_ab (ele, .false., ix_pole_max, an,      bn,      magnetic$, include_kicks = .true.)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

call offset_particle (ele, param, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, orbit, magnetic$, 1.0_rp/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, orbit, electric$, ele%value(l$)/2, mat6, make_matrix)

! The RF phase is defined with respect to the time at the beginning of the element.
! So if dealing with a slave element and absolute time tracking then need to correct.

phase = twopi * (ele%value(phi0_err$) + ele%value(phi0_autoscale$) + &
           ele%value(phi0$) + ele%value(phi0_multipass$) + &
           (particle_rf_time (orbit, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))

call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit, mat6, make_matrix)

!

rel_p = 1 + orbit%vec(6)
E_start_ref  = ele%value(E_tot_start$)
E_end_ref    = ele%value(E_tot$)
gradient_ref = (E_end_ref - E_start_ref) / length
pc_start_ref = ele%value(p0c_start$)
pc_end_ref   = ele%value(p0c$)
beta_start_ref = pc_start_ref / E_start_ref
beta_end_ref   = pc_end_ref / E_end_ref
mc2 = mass_of(orbit%species)

pc_start = pc_start_ref * rel_p
beta_start = orbit%beta
E_start = pc_start / beta_start 


gradient_max = e_accel_field(ele, gradient$)

cos_phi = cos(phase)
sin_phi = sin(phase)
gradient_net = gradient_max * cos_phi + gradient_shift_sr_wake(ele, param)

dE = gradient_net * length
E_end = E_start + dE
if (E_end <= mass_of(orbit%species)) then
  orbit%state = lost_pz_aperture$
  orbit%vec(6) = -1.01  ! Something less than -1
  if (present(mat6)) mat6 = 0
  return
endif

call convert_total_energy_to (E_end, orbit%species, pc = pc_end, beta = beta_end)
E_ratio = E_end / E_start
mc2 = mass_of(orbit%species)

! Body tracking longitudinal

if (abs(dE) <  1d-4*(pc_end+pc_start)) then
  dp_dg = length * (1 / beta_start - mc2**2 * dE / (2 * pc_start**3) + (mc2 * dE)**2 * E_start / (2 * pc_start**5))
else
  dp_dg = (pc_end - pc_start) / gradient_net
endif

if (logic_option(.false., make_matrix)) then
  om = twopi * ele%value(rf_frequency$) / c_light
  om_g = om * gradient_max * length
  dbeta1_dE1 = mc2**2 / (pc_start * E_start**2)
  dbeta2_dE2 = mc2**2 / (pc_end * E_end**2)

  ! Convert from (x, px, y, py, z, pz) to (x, x', y, y', c(t_ref-t), E) coords 
  mat6(2,:) = mat6(2,:) / rel_p - orbit%vec(2) * mat6(6,:) / rel_p**2
  mat6(4,:) = mat6(4,:) / rel_p - orbit%vec(4) * mat6(6,:) / rel_p**2

  m2(1,:) = [1/orbit%beta, -orbit%vec(5) * mc2**2 * orbit%p0c / (pc_start**2 * E_start)]
  m2(2,:) = [0.0_rp, orbit%p0c * orbit%beta]
  mat6(5:6,:) = matmul(m2, mat6(5:6,:))

  ! longitudinal track
  kmat = 0
  kmat(6,5) = om_g * sin_phi
  kmat(6,6) = 1

  if (abs(dE) <  1d-4*(pc_end+pc_start)) then
    kmat(5,5) = 1 - length * (-mc2**2 * kmat(6,5) / (2 * pc_start**3) + mc2**2 * dE * kmat(6,5) * E_start / pc_start**5)
    kmat(5,6) = -length * (-dbeta1_dE1 / beta_start**2 + 2 * mc2**2 * dE / pc_start**4 + &
                    (mc2 * dE)**2 / (2 * pc_start**5) - 5 * (mc2 * dE)**2 / (2 * pc_start**5))
  else
    kmat(5,5) = 1 - kmat(6,5) / (beta_end * gradient_net) + kmat(6,5) * (pc_end - pc_start) / (gradient_net**2 * length)
    kmat(5,6) = -1 / (beta_end * gradient_net) + 1 / (beta_start * gradient_net)
  endif
endif


! Convert to (x', y', c(t_ref-t), E) coords

orbit%vec(2) = orbit%vec(2) / rel_p    ! Convert to x'
orbit%vec(4) = orbit%vec(4) / rel_p    ! Convert to y'
orbit%vec(5) = orbit%vec(5) / beta_start
orbit%vec(6) = rel_p * orbit%p0c / orbit%beta - 1
orbit%t = orbit%t + dp_dg / c_light

! Body tracking. Transverse kick only happens with standing wave cavities.

! Traveling wave
if (nint(ele%value(cavity_type$)) == traveling_wave$) then

  if (logic_option(.false., make_matrix)) then
    kmat(1,1) = 1
    kmat(1,2) = length
    kmat(2,2) = 1

    kmat(3,3) = 1
    kmat(3,4) = length
    kmat(4,4) = 1

    kmat(5,2) = -length * orbit%vec(4)
    kmat(5,4) = -length * orbit%vec(4)

    mat6 = matmul(kmat, mat6)
  endif

  orbit%vec(1) = orbit%vec(1) + orbit%vec(2) * length
  orbit%vec(3) = orbit%vec(3) + orbit%vec(4) * length

  dz_factor = (orbit%vec(2)**2 + orbit%vec(4)**2) * dp_dg / 2
  orbit%vec(5) = orbit%vec(5) - dz_factor
  orbit%t = orbit%t + dz_factor / c_light

! Standing wave
else
  sqrt_8 = 2 * sqrt_2
  voltage_max = gradient_max * length

  if (abs(voltage_max * cos_phi) < 1d-5 * E_start) then
    f = voltage_max / E_start
    alpha = f * (1 + f * cos_phi / 2)  / sqrt_8
    coef = length * beta_start * (1 - voltage_max * cos_phi / (2 * E_start))
  else
    alpha = log(E_ratio) / (sqrt_8 * cos_phi)
    coef = sqrt_8 * pc_start * sin(alpha) / gradient_max
  endif

  cos_a = cos(alpha)
  sin_a = sin(alpha)

  r_mat(1,1) =  cos_a
  r_mat(1,2) =  coef 
  r_mat(2,1) = -sin_a * gradient_max / (sqrt_8 * pc_end)
  r_mat(2,2) =  cos_a * pc_start / pc_end

  if (logic_option(.false., make_matrix)) then
    if (abs(voltage_max * cos_phi) < 1d-5 * E_start) then
      dalpha_dt1 = f * f * om * sin_phi / (2 * sqrt_8) 
      dalpha_dE1 = -(voltage_max / E_start**2 + voltage_max**2 * cos_phi / E_start**3) / sqrt_8
      dcoef_dt1 = -length * beta_start * sin_phi * om_g / (2 * E_start)
      dcoef_dE1 = length * beta_start * voltage_max * cos_phi / (2 * E_start**2) + coef * dbeta1_dE1 / beta_start
    else
      dalpha_dt1 = kmat(6,5) / (E_end * sqrt_8 * cos_phi) - log(E_ratio) * om * sin_phi / (sqrt_8 * cos_phi**2)
      dalpha_dE1 = 1 / (E_end * sqrt_8 * cos_phi) - 1 / (E_start * sqrt_8 * cos_phi)
      dcoef_dt1 = sqrt_8 * pc_start * cos(alpha) * dalpha_dt1 / gradient_max
      dcoef_dE1 = coef / (beta_start * pc_start) + sqrt_8 * pc_start * cos(alpha) * dalpha_dE1 / gradient_max
    endif

    z21 = -gradient_max / (sqrt_8 * pc_end)
    z22 = pc_start / pc_end  

    c_min = cos_a - sqrt_2 * beta_start * sin_a * cos_phi
    c_plu = cos_a + sqrt_2 * beta_end * sin_a * cos_phi
    dc_min = -sin_a - sqrt_2 * beta_start * cos_a * cos_phi 
    dc_plu = -sin_a + sqrt_2 * beta_end * cos_a * cos_phi 

    cos_term = 1 + 2 * beta_start * beta_end * cos_phi**2
    dcos_phi = om * sin_phi

    kmat(1,1) =  c_min
    kmat(1,2) =  coef 
    kmat(2,1) =  z21 * (sqrt_2 * (beta_start - beta_end) * cos_phi * cos_a + cos_term * sin_a)
    kmat(2,2) =  c_plu * z22

    kmat(1,5) = orbit%vec(1) * (dc_min * dalpha_dt1 - sqrt_2 * beta_start * sin_a * dcos_phi) + & 
                 orbit%vec(2) * (dcoef_dt1)

    kmat(1,6) = orbit%vec(1) * (dc_min * dalpha_dE1 - sqrt_2 * dbeta1_dE1 * sin_a * cos_phi) + &
                 orbit%vec(2) * (dcoef_dE1)

    kmat(2,5) = orbit%vec(1) * z21 * (sqrt_2 * (beta_start - beta_end) * (dcos_phi * cos_a - cos_phi * sin_a * dalpha_dt1)) + &
                 orbit%vec(1) * z21 * sqrt_2 * (-dbeta2_dE2 * kmat(6,5)) * cos_phi * cos_a + &
                 orbit%vec(1) * z21 * (4 * beta_start * beta_end *cos_phi * dcos_phi * sin_a + cos_term * cos_a * dalpha_dt1) + &
                 orbit%vec(1) * z21 * (2 * beta_start * dbeta2_dE2 * kmat(6,5) * sin_a) + &
                 orbit%vec(1) * (-kmat(2,1) * kmat(6,5) / (beta_end * pc_end)) + &
                 orbit%vec(2) * z22 * (dc_plu * dalpha_dt1 + sqrt_2 * sin_a * (beta_end * dcos_phi + dbeta2_dE2 * kmat(6,5) * cos_phi)) + &
                 orbit%vec(2) * z22 * (-c_plu * kmat(6,5) / (beta_end * pc_end))

    kmat(2,6) = orbit%vec(1) * z21 * (sqrt_2 * cos_phi * ((dbeta1_dE1 - dbeta2_dE2) * cos_a - (beta_start - beta_end) * sin_a * dalpha_dE1)) + &
                 orbit%vec(1) * z21 * (2 * cos_phi**2 * (dbeta1_dE1 * beta_end + beta_start * dbeta2_dE2) * sin_a + cos_term * cos_a * dalpha_dE1) + &
                 orbit%vec(1) * (-kmat(2,1) / (beta_end * pc_end)) + &
                 orbit%vec(2) * z22 * (dc_plu * dalpha_dE1 + sqrt_2 * dbeta2_dE2 * sin_a * cos_phi) + &
                 orbit%vec(2) * c_plu * (1 / (beta_start * pc_end) - pc_start / (beta_end * pc_end**2))

    kmat(3:4,3:4) = kmat(1:2,1:2)

    kmat(3,5) = orbit%vec(3) * (dc_min * dalpha_dt1 - sqrt_2 * beta_start * sin_a * dcos_phi) + & 
                 orbit%vec(4) * (dcoef_dt1)

    kmat(3,6) = orbit%vec(3) * (dc_min * dalpha_dE1 - sqrt_2 * dbeta1_dE1 * sin_a * cos_phi) + &
                 orbit%vec(4) * (dcoef_dE1)

    kmat(4,5) = orbit%vec(3) * z21 * (sqrt_2 * (beta_start - beta_end) * (dcos_phi * cos_a - cos_phi * sin_a * dalpha_dt1)) + &
                 orbit%vec(3) * z21 * sqrt_2 * (-dbeta2_dE2 * kmat(6,5)) * cos_phi * cos_a + &
                 orbit%vec(3) * z21 * (4 * beta_start * beta_end *cos_phi * dcos_phi * sin_a + cos_term * cos_a * dalpha_dt1) + &
                 orbit%vec(3) * z21 * (2 * beta_start * dbeta2_dE2 * kmat(6,5) * sin_a) + &
                 orbit%vec(3) * (-kmat(2,1) * kmat(6,5) / (beta_end * pc_end)) + &
                 orbit%vec(4) * z22 * (dc_plu * dalpha_dt1 + sqrt_2 * sin_a * (beta_end * dcos_phi + dbeta2_dE2 * kmat(6,5) * cos_phi)) + &
                 orbit%vec(4) * z22 * (-c_plu * kmat(6,5) / (beta_end * pc_end))

    kmat(4,6) = orbit%vec(3) * z21 * (sqrt_2 * cos_phi * ((dbeta1_dE1 - dbeta2_dE2) * cos_a - (beta_start - beta_end) * sin_a * dalpha_dE1)) + &
                 orbit%vec(3) * z21 * (2 * cos_phi**2 * (dbeta1_dE1 * beta_end + beta_start * dbeta2_dE2) * sin_a + cos_term * cos_a * dalpha_dE1) + &
                 orbit%vec(3) * (-kmat(2,1) / (beta_end * pc_end)) + &
                 orbit%vec(4) * z22 * (dc_plu * dalpha_dE1 + sqrt_2 * dbeta2_dE2 * sin_a * cos_phi) + &
                 orbit%vec(4) * c_plu * (1 / (beta_start * pc_end) - pc_start / (beta_end * pc_end**2))


    ! Correction to z for finite x', y'
    ! Note: Corrections to kmat(5,5) and kmat(5,6) are ignored since these are small (quadratic
    ! in the transvers coords).

    c_plu = sqrt_2 * cos_phi * cos_a + sin_a

    drp1_dr0  = -gradient_net / (2 * E_start)
    drp1_drp0 = 1

    xp1 = drp1_dr0 * orbit%vec(1) + drp1_drp0 * orbit%vec(2)
    yp1 = drp1_dr0 * orbit%vec(3) + drp1_drp0 * orbit%vec(4)

    drp2_dr0  = (c_plu * z21)
    drp2_drp0 = (cos_a * z22)

    xp2 = drp2_dr0 * orbit%vec(1) + drp2_drp0 * orbit%vec(2)
    yp2 = drp2_dr0 * orbit%vec(3) + drp2_drp0 * orbit%vec(4)

    kmat(5,1) = -(orbit%vec(1) * (drp1_dr0**2 + drp1_dr0*drp2_dr0 + drp2_dr0**2) + &
                   orbit%vec(2) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * dp_dg / 3

    kmat(5,2) = -(orbit%vec(2) * (drp1_drp0**2 + drp1_drp0*drp2_drp0 + drp2_drp0**2) + &
                   orbit%vec(1) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * dp_dg / 3

    kmat(5,3) = -(orbit%vec(3) * (drp1_dr0**2 + drp1_dr0*drp2_dr0 + drp2_dr0**2) + &
                   orbit%vec(4) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * dp_dg / 3

    kmat(5,4) = -(orbit%vec(4) * (drp1_drp0**2 + drp1_drp0*drp2_drp0 + drp2_drp0**2) + &
                   orbit%vec(3) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * dp_dg / 3

    !

    mat6 = matmul(kmat, mat6)
  endif

  k1 = -gradient_net / (2 * E_start)
  orbit%vec(2) = orbit%vec(2) + k1 * orbit%vec(1)    ! Entrance kick
  orbit%vec(4) = orbit%vec(4) + k1 * orbit%vec(3)    ! Entrance kick

  xp1 = orbit%vec(2)
  yp1 = orbit%vec(4)

  orbit%vec(1:2) = matmul(r_mat, orbit%vec(1:2))   ! R&S Eq 9.
  orbit%vec(3:4) = matmul(r_mat, orbit%vec(3:4))

  xp2 = orbit%vec(2)
  yp2 = orbit%vec(4)

  ! Correction of z for finite transverse velocity assumes a uniform change in slope.

  dz_factor = (xp1**2 + xp2**2 + xp1*xp2 + yp1**2 + yp2**2 + yp1*yp2) * dp_dg / 6
  orbit%vec(5) = orbit%vec(5) - dz_factor
  orbit%t = orbit%t + dz_factor / c_light

  !

  k2 = gradient_net / (2 * E_end) 
  orbit%vec(2) = orbit%vec(2) + k2 * orbit%vec(1)         ! Exit kick
  orbit%vec(4) = orbit%vec(4) + k2 * orbit%vec(3)         ! Exit kick

endif

!

orbit%vec(5) = orbit%vec(5) - (dp_dg - c_light * ele%value(delta_ref_time$))
orbit%vec(6) = (pc_end - pc_end_ref) / pc_end_ref 
orbit%p0c = pc_end_ref

! Convert back from (x', y', c(t_ref-t), E) coords

if (logic_option(.false., make_matrix)) then
  rel_p = pc_end / pc_end_ref
  mat6(2,:) = rel_p * mat6(2,:) + orbit%vec(2) * mat6(6,:) / (pc_end_ref * beta_end)
  mat6(4,:) = rel_p * mat6(4,:) + orbit%vec(4) * mat6(6,:) / (pc_end_ref * beta_end)

  m2(1,:) = [beta_end, orbit%vec(5) * mc2**2 / (pc_end * E_end**2)]
  m2(2,:) = [0.0_rp, 1 / (pc_end_ref * beta_end)]

  mat6(5:6,:) = matmul(m2, mat6(5:6,:))
endif

orbit%vec(2) = orbit%vec(2) * (1 + orbit%vec(6))  ! Convert back to px
orbit%vec(4) = orbit%vec(4) * (1 + orbit%vec(6))  ! Convert back to py
orbit%vec(5) = orbit%vec(5) * beta_end

orbit%beta = beta_end

! Coupler kick

call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit, mat6, make_matrix)

if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, orbit, magnetic$, 1.0_rp/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, orbit, electric$, ele%value(l$)/2, mat6, make_matrix)

call offset_particle (ele, param, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

end subroutine
