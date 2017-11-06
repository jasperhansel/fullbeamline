!+
! Subroutine track_a_rfcavity (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an rfcavity element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: rfcavity element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_rfcavity (orbit, ele, param, mat6, make_matrix)

use track1_mod, except_dummy => track_a_rfcavity

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) beta_ref, dt_ref, voltage, phase0, phase, dE, z, charge_dir, length
real(rp) mc2, p0c, m2(2,2), t0, factor, pc, E, ff, new_pc, new_beta
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer i, n_slice, orientation
integer ix_pole_max, ix_elec_max

logical, optional :: make_matrix
logical drifting

!

call multipole_ele_to_ab (ele, .false., ix_pole_max, an,      bn,      magnetic$, include_kicks = .true.)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

call offset_particle (ele, param, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, orbit, magnetic$, 1.0_rp/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, orbit, electric$, ele%value(l$)/2, mat6, make_matrix)

!

length = ele%value(l$)
beta_ref = ele%value(p0c$) / ele%value(e_tot$)
n_slice = max(1, nint(length / ele%value(ds_step$))) 
dt_ref = length / (c_light * beta_ref)
charge_dir = rel_tracking_charge_to_mass(orbit, param) * ele%orientation
mc2 = mass_of(orbit%species)
p0c = orbit%p0c

! The cavity field is modeled as a standing wave antisymmetric wrt the center.
! Thus if the cavity is flipped (orientation = -1), the wave of interest, which is 
! always the accelerating wave, is the "backward" wave. And the phase of the backward 
! wave is different from the phase of the forward wave by a constant dt_ref * freq.

voltage = e_accel_field(ele, voltage$) * charge_dir

phase0 = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_autoscale$) - &
        (particle_rf_time (orbit, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))
if (ele%orientation == -1) phase0 = phase0 + twopi * ele%value(rf_frequency$) * dt_ref
phase = phase0

call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit, mat6, make_matrix)

! Track through slices.
! The phase of the accelerating wave traveling in the same direction as the particle is
! assumed to be traveling with a phase velocity the same speed as the reference velocity.

t0 = orbit%t

do i = 0, n_slice

  if (logic_option(.false., make_matrix)) then
    factor = voltage / n_slice
    if (i == 0 .or. i == n_slice) factor = factor / 2

    dE = factor * sin(phase)
    pc = (1 + orbit%vec(6)) * p0c 
    E = pc / orbit%beta
    call convert_total_energy_to (E + dE, orbit%species, pc = new_pc, beta = new_beta)
    ff = twopi * factor * ele%value(rf_frequency$) * cos(phase) / (p0c * new_beta * c_light)

    m2(2,1) = ff / orbit%beta
    m2(2,2) = orbit%beta / new_beta - ff * orbit%vec(5) *mc2**2 * p0c / (E * pc**2) 
    m2(1,1) = new_beta / orbit%beta + orbit%vec(5) * (mc2**2 * p0c * m2(2,1) / (E+dE)**3) / orbit%beta
    m2(1,2) = orbit%vec(5) * mc2**2 * p0c * (m2(2,2) / ((E+dE)**3 * orbit%beta) - new_beta / (pc**2 * E))

    mat6(5:6, :) = matmul(m2, mat6(5:6, :))
  endif

  dE = voltage * sin(phase) / n_slice
  if (i == 0 .or. i == n_slice) dE = dE / 2

  call apply_energy_kick (dE, orbit)
  
  if (orbit%vec(6) == -1) then
    orbit%state = lost_pz_aperture$
    return
  endif

  if (i /= n_slice) then
    z = orbit%vec(5)
    call track_a_drift (orbit, length/n_slice, mat6, make_matrix)
    phase = phase + twopi * ele%value(rf_frequency$) * (orbit%vec(5) - z) / (c_light * orbit%beta)
  endif

enddo

! coupler kick, multipoles, back to lab coords.

call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit, mat6, make_matrix)

if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, orbit, magnetic$, 1.0_rp/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, orbit, electric$, ele%value(l$)/2, mat6, make_matrix)

call offset_particle (ele, param, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

end subroutine
