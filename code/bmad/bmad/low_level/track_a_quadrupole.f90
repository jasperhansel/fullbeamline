!+
! Subroutine track_a_quadrupole (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through a quadrupole element. 
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Quadrupole element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_quadrupole (orbit, ele, param, mat6, make_matrix)

use track1_mod, except_dummy => track_a_quadrupole

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (fringe_edge_info_struct) fringe_info

real(rp), optional :: mat6(6,6)
real(rp) kmat(6,6), mat2(2,2), rel_p, dz_x(3), dz_y(3), ddz_x(3), ddz_y(3)
real(rp) k1, rel_tracking_charge, charge_dir, r_step, step_len, s_off, mass, e_tot
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer i, n_step, orientation, ix_pole_max, ix_elec_max

logical, optional :: make_matrix
logical drifting

!

start_orb = orbit
orientation = ele%orientation * start_orb%direction
rel_tracking_charge = rel_tracking_charge_to_mass(start_orb, param)
charge_dir = rel_tracking_charge * orientation

call multipole_ele_to_ab (ele, .false., ix_pole_max, an,      bn,      magnetic$, include_kicks = .true.)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

n_step = 1
if (ix_pole_max > -1 .or. ix_elec_max > -1) n_step = max(nint(ele%value(l$) / ele%value(ds_step$)), 1)

r_step = 1.0_rp / n_step
step_len = ele%value(l$) * r_step

! Entrance edge

call offset_particle (ele, param, set$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

nullify(fringe_info%hard_ele)
fringe_info%particle_at = first_track_edge$
call apply_element_edge_kick(orbit, fringe_info, ele, param, .false., mat6, make_matrix)
if (orbit%state /= alive$) return

! Multipole kicks. Notice that the magnetic multipoles have already been normalized by the length.

if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, orbit, magnetic$, r_step/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, orbit, electric$, step_len/2, mat6, make_matrix)

! Body

do i = 1, n_step

  if (logic_option(.false., make_matrix)) call mat_make_unit (kmat)

  rel_p = 1 + orbit%vec(6)  ! Can change when there are electric fields
  k1 = charge_dir * ele%value(k1$) / rel_p

  call quad_mat2_calc (-k1, step_len, rel_p, kmat(1:2,1:2), dz_x, ddz_x)
  call quad_mat2_calc ( k1, step_len, rel_p, kmat(3:4,3:4), dz_y, ddz_y)

  ! The mat6(i,6) terms are constructed so that mat6 is sympelctic

  if (logic_option(.false., make_matrix)) then
    if (any(orbit%vec(1:4) /= 0)) then
      kmat(5,1) = 2 * orbit%vec(1) * dz_x(1) +     orbit%vec(2) * dz_x(2)
      kmat(5,2) =     orbit%vec(1) * dz_x(2) + 2 * orbit%vec(2) * dz_x(3)
      kmat(5,3) = 2 * orbit%vec(3) * dz_y(1) +     orbit%vec(4) * dz_y(2)
      kmat(5,4) =     orbit%vec(3) * dz_y(2) + 2 * orbit%vec(4) * dz_y(3)
      kmat(5,6) = orbit%vec(1)**2 * ddz_x(1) + orbit%vec(1)*orbit%vec(2) * ddz_x(2) + orbit%vec(2)**2 * ddz_x(3) + &
                   orbit%vec(3)**2 * ddz_y(1) + orbit%vec(3)*orbit%vec(4) * ddz_y(2) + orbit%vec(4)**2 * ddz_y(3)  
    endif

    if (any(kmat(5,1:4) /= 0)) then
      kmat(1,6) = kmat(5,2) * kmat(1,1) - kmat(5,1) * kmat(1,2)
      kmat(2,6) = kmat(5,2) * kmat(2,1) - kmat(5,1) * kmat(2,2)
      kmat(3,6) = kmat(5,4) * kmat(3,3) - kmat(5,3) * kmat(3,4)

      kmat(4,6) = kmat(5,4) * kmat(4,3) - kmat(5,3) * kmat(4,4)
    endif

    mat6 = matmul(kmat, mat6)
  endif

  !

  orbit%vec(5) = orbit%vec(5) + &
                  dz_x(1) * orbit%vec(1)**2 + dz_x(2) * orbit%vec(1) * orbit%vec(2) + dz_x(3) * orbit%vec(2)**2 + &
                  dz_y(1) * orbit%vec(3)**2 + dz_y(2) * orbit%vec(3) * orbit%vec(4) + dz_y(3) * orbit%vec(4)**2 

  orbit%vec(1:2) = matmul(kmat(1:2,1:2), orbit%vec(1:2))
  orbit%vec(3:4) = matmul(kmat(3:4,3:4), orbit%vec(3:4))

  orbit%vec(5) = orbit%vec(5) + low_energy_z_correction (orbit, ele, step_len, mat6, make_matrix)

  if (i == n_step) then
    if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, orbit, magnetic$, r_step/2,   mat6, make_matrix)
    if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, orbit, electric$, step_len/2, mat6, make_matrix)
  else
    if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, orbit, magnetic$, r_step,   mat6, make_matrix)
    if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, orbit, electric$, step_len, mat6, make_matrix)
  endif

enddo

! Exit edge

fringe_info%particle_at = second_track_edge$
call apply_element_edge_kick(orbit, fringe_info, ele, param, .false., mat6, make_matrix)
if (orbit%state /= alive$) return

call offset_particle (ele, param, unset$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

orbit%t = start_orb%t + ele%value(delta_ref_time$) + (start_orb%vec(5) - orbit%vec(5)) / (orbit%beta * c_light)

!

end subroutine
