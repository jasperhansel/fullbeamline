!+
! Subroutine track1_bmad (start_orb, ele, param, end_orb, err_flag, mat6, make_matrix)
!
! Particle tracking through a single element BMAD_standard style.
!
! Input:
!   start_orb   -- Coord_struct: Starting position
!   ele         -- Ele_struct: Element
!   param       -- lat_param_struct:
!     %particle     -- Particle type
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before the element.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   end_orb     -- Coord_struct: End position
!   err_flag    -- Logical, optional: Set true if there is an error. False otherwise.
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix propagated through the element.
!-

subroutine track1_bmad (start_orb, ele, param, end_orb, err_flag, mat6, make_matrix)

use track1_mod, dummy2 => track1_bmad
use geometry_mod, dummy4 => track1_bmad

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp)  knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

integer, parameter :: do_nothing$ = 9999
integer key, ix_pole_max

logical, optional :: err_flag
logical, optional :: make_matrix
logical err

character(*), parameter :: r_name = 'track1_bmad'

! initially set end_orb = start_orb

if (present(err_flag)) err_flag = .false.

end_orb = start_orb     ! transfer start to end

!-----------------------------------------------
! If element is off... 

key = ele%key

if (.not. ele%is_on) then
  select case (key)
  case (ab_multipole$, multipole$, taylor$, match$, fiducial$, floor_shift$)
    key = do_nothing$
  case (lcavity$, sbend$, patch$)
    ! Note: LCavities will still do wakefields.
  case default
    key = drift$
  end select
endif

! Select.

select case (key)

!-----------------------------------------------
! beambeam
                        
case (beambeam$)
  call track_a_beambeam(end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! Thick multipoles

case (rcollimator$, ecollimator$, monitor$, instrument$, pipe$, ac_kicker$, kicker$, hkicker$, vkicker$) 
  call track_a_thick_multipole (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! drift
 
case (drift$) 
  call track_a_drift (end_orb, ele%value(l$), mat6, make_matrix)

!-----------------------------------------------
! elseparator

case (elseparator$)
  call track_a_thick_multipole (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! LCavity: Linac rf cavity.

case (lcavity$)
  call track_a_lcavity (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! marker, etc.
! Note: floor_shift elements can have finite length in the case where it is a slice_slave of a taylor 
! element (the first slice is a taylor element and all other slices are floor_shifts).

case (marker$, fork$, photon_fork$, floor_shift$, fiducial$, detector$)

  end_orb%t = end_orb%t + ele%value(delta_ref_time$)

!-----------------------------------------------
! mask

case (mask$)
  call track_a_mask (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! match

case (match$)
  call track_a_match (end_orb, ele, param, err_flag, mat6, make_matrix)

!-----------------------------------------------
! multipole, ab_multipole

case (multipole$, ab_multipole$) 

  call offset_particle (ele, param, set$, end_orb, set_tilt = .false.)

  call multipole_ele_to_kt(ele, .true., ix_pole_max, knl, tilt)

  if (ix_pole_max > -1) then
    call multipole_kicks (knl, tilt, param%particle, ele, end_orb, ref_orb_offset = (ele%key == multipole$))

    if (logic_option(.false., make_matrix)) then
      call multipole_kick_mat (knl, tilt, param%particle, ele, end_orb, 1.0_rp, ele%mat6)

      ! if knl(0) is non-zero then the reference orbit itself is bent
      ! and we need to account for this.

      if (knl(0) /= 0 .and. ele%key == multipole$) then
        ele%mat6(2,6) = knl(0) * cos(tilt(0))
        ele%mat6(4,6) = knl(0) * sin(tilt(0))
        ele%mat6(5,1) = -ele%mat6(2,6)
        ele%mat6(5,3) = -ele%mat6(4,6)
      endif
    endif
  endif

  call offset_particle (ele, param, unset$, end_orb, set_tilt = .false.)

!-----------------------------------------------
! octupole
! The octupole is modeled using kick-drift.

case (octupole$)
  call track_a_thick_multipole (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! patch

case (patch$)
  call track_a_patch (ele, end_orb, mat6 = mat6, make_matrix = make_matrix)

!-----------------------------------------------
! quadrupole

case (quadrupole$)
  call track_a_quadrupole (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! rfcavity

case (rfcavity$)
  call track_a_rfcavity (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! sad_multipole

case (sad_mult$)
  call track_a_sad_mult (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! sbend

case (sbend$)
  call track_a_bend (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! sextupole
! The sextupole is modeled using kick-drift.

case (sextupole$)
  call track_a_thick_multipole (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! Solenoid

case (sol_quad$, solenoid$)
  call track_a_sol_quad (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! Taylor

case (taylor$)
  call track1_taylor (end_orb, ele, param, end_orb, mat6 = mat6, make_matrix = make_matrix)

!-----------------------------------------------
! wiggler:

case (wiggler$, undulator$)
  call track_a_wiggler (end_orb, ele, param, mat6, make_matrix)

!-----------------------------------------------
! do nothing case

case (do_nothing$)

!-----------------------------------------------
! unknown

case default

  if (present(err_flag)) err_flag = .true.
  call out_io (s_fatal$, r_name, &
          'BMAD_STANDARD TRACKING_METHOD NOT IMPLMENTED FOR: ' // key_name(ele%key), &
          'FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return

end select

!-----------------------------------------------------------------------------------
! Set s-position

if (end_orb%direction == 1) then
  end_orb%s = ele%s
else
  end_orb%s = ele%s_start
endif

end subroutine track1_bmad
