module custom_bmad_interface

use bmad_struct

interface 

subroutine apply_element_edge_kick_hook (orb, fringe_info, track_ele, param, finished, mat6, make_matrix, rf_time)
  import
  implicit none
  type (coord_struct) orb
  type (fringe_edge_info_struct) fringe_info
  type (ele_struct) track_ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6), rf_time
  logical, optional :: make_matrix
  logical finished
end subroutine

subroutine check_aperture_limit_custom (orb, ele, particle_at, param, err_flag)
  import
  implicit none
  type (coord_struct) :: orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  integer particle_at
  logical err_flag
end subroutine

function distance_to_aperture_custom (orbit, particle_at, ele, no_aperture_here) result (dist)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp) dist
  integer particle_at
  logical no_aperture_here
end function

subroutine ele_geometry_hook (floor0, ele, floor, finished, len_scale)
  import
  implicit none
  type (ele_struct) ele
  type (floor_position_struct) floor0, floor
  real(rp) len_scale
  logical finished
end subroutine ele_geometry_hook

subroutine wall_hit_handler_custom (orb, ele, s)
  import
  implicit none
  type (coord_struct) :: orb
  type (ele_struct) :: ele
  real(rp) s
end subroutine

subroutine em_field_custom (ele, param, s_rel, orb, local_ref_frame, field, calc_dfield, err_flag)
  import
  implicit none
  type (ele_struct) :: ele
  type (lat_param_struct) param
  type (coord_struct), intent(in) :: orb
  real(rp), intent(in) :: s_rel
  logical local_ref_frame
  type (em_field_struct) :: field
  logical, optional :: err_flag
  logical, optional :: calc_dfield
end subroutine

subroutine ele_to_fibre_hook (ele, ptc_fibre, param)
  import
  implicit none
  type (ele_struct) ele
  type (fibre) ptc_fibre
  type (lat_param_struct) param
end subroutine

subroutine radiation_integrals_custom (lat, ir, orb, err_flag)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct) orb(0:)
  integer ir
  logical err_flag
end subroutine

subroutine init_custom (ele, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  logical err_flag
end subroutine

subroutine make_mat6_custom (ele, param, start_orb, end_orb, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical err_flag, finished
end subroutine

subroutine time_runge_kutta_periodic_kick_hook (orbit, z_phase_space, ele, param, stop_time, init_needed)
  import
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp) z_phase_space, stop_time
  integer :: init_needed
end subroutine

subroutine track1_beam_hook (beam_start, lat, ele, beam_end, err, centroid, direction, finished)
  import
  implicit none
  type (beam_struct) beam_start
  type (beam_struct) :: beam_end
  type (lat_struct) :: lat
  type (ele_struct) ele
  type (coord_struct), optional :: centroid(0:)
  integer, optional :: direction
  logical err, finished
end subroutine

subroutine track1_custom (start_orb, ele, param, end_orb, err_flag, finished, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
  logical err_flag, finished, radiation_included
end subroutine

subroutine track1_postprocess (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_preprocess (start_orb, ele, param, err_flag, finished, radiation_included, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
  logical err_flag, finished, radiation_included
end subroutine

subroutine track1_spin_custom (start_orb, ele, param, end_orb, err_flag, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
  logical err_flag
end subroutine

subroutine track1_wake_hook (bunch, ele, finished)
  import
  implicit none
  type (bunch_struct) bunch
  type (ele_struct) ele
  logical finished
end subroutine

end interface

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
