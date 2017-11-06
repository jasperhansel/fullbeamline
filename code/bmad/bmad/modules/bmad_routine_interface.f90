!+
! This file defines the interfaces for the BMAD subroutines
!-

module bmad_routine_interface

use bmad_struct

interface

function ac_kicker_amp(ele, time) result (ac_amp)
  import
  implicit none
  type (ele_struct) ele
  type (ac_kicker_struct), pointer :: ac
  real(rp) time, ac_amp
end function

subroutine add_lattice_control_structs (ele, n_add_slave, n_add_lord, n_add_slave_field, n_add_lord_field, add_at_end)
  import
  implicit none
  type (ele_struct) ele
  integer, optional :: n_add_slave, n_add_lord, n_add_slave_field, n_add_lord_field
  logical, optional :: add_at_end
end subroutine

subroutine aml_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) lat_file
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  character(*), optional :: use_line
end subroutine

subroutine apply_element_edge_kick (orb, fringe_info, track_ele, param, track_spin, mat6, make_matrix, rf_time, apply_sol_fringe)
  import
  implicit none
  type (coord_struct) orb
  type (fringe_edge_info_struct) fringe_info
  type (ele_struct) hard_ele, track_ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6), rf_time
  logical, optional :: make_matrix, apply_sol_fringe
  logical track_spin
end subroutine

subroutine apply_energy_kick (dE, orbit, ddE_dr, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  real(rp) dE
  real(rp), optional :: ddE_dr(2), mat6(6,6)
  logical, optional :: make_matrix
end subroutine

function at_this_ele_end (now_at, where_at) result (is_at_this_end)
  import
  implicit none
  integer now_at, where_at
  logical is_at_this_end
end function

subroutine attribute_bookkeeper (ele, param, force_bookkeeping)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  logical, optional :: force_bookkeeping
end subroutine

subroutine bbi_kick (x, y, r, kx, ky)
  import
  implicit none
  real(rp) x, y, r, kx, ky
end subroutine

subroutine bend_exact_multipole_field (ele, param, orbit, local_ref_frame, field, calc_dfield, potential)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  type (em_field_struct) field
  logical local_ref_frame
  logical, optional :: calc_dfield
  type (em_potential_struct), optional :: potential
end subroutine

subroutine bmad_and_xsif_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) lat_file
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  character(*), optional :: use_line
end subroutine

subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) lat_file
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  character(*), optional :: use_line
end subroutine

subroutine bmad_parser2 (in_file, lat, orbit, make_mats6, err_flag)
  import
  implicit none
  character(*) in_file
  type (lat_struct), target :: lat
  type (coord_struct), optional :: orbit(0:)
  logical, optional :: make_mats6, err_flag
end subroutine

function branch_name(branch) result (name)
  import
  implicit none
  type (branch_struct), target :: branch
  character(40) name
end function

function c_multi (n, m, no_n_fact, c_full) result (c_out)
  import
  implicit none
  integer, intent(in) :: n, m
  real(rp) c_out
  real(rp), optional :: c_full(0:n_pole_maxx, 0:n_pole_maxx)
  logical, optional :: no_n_fact
end function

subroutine c_to_cbar (ele, cbar_mat)
  import
  implicit none
  type (ele_struct) ele
  real(rp) cbar_mat(2,2)
end subroutine

subroutine calc_z_tune (lat, ix_branch)
  import
  implicit none
  type (lat_struct) lat
  integer, optional :: ix_branch
end subroutine

subroutine cbar_to_c (cbar_mat, a, b, c_mat)
  import
  implicit none
  real(rp) cbar_mat(2,2), c_mat(2,2)
  type (twiss_struct) a, b
end subroutine

recursive subroutine check_aperture_limit (orb, ele, particle_at, param, old_orb, check_momentum)
  import
  implicit none
  type (coord_struct) :: orb
  type (ele_struct) :: ele
  type (lat_param_struct), intent(inout) :: param
  integer particle_at
  type (coord_struct), optional :: old_orb
  logical, optional :: check_momentum
end subroutine

subroutine compute_even_steps (ds_in, length, ds_default, ds_out, n_step)
  import
  implicit none
  real(rp) ds_in, length, ds_default, ds_out
  integer n_step
end subroutine

subroutine convert_total_energy_to (E_tot, particle, gamma, kinetic, beta, pc, brho, beta1, err_flag)
  import
  implicit none
  real(rp), intent(in) :: E_tot
  real(rp), intent(out), optional :: pc, kinetic, beta, brho, gamma, beta1
  integer, intent(in) :: particle
  logical, optional :: err_flag
end subroutine

subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho, beta1, err_flag)
  import
  implicit none
  real(rp), intent(in) :: pc
  real(rp), intent(out), optional :: E_tot, kinetic, beta, brho, gamma, beta1
  integer, intent(in) :: particle
  logical, optional :: err_flag
end subroutine

subroutine crystal_attribute_bookkeeper (ele)
  import
  type (ele_struct) ele
end subroutine

subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, err_flag, &
                       pz, low_E_lat, high_E_lat, low_E_orb, high_E_orb, ix_branch)
  import
  implicit none
  type (lat_struct) lat
  type (lat_struct), optional, target :: low_E_lat, high_E_lat
  type (coord_struct), allocatable, optional, target :: low_E_orb(:), high_E_orb(:)
  real(rp) delta_e
  real(rp) chrom_x
  real(rp) chrom_y
  real(rp), optional :: pz
  logical, optional, intent(out) :: err_flag
  integer, optional :: ix_branch
end subroutine

subroutine chrom_tune (lat, delta_e, chrom_x, chrom_y, err_tol, err_flag)
  import
  implicit none
  type (lat_struct) lat
  real(rp) delta_e
  real(rp) chrom_x
  real(rp) chrom_y
  real(rp) err_tol
  logical err_flag
end subroutine

subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag, print_err)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), allocatable, target :: closed_orb(:)
  integer, optional :: direction, ix_branch, i_dim
  logical, optional, intent(out) :: err_flag
  logical, optional, intent(in) :: print_err
end subroutine

subroutine closed_orbit_from_tracking (lat, closed_orb, i_dim, &
                                     eps_rel, eps_abs, init_guess, err_flag)
  import
  type (lat_struct) lat
  type (coord_struct), allocatable :: closed_orb(:)
  type (coord_struct), optional :: init_guess
  real(rp), intent(in), optional :: eps_rel(:), eps_abs(:)
  integer i_dim
  logical, optional :: err_flag
end subroutine

subroutine combine_consecutive_elements (lat)
  import
  type (lat_struct), target :: lat
end subroutine

subroutine convert_bend_exact_multipole (bend_in, bend_out, out_type)
  import
  implicit none
  type (ele_struct) bend_in, bend_out
  integer out_type
end subroutine

subroutine create_field_overlap (lat, lord_name, slave_name, err_flag)
  import
  implicit none
  type (lat_struct) lat
  character(*) lord_name, slave_name
  logical err_flag
end subroutine

subroutine create_girder (lat, ix_ele, con, init_ele, err_flag)
  import
  implicit none
  type (lat_struct) lat
  type (ele_struct) :: init_ele
  type (control_struct) con(:)
  integer, intent(in) :: ix_ele
  logical err_flag
end subroutine

subroutine create_group (lord, con, err, err_print_flag)
  import
  implicit none
  type (ele_struct) lord
  type (control_struct) con(:)
  logical err
  logical, optional :: err_print_flag
end subroutine

subroutine create_overlay (lord, contl, err, err_print_flag)
  import
  implicit none
  type (ele_struct) lord
  type (control_struct) contl(:)
  logical err
  logical, optional :: err_print_flag
end subroutine

subroutine create_uniform_element_slice (ele, param, i_slice, n_slice_tot, sliced_ele, s_start, s_end)
  import
  implicit none
  type (ele_struct) ele, sliced_ele
  type (lat_param_struct) param
  integer i_slice, n_slice_tot
  real(rp), optional :: s_start, s_end
end subroutine

subroutine create_unique_ele_names (lat, key, suffix)
  import
  type (lat_struct), target :: lat
  integer key
  character(*) suffix
end subroutine

subroutine convert_coords (in_type_str, coord_in, ele, out_type_str, coord_out, err_flag)
  import
  implicit none
  character(*) in_type_str
  character(*) out_type_str
  type (coord_struct) coord_in
  type (coord_struct) coord_out
  type (ele_struct) ele
  logical, optional :: err_flag
end subroutine

function default_tracking_species (param) result (species)
  import
  implicit none
  type (lat_param_struct) param
  integer species
end function

function diffraction_plate_or_mask_hit_spot (ele, orbit) result (ix_section)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) orbit
  integer :: ix_section
end function

recursive function distance_to_aperture (orbit, particle_at, ele, no_aperture_here) result (dist)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp) dist
  integer particle_at
  logical no_aperture_here
end function

subroutine do_mode_flip (ele, err_flag)
  import
  implicit none
  type (ele_struct) ele
  logical, optional :: err_flag
end subroutine

function e_accel_field (ele, voltage_or_gradient) result (field)
  import
  implicit none
  type (ele_struct) ele
  real(rp) field
  integer voltage_or_gradient 
end function

recursive subroutine ele_compute_ref_energy_and_time (ele0, ele, param, err_flag)
  import
  type (ele_struct) ele0, ele
  type (lat_param_struct) param
  real(rp) e_tot_start, p0c_start, ref_time_start
  logical err_flag
end subroutine

function ele_has_constant_ds_dt_ref (ele) result (is_const)
  import
  implicit none
  type (ele_struct) ele
  logical is_const
end function

subroutine fibre_to_ele (ptc_fibre, branch, ix_ele, err_flag, from_mad)
  import
  implicit none
  type (fibre) ptc_fibre
  type (branch_struct) branch
  integer ix_ele
  logical err_flag
  logical, optional :: from_mad
end subroutine

subroutine find_element_ends (ele, ele1, ele2, ix_multipass)
  import
  implicit none
  type (ele_struct) ele
  type (ele_struct), pointer :: ele1, ele2
  integer, optional :: ix_multipass
end subroutine

subroutine find_matching_fieldmap (file_name, ele, t_type, match_ele, ix_field)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: match_ele
  integer t_type, ix_field
  character(*) file_name
end subroutine

function hard_edge_model_length (ele) result (l_hard)
  import
  implicit none
  type (ele_struct) ele
  real(rp) l_hard
end function

subroutine init_a_photon_from_a_photon_init_ele (ele, param, orbit)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orbit
end subroutine

subroutine init_wake (wake, n_sr_long, n_sr_trans, n_lr_mode, n_lr_spline, always_allocate)
  import
  implicit none
  type (wake_struct), pointer :: wake
  integer n_sr_long, n_sr_trans, n_lr_mode, n_lr_spline
  logical, optional :: always_allocate
end subroutine

subroutine insert_element (lat, insert_ele, insert_index, ix_branch, orbit)
  import
  implicit none
  type (lat_struct) lat
  type (ele_struct) insert_ele
  integer insert_index
  integer, optional :: ix_branch
  type (coord_struct), optional, allocatable :: orbit(:)
end subroutine

subroutine ion_kick (orbit, r_beam, n_beam_part, a_twiss, b_twiss, sig_ee, kick)
  import
  implicit none
  type (coord_struct) orbit
  type (twiss_struct) a_twiss, b_twiss
  real(rp) r_beam(2), n_beam_part, sig_ee, kick(3)
end subroutine

function key_name_to_key_index (key_str, abbrev_allowed) result (key_index)
  import
  implicit none
  character(*) key_str
  logical, optional :: abbrev_allowed
  integer key_index
end function

subroutine kill_ptc_layouts (lat)
  import
  implicit none
  type (lat_struct) lat
end subroutine

subroutine lat_compute_ref_energy_and_time (lat, err_flag)
  import
  type (lat_struct) lat
  logical err_flag
end subroutine

recursive subroutine lat_make_mat6 (lat, ix_ele, ref_orb, ix_branch, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), optional :: ref_orb(0:)
  integer, optional :: ix_ele, ix_branch
  logical, optional :: err_flag
end subroutine

subroutine lat_sanity_check (lat, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  logical, intent(out) :: err_flag
end subroutine

function low_energy_z_correction (orbit, ele, ds, mat6, make_matrix) result (dz)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp) ds, dz
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end function

subroutine make_g_mats (ele, g_mat, g_inv_mat)
  import
  implicit none
  type (ele_struct) ele
  real(rp) g_mat(4,4)
  real(rp) g_inv_mat(4,4)
end subroutine

subroutine make_hybrid_lat (r_in, r_out, use_taylor, orb0)
  import
  implicit none
  type (lat_struct), target :: r_in
  type (lat_struct), target :: r_out
  logical, optional :: use_taylor
  type (coord_array_struct), optional :: orb0(0:)
end subroutine

recursive subroutine make_mat6 (ele, param, start_orb, end_orb, err_flag)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct), optional :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: err_flag
end subroutine

subroutine make_mat6_taylor (ele, param, start_orb, end_orb, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: err_flag
end subroutine

subroutine make_mat6_bmad (ele, param, start_orb, end_orb, err)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: err
end subroutine

subroutine make_mat6_bmad_photon (ele, param, start_orb, end_orb, err)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: err
end subroutine

subroutine make_mat6_runge_kutta (ele, param, start_orb, end_orb)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
end subroutine

subroutine make_mat6_symp_lie_ptc (ele, param, start_orb, end_orb)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
end subroutine

subroutine make_mat6_tracking (ele, param, start_orb, end_orb)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
end subroutine

subroutine make_v_mats (ele, v_mat, v_inv_mat)
  import
  implicit none
  type (ele_struct) ele
  real(rp), optional :: v_mat(4,4)
  real(rp), optional :: v_inv_mat(4,4)
end subroutine

subroutine mat6_add_offsets (ele, param)
  import
  type (ele_struct) ele
  type (lat_param_struct) param
end subroutine

subroutine match_ele_to_mat6 (ele, start_orb, mat6, vec0, err_flag, twiss_ele, include_delta_time)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) start_orb
  real(rp) mat6(6,6), vec0(6)
  type (ele_struct), optional, target :: twiss_ele
  logical :: err_flag
  logical, optional :: include_delta_time
end subroutine

function mexp (x, m) result (this_exp)
  import
  implicit none
  real(rp) x, this_exp
  integer m
end function

subroutine multi_turn_tracking_analysis (track, i_dim, track0, ele, stable, growth_rate, chi, err_flag)
  import
  implicit none
  type (coord_struct), intent(in) :: track(:)
  type (coord_struct), intent(out) :: track0
  type (ele_struct) :: ele
  real(rp), intent(out) :: growth_rate, chi
  integer, intent(in) :: i_dim
  logical, intent(out) :: stable, err_flag
end subroutine

subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, map0, track0, chi)
  import
  implicit none
  type (coord_struct), intent(in), target :: track(:)
  type (coord_struct), intent(out) :: track0
  real(rp), intent(out) :: mat1(:,:), map0(:)
  real(rp), intent(out) :: chi
  integer, intent(in) :: i_dim
end subroutine

subroutine multipole_init (ele, who, zero)
  import
  implicit none
  type (ele_struct) ele
  integer who
  logical, optional :: zero
end subroutine

subroutine name_to_list (lat, ele_names)
  import
  implicit none
  type (lat_struct) lat
  character(*) ele_names(:)
end subroutine

subroutine new_control (lat, ix_ele)
  import
  implicit none
  type (lat_struct) lat
  integer ix_ele
end subroutine

function num_lords (slave, lord_type) result (num)
  import
  implicit none
  type (ele_struct) slave
integer lord_type, num
end function

subroutine offset_particle (ele, param, set, coord, set_tilt, set_hvkicks, drift_to_edge, s_pos, s_out, set_spin, mat6, make_matrix)
  import
  implicit none
  type (ele_struct) :: ele
  type (lat_param_struct) param
  type (coord_struct), intent(inout) :: coord
  integer particle
  logical, intent(in) :: set
  logical, optional, intent(in) :: set_tilt, set_hvkicks, drift_to_edge, set_spin
  real(rp), optional :: s_pos, mat6(6,6), s_out
  logical, optional :: make_matrix
end subroutine

subroutine offset_photon (ele, coord, set, offset_position_only, rot_mat)
  import
  implicit none
  type (ele_struct) :: ele
  type (coord_struct) :: coord
  logical :: set
  logical, optional :: offset_position_only
  real(rp), optional :: rot_mat(3,3)
end subroutine

subroutine one_turn_mat_at_ele (ele, phi_a, phi_b, mat4)
  import
  type (ele_struct) ele
  real(rp) phi_a
  real(rp) phi_b
  real(rp) mat4(4,4)
end subroutine

subroutine orbit_amplitude_calc (ele, orb, amp_a, amp_b, amp_na, amp_nb)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct) orb
  real(rp), optional :: amp_a, amp_b, amp_na, amp_nb
end subroutine

function orbit_too_large (orbit, param) result (is_too_large)
  import
  implicit none
  type (coord_struct) orbit
  type (lat_param_struct), optional :: param
  logical is_too_large
  real(rp) rel_p
end function

subroutine order_super_lord_slaves (lat, ix_lord)
  import
  implicit none
  type (lat_struct), target :: lat
  integer ix_lord
end subroutine

function particle_rf_time (orbit, ele, apply_hard_edge_offset, s_rel) result (time)
  import
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp), optional :: s_rel
  real(rp) time
  logical apply_hard_edge_offset
end function

subroutine phase_space_fit (x, xp, twiss, tune, emit, x_0, xp_0, chi, tol)
  import
  implicit none
  type (twiss_struct) twiss
  real(rp), optional :: tol
  real(rp) x(:), xp(:)
  real(rp) tune, emit
  real(rp) x_0, xp_0, chi
end subroutine

function physical_ele_end (track_end, track_direction, ele_orientation, return_stream_end) result (physical_end)
  import
  integer track_end, track_direction, ele_orientation, physical_end
  logical, optional :: return_stream_end
end function

subroutine pointer_to_attribute (ele, attrib_name, do_allocation, a_ptr, err_flag, err_print_flag, ix_attrib)
  import
  implicit none
  type (ele_struct), target :: ele
  type (all_pointer_struct) a_ptr
  character(*) attrib_name
  logical err_flag
  logical do_allocation
  logical, optional :: err_print_flag
  integer, optional :: ix_attrib
end subroutine

function pointer_to_lord (slave, ix_lord, control, ix_slave, field_overlap_ptr) result (lord_ptr)
  import
  implicit none
  type (ele_struct), target :: slave
  type (control_struct), pointer, optional :: control
  type (ele_struct), pointer :: lord_ptr
  integer, optional :: ix_slave
  integer ix_lord
  logical, optional :: field_overlap_ptr
end function

function pointer_to_next_ele (this_ele, offset, skip_beginning, follow_fork) result (next_ele)
  import
  implicit none
  type (ele_struct), target :: this_ele
  type (ele_struct), pointer :: next_ele
  integer, optional :: offset
  logical, optional :: skip_beginning, follow_fork
end function

function pointer_to_slave (lord, ix_slave, control, field_overlap_ptr) result (slave_ptr)
  import
  implicit none
  type (ele_struct), target :: lord
  type (control_struct), pointer, optional :: control
  type (ele_struct), pointer :: slave_ptr
  integer ix_slave
  logical, optional :: field_overlap_ptr
end function

function pointer_to_wake_ele (ele, delta_s) result (wake_ele)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: wake_ele
  real(rp), optional :: delta_s
end function

subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation, &
                  ptr_array, err_flag, err_print_flag, eles, ix_attrib)
  import
  implicit none
  type (lat_struct) lat
  type (all_pointer_struct), allocatable :: ptr_array(:)
  character(*) ele_name, attrib_name
  logical err_flag
  logical do_allocation
  logical, optional :: err_print_flag
  type (ele_pointer_struct), optional, allocatable :: eles(:)
  integer, optional :: ix_attrib
end subroutine

subroutine ptc_bookkeeper (lat)
  import
  implicit none
  type (lat_struct) lat
end subroutine

subroutine ptc_read_flat_file (flat_file, err_flag, lat, create_end_marker, from_mad)
  import
  implicit none
  type (lat_struct), optional :: lat
  character(*) flat_file(:)
  logical err_flag
  logical, optional :: create_end_marker, from_mad
end subroutine

subroutine quad_beta_ave (ele, beta_a_ave, beta_b_ave)
  import
  implicit none
  type (ele_struct) ele
  real(rp) beta_a_ave
  real(rp) beta_b_ave
end subroutine

subroutine radiation_integrals (lat, orb, mode, ix_cache, ix_branch, rad_int_by_ele)
  import
  implicit none
  type (lat_struct), target :: lat
  type (rad_int_all_ele_struct), optional :: rad_int_by_ele
  type (coord_struct), target :: orb(0:)
  type (normal_modes_struct) mode
  integer, optional :: ix_cache, ix_branch
end subroutine

subroutine reallocate_expression_stack (stack, n, exact)
  import
  implicit none
  type (expression_atom_struct), allocatable :: stack(:)
  integer n
  logical, optional :: exact
end subroutine

subroutine reference_energy_correction (ele, orbit, particle_at)
  import
  implicit none
  type (ele_struct) :: ele
  type (coord_struct) :: orbit
  integer particle_at
end subroutine

subroutine remove_eles_from_lat (lat, check_sanity)
  import
  implicit none
  type (lat_struct) lat
  logical, optional :: check_sanity
end subroutine

subroutine read_digested_bmad_file (in_file_name, lat, version, err_flag)
  import
  implicit none
  type (lat_struct), target, intent(inout) :: lat
  integer version
  character(*) in_file_name
  logical, optional :: err_flag
end subroutine

subroutine reallocate_beam (beam, n_bunch, n_particle)
  import
  type (beam_struct) beam
  integer n_bunch, n_particle
end subroutine

subroutine reallocate_bunch (bunch, n_particle)
  import
  type (bunch_struct) bunch
  integer n_particle
end subroutine

function relative_mode_flip (ele1, ele2)
  import
  implicit none
  logical relative_mode_flip
  type (ele_struct) ele1
  type (ele_struct) ele2
end function

subroutine reverse_lat (lat_in, lat_rev, track_antiparticle)
  import
  implicit none
  type (lat_struct), target :: lat_in, lat_rev
  logical, optional :: track_antiparticle
end subroutine

subroutine rf_coupler_kick (ele, param, particle_at, phase, orbit, mat6, make_matrix)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  real(rp) phase
  real(rp), optional :: mat6(6,6)
  integer particle_at
  logical, optional :: make_matrix
end subroutine

subroutine s_calc (lat)
  import
  implicit none
  type (lat_struct) lat
end subroutine

subroutine save_a_step (track, ele, param, local_ref_frame, orb, s_rel, save_field, rf_time)
  import
  implicit none
  type (track_struct) track
  type (ele_struct), target :: ele
  type (lat_param_struct), intent(in) :: param
  type (coord_struct) orb
  real(rp) s_rel
  real(rp), optional :: rf_time
  logical local_ref_frame
  logical, optional :: save_field
end subroutine

recursive subroutine set_lords_status_stale (ele, stat_group, control_bookkeeping, flag)
  import
  implicit none
  type (ele_struct) ele
  integer stat_group
  logical, optional :: control_bookkeeping
  integer, optional :: flag
end subroutine

subroutine set_particle_from_rf_time (rf_time, ele, apply_hard_edge_offset, orbit)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct) orbit
  real(rp) rf_time
  logical apply_hard_edge_offset
end subroutine

subroutine set_status_flags (bookkeeping_state, stat)
  import
  implicit none
  type (bookkeeping_state_struct) bookkeeping_state
  integer stat
end subroutine

subroutine save_bunch_track (bunch, ele, s_travel)
  import
  implicit none
  type (bunch_struct) bunch
  type (ele_struct) ele
  real(rp) s_travel
end subroutine

subroutine sbend_body_with_k1_map (ele, param, n_step, orbit, mat6, make_matrix)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  integer n_step
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine set_on (key, lat, on_switch, orb)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), optional :: orb(0:)
  integer key
  logical on_switch
end subroutine

subroutine set_ele_attribute (ele, set_string, lat, err_flag, err_print_flag)
  import
  implicit none
  type (ele_struct) ele
  type (lat_struct) lat
  character(*) set_string
  logical err_flag
  logical, optional :: err_print_flag
end subroutine

subroutine set_ele_defaults (ele, do_allocate)
  import
  implicit none
  type (ele_struct) ele
  logical, optional :: do_allocate
end subroutine

subroutine set_tune (phi_a_set, phi_b_set, dk1, lat, orb, ok)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:)
  real(rp) phi_a_set
  real(rp) phi_b_set
  real(rp) dk1(:)
  logical ok
end subroutine

subroutine sol_quad_mat6_calc (ks, k1, s_len, ele, orbit, mat6, make_matrix)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct) orbit
  real(rp) ks, k1, s_len
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine solenoid_track_and_mat (ele, length, param, start_orb, end_orb, mat6)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) start_orb, end_orb
  real(rp) length
  real(rp), optional :: mat6(:,:)
end subroutine

subroutine split_lat (lat, s_split, ix_branch, ix_split, split_done, add_suffix, check_sanity, save_null_drift, err_flag, choose_max)
  import
  implicit none
  type (lat_struct), target :: lat
  real(rp) s_split
  integer ix_branch
  integer ix_split
  logical split_done
  logical, optional :: add_suffix, check_sanity, save_null_drift, err_flag, choose_max
end subroutine

subroutine tilt_coords (tilt_val, coord, mat6, make_matrix)
  import
  implicit none
  real(rp) tilt_val
  real(rp) coord(:)
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_beambeam (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_bend (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_drift (orb, length, mat6, make_matrix, include_ref_motion)
  import
  implicit none
  type (coord_struct) orb
  real(rp) length
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix, include_ref_motion
end subroutine

subroutine track_a_drift_photon (orb, length, phase_relative_to_ref)
  import
  implicit none
  type (coord_struct) orb
  real(rp) length
  logical phase_relative_to_ref
end subroutine

subroutine track_a_lcavity (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_mask (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_match (orbit, ele, param, err_flag, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix, err_flag
end subroutine

subroutine track_a_patch (ele, orbit, drift_to_exit, s_ent, ds_ref, track_spin, mat6, make_matrix)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) orbit
  real(rp), optional :: mat6(6,6), s_ent, ds_ref
  logical, optional :: drift_to_exit, track_spin, make_matrix
end subroutine

subroutine track_a_quadrupole (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_rfcavity (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_sad_mult (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_sol_quad (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_thick_multipole (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_wiggler (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_zero_length_element (start_orb, ele, param, end_orb, err_flag, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct), target :: ele
  type (lat_param_struct), target :: param
  logical err_flag
  type (track_struct), optional :: track
end subroutine

subroutine track_all (lat, orbit, ix_branch, track_state, err_flag, orbit0)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), allocatable, target :: orbit(:)
  type (coord_struct), optional, allocatable, target :: orbit0(:)
  integer, optional :: ix_branch, track_state
  logical, optional :: err_flag
end subroutine

subroutine track_from_s_to_s (lat, s_start, s_end, orbit_start, orbit_end, all_orb, ix_branch, track_state)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct) orbit_start, orbit_end
  type (coord_struct), optional, allocatable :: all_orb(:)
  real(rp) s_start, s_end
  integer, optional :: ix_branch, track_state
end subroutine

subroutine track_many (lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
  import
  implicit none
  type (lat_struct)  lat
  type (coord_struct)  orbit(0:)
  integer ix_start
  integer ix_end
  integer direction
  integer, optional :: ix_branch, track_state
end subroutine

subroutine track_backwards_time (lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
  import
  implicit none
  type (lat_struct)  lat
  type (coord_struct)  orbit(0:)
  integer ix_start
  integer ix_end
  integer direction
  integer, optional :: ix_branch, track_state
end subroutine

recursive subroutine track1 (start_orb, ele, param, end_orb, track, err_flag, ignore_radiation, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct)   :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
  logical, optional :: err_flag, ignore_radiation
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track1_backwards_time (end_orb, ele, param, start_orb, err_flag)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  logical, optional :: err_flag
end subroutine

subroutine track1_bmad (start_orb, ele, param, end_orb, err_flag, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  logical, optional :: err_flag
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track1_bmad_photon (start_orb, ele, param, end_orb, err_flag)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  logical, optional :: err_flag
end subroutine

subroutine track1_linear (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct), target :: ele
  type (lat_param_struct), target :: param
  logical err_flag
  type (track_struct), optional :: track
end subroutine

subroutine track1_symp_lie_ptc (start_orb, ele, param, end_orb, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
end subroutine

subroutine track1_symp_map (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_taylor (start_orb, ele, param, end_orb, taylor, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
  type (taylor_struct), optional, target :: taylor(6)
end subroutine

subroutine track1_time_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct), target :: ele
  type (lat_param_struct), target :: param
  logical err_flag
  type (track_struct), optional :: track
end subroutine

subroutine transfer_ac_kick (ac_kick_in, ac_kick_out)
  import
  implicit none
  type (ac_kicker_struct), pointer :: ac_kick_in, ac_kick_out
end subroutine transfer_ac_kick

subroutine transfer_map_calc (lat, t_map, err_flag, ix1, ix2, ref_orb, ix_branch, one_turn, unit_start)
  import
  implicit none
  type (lat_struct) lat
  type (taylor_struct) :: t_map(:)
  type (coord_struct), optional :: ref_orb
  integer, intent(in), optional :: ix1, ix2, ix_branch
  logical err_flag
  logical, optional :: one_turn, unit_start
end subroutine

subroutine transfer_mat_from_twiss (ele1, ele2, orb1, orb2, m)
  import
  implicit none
  type (ele_struct) ele1, ele2
  real(rp) orb1(6), orb2(6)
  real(rp) m(6,6)
end subroutine

subroutine transfer_matrix_calc (lat, xfer_mat, xfer_vec, ix1, ix2, ix_branch, one_turn)
  import
  implicit none
  type (lat_struct) lat
  real(rp) :: xfer_mat(:,:)
  real(rp), optional :: xfer_vec(:)
  integer, optional :: ix1, ix2, ix_branch
  logical, optional :: one_turn
end subroutine

subroutine transfer_wake (wake_in, wake_out)
  import
  implicit none
  type (wake_struct), pointer :: wake_in, wake_out
end subroutine

subroutine twiss_and_track_from_s_to_s (branch, orbit_start, s_end, orbit_end, &
                                                               ele_start, ele_end, err, compute_floor_coords)
  import
  implicit none
  type (coord_struct) :: orbit_start, orbit_end
  type (ele_struct), optional :: ele_start, ele_end
  type (branch_struct) branch
  real(rp) s_end
  logical, optional, intent(inout) :: err
  logical, optional :: compute_floor_coords
end subroutine

recursive subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_upstream_end, &
                 track_downstream_end, orbit_start, orbit_end, ele_start, ele_end, err, compute_floor_coords)
  import
  implicit none
  type (coord_struct), optional :: orbit_start, orbit_end
  type (ele_struct), optional :: ele_start, ele_end
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp) l_start, l_end
  logical track_upstream_end, track_downstream_end
  logical, optional :: err, compute_floor_coords
end subroutine

recursive subroutine twiss_at_element (ele, start_ele, end_ele, average)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), optional :: start_ele
  type (ele_struct), optional :: end_ele
  type (ele_struct), optional :: average
  integer ix_ele
end subroutine

subroutine twiss_at_start (lat, status, ix_branch)
  import
  implicit none
  type (lat_struct) lat
  integer, optional, intent(in) :: ix_branch
  integer, optional, intent(out) :: status
end subroutine

subroutine twiss1_propagate (twiss1, mat2, ele_type, length, twiss2, err)
  import
  implicit none
  type (twiss_struct) twiss1, twiss2
  integer ele_type
  real(rp) mat2(2,2), length
  logical err
end subroutine

subroutine twiss_from_mat6 (mat6, map0, ele, stable, growth_rate, status, type_out)
  import
  implicit none
  type (ele_struct) :: ele
  real(rp) :: mat6(:,:), map0(:)
  real(rp) :: growth_rate
  logical :: stable, type_out
  integer :: status
end subroutine

subroutine twiss_from_tracking (lat, ref_orb0, symp_err, err_flag, d_orb) 
  import
  type (lat_struct), intent(inout) :: lat
  type (coord_struct), intent(in) :: ref_orb0
  real(rp), intent(in), optional :: d_orb(:)   
  real(rp), intent(out) :: symp_err
  logical err_flag
end subroutine

subroutine twiss_propagate1 (ele1, ele2, err)
  import
  implicit none
  type (ele_struct) ele1
  type (ele_struct) ele2
  logical, optional :: err
end subroutine

subroutine twiss_propagate_all (lat, ix_branch, err_flag, ie_start, ie_end, zero_uncalculated)
  import
  implicit none
  type (lat_struct) lat
  integer, optional :: ix_branch, ie_start, ie_end
  logical, optional :: err_flag, zero_uncalculated
end subroutine

subroutine type_coord (coord)
  import
  implicit none
  type (coord_struct) coord
end subroutine

subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, twiss_out, &
      type_control, type_wake, type_floor_coords, type_field, type_wall, lines, n_lines)
  import
  implicit none
  type (ele_struct), target :: ele
  integer, optional, intent(in) :: type_mat6
  integer, optional, intent(out) :: n_lines
  integer, optional, intent(in) :: twiss_out
  logical, optional, intent(in) :: type_control, type_taylor, type_floor_coords
  logical, optional, intent(in) :: type_zero_attrib, type_wake
  logical, optional :: type_field, type_wall
  character(*), optional, allocatable :: lines(:)
end subroutine

subroutine type_twiss (ele, frequency_units, compact_format, lines, n_lines)
  import
  implicit none
  type (ele_struct) ele
  integer, optional :: frequency_units
  integer, optional :: n_lines
  character(*), optional :: lines(:)
  logical, optional :: compact_format
end subroutine

function value_of_attribute (ele, attrib_name, err_flag, err_print_flag) result (value)
  import
  implicit none
  type (ele_struct), target :: ele
  type (all_pointer_struct) a_ptr
  real(rp) value
  character(*) attrib_name
  logical, optional :: err_print_flag, err_flag
end function

subroutine write_digested_bmad_file (digested_name, lat,  n_files, file_names, extra, err_flag)
  import
  implicit none
  type (lat_struct), target, intent(in) :: lat
  integer, intent(in), optional :: n_files
  character(*) digested_name
  character(*), optional :: file_names(:)
  type (extra_parsing_info_struct), optional :: extra
  logical, optional :: err_flag
end subroutine

subroutine xsif_parser (xsif_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) xsif_file
  character(*), optional :: use_line
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  end subroutine

end interface

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
