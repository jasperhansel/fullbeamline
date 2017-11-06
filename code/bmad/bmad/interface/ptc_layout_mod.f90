!+
! Module ptc_layout_mod
!
! Module of PTC layout interface routines.
! Also see: ptc_interface_mod
!-

module ptc_layout_mod

use ptc_interface_mod
use multipass_mod
use track1_mod

implicit none

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine type_ptc_layout (lay)
!
! Subroutine to print the global information in a layout
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   lay - layout: layout to use.
!+

subroutine type_ptc_layout (lay)

use s_def_all_kinds, only: layout

type (layout) lay

!

if (.not. associated(lay%start)) then
  print *, 'Warning from TYPE_LAYOUT: Layout NOT Associated'
  return
endif

print *, 'Name:         ', lay%name
print *, 'N:            ', lay%N,        '  ! Number of Elements'
print *, 'LatPos:       ', lay%lastpos,  '  ! Last position'

end subroutine type_ptc_layout

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine lat_to_ptc_layout (lat, use_hard_edge_drifts)
!
! Subroutine to create a PTC layout from a Bmad lattice branch.
! Note: If lat_to_ptc_layout has already been setup, you should first do a 
!           call kill_ptc_layouts(lat)
! This deallocates the pointers in PTC
!
! Note: If not already done, before you call this routine you need to first call:
!    call set_ptc (...)
! [This is normally done in bmad_parser.]
!
! Note: If a Bmad element is using a hard edge model (EG: RFcavity element), there 
! will be three corresponding PTC fibre elements: (drift, RF. drift) for example.
! In this case, ele%ptc_fibre will be set to point to the last PTC fibre. That is the 
! exit end of ele will correspond to the exit end of ele%ptc_fibre.
!
! Note: Photon branches will not be included in the layout.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   lat                   -- lat_struct: Input lattice
!   use_hard_edge_drifts  -- logical, optional: Use bmad hard edge model for cavities, etc?
!                              default is set by bmad_com%use_hard_edge_drifts.
!
! Output:
!   lat%branch(:)%ptc              -- Pointers to generated layouts.
!   lat%branch(:)%ele(:)%ptc_fibre -- Pointer to PTC fibres
!-

subroutine lat_to_ptc_layout (lat, use_hard_edge_drifts)

use madx_ptc_module, only: m_u, m_t, fibre, append_empty_layout, survey, make_node_layout, &
                           append_point, set_up, ring_l

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, ele0
type (fibre), pointer :: save_fib
type (layout), pointer :: lay

real(rp) ptc_orientation(3,3), ang(3)

integer i, j, ix_pass
logical logic
logical, optional :: use_hard_edge_drifts

character(20), parameter :: r_name = 'lat_to_ptc_layout'

! Setup m_u

lat%ptc_uses_hard_edge_drifts = logic_option(bmad_com%use_hard_edge_drifts, use_hard_edge_drifts) 

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (branch%param%particle == photon$) cycle
  call branch_to_ptc_m_u (branch, use_hard_edge_drifts)
enddo

! setup m_t

do i = 0, ubound(lat%branch, 1)

  branch => lat%branch(i)
  if (branch%param%particle == photon$) cycle

  call append_empty_layout(m_t)
  call set_up (m_t%end)
  write (m_t%end%name, '(a, i4)') 'm_t Bmad branch:', i  ! For diagnostic purposes

  branch%ptc%m_t_layout => m_t%end   ! Save layout

  lay => m_t%end

  do j = 0, branch%n_ele_track
    ele => branch%ele(j)

    if (.not. associated(ele%ptc_fibre)) then
      call out_io (s_fatal$, r_name, 'NO FIBRE ASSOCIATED WITH ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
    endif

    if (tracking_uses_end_drifts(ele, use_hard_edge_drifts)) then
      call append_this_fibre(ele%ptc_fibre%previous)
    endif

    save_fib => ele%ptc_fibre%next
    call append_this_fibre(ele%ptc_fibre, .true.)

    if (tracking_uses_end_drifts(ele, use_hard_edge_drifts)) then
      call append_this_fibre(save_fib)
    endif

    ! Must add an energy patch if the reference energy shifts.

    if (j > 0) then
      ele0 => branch%ele(j-1)
      if (ele0%ptc_fibre%mag%p%p0c /= ele%ptc_fibre%mag%p%p0c) then
        save_fib => ele%ptc_fibre
        if (tracking_uses_end_drifts(ele, use_hard_edge_drifts))  save_fib => ele%ptc_fibre%previous
        save_fib%patch%energy = 1
      endif
    endif

  enddo

  lay%closed = .true.

  logic = .true.
  call ring_l(lay, logic)

  ele => branch%ele(0)
  ang = [ele%floor%phi, ele%floor%theta, ele%floor%psi]
  call bmad_patch_parameters_to_ptc (ang, ptc_orientation)

  call survey (m_t%end, ptc_orientation, ele%floor%r)

  !

  call make_node_layout(m_t%end)

enddo

!-----------------------------------------------------------------------------
contains

subroutine append_this_fibre(ele_fib, do_point)

type (fibre), pointer :: ele_fib
type (fibre), pointer :: this_fib
logical, optional :: do_point

!

call append_point(lay, ele_fib)
this_fib => lay%end

if (ele%key == patch$ .or. ele%key == floor_shift$) then
  this_fib%dir = ele%value(upstream_ele_dir$)
else
  this_fib%dir = ele%orientation
endif

this_fib%charge = charge_of(default_tracking_species(ele%branch%param)) / charge_of(ele%branch%param%particle)

if (logic_option(.false., do_point)) ele%ptc_fibre => this_fib

end subroutine append_this_fibre

end subroutine lat_to_ptc_layout 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine branch_to_ptc_m_u (branch, use_hard_edge_drifts)
!
! Subroutine to create a PTC layout from a Bmad lattice branch.
! Note: If lat_to_ptc_layout has already been setup, you should first do a 
!           call kill_ptc_layouts(lat)
! This deallocates the pointers in PTC
!
! Note: If not already done, before you call this routine you need to first call:
!    call set_ptc (...)
! [This is normally done in bmad_parser.]
!
! Note: If a Bmad element is using a hard edge model (EG: RFcavity element), there 
! will be three corresponding PTC fibre elements: (drift, RF. drift) for example.
! In this case, ele%ptc_fibre will be set to point to the last PTC fibre. That is the 
! exit end of ele will correspond to the exit end of ele%ptc_fibre.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   branch                -- branch_struct: Input branch.
!   use_hard_edge_drifts  -- logical, optional: Use bmad hard edge model for cavities, etc?
!                              default is set by bmad_com%use_hard_edge_drifts.
!
! Output:
!   branch(:)%ptc              -- Pointers to generated layouts.
!   branch(:)%ele(:)%ptc_fibre -- Pointer to PTC fibres
!-

subroutine branch_to_ptc_m_u (branch, use_hard_edge_drifts)

use s_fibre_bundle, only: ring_l, append, lp, layout, fibre
use mad_like, only: set_up, kill, lielib_print
use madx_ptc_module, only: m_u, m_t, append_empty_layout, survey, make_node_layout

type (branch_struct) :: branch
type (ele_struct) drift_ele
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: chain_ele(:)

integer n, ib, ie, ix_pass, ix_pass0, n_links
logical doneit, ele_inserted_in_layout
logical, optional :: use_hard_edge_drifts

! Transfer elements.

ix_pass0 = 0
ele_inserted_in_layout = .false.

call layout_init_stuff
call add_ptc_layout_to_list (branch%ptc,  m_u%end)   ! Save layout

do ie = 0, branch%n_ele_track
  ele => branch%ele(ie)

  call multipass_chain(ele, ix_pass, n_links, chain_ele)
  if (ix_pass > 1) then
    ele%ptc_fibre => chain_ele(1)%ele%ptc_fibre
    ix_pass0 = ix_pass
    cycle
  endif

  ! Use new layout for multipass regions.

  if (ix_pass0 /= ix_pass .and. ele_inserted_in_layout) then
    call layout_end_stuff
    call layout_init_stuff
    call add_ptc_layout_to_list (branch%ptc, m_u%end)   ! Save layout
    ele_inserted_in_layout = .false.
  endif

  ! If there are end drifts then ele%ptc_fibre points to middle element.

  if (tracking_uses_end_drifts(ele, use_hard_edge_drifts)) then
    call create_hard_edge_drift (ele, upstream_end$, drift_ele)
    call ele_to_fibre (drift_ele, drift_ele%ptc_fibre, branch%param, .true., &
                                    for_layout = .true., use_hard_edge_drifts = use_hard_edge_drifts)
  endif

  call ele_to_fibre (ele, ele%ptc_fibre, branch%param, .true., &
                                    for_layout = .true., use_hard_edge_drifts = use_hard_edge_drifts)

  if (tracking_uses_end_drifts(ele, use_hard_edge_drifts)) then
    call create_hard_edge_drift (ele, downstream_end$, drift_ele)
    call ele_to_fibre (drift_ele, drift_ele%ptc_fibre, branch%param, .true., &
                                    for_layout = .true., use_hard_edge_drifts = use_hard_edge_drifts)
  endif

  ele_inserted_in_layout = .true.
  ix_pass0 = ix_pass

enddo

call layout_end_stuff

! Set bookkeeping state

do ie = 0, branch%n_ele_max
  branch%ele(ie)%bookkeeping_state%ptc = ok$
enddo
branch%param%bookkeeping_state%ptc = ok$

!-----------------------------------------------------------------------------
contains

subroutine layout_init_stuff ()

call append_empty_layout(m_u)
call set_up(m_u%end)

write (m_u%end%name, '(a, i4)') 'm_u Bmad branch:', branch%ix_branch ! Diagnostic

end subroutine layout_init_stuff

!-----------------------------------------------------------------------------
! contains

subroutine layout_end_stuff ()

if (branch%param%geometry == closed$ .and. ie == branch%n_ele_track) then
  m_u%end%closed = .true.
else
  m_u%end%closed = .false.
endif

n = lielib_print(12)
lielib_print(12) = 0  ! No printing info messages

m_u%end%closed = .true.
doneit = .true.
call ring_l (m_u%end, doneit)
call survey (m_u%end)
call make_node_layout (m_u%end)

lielib_print(12) = n

end subroutine layout_end_stuff

end subroutine branch_to_ptc_m_u

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine add_ptc_layout_to_list (branch_ptc_info, layout_end)   
! 
! Routine to add a layout the a list of layouts.
!
! Module needed:
!   use ptc_layout_mod
!
! Input:
!   branch_ptc_info  -- ptc_branch1_info_struct: List of layouts
!   layout_end       -- layout: ptc layout
!
! Output:
!   branch_ptc_info  -- ptc_branch1_info_struct: Updated list.
!-

subroutine add_ptc_layout_to_list (branch_ptc_info, layout_end)   

type (ptc_branch1_info_struct) branch_ptc_info, temp_info
type (layout), target :: layout_end
integer n

!

if (allocated(branch_ptc_info%m_u_layout)) then
  n = size(branch_ptc_info%m_u_layout)
  call move_alloc(branch_ptc_info%m_u_layout, temp_info%m_u_layout)
  allocate(branch_ptc_info%m_u_layout(n+1))
  branch_ptc_info%m_u_layout(1:n) = temp_info%m_u_layout
  branch_ptc_info%m_u_layout(n+1)%ptr => layout_end
  deallocate(temp_info%m_u_layout)

else
  allocate(branch_ptc_info%m_u_layout(1))
  branch_ptc_info%m_u_layout(1)%ptr => layout_end  
endif 

end subroutine add_ptc_layout_to_list

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_one_turn_mat_and_closed_orbit_calc (branch, to_bmad_coords, pz)
!
! Routine to compute the one turn matrices and closed orbit for a lattice branch with closed geometry.
!
! Note: PTC itself does not compute Twiss parameters. Use twiss_from_mat6 to compute this.
!
! Input:
!   branch            -- branch_struct: Lattice branch.
!   to_bmad_coords    -- logical: If False then leave results in PTC phase space units.
!                                 If True, convert to Bmad phase space units
!   pz                -- real(rp), optional: energy offset around which to 
!                          calculate the matrices if there is no RF.
! Output:
!   branch            -- branch_struct: Lattice branch containing the matrices.
!     %ele(i)%fibre%i%m(6,6)    -- matrices.
!     %ele(i)%fibre%i%fix0(6)   -- Closed orbit.
!-

subroutine ptc_one_turn_mat_and_closed_orbit_calc (branch, to_bmad_coords, pz)

use s_fitting_new, only: default, compute_linear_one_turn_maps, fibre, layout

type (branch_struct), target :: branch
type (fibre), pointer :: ptc_fibre
type (layout), pointer :: ptc_layout

real(rp), optional :: pz
real(rp) c_mat(6,6), c_mat_inv(6,6)
integer i
logical to_bmad_coords

!

ptc_fibre => branch%ele(1)%ptc_fibre
call compute_linear_one_turn_maps (ptc_fibre, DEFAULT, pz)

if (to_bmad_coords) then
  ptc_layout => ptc_fibre%parent_layout
  do i = 1, ptc_layout%n
    call vec_ptc_to_bmad (ptc_fibre%i%fix0, ptc_fibre%beta0, ptc_fibre%i%fix0, c_mat)
    call mat_inverse(c_mat, c_mat_inv)
    ptc_fibre%i%m = matmul(matmul(c_mat, ptc_fibre%i%m), c_mat_inv)
    ptc_fibre => ptc_fibre%next
  enddo
endif

end subroutine ptc_one_turn_mat_and_closed_orbit_calc

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_transfer_mat_and_closed_orbit_calc (branch, to_bmad_coords, pz)
!
! Routine to compute the transfer matrices for the individual elements and closed orbit 
! for a lattice branch with closed geometry.
!
! Note: PTC itself does not compute Twiss parameters. Use twiss_from_mat6 to compute this.
!
! Input:
!   branch            -- branch_struct: Lattice branch.
!   to_bmad_coords    -- logical: If False then leave results in PTC phase space units.
!                                 If True, convert to Bmad phase space units
!   pz                -- real(rp), optional: energy offset around which to 
!                          calculate the matrices if there is no RF.
! Output:
!   branch            -- branch_struct: Lattice branch containing the matrices.
!     %ele(i)%fibre%i%m(6,6)    -- matrices.
!     %ele(i)%fibre%i%fix0(6)   -- Closed orbit at entrance.
!     %ele(i)%fibre%i%fix(6)    -- Closed orbit at exit.
!-

subroutine ptc_transfer_mat_and_closed_orbit_calc (branch, to_bmad_coords, pz)

use s_fitting_new, only: default, compute_linear_one_magnet_maps, fibre, layout

type (branch_struct) branch
type (fibre), pointer :: ptc_fibre
type (layout), pointer :: ptc_layout

real(rp), optional :: pz
real(rp) c_mat(6,6), c_mat_inv(6,6)
integer i
logical to_bmad_coords

!

ptc_fibre => branch%ele(1)%ptc_fibre
call compute_linear_one_magnet_maps (ptc_fibre, DEFAULT, pz)

if (to_bmad_coords) then
  ptc_layout => ptc_fibre%parent_layout
  do i = 1, ptc_layout%n
    call vec_ptc_to_bmad (ptc_fibre%i%fix0, ptc_fibre%beta0, ptc_fibre%i%fix0, c_mat)
    call vec_ptc_to_bmad (ptc_fibre%i%fix, ptc_fibre%beta0, ptc_fibre%i%fix)
    call mat_inverse(c_mat, c_mat_inv)
    ptc_fibre%i%m = matmul(matmul(c_mat, ptc_fibre%i%m), c_mat_inv)
    ptc_fibre => ptc_fibre%next
  enddo
endif

end subroutine ptc_transfer_mat_and_closed_orbit_calc

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_emit_calc (ele, norm_mode, sigma_mat, closed_orb)
!
! Routine to calculate emittances, etc.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input: 
!   ele -- ele_struct: Element at which to evaluate the parameters.
!
! Output:
!   norm_mode       -- Normal_modes_struct
!     %a%tune, %b%tune, %z%tune
!     %a%alpha_damp, etc.
!     %a%emittance, etc.
!   sigma_map(6,6)  -- real(rp): Sigma matrix (Bmad coordinates).
!   closed_orb      -- coord_struct: Closed orbit at ele (Bmad coordinates).
!                        Notice: This closed orbit is the closed orbit with radiation on.
!-

subroutine ptc_emit_calc (ele, norm_mode, sigma_mat, closed_orb)

use pointer_lattice

type (ele_struct) ele
type (normal_modes_struct) norm_mode
type (coord_struct) closed_orb
type (fibre), pointer :: ptc_fibre

real(rp) sigma_mat(6,6), emit(3), ptc_sigma_mat(6,6), tune(3), damp(3)
complex(rp) cmplx_sigma_mat(6,6)


ptc_fibre => ptc_reference_fibre(ele)
call radia_new (ele%branch%ptc%m_t_layout, ptc_fibre%pos, DEFAULT, fix = closed_orb%vec, em = emit, &
                sij = ptc_sigma_mat, sijr = cmplx_sigma_mat, tune = tune, damping = damp)

call init_coord(closed_orb, closed_orb, ele, downstream_end$)

norm_mode%a%tune = tune(1)   ! Fractional tune with damping
norm_mode%b%tune = tune(2)
norm_mode%z%tune = tune(3)

norm_mode%a%alpha_damp = damp(1)
norm_mode%b%alpha_damp = damp(2)
norm_mode%z%alpha_damp = damp(3)

norm_mode%a%emittance = emit(1)
norm_mode%b%emittance = emit(2)
norm_mode%z%emittance = emit(3)

call sigma_mat_ptc_to_bmad (ptc_sigma_mat, ele%ptc_fibre%beta0, sigma_mat)

call init (DEFAULT, ptc_com%taylor_order_ptc, 0)
end subroutine ptc_emit_calc 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function ptc_reference_fibre (ele) result (ref_fibre)
!
! Routine to return the reference fibre for a bmad element.
!
! The reference fibre is the fibre whose upstream edge corresponds
! to the downstream end fo the bmad element. 
!
! The reference fibre is so choisen since the referece edge of a
! Bmad element where such things as Twiss parameters are computed
! is the downstream edge while a PTC fibre uses the upstream edge 
! for the reference.
!
! Module Needed:
!   use patc_layout_mod
!
! Input:
!   ele   -- ele_struct: Bmad element
!
! Output:
!   ref_fibre -- fibre, pointer: Pointer to the corresponding reference fibre.
!-

function ptc_reference_fibre (ele) result (ref_fibre)

type (ele_struct), target :: ele
type (fibre), pointer :: ref_fibre

!

if (tracking_uses_end_drifts(ele, ele%branch%lat%ptc_uses_hard_edge_drifts)) then
  ref_fibre => ele%ptc_fibre%next%next
else
  ref_fibre => ele%ptc_fibre%next
endif

end function ptc_reference_fibre

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_track_all (branch, orbit, track_state, err_flag)
!
! Routine to track from the start to the end of a lattice branch. 
! 
! Modules Needed:
!   use ptc_layout_mod
!
! Input:
!   branch      -- lat_struct: Lat to track through.
!   orbit(0)    -- Coord_struct: Coordinates at beginning of branch.
!
! Output:
!   orbit(0:*)   -- Coord_struct: Orbit array.
!   track_state  -- Integer, optional: Set to moving_forward$ if everything is OK.
!                     Otherwise: set to index of element where particle was lost.
!   err_flag     -- Logical, optional: Set true if particle lost or error. False otherwise
!-

subroutine ptc_track_all (branch, orbit, track_state, err_flag)

use madx_ptc_module

type (branch_struct), target :: branch
type (coord_struct), allocatable :: orbit(:)
type (ele_struct), pointer :: ele
type (fibre), pointer :: fib

real(rp) vec(6)
real(dp) x(6)

integer i
integer, optional ::track_state

logical, optional :: err_flag

! Init

if (present(err_flag)) err_flag = .true.
if (present(track_state)) track_state = moving_forward$

!

if (orbit(0)%state == not_set$) call init_coord(orbit(0), orbit(0)%vec, branch%ele(0), downstream_end$) 
call vec_bmad_to_ptc (orbit(0)%vec, branch%ele(0)%value(p0c$) / branch%ele(0)%value(E_tot$), x)

do i = 1, branch%n_ele_track
  ele => branch%ele(i)

  orbit(i) = orbit(i-1)
  call check_aperture_limit (orbit(i), ele, first_track_edge$, branch%param)
  if (orbit(i)%state /= alive$) then
    if (present(err_flag)) err_flag = .false.
    return
  endif

  fib => branch%ele(i)%ptc_fibre%next
  call track_probe_x (x, DEFAULT, branch%ele(i-1)%ptc_fibre%next, fib)
  call vec_ptc_to_bmad (x, fib%beta0, vec)
  call init_coord (orbit(i), vec, ele, downstream_end$)

  call check_aperture_limit (orbit(i), ele, second_track_edge$, branch%param)
  if (orbit(i)%state /= alive$) then
    if (present(err_flag)) err_flag = .false.
    return
  endif

enddo

end subroutine ptc_track_all

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_closed_orbit_calc (branch, closed_orbit, radiation_damping_on)
!
! Routine to calculate the closed orbit of a lattice branch using PTC.
! This routine assumes the associated PTC layout has been crated 
! with lat_to_ptc_layout.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   branch          -- branch_struct: Branch of a lattice.
!   radiation_damping_on
!                   -- logical, optional: If True, radiation dampling is included in
!                        the calculation. 
!                        Default is the setting of bmad_com%%radiation_damping_on.
!
! Output:
!   closed_orbit(:) -- coord_struct, allocatable: closed_orbit
!-

subroutine ptc_closed_orbit_calc (branch, closed_orbit, radiation_damping_on)

use madx_ptc_module

type (branch_struct), target :: branch
type (coord_struct), allocatable :: closed_orbit(:)
type (fibre), pointer :: fib
type (internal_state) ptc_state

real(dp) x(6)
real(rp) vec(6)

integer i

logical, optional :: radiation_damping_on

!

if (logic_option(bmad_com%radiation_damping_on, radiation_damping_on)) then
  ptc_state = DEFAULT + radiation0
else
  ptc_state = DEFAULT - radiation0
endif

call reallocate_coord(closed_orbit, branch%n_ele_max)

x = 0
fib => branch%ele(0)%ptc_fibre%next
call find_orbit_x (x, ptc_state, 1.0d-7, fibre1 = fib)  ! find closed orbit
call vec_ptc_to_bmad (x, fib%beta0, vec)
call init_coord (closed_orbit(0), vec, branch%ele(0), downstream_end$)

do i = 1, branch%n_ele_track
  fib => branch%ele(i)%ptc_fibre%next
  call track_probe_x (x, ptc_state, branch%ele(i-1)%ptc_fibre%next, fib)
  call vec_ptc_to_bmad (x, fib%beta0, vec)
  call init_coord (closed_orbit(i), vec, branch%ele(i), downstream_end$)
enddo

end subroutine ptc_closed_orbit_calc 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ptc_one_turn_map_at_ele (ele, map, rf_on, pz, order)
!
! Routine to calculate the one turn map for a ring.
! Note: Use set_ptc(no_cavity = True/False) set turn on/off the RF cavities.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   ele     -- ele_struct: Element determining start/end position for one turn map.
!   rf_on   -- logical: calculate with RF on
!   pz      -- real(rp), optional: momentum deviation of closed orbit. 
!                                  Default = 0
!   order   -- integer, optional: Order of the map. If not given then default order is used.
!
! Output:
!   map(6)  -- taylor_struct: Bmad taylor map
!-

subroutine ptc_one_turn_map_at_ele (ele, map, rf_on, pz, order)

use madx_ptc_module

type (ele_struct), target :: ele
type (taylor_struct) map(6)
type (internal_state) ptc_state
type (fibre), pointer :: fib
type (damap) da_map
type (real_8) ray(6)

real(rp), optional :: pz
real(rp) x_bmad(6)
real(dp) x(6)

integer :: map_order
integer, optional :: order

logical rf_on 

!

map_order = integer_option(order, ptc_com%taylor_order_ptc)

if (rf_on) then
  ptc_state = default - nocavity0
else
  ptc_state = default + nocavity0
endif

call init(ptc_state, map_order, 0) ! The third argument is the number of parametric variables

! Find closed orbit

x = 0
if (present(pz)) x(5) = pz
fib => ele%ptc_fibre%next
call find_orbit_x (x, ptc_state, 1.0d-5, fibre1 = fib)  ! find closed orbit

! Construct map.

call alloc(da_map)
call alloc(ray)
da_map = 1   ! Identity
ray = da_map + x
call track_probe_x (ray, ptc_state, fibre1 = fib)

call vec_ptc_to_bmad (x, fib%beta0, x_bmad)
map%ref = x_bmad
call real_8_to_taylor(ray, fib%beta0, fib%beta0, map)

call kill(ray)
call kill(da_map)
call init (DEFAULT, ptc_com%taylor_order_ptc, 0)

end Subroutine ptc_one_turn_map_at_ele

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine normal_form_taylors(one_turn_taylor, rf_on, dhdj, A, A_inverse)
!
! Do a normal form decomposition on a one-turn taylor map M:
!   M = A o R o A_inverse
! where A maps Floquet (fully normalized) coordinates to lab coordinates. 
! In Floquet coordinates, the amplitudes are defined as J_i = (1/2) (x_i^2 + p_i^2).
! The map R = exp(:h:) is a pure rotation with h = h(J) is a function of the amplitudes only.
! The angles (phase advances) are given by phi_i = 2pi*dh/dJ_i.
! The taylor terms of dhdj are therefore the tunes, chromaticities, amplitude dependent tune shifts, etc.
!
! The mapping procedure for one turn is:
!  z_Floquet_in = A_inverse o z_Lab_in
!  [phi_a, phi_b, phi_c] = 2 pi * dhdj o z_Floquet_in
!  z_Floquet_out = RotationMatrix(phi_a, phi_b, phi_c) . z_Floquet_in
!  z_Lab_out = A o z_Floquet_out
!
!
! Module Needed:
!   use ptc_layout_mod
!
! Input: 
!   one_turn_taylor -- taylor_struct      : one turn taylor map
!   rf_on           -- logical            : Was the map calculated with RF on?
!
! Output:
!   A             -- taylor_struct, optional: Map from Floquet coordinates to Lab coordinates
!   A_inverse     -- taylor_struct, optional: Map from Lab coordinates to Floquet coordinates
!   dhdj            -- taylor_struct, optional: Map from Floquet coordinates to phase advances
!-
subroutine normal_form_taylors (one_turn_taylor, rf_on, dhdj, A, A_inverse)

use madx_ptc_module

type (taylor_struct) :: one_turn_taylor(6)
type (taylor_struct), optional :: A_inverse(6), dhdj(6), A(6)
type (damap) :: da_map
type (real_8) :: map8(6)
type (normalform) :: normal
type (internal_state) :: state
logical :: rf_on
!

if (rf_on) then
  state = default - nocavity0
else
  state = default + nocavity0
  ndpt_bmad = 1 ! Indicates that delta is in position 6 and not 5
endif

call init (state, ptc_com%taylor_order_ptc, 0) 
call alloc(map8)
call alloc(da_map)
call alloc(normal)

! Convert to real_8, then a damap
map8 = one_turn_taylor

! Get the normal form
da_map = map8
normal = da_map

! Convert to taylor_structs

! A
if (present(A)) then
  A(1:6) = normal%A_t%v(1:6)
endif

! A_inverse
if (present(A_inverse)) then
  map8 = normal%A_t**(-1)
  A_inverse = map8
endif

! dhdj 
if (present(dhdj)) then
  dhdj(1:6) = normal%dhdj%v(1:6)
endif

! Cleanup

call kill(map8)
call kill(da_map)
call kill(normal)

ndpt_bmad = 0
call init (DEFAULT, ptc_com%taylor_order_ptc, 0)

end subroutine normal_form_taylors

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine normal_form_complex_taylors
!
! UNDER DEVELOPMENT
!-
subroutine normal_form_complex_taylors(one_turn_taylor, rf_on, F, L, A, A_inverse, order)

use madx_ptc_module

type (taylor_struct) :: one_turn_taylor(6)
type (complex_taylor_struct), optional :: F(6), L(6)
type (taylor_struct), optional ::  A(6), A_inverse(6)
type (c_damap) :: cda, cdaLinear
type (damap) :: da
type (real_8) :: map8(6)
type (c_normal_form) :: complex_normal_form
type(c_vector_field) :: fvecfield
type (internal_state) :: state
integer :: i, order_for_normal_form
integer, optional :: order
logical :: rf_on, c_verbose_save
!

order_for_normal_form = integer_option(1, order)

if (rf_on) then
  state = default - nocavity0
else
  state = default + nocavity0
  ndpt_bmad = 1  ! Indicates that delta is in position 6 and not 5
endif

! Set PTC state
!no longer needed: use_complex_in_ptc=my_true
c_verbose_save = c_verbose
c_verbose = .false.
call init_all (state, ptc_com%taylor_order_ptc, 0) 
call alloc(map8)
call alloc(da)
call alloc(cda)
call alloc(cdaLinear)
call alloc(fvecfield)
call alloc(complex_normal_form)

! Convert to real_8, then a da map, then complex da map
map8 = one_turn_taylor
da = map8
cda = da

! Complex normal form in phasor basis
! See: fpp-ptc-read-only/build_book_example_g95/the_fpp_on_line_glossary/complex_normal.htm
! M = A o N o A_inverse.
call c_normal(cda, complex_normal_form, dospin=my_false, no_used=order_for_normal_form) 

if (present(F) .or. present(L)) then
  cda = complex_normal_form%N
  
  ! Move to the phasor basis
  cda=from_phasor(-1)*cda*from_phasor(1)

  !   N = L exp(F.grad), where L is the linear (rotation) map. 
  call c_factor_map(cda, cdaLinear, fvecfield, 1) ! 1 => Dragt-Finn direction

  ! Zero out small coefficients 
  call c_clean_vector_field(fvecfield, fvecfield, 1.d-8 )
  ! Output
endif

if (present(L)) then
  L(1:6) = cdaLinear%v(1:6)
endif

if (present(F)) then
  F(1:6) = fvecfield%v(1:6)
endif

if(present(A)) then
  da = complex_normal_form%a_t
  A(1:6) = da%v(1:6)
endif

if(present(A_inverse)) then
  da = complex_normal_form%a_t**(-1)
  A_inverse(1:6) = da%v(1:6)
endif

! Cleanup
call kill(map8)
call kill(da)
call kill(cda)
call kill(cdaLinear)
call kill(fvecfield)
call kill(complex_normal_form)

! Reset PTC state
use_complex_in_ptc=my_false
c_verbose = c_verbose_save
ndpt_bmad = 0
call init (DEFAULT, ptc_com%taylor_order_ptc, 0)

end subroutine normal_form_complex_taylors

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+

subroutine set_ptc_verbose(on)

use madx_ptc_module

logical :: on

!

c_verbose = on

end subroutine set_ptc_verbose

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine write_ptc_flat_file_lattice (file_name, branch)
!
! Routine to create a PTC flat file lattice from a Bmad branch.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   file_name     -- character(*): Flat file name.
!   branch        -- branch_struct: Branch containing a layout.
!-

subroutine write_ptc_flat_file_lattice (file_name, branch)

use pointer_lattice

type (branch_struct) branch
character(*) file_name

!

call print_complex_single_structure (branch%ptc%m_t_layout, file_name)

end subroutine write_ptc_flat_file_lattice

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine update_ptc_fibre_from_bmad (ele)
!
! Routine to update a fibre when the associated Bmad ele has been modified.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   ele           -- ele_struct: Element with corresponding PTC fibre.
!
! Output:
!   ele%ptc_fibre -- PTC fibre.
!-

subroutine update_ptc_fibre_from_bmad (ele)

use madx_ptc_module

type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (keywords) ptc_key
type (element), pointer :: mag
type (elementp), pointer :: magp
type (magnet_chart), pointer :: p, pp

real(rp) value, hk, vk, phi_tot
real(rp), pointer :: val(:)

integer i, ix

character(*), parameter :: r_name = 'update_ptc_fibre_from_bmad'

! "0" argument in add routine means set k/ks to value given.
! As opposed to "1" which means add to existing value.

branch => pointer_to_branch(ele)
call ele_to_an_bn (ele, branch%param, .false., ptc_key%list%k, ptc_key%list%ks, ptc_key%list%nmul)

! Must set all poles even if zero since they might have been non-zero beforehand.

do i = ptc_key%list%nmul, 1, -1
  call add (ele%ptc_fibre,  i, 0, ptc_key%list%k(i))
  call add (ele%ptc_fibre, -i, 0, ptc_key%list%ks(i))
enddo

! Note: ele_to_an_bn takes care of such stuff as sextupole strength conversion so
! only have to worry about non-multipole components here.

val => ele%value

mag  => ele%ptc_fibre%mag
magp => ele%ptc_fibre%magp

p => mag%p
pp => magp%p

! call set_integer (p%method, pp%method, nint(ele%value(integrator_order$)))
! call set_integer (p%nst, pp%nst, nint(ele%value(num_steps$)))

select case (ele%key)

case (elseparator$)
  hk = val(hkick$) / val(l$)
  vk = val(vkick$) / val(l$)
  if (hk == 0 .and. vk == 0) then
    ptc_key%tiltd = 0
  else
    if (branch%param%particle < 0) then
      hk = -hk
      vk = -vk
    endif
    ptc_key%tiltd = -atan2 (hk, vk) + val(tilt_tot$)
  endif
  call set_real (mag%volt, magp%volt, 1d-6 * val(e_tot$) * sqrt(hk**2 + vk**2))

case (solenoid$)
  call set_real (mag%b_sol, magp%b_sol, val(ks$))

case (sol_quad$)
  call set_real (mag%b_sol, magp%b_sol, val(ks$))

case (bend_sol_quad$)
  call set_real (mag%b_sol, magp%b_sol, val(ks$))

case (rfcavity$, lcavity$)
  phi_tot = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + ele%value(phi0_autoscale$))
  if (ele%key == lcavity$) phi_tot = pi / 2 - twopi * phi_tot
  call set_real (mag%phas, magp%phas, -mag%lag)

  call set_real (mag%volt, magp%volt, 2d-6 * e_accel_field(ele, voltage$))

case (sad_mult$)

  if (ele%value(l$) /= 0) then
    call set_real (mag%b_sol, magp%b_sol, val(ks$))
    call set_real (mag%va, magp%va, -sign(sqrt(24 * abs(ele%value(fq1$))), ele%value(fq1$)))
    call set_real (mag%vs, magp%vs, ele%value(fq2$))
  endif

  call set_real (mag%b_sol, magp%b_sol, val(ks$))

case (sbend$)
  call set_real (mag%hgap, magp%hgap, val(hgap$))
  call set_real (mag%fint, magp%fint, val(fint$))
  call set_real_all (mag%p%edge(1), magp%p%edge(1), val(e1$))
  call set_real_all (mag%p%edge(2), magp%p%edge(2), val(e2$))

end select

! Fringe

if (ele%key == sbend$) then
  ix = nint(val(ptc_fringe_geometry$))
  call set_logic (p%bend_fringe, pp%bend_fringe, (ix == x_invariant$))

  ix = nint(val(fringe_type$))
  select case (ix)
  case (none$)
    call set_integer (p%permfringe, pp%permfringe, 0)
  case (basic_bend$, linear_edge$)
    call set_integer (p%permfringe, pp%permfringe, 0)
  case (full$)
    call set_integer (p%permfringe, pp%permfringe, 1)
  case (hard_edge_only$)
    call set_integer (p%permfringe, pp%permfringe, 1)
  case (soft_edge_only$)
    call set_integer (p%permfringe, pp%permfringe, 2)
  case (sad_full$)
    call set_integer (p%permfringe, pp%permfringe, 3)
  end select

elseif (attribute_index(ele, 'FRINGE_TYPE') > 0) then  ! If fringe_type is a valid attribute
  call set_logic (p%bend_fringe, pp%bend_fringe, .false.)

  ix = nint(val(fringe_type$))
  select case (ix)
  case (none$)
    call set_integer (p%permfringe, pp%permfringe, 0)
  case (hard_edge_only$)
    call set_integer (p%permfringe, pp%permfringe, 1)
  case (soft_edge_only$)
    call set_integer (p%permfringe, pp%permfringe, 2)
  case (full$)
    call set_integer (p%permfringe, pp%permfringe, 3)
  end select

  if (ele%key == sad_mult$ .and. ele%value(l$) == 0) call set_integer (p%permfringe, pp%permfringe, 0)
endif

if (attribute_index(ele, 'FRINGE_AT') > 0) then  ! If fringe_at is a valid attribute
  ix = nint(val(fringe_at$))
  call set_logic (p%kill_ent_fringe, pp%kill_ent_fringe, (ix == downstream_end$ .or. ix == no_end$))
  call set_logic (p%kill_exi_fringe, pp%kill_exi_fringe, (ix == upstream_end$ .or. ix == no_end$))
endif

! misalign

call misalign_ele_to_fibre (ele, .true., ele%ptc_fibre)

ele%bookkeeping_state%ptc = ok$

!-------------------------------------------------------------------------
contains

subroutine set_real (to1, to2, value)
real(rp) value, to1
type(real_8) to2

to1 = value
to2 = value

end subroutine set_real

!-------------------------------------------------------------------------
! contains

subroutine set_real_all (to1, to2, value)
real(rp) value, to1, to2

to1 = value
to2 = value

end subroutine set_real_all

!-------------------------------------------------------------------------
! contains

subroutine set_integer (to1, to2, value)
integer to1, to2, value

to1 = value
to2 = value

end subroutine set_integer

!-------------------------------------------------------------------------
! contains

subroutine set_logic (to1, to2, value)
logical to1, to2, value

to1 = value
to2 = value

end subroutine set_logic

end subroutine update_ptc_fibre_from_bmad 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine update_bmad_ele_from_ptc (ele)
!
! Routine to update a bmad lattice element when the associated PTC fibre has been modified.
! Remember to call lattice_bookkeeper after calling this routine.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   ele           -- ele_struct: Element with corresponding ele%ptc_fibre fibre.
!
! Output:
!   ele       -- ele_struct: Modified element. 
!-

subroutine update_bmad_ele_from_ptc (ele)

use precision_constants, only: volt_c

type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (fibre), pointer :: fib

real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), tn(0:n_pole_maxx), tilt, kick
integer nmul, i, ix
character(40) name

!

fib => ele%ptc_fibre
branch => pointer_to_branch(ele)

call update_this_real (ele%value(l$), fib%mag%p%ld)

if (attribute_name(ele, num_steps$) == 'NUM_STEPS') then
  call update_this_real (ele%value(num_steps$), real(fib%mag%p%nst, rp))
  if (ele%value(num_steps$) == 0) then
    ele%value(ds_step$) = 0
  else
    ele%value(ds_step$) = ele%value(l$) / ele%value(num_steps$)
  endif
endif

! If integrator_order is defined for this element then update

name = attribute_name(ele, integrator_order$)
if (name(1:1) /= '!') call update_this_real (ele%value(integrator_order$), real(fib%mag%p%method, rp))

! Multipole

a_pole = 0
b_pole = 0
nmul = min(fib%mag%p%nmul, n_pole_maxx+1)
a_pole(0:nmul-1) = fib%mag%an(1:nmul)
b_pole(0:nmul-1) = fib%mag%bn(1:nmul)

call multipole_ab_to_kt (a_pole, b_pole, knl, tn)

! Electric Multipole

if (associated(fib%mag%tp10)) then
  if (associated(fib%mag%tp10%ae)) then
    if (.not. associated(ele%a_pole_elec)) call multipole_init(ele, electric$)
    nmul = min(size(fib%mag%tp10%ae), n_pole_maxx+1)
    ele%a_pole_elec(0:nmul-1) = 1d9 * VOLT_C * fib%mag%tp10%ae(1:nmul)
    ele%b_pole_elec(0:nmul-1) = 1d9 * VOLT_C * fib%mag%tp10%be(1:nmul)
  endif
endif

!

if (ele%key == sbend$) then
  call update_this_real (ele%value(ref_tilt_tot$), fib%mag%p%tiltd)
else
  call update_this_real (ele%value(tilt_tot$), fib%mag%p%tiltd)
endif

!

select case (ele%key)
case (ab_multipole$)
  ele%a_pole = a_pole
  ele%b_pole = b_pole

case (drift$)

! Use dsin & dcos due to bug in ifort 13.1 compiler. 
! Can rename to sin & cos when bug is fixed.

case (elseparator$)
  kick = fib%mag%volt * 1d6 / ele%value(e_tot$)
  if (branch%param%particle < 0) kick = -kick
  tilt = fib%mag%p%tiltd - ele%value(tilt_tot$)
  call update_this_real (ele%value(hkick$), -kick * dsin(tilt))
  call update_this_real (ele%value(vkick$), -kick * dcos(tilt))

case (hkicker$)
  call update_this_real (ele%value(kick$), knl(1))

case (lcavity$, rfcavity$)
  call update_this_real (ele%value(rf_frequency$), fib%mag%freq)

  select case (nint(ele%value(cavity_type$)))
  case (traveling_wave$);     call update_this_real (ele%value(voltage$), fib%mag%volt*1d6)
  case (standing_wave$);      call update_this_real (ele%value(voltage$), fib%mag%volt*0.5d6)
  case (ptc_standard$);       call update_this_real (ele%value(voltage$), fib%mag%volt*1d6)
  end select

  if (ele%key == lcavity$) then
    call update_this_real (ele%value(phi0$), pi/2 - fib%mag%lag/twopi)
  else
    call update_this_real (ele%value(phi0$), fib%mag%lag/twopi)
  endif

case (multipole$)
  ele%a_pole = knl
  ele%b_pole = tn

case (octupole$)
  call update_this_real (ele%value(k3$), knl(3))
  call update_this_real (ele%value(tilt$), tn(3))
  knl(3) = 0
  tn(3) = 0

case (quadrupole$)
  call update_this_real (ele%value(k1$), knl(1))
  call update_this_real (ele%value(tilt$), tn(1))
  knl(1) = 0
  tn(1) = 0

case (sbend$)
  call update_this_real (ele%value(g$), fib%mag%p%b0)
  call update_this_real (ele%value(angle$), ele%value(g$) * ele%value(l$))
  ix = nint(ele%value(ptc_field_geometry$))
  if (ix == straight$ .or. ix == true_rbend$) then
    call update_this_real (ele%value(e1$), fib%mag%p%edge(1) + ele%value(angle$)/2)
    call update_this_real (ele%value(e2$), fib%mag%p%edge(2) + ele%value(angle$)/2)
  else
    call update_this_real (ele%value(e1$), fib%mag%p%edge(1))
    call update_this_real (ele%value(e2$), fib%mag%p%edge(2))
  endif
  call update_this_real (ele%value(hgap$), fib%mag%hgap)
  call update_this_real (ele%value(fint$), fib%mag%fint)

case (sextupole$)
  call update_this_real (ele%value(k2$), knl(2))
  call update_this_real (ele%value(tilt$), tn(2))
  knl(2) = 0
  tn(2) = 0

case (solenoid$)
  call update_this_real (ele%value(ks$), fib%mag%b_sol)

case (sol_quad$)
  call update_this_real (ele%value(ks$), fib%mag%b_sol)
  call update_this_real (ele%value(k1$), knl(1))
  call update_this_real (ele%value(tilt$), tn(1))
  knl(1) = 0
  tn(1) = 0

case (wiggler$, undulator$)


case default
end select

! multipoles

if (any(knl /= 0)) then
endif


! kicks

! Fringes

!------------------------------------------------------------------------
contains

subroutine update_this_real (var, value)

real(rp) var, value

if (var == value) return

var = value
call set_flags_for_changed_attribute (ele, var)

end subroutine update_this_real

!------------------------------------------------------------------------
! contains

subroutine update_this_integer (var, value)

integer var, value

if (var == value) return

var = value
call set_flags_for_changed_attribute (ele, var)

end subroutine update_this_integer

end subroutine update_bmad_ele_from_ptc

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ptc_calculate_tracking_step_size (ptc_layout, kl_max, ds_max, 
!                                 even_steps, r_typical, dx_tol_bend, use_2nd_order)
!
! Routine to calculate the optimum number of tracking steps and order
! of the integrator (2, 4, or 6) for each fibre in a layout.
!
! See the Bmad manual chapter on PTC for more details.
!
! Module needed:
!   use ptc_interface_mod
!
! Input:
!   ptc_layout    -- layout:
!   kl_max        -- real(rp): Maximum K1*L per tracking step.
!   ds_max        -- real(rp): Maximum ds for any step. 
!                      Useful when including other physicas like space charge.
!   even_steps    -- logical, optional: Always use an even number of steps for a fibre?
!                      Useful if need to evaluate at the center of fibres.
!   r_typical     -- real(rp), optional: Typical transverse offset. Used for computing the
!                      effective contribution of K1*L due to sextupoles.
!   dx_tol_bend   -- real(rp): Tolerable residual orbit in a bend.
!   use_2nd_order -- logical, optional: If present and True then force the use of 2nd order
!                       integrator.
!   crossover(2)  -- integer, optional: crossover points between orders for all
!                       elements except wigglers. Default is [4, 18].
!   crossover_wiggler(2)
!                 -- integer, optional: crossover points for wigglers. Default is [30, 60].
!
! Output:
!   ptc_layout -- layout: Lattice with the optimum number of tracking steps and integrator order.
!-

subroutine ptc_calculate_tracking_step_size (ptc_layout, kl_max, ds_max, &
                    even_steps, r_typical, dx_tol_bend, use_2nd_order, crossover, crossover_wiggler)

use madx_ptc_module

type (layout) ptc_layout

real(rp) kl_max
real(rp), optional :: ds_max, dx_tol_bend, r_typical

integer, optional :: crossover(2), crossover_wiggler(2)
integer :: limit_int(2), lp14

logical, optional :: even_steps(2), use_2nd_order

!

resplit_cutting = 0   ! Ignore ds_max
if (present(ds_max)) then
  if (ds_max > 0) resplit_cutting = 2
endif

limit_int = [4, 18]

if (logic_option(.false., use_2nd_order)) limit_int = [10000, 10001] ! Something big.

lp14 = lielib_print(14)
lielib_print(14) = 0

call thin_lens_resplit (ptc_layout, kl_max, lim = crossover, &
            limit_wiggler = crossover_wiggler, lmax0 = ds_max, sexr = r_typical, xbend = dx_tol_bend)

lielib_print(14) = lp14

end subroutine ptc_calculate_tracking_step_size

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ptc_layouts_resplit (dKL_max, l_max, l_max_drift_only, bend_dorb, sex_dx, 
!                                                           even, crossover, crossover_wiggler)
!
! Routine to resplit (that is, recalculate the number of integration steps for an element)
! For the fibres in all layouts. After doing a resplit, the tune (and any other relavent
! "adjustable" parameters) should be adjusted to the correct values.
!
! Module needed:
!   use ptc_layout_mod
!
! Input:
!   dKL_max      -- real(rp): Maximum K1 * L quadrupole strength allowed for an integration step.
!                     Reasonable value would be something like 0.04.
!   l_max        -- real(rp): Maximum step length. Ignored if set to 0.
!   l_max_drift_only
!                -- logical: If True then l_max is only used for splitting drifts.
!   bend_dorb    -- real(rp): Residual bend orbit error. With some integration methods a zero
!                     orbit at the start of the bend will not be zero at the end. In this case,
!                     bend_dorb sets a maximum allowable orbit deviation. If set to zero, 
!                     this argument will be ignored. A resonable value is 10d-7. Note that the 
!                     actual orbit deviation is not simply related to bend_dorb and can be larger.
!                     In any case, lowering bend_dorb (without making it zero) will lower the 
!                     orbit deviation.
!   sex_dx       -- real(rp): To split sextupoles, sex_dx is used as the reference position 
!                     about which the quadrupole strength is calculated. This quadrupole strength
!                     is then used with dKL_max to calculate the number of integration steps.
!                     Set to zero to ignore.
!   even         -- logical, optional: If True then each fibre  will have an even number of steps.
!                     If False then the number of steps will be odd. If not present then number
!                     of steps is not constrained to be even or odd.
!   crossover(2) -- integer, optional: crossover(1) sets the maximum number of 2nd order integration
!                     steps to use. If the number of steps would exceed crossover(1) then integration
!                     is switched to 4th order. crossover(2) sets the maximum number of 4th order
!                     integration steps. If this number is exceeded, 6th order integration is used.
!                     Currently the default in PTC is [4, 18].
!   crossover_wiggler(2)
!                 -- integer, optional: crossover for wiggler elements.
!-


subroutine ptc_layouts_resplit (dKL_max, l_max, l_max_drift_only, bend_dorb, sex_dx, &
                                                            even, crossover, crossover_wiggler)

use s_fitting, only: thin_lens_restart, thin_lens_resplit
use madx_ptc_module, only: m_u, lielib_print

type(layout), pointer :: r

real(rp) dKL_max, bend_dorb, sex_dx, l_max
integer, optional :: crossover(2), crossover_wiggler(2)
integer lp14
logical l_max_drift_only
logical, optional :: even

!

r => m_u%start

lp14 = lielib_print(14)
lielib_print(14) = 0

call thin_lens_restart(r, universe=.true.)
call thin_lens_resplit(r, dKL_max, even, crossover, crossover_wiggler, l_max, bend_dorb, sex_dx, universe=.true.)

lielib_print(14) = lp14

end subroutine ptc_layouts_resplit

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ptc_check_for_lost_particle (state, ptc_fibre, do_reset)
!
! Routine to check if a particle has been lost when tracking with PTC.
!
! Input:
!   do_reset -- logical: If True then reset ptc flags.
!
! Output:
!   state     -- integer: Same as coord_struct%state. alive$, lost$, lost_neg_x_aperture$, etc.
!   ptc_fibre -- fibre, pointer: Pointer to fibre where particle lost. Nullified if particle alive.
!-

subroutine ptc_check_for_lost_particle (state, ptc_fibre, do_reset)

use definition, only: check_stable, lost_node, lost_fibre, xlost, reset_aperture_flag

type (fibre), pointer :: ptc_fibre
integer state
logical do_reset

!

if (check_stable) then
  state = alive$
  return
endif

!

ptc_fibre => lost_fibre
state = lost$

!! ptc_fibre%cylindrical_map%p%aperture%...
!! xlost(1:6) is lost position

!

if (do_reset) call reset_aperture_flag()

end subroutine ptc_check_for_lost_particle

end module
