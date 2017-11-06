!+
! Bmad_struct holds the structure definitions for Bmad routines.
!-

module bmad_struct

use taylor_mod
use complex_taylor_mod
use random_mod
use spline_mod

use definition, only: genfield, fibre, layout

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! IF YOU CHANGE THE LAT_STRUCT OR ANY ASSOCIATED STRUCTURES YOU MUST INCREASE THE VERSION NUMBER !!!
! THIS IS USED BY BMAD_PARSER TO MAKE SURE DIGESTED FILES ARE OK.

integer, parameter :: bmad_inc_version$ = 201

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   hard_ele     -- ele_struct, pointer: Points to element with the hard edge.
!                     Will be nullified if there is no hard edge.
!                     This will be track_ele (if it has a hard edge) unless track_ele is a super_slave.
!   s_edge_hard  -- Real(rp): S-position of next hard egde in hard_ele frame.
!   particle_at  -- Edge of hard_ele being tracked through. first_track_edge$, second_track_edge$, or none$

type fringe_edge_info_struct
  type (ele_struct), pointer :: hard_ele
  real(rp) s_edge_hard
  real(rp) ds_edge                     ! Distance from particle to edge
  integer particle_at                  !
  integer, pointer :: hard_location    ! Particle location wrt hard_ele.
  integer, allocatable :: location(:)  ! Particle location in an element.
                                       ! elements in list are the tracking element or its lords.
end type

!

integer, parameter :: n_pole_maxx = 21  ! maximum multipole order

integer, parameter :: old_control_var_offset$ = 1000  ! For indexing into ele%control_var(:) array
integer, parameter :: var_offset$ = 2000              ! Important: var_offset$ > old_control_var_offset$
integer, parameter :: taylor_offset$ = 1000000000     ! Taylor term index offset.

type expression_atom_struct
  character(40) :: name = ''
  integer :: type = 0
  real(rp) :: value = 0
end type

type controller_var_struct
  character(40) :: name = ''
  real(rp) :: value = 0
  real(rp) :: old_value = 0
end type

!-------------------------------------------------------------------------
! Note: custom$ = 7, and taylor$ = 8 are taken from the element key list.

integer, parameter :: bmad_standard$ = 1, symp_lie_ptc$ = 2
integer, parameter :: runge_kutta$ = 3 
integer, parameter :: linear$ = 4, tracking$ = 5, symp_map$ = 6
integer, parameter :: hard_edge_model$ = 9, symp_lie_bmad$ = 10, static$ = 11
integer, parameter :: boris$ = 12, mad$ = 14
integer, parameter :: time_runge_kutta$ = 15
integer, parameter :: n_methods$ = 15

character(16), parameter :: tracking_method_name(0:n_methods$) = [ &
      'GARBAGE!        ', 'Bmad_Standard   ', 'Symp_Lie_PTC    ', 'Runge_Kutta     ', &
      'Linear          ', 'GARBAGE!        ', 'Symp_Map        ', 'Custom          ', &
      'Taylor          ', 'GARBAGE!        ', 'Symp_Lie_Bmad   ', 'Static          ', &
      'Boris           ', 'GARBAGE!        ', 'MAD             ', 'Time_Runge_Kutta']

character(16), parameter :: spin_tracking_method_name(0:n_methods$) = [ &
      'GARBAGE!        ', 'GARBAGE!        ', 'Symp_Lie_PTC    ', 'GARBAGE!        ', &
      'GARBAGE!        ', 'Tracking        ', 'GARBAGE!        ', 'Custom          ', &
      'GARBAGE!        ', 'GARBAGE!        ', 'GARBAGE!        ', 'GARBAGE!        ', &
      'GARBAGE!        ', 'GARBAGE!        ', 'GARBAGE!        ', 'GARBAGE!        ']

character(16), parameter :: mat6_calc_method_name(0:n_methods$) = [ &
      'GARBAGE!        ', 'Bmad_Standard   ', 'Symp_Lie_PTC    ', 'GARBAGE!        ', &
      'Linear          ', 'Tracking        ', 'Symp_Map        ', 'Custom          ', &
      'Taylor          ', 'GARBAGE!        ', 'Symp_Lie_Bmad   ', 'Static          ', &
      'GARBAGE!        ', 'GARBAGE!        ', 'MAD             ', 'GARBAGE!       a']

integer, parameter :: drift_kick$ = 1, matrix_kick$ = 2, ripken_kick$ = 3
character(16), parameter :: ptc_integration_type_name(0:3) = [&
         'GARBAGE!   ', 'Drift_Kick ', 'Matrix_Kick', 'Ripken_Kick']

! sbend$ and rbend$ are from key definitions.

integer, parameter :: map_type$ = 1, periodic_type$ = 3
character(16), parameter :: sub_key_name(0:18) = ['GARBAGE!     ', 'Map          ', &
    'SBend        ', 'Periodic     ', 'GARBAGE!     ', 'GARBAGE!     ', &
    'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', &
    'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', &
    'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', &
    'RBend        ']

! field_calc names.
! Note: refer_to_lords is an "internal" value which is not valid for use in a lattice file.
!   The period in "Refer_to_Lords." is used to prevent sets in the lattice file.

integer, parameter :: fieldmap$ = 2, Refer_to_lords$ = 4, no_field$ = 5
character(16), parameter :: field_calc_name(0:7) = &
    ['GARBAGE!       ', 'Bmad_Standard  ', 'FieldMap       ', 'GARBAGE!       ', &
     'Refer_to_Lords.', 'No_Field       ', 'GARBAGE!       ', 'Custom         ']

! Crystal sub_key values.

integer, parameter :: bragg$ = 1, laue$ = 2
character(8), parameter :: diffraction_type_name(0:2) = ['GARBAGE!', 'Bragg   ', 'Laue    ']

! Distribution

integer, parameter :: uniform$ = 1, gaussian$ = 2, spherical$ = 3
character(12), parameter :: distribution_name(0:3) = ['GARBAGE! ', 'Uniform  ', 'Gaussian ', 'Spherical']

! Control element logicals
! Idea: Combine girder_lord, overlay_lord and group_lord -> control_lord

integer, parameter :: minor_slave$ = 1, super_slave$ = 2, free$ = 3
integer, parameter :: group_lord$ = 4, super_lord$ = 5, overlay_lord$ = 6
integer, parameter :: girder_lord$ = 7, multipass_lord$ = 8, multipass_slave$ = 9
integer, parameter :: not_a_lord$ = 10, slice_slave$ = 11, control_lord$ = 12

character(20), parameter :: control_name(12) = [character(20):: &
            'Minor_Slave', 'Super_Slave', 'Free', 'Group_Lord', &
            'Super_Lord', 'Overlay_Lord', 'Girder_Lord', 'Multipass_Lord ', &
            'Multipass_Slave', 'Not_a_Lord', 'Slice_Slave', 'Control_Lord']

logical, parameter :: set$ = .true., unset$ = .false.

!

integer, parameter :: grided$ = 123
integer, parameter :: auto_aperture$ = 1, rectangular$ = 2, elliptical$ = 3, wall3d$ = 5, custom_aperture$ = 7
character(16), parameter :: aperture_type_name(0:7) = &
                                    ['garbage!   ', 'Auto       ', 'Rectangular', 'Elliptical ', &
                                     'Surface    ', 'Wall3D     ', 'garbage!   ', 'Custom     ']

! fringe_type
! non-bend fringe type names are in the range fringe_type(1:n_non_bend_fringe_type$)

integer, parameter :: soft_edge_only$ = 2, hard_edge_only$ = 3, full$ = 4
integer, parameter :: sad_full$ = 5, linear_edge$ = 6, basic_bend$ = 7, test_edge$ = 8
integer, parameter :: n_non_bend_fringe_type$ = 4

character(16), parameter :: fringe_type_name(0:8) = ['Garbage!      ', &
                                   'None          ', 'Soft_Edge_Only', 'Hard_edge_only', 'Full          ', &
                                   'SAD_Full      ', 'Linear_Edge   ', 'Basic_Bend    ', 'Test          ']

character(16), parameter :: higher_order_fringe_type_name(0:4) = fringe_type_name(0:4)

integer, parameter :: standing_wave$ = 1, traveling_wave$ = 2, ptc_standard$ = 3
character(16), parameter :: cavity_type_name(0:3) = ['Garbage!      ', 'Standing_Wave ', 'Traveling_Wave', 'PTC_Standard  ']


integer, parameter :: x_invariant$ = 1, multipole_symmetry$ = 2
character(16), parameter :: ptc_fringe_geometry_name(0:2) = ['Garbage!          ', 'x_invariant       ', 'multipole_symmetry']

integer, parameter :: control_var$ = 1, old_control_var$ = 2, all_control_var$ = 3, elec_multipole$ = 4

!

integer, parameter :: ok$ = 1, in_stop_band$ = 2, non_symplectic$ = 3, unstable$ = 4
integer, parameter :: unstable_a$ = 5, unstable_b$ = 6
integer, parameter :: xfer_mat_calc_failure$ = 7, twiss_propagate_failure$ = 8, no_closed_orbit$ = 9

character(24) :: matrix_status_name(9) = [character(24) :: 'OK', 'IN_STOP_BAND', 'NON_SYMPLECTIC', &
                       'UNSTABLE', 'UNSTABLE A-MODE', 'UNSTABLE B-MODE', 'XFER_MAT_CALC_FAILURE', &
                       'TWISS_PROPAGATE_FAILURE', 'NO_CLOSED_ORBIT']


type twiss_struct
  real(rp) :: beta = 0, alpha = 0, gamma = 0, phi = 0, eta = 0, etap = 0
  real(rp) :: sigma = 0, sigma_p = 0, emit = 0, norm_emit = 0
end type

!-------------------------------------------------------------------------
! Structure for holding the photon reflection probability tables.

type interval1_coef_struct
  real(rp) c0, c1, n_exp
end type

type photon_reflect_table_struct
  real(rp), allocatable :: angle(:)              ! Vector of angle values for %p_reflect
  real(rp), allocatable :: energy(:)             ! Vector of energy values for %p_reflect
  type (interval1_coef_struct), allocatable :: int1(:)
  real(rp), allocatable :: p_reflect(:,:)        ! (angle, ev) Logarithm of smooth surface reflection probability
  real(rp) max_energy                            ! maximum energy for this table
  real(rp), allocatable :: p_reflect_scratch(:)       ! Scratch space
end type

! Each photon_reflect_reflect_table_array(:) represents a different surface type.
! photon_reflect_reflect_table_array(1) is initialized by the photon_reflection_init routine
! All others can be set by an outside programmer. 

type photon_reflect_surface_struct
  character(40) :: name = ''
  character(80) :: description = ''                    ! Descriptive name
  character(200) :: reflectivity_file = ''
  type (photon_reflect_table_struct), allocatable :: table(:)
  real(rp) :: surface_roughness_rms = 0       ! sigma in Dugan's notation
  real(rp) :: roughness_correlation_len = 0   ! T in Dugan's notation
  logical :: initialized = .false.
  integer :: ix_surface = -1
end type

!-------------------------------------------------------------------------

integer, parameter :: ascii$ = 1, binary$ = 2, hdf5$ = 3

! num_ele_attrib$ is size of ele%value(:) array.

integer, parameter :: num_ele_attrib$ = 80

integer, parameter :: x_polarization$ = 2, y_polarization$ = 3
character(8), parameter :: sr_polarization_name(3) = ['None  ', 'X_Axis', 'Y_Axis']

integer, parameter :: off$ = 1, on$ = 2
integer, parameter :: none$ = 1

integer, parameter :: horizontally_pure$ = 2, vertically_pure$ = 3
character(20), parameter :: exact_multipoles_name(3) = [character(20):: 'Off', 'Horizontally_Pure', 'Vertically_Pure']

! Polarization is not 1 when the spin_polar struct represents an ensamble of spins.

type spin_polar_struct
  real(rp) :: polarization = 1
  real(rp) :: theta = 0
  real(rp) :: phi   = 0
  real(rp) :: xi    = 0
end type

real(rp), parameter :: x_unit_vec(3) = [1, 0, 0], y_unit_vec(3) = [0, 1, 0], z_unit_vec(3) = [0, 0, 1]

integer, parameter :: magnetic$ = 1, electric$ = 2, mixed$ = 3
character(8), parameter :: em_field_type_name(3) = ['Magnetic', 'Electric', 'Mixed   ']

! Diffraction

integer, parameter :: bragg_diffracted$ = 1, forward_diffracted$ = 2, undiffracted$ = 3
character(20), parameter :: ref_orbit_follows_name(0:3) = [character(20) :: 'GARBAGE!', &
                                             'Bragg_Diffracted', 'Forward_Diffracted', 'Undefracted']

integer, parameter :: reflection$ = 1, transmission$ = 2
character(16), parameter :: mode_name(0:2) = [character(16) :: 'GARBAGE!', 'Reflection', 'Transmission']

! wall3d definitions.

integer, parameter :: anchor_beginning$ = 1, anchor_center$ = 2, anchor_end$ = 3
character(12), parameter :: anchor_pt_name(0:3) = ['GARBAGE! ', 'Beginning', 'Center   ', 'End      ']

! Note: upstream_end$ = entrance_end$ & downstream_end$ = exit_end$

! first_track_edge$ is the edge a particle enters the element at. 
! This edge will depend upon whether a particle is moving in +s or -s direction.
! Similarly, second_track_edge$ is the edge a particle leaves the element at.

integer, parameter :: entrance_end$ = 1, exit_end$ = 2, both_ends$ = 3, no_end$ = 4, no_aperture$ = 4
integer, parameter :: continuous$ = 5, surface$ = 6, wall_transition$ = 7
integer, parameter :: first_track_edge$ = 11, second_track_edge$ = 12, in_between$ = 13

character(16), parameter :: aperture_at_name(0:7) = [ &
      'GARBAGE!       ', 'Entrance_End   ', 'Exit_End       ', 'Both_Ends      ', &
      'No_Aperture    ', 'Continuous     ', 'Surface        ', 'Wall_Transition']

character(16), parameter :: end_at_name(0:4) = [ &
      'GARBAGE!     ', 'Entrance_End ', 'Exit_End     ', 'Both_Ends    ', &
      'No_End       ']

integer, parameter :: upstream_end$ = 1, downstream_end$ = 2
integer, parameter :: inside$ = 3, center_pt$ = 3, start_end$ = 99

character(16), parameter :: ref_pt_name(0:3) = [ &
      'GARBAGE!      ', 'Entrance_End  ', 'Exit_End      ', 'Center        ']

character(16), parameter :: location_name(0:3) = [ &
      'GARBAGE!      ', 'Upstream_End  ', 'Downstream_End', 'Inside        ']

! Structures for defining cross-sections of beam pipes and capillaries
! A cross-section is defined by an array v(:) of wall3d_section_vertex_structs.
! Each vertex v(i) defines a point on the pipe/capillary.
! Vertices are connected by straight lines, circular arcs, or ellipses.
! The radius and tilt values are for the arc from the preceding vertex to this one.
! For v(1), the radius and tilt values are for the arc between v(n) and v(1) where
!   n = upper bound of v(:) array.

integer, parameter :: normal$ = 1, clear$ = 2, opaque$ = 3, trunk$ = 4, trunk1$ = 5, trunk2$ = 6
integer, parameter :: leg1$ = 7, leg2$ = 8, wall_start$ = 9, wall_end$ = 10
character(16), parameter :: wall3d_section_type_name(10) = [ &
                     'Normal     ', 'Clear      ', 'Opaque     ', 'Trunk      ', &
                     'Trunk1     ', 'Trunk2     ', 'Leg1       ', 'Leg2       ', &
                     'Wall_Start ', 'Wall_End   ']

! Note: Component order in wall3d_vertex_struct is important since sr3d_read_wall_file uses
! a namelist read to input vertex points

type wall3d_vertex_struct
  real(rp) :: x = 0, y = 0      ! Coordinates of the vertex.
  real(rp) :: radius_x = 0      ! Radius of arc or ellipse x-axis half width. 0 => Straight line.
  real(rp) :: radius_y = 0      ! Ellipse y-axis half height. 
  real(rp) :: tilt = 0          ! Tilt of ellipse
  real(rp) :: angle = 0         ! Angle of (x, y) point.
  real(rp) :: x0 = 0, y0 = 0    ! Center of ellipse
  integer :: type = normal$     ! No longer used.
end type

! A beam pipe or capillary cross-section is a collection of vertexes.
! Vertices are always ordered in increasing angle.
! Note: %surface is not saved in digested files.

type wall3d_section_struct
  character(40) :: name = ''                ! Identifying name
  character(20) :: material = ''            ! Material.
  type (wall3d_vertex_struct), allocatable :: v(:)  ! Array of vertices. Always stored relative.
  type (photon_reflect_surface_struct), pointer :: surface => null()
                                            ! Surface reflectivity tables.
  integer :: type = normal$                 ! normal$, clear$, opaque$, trunk$, trunk1$, leg2$, ...
  integer :: n_vertex_input = 0             ! Number of vertices specified by the user.
  integer :: ix_ele = 0                     ! index of lattice element containing section
  integer :: ix_branch = 0                  ! Index of branch lattice element is in.
  logical :: patch_in_region = .false.      ! Patch element exists between this section and previous one?
  logical :: absolute_vertices_input = .false. ! Input vertex numbers absolute or relative to r0?
  real(rp) :: thickness = -1                ! Material thickness.
  real(rp) :: s = 0                         ! Longitudinal position
  real(rp) :: r0(2) = 0                     ! Center of section
  ! Section-to-section spline interpolation of the center of the section
  real(rp) :: dx0_ds = 0                    ! Center of wall derivative
  real(rp) :: dy0_ds = 0                    ! Center of wall derivative
  real(rp) :: x0_coef(0:3) = 0              ! Spline coefs for x-center
  real(rp) :: y0_coef(0:3) = 0              ! Spline coefs for y-center
  ! Section-to_section spline interpolation of the wall.
  real(rp) :: dr_ds = real_garbage$         ! derivative of wall radius 
  real(rp) :: p1_coef(3) = 0                ! Spline coefs for p0 function
  real(rp) :: p2_coef(3) = 0                ! Spline coefs for p1 function
end type

! If, say, %ele_anchor_pt = center$ then center of wall is at the center of the element.

integer, parameter :: chamber_wall$ = 1, mask_plate$ = 2
character(12), parameter :: wall3d_name(2) = [character(12) :: 'Chamber_Wall', 'Mask_Plate']

type wall3d_struct
  character(40) :: name = ''
  integer :: type = chamber_wall$                 ! or mask_plate$
  integer :: ix_wall3d = 0                        ! Index in branch%wall3d(:) array.
  integer :: n_link = 1                           ! For memory management of ele%wall3d
  real(rp) :: thickness = -1                      ! For diffraction_plate elements
  character(20) :: clear_material = ''            !
  character(20) :: opaque_material = ''           !
  logical :: superimpose = .false.                ! Can overlap another wall
  integer :: ele_anchor_pt = anchor_beginning$    ! anchor_beginning$, anchor_center$, or anchor_end$
  type (wall3d_section_struct), allocatable :: section(:) ! Indexed from 1.
end type  

! plane list, etc

integer, parameter :: x_plane$ = 1, y_plane$ = 2
integer, parameter :: z_plane$ = 3, n_plane$ = 4, s_plane$ = 5

character(1), parameter :: plane_name(6) = ['X', 'Y', 'Z', 'N', 'S', ' ']
character(2), parameter :: field_plane_name(3) = ['Bx', 'By', 'Bz']

! coordinate def
! Use coord_state_name for getting the string representation of coord%state

integer, parameter :: moving_forward$ = -9

integer, parameter :: alive$ = 1, lost$ = 2
integer, parameter :: lost_neg_x_aperture$ = 3, lost_pos_x_aperture$ = 4 
integer, parameter :: lost_neg_y_aperture$ = 5, lost_pos_y_aperture$ = 6
integer, parameter :: lost_pz_aperture$ = 7  ! Particle "turned around" when not tracking with time_runge_kutta.

type coord_struct                 ! Particle coordinates at a single point
  real(rp) :: vec(6) = 0          ! (x, px, y, py, z, pz)
  real(rp) :: s = 0               ! Longitudinal position 
  real(rp) :: t = 0               ! Absolute time (not relative to reference).
  real(rp) :: spin(3) = 0         ! Spin.
  real(rp) :: field(2) = 0        ! Photon E-field intensity (x,y).
  real(rp) :: phase(2) = 0        ! Photon E-field phase (x,y)
  real(rp) :: charge = 0          ! Macro charge (Coul) of particle. 
  real(rp) :: path_len = 0        ! path length (used by coherent photons).
  real(rp) :: p0c = 0             ! For non-photons: Reference momentum.
                                  !     For photons: Photon momentum (not reference).
  real(rp) :: beta = -1           ! Velocity / c_light.
  integer :: ix_ele = -1          ! Index of element particle was tracked through.
                                  !   May be -1 if element is not associated with a lattice.
  integer :: ix_user = -1         ! For general use, not used by Bmad.
  integer :: state = not_set$     ! alive$, lost$, lost_neg_x_aperture$, etc.
  integer :: direction = 1        ! Sign of longitudinal direction of motion (ds/dt). = +/- 1.
                                  !  This is independent of the element orientation. See manual.
  integer :: species = not_set$   ! positron$, proton$, etc.  
  integer :: location = upstream_end$  ! upstream_end$, inside$, or downstream_end$
end type

type coord_array_struct
  type (coord_struct), allocatable :: orbit(:)
end type

! Coupling structure

type bpm_phase_coupling_struct
  real(rp) K_22a  ! In-phase y/x for a-mode oscillations.
  real(rp) K_12a  ! Out-of-phase y/x for a-mode oscillations.
  real(rp) K_11b  ! In-phase x/y for b-mode oscillations.
  real(rp) K_12b  ! Out-of-phase x/y for b-mode oscillations.
  real(rp) Cbar22_a ! Cbar22 as calculated from K_22a.
  real(rp) Cbar12_a ! Cbar12 as calculated from K_12a.
  real(rp) Cbar11_b ! Cbar11 as calculated from K_11b.
  real(rp) Cbar12_b ! Cbar12 as calculated from K_12b.
  real(rp) phi_a    ! a-mode betatron phase.
  real(rp) phi_b    ! b-mode betatron phase.
end type

! Wakefield structs...

integer, parameter :: linear_leading$ = 2, linear_trailing$ = 3
character(16), parameter :: sr_transverse_dependence_name(3) = ['none           ', 'linear_leading ', 'linear_trailing']

type wake_sr_mode_struct    ! Psudo-mode Short-range wake struct 
  real(rp) :: amp = 0       ! Amplitude
  real(rp) :: damp = 0      ! Dampling factor.
  real(rp) :: k = 0         ! k factor
  real(rp) :: phi = 0       ! Phase in radians
  real(rp) :: b_sin = 0     ! non-skew (x) sin-like component of the wake
  real(rp) :: b_cos = 0     ! non-skew (x) cos-like component of the wake
  real(rp) :: a_sin = 0     ! skew (y) sin-like component of the wake
  real(rp) :: a_cos = 0     ! skew (y) cos-like component of the wake
  integer :: polarization = none$                 ! none$, x_axis$ or y_axis$
  integer :: transverse_dependence = not_set$     ! linear_leading$, linear_trailing, none$
end type

type wake_sr_struct  ! Psudo-mode short-Range Wake struct 
  type (wake_sr_mode_struct), allocatable :: mode(:)
  real(rp) :: z_ref = 0    ! z reference value for computing the wake amplitude.
                           !  This is used to prevent value overflow with long bunches.
end type

! Each wake_lr_struct represents a different mode.
! A non-zero lr_freq_spread attribute value will make freq different from freq_in.

type wake_lr_mode_struct    ! Long-Range Wake struct.
  real(rp) :: freq = 0       ! Actual Frequency in Hz.
  real(rp) :: freq_in = 0    ! Input frequency in Hz.
  real(rp) :: R_over_Q = 0   ! Strength in V/C/m^2.
  real(rp) :: Q = 0          ! Quality factor.
  real(rp) :: angle = 0      ! polarization angle (radians/2pi).
  real(rp) :: b_sin = 0      ! non-skew sin-like component of the wake.
  real(rp) :: b_cos = 0      ! non-skew cos-like component of the wake.
  real(rp) :: a_sin = 0      ! skew sin-like component of the wake.
  real(rp) :: a_cos = 0      ! skew cos-like component of the wake.
  real(rp) :: t_ref = 0      ! time reference value for computing the wake amplitude.
                             !  This is used to prevent value overflow with long trains.
  integer :: m = 0           ! Order (1 = dipole, 2 = quad, etc.)
  logical :: polarized = .false.   ! Polaraized mode?
end type

type wake_lr_position1_struct
  real(rp) :: vec(6) = 0   ! (x, px, y, py, z, pz) of bunch centroid
  real(rp) :: charge = 0   ! Bunch charge (Coul).
  real(rp) :: t = 0        ! Time of bunch
end type

type wake_lr_spline_struct
  type (spline_struct), allocatable :: spline(:)
  type (wake_lr_position1_struct), allocatable :: bunch(:) ! Array of past bunch positions
  real(rp) :: t_max = 0         ! Maximum time beyond which the wake is taken to be zero.
  real(rp) :: polarization_angle = 0      ! polarization angle (radians/2pi).
  logical :: polarized = .false.   ! Polaraized mode?
  integer :: transverse_dependence = not_set$     ! linear_leading$, linear_trailing, none$
end type

!

type wake_struct
  character(200) :: sr_file = ''
  character(200) :: lr_file = ''
  type (wake_sr_struct) :: sr_long
  type (wake_sr_struct) :: sr_trans
  type (wake_lr_mode_struct), allocatable :: lr_mode(:)
  type (wake_lr_spline_struct), allocatable :: lr_spline(:)
  real(rp) :: z_sr_max = 0              ! Max allowable z value sr_mode. 
  real(rp) :: lr_freq_spread = 0        ! Random frequency spread of long range modes.
  logical :: lr_self_wake_on = .true.   ! Long range self-wake used in tracking?
end type

! ac_kicker structure

type ac_kicker_time_struct
  real(rp) :: amp = 0
  real(rp) :: time = 0
  type (spline_struct) :: spline = spline_struct()
end type

type ac_kicker_freq_struct
  real(rp) :: f = 0
  real(rp) :: amp = 0
  real(rp) :: phi = 0
end type

type ac_kicker_struct
  type (ac_kicker_time_struct), allocatable :: amp_vs_time(:)
  type (ac_kicker_freq_struct), allocatable :: frequencies(:)
end type

! Cartesian field decomposition

integer, parameter :: hyper_y_family_y$ = 1, hyper_xy_family_y$ = 2, hyper_x_family_y$ = 3
integer, parameter :: hyper_y_family_x$ = 4, hyper_xy_family_x$ = 5, hyper_x_family_x$ = 6
integer, parameter :: hyper_y_family_qu$ = 7, hyper_xy_family_qu$ = 8, hyper_x_family_qu$ = 9
integer, parameter :: hyper_y_family_sq$ = 10, hyper_xy_family_sq$ = 11, hyper_x_family_sq$ = 12
integer, parameter :: x_family$ = 1, y_family$ = 2, qu_family$ = 3, sq_family$ = 4

character(20), parameter :: cartesian_map_term_name(0:12) = [ 'GARBAGE!          ', &
                                        'Hyper_Y_Family_Y  ', 'Hyper_XY_Family_Y ', 'Hyper_X_Family_Y  ', &
                                        'Hyper_Y_Family_X  ', 'Hyper_XY_Family_X ', 'Hyper_X_Family_X  ', &
                                        'Hyper_Y_Family_QU ', 'Hyper_XY_Family_QU', 'Hyper_X_Family_QU ', &
                                        'Hyper_Y_Family_SQ ', 'Hyper_XY_Family_SQ', 'Hyper_X_Family_SQ ']

type cartesian_map_term1_struct
  real(rp) :: coef = 0
  real(rp) :: kx = 0, ky = 0, kz = 0
  real(rp) :: x0 = 0, y0 = 0, phi_z = 0
  integer :: type = 0        ! hyper_y_family_x$, etc.
end type

type cartesian_map_term_struct
  character(200) :: file = ''      ! Input file name. Used also as ID for instances. 
  integer :: n_link = 1            ! For memory management of %term
  type (cartesian_map_term1_struct), allocatable :: term(:)  
end type

type cartesian_map_struct
  real(rp) :: field_scale = 1      ! Factor to scale the fields by
  real(rp) :: r0(3) = 0            ! Field origin offset.
  integer :: master_parameter = 0  ! Master parameter in ele%value(:) array to use for scaling the field.
  integer :: ele_anchor_pt = anchor_beginning$  ! anchor_beginning$, anchor_center$, or anchor_end$
  integer :: field_type = magnetic$
  type (cartesian_map_term_struct), pointer :: ptr
end type

! Cylindrical map field 

type cylindrical_map_term1_struct
  complex(rp) :: e_coef = 0
  complex(rp) :: b_coef = 0
end type

type cylindrical_map_term_struct
  character(200) :: file = ''   ! Input file name. Used also as ID for instances. 
  integer :: n_link = 1         ! For memory management of this structure
  type (cylindrical_map_term1_struct), allocatable :: term(:)
end type

type cylindrical_map_struct
  integer m                        ! Azimuthal Mode: varies as cos(m*phi - theta0_azimuth)
  integer :: harmonic = 0          ! Harmonic of fundamental
  real(rp) :: phi0_fieldmap = 0    ! Mode oscillates as: twopi * (f * t + phi0_fieldmap)
  real(rp) :: theta0_azimuth = 0    ! Azimuthal ((x, y) plane) orientation of mode.
  real(rp) :: field_scale = 1      ! Factor to scale the fields by
  integer :: master_parameter = 0  ! Master parameter in ele%value(:) array to use for scaling the field.
  integer :: ele_anchor_pt = anchor_beginning$  ! anchor_beginning$, anchor_center$, or anchor_end$
  real(rp) :: dz = 0               ! Distance between sampled field points.
  real(rp) :: r0(3) = 0            ! Field origin offset.
  type (cylindrical_map_term_struct), pointer :: ptr
end type

! Grid field

type grid_field_pt1_struct
  complex(rp) :: E(3) = 0
  complex(rp) :: B(3) = 0
end type

type grid_field_pt_struct
  character(200) :: file = ''   ! Input file name. Used also as ID for instances. 
  integer :: n_link = 1         ! For memory management of this structure
  type (grid_field_pt1_struct), allocatable :: pt(:,:,:)
end type

type grid_field_struct
  integer :: geometry = 0          ! Type of grid structure
  integer :: harmonic = 0          ! Harmonic of fundamental
  real(rp) :: phi0_fieldmap = 0    ! Mode oscillates as: twopi * (f * t + phi0_fieldmap)
  real(rp) :: field_scale = 1      ! Factor to scale the fields by
  integer :: field_type = mixed$
  integer :: master_parameter = 0  ! Master parameter in ele%value(:) array to use for scaling the field.
  integer :: ele_anchor_pt = anchor_beginning$  ! anchor_beginning$, anchor_center$, or anchor_end$
  real(rp) :: dr(3) = 0   ! Grid spacing.
  real(rp) :: r0(3) = 0   ! Field origin relative to ele_anchor_pt.
  logical :: curved_ref_frame = .false.
  type (grid_field_pt_struct), pointer :: ptr
end type

! The Taylor field is a set of evenly spaced planes.
! For each plane a taylor series is used to calculate the field components.

type taylor_field_plane1_struct
  type (em_taylor_struct) :: field(3)    ! [Bx, By, Bz] or [Ex, Ey, Ez]
end type

type taylor_field_plane_struct
  character(200) :: file = ''   ! Input file name. Used also as ID for instances. 
  integer :: n_link = 1         ! For memory management of this structure
  type (taylor_field_plane1_struct), allocatable :: plane(:)
end type

type taylor_field_struct
  integer :: ele_anchor_pt = anchor_beginning$  ! anchor_beginning$, anchor_center$, or anchor_end$
  integer :: field_type = magnetic$
  real(rp) :: dz = 0               ! Plane spacing.
  real(rp) :: r0(3) = 0            ! field origin relative to ele_anchor_pt.
  real(rp) :: field_scale = 1      ! Factor to scale the fields by
  integer :: master_parameter = 0  ! Master parameter in ele%value(:) array to use for scaling the field.
  logical :: curved_ref_frame = .false.
  logical :: canonical_tracking = .false.
  type (taylor_field_plane_struct), pointer :: ptr
end type

! Local reference frame position with respect to the global (floor) coordinates

real(rp), parameter :: r0_vec(3) = 0
real(rp), parameter :: w_unit(3,3) = reshape( [1, 0, 0, 0, 1, 0, 0, 0, 1], [3,3])

type floor_position_struct
  real(rp) :: r(3) = 0                      ! (x, y, z) offset from origin
  real(rp) :: w(3,3) =  w_unit              ! W matrix. Columns are unit vectors of the frame axes.
  real(rp) :: theta = 0, phi = 0, psi = 0   ! angular orientation consistent with W matrix
end type

! Space charge structure. This structure contains information about the beam as a whole.

type space_charge_struct
  type (coord_struct) closed_orb   ! beam orbit
  real(rp) kick_const
  real(rp) sig_x
  real(rp) sig_y
  real(rp) phi      ! Rotation angle to go from lab frame to rotated frame.
  real(rp) sin_phi
  real(rp) cos_phi
  real(rp) sig_z
end type    

type xy_disp_struct
  real(rp) :: eta = 0, etap = 0
end type

! Structure to hold the information of where an individual element is in the lattice.

type lat_ele_loc_struct
  integer :: ix_ele = -1
  integer :: ix_branch = 0
end type

! Structure to be used for an array of pointers to elements.
! The id component is not set by any Bmad routines and can be used, for example, by 
! programs that handle multiple lattices to indicate which lattice 
! the element pointer is pointing to.
! A pointer to an element in a lattice is not usable if the number of elements in the 
! lattice is modified (so that the lattice element array is reallocated). In this case,
! the %loc component is potentially useful (as long as the element pointed to does not move).

type ele_pointer_struct
  type (ele_struct), pointer :: ele => null()
  type (lat_ele_loc_struct) :: loc = lat_ele_loc_struct()
  integer :: id = -1                    ! For general use. Not used by Bmad.
end type

! Structure for ptc genfield

type ptc_genfield_struct
  type (genfield), pointer :: field => null()           ! For symp_map$
  real(rp) :: vec0(6) = 0                               ! constant part of the genfield map.
end type

! The mode3_struct is used for normal mode analysis of the full 6x6 transfer matrix.

type mode3_struct
  real(rp) v(6,6)
  type (twiss_struct) a, b, c
  type (twiss_struct) x, y
end type

integer, parameter :: super_ok$ = 0, stale$ = 2
integer, parameter :: attribute_group$ = 1, control_group$ = 2, floor_position_group$ = 3
integer, parameter :: s_position_group$ = 4, ref_energy_group$ = 5, mat6_group$ = 6
integer, parameter :: rad_int_group$ = 7, all_groups$ = 8, s_and_floor_position_group$ = 9

! The bookkeeping_state_struct is used for keeping track of what bookkeeping has
! been done on an element. NOTE: The information in this structure is ignored if 
! bmad_com%auto_bookkeeper = True which is the default.
! See the Bmad manual for more details.

type bookkeeping_state_struct
  integer :: attributes = stale$      ! Element dependent attributes: super_ok$, ok$ or stale$
  integer :: control = stale$         ! Lord/slave bookkeeping status: super_ok$, ok$ or stale$ 
  integer :: floor_position = stale$  ! Global (floor) geometry: super_ok$, ok$ or stale$
  integer :: s_position = stale$      ! Longitudinal position & element length: super_ok$, ok$ or stale$
  integer :: ref_energy = stale$      ! Reference energy and ref time: super_ok$, ok$ or stale$
  integer :: mat6 = stale$            ! Linear transfer map status: super_ok$, ok$ or stale$
  integer :: rad_int = stale$         ! Radiation integrals cache status
  integer :: ptc = stale$             ! Associated PTC fibre (or layout) status.
end type

! radiation integral data cache

type rad_int_ele_cache_struct
  real(rp) :: orb0(6) = 0       ! Reference orbit for the calculation
  real(rp) :: g2_0 = 0          ! g2 factor when orbit = %vec0
  real(rp) :: g3_0 = 0          ! g3 factor when orbit = %vec0
  real(rp) :: dg2_dorb(6) = 0   ! variation of g2 with respect to orbit.
  real(rp) :: dg3_dorb(6) = 0   ! Variation of g3 with respect to orbit.
  logical :: stale = .true.
end type 

! Structure for surfaces of detectors, mirrors, crystals, etc.
! Rule: This structure is always allocated in the ele_struct for elements that can utilize it.

type surface_orientation_struct
  real(rp) :: x_pitch = 0, y_pitch = 0, x_pitch_rms = 0, y_pitch_rms = 0
end type

type surface_grid_pt_struct
  type (surface_orientation_struct) :: orientation = surface_orientation_struct()
  ! Photon statistics...
  integer :: n_photon = 0
  complex(rp) :: E_x  = 0, E_y = 0
  real(rp) ::  intensity_x = 0, intensity_y = 0, intensity = 0
  real(rp) :: orbit(6) = 0            ! x, Vx/c, y, Vy/c, dummy, E - E_ref.
  real(rp) :: orbit_rms(6) = 0        ! RMS statistics.
  real(rp) :: init_orbit(6) = 0       ! Initial orbit at start of lattice statistics.
  real(rp) :: init_orbit_rms(6) = 0   ! Initial orbit at start of lattice RMS statistics.
end type

integer, parameter :: segmented$ = 2, h_misalign$ = 3, diffract_target$ = 4
character(16), parameter :: surface_grid_type_name(0:4) = ['GARBAGE!       ', 'Off            ',  &
                                        'Segmented      ', 'H_Misalign     ', 'Diffract_Target']

type surface_grid_struct
  character(200) :: file = ''
  integer :: type = off$   ! or segmented$, or h_misalign$, or diffract_target$
  real(rp) :: dr(2) = 0, r0(2) = 0
  type (surface_grid_pt_struct), allocatable :: pt(:,:) 
end type

! Scratch space for segmented surface calculations

type segmented_surface_struct
  integer :: ix = int_garbage$, iy = int_garbage$    ! Index of segment
  real(rp) :: x0 = 0, y0 = 0, z0 = 0         ! Center of segment
  real(rp) :: slope_x = 0, slope_y = 0       ! Slopes of segment
end type

! Surface container structure

type photon_surface_struct
  type (surface_grid_struct) :: grid = surface_grid_struct('', off$, 0, 0, null())
  type (segmented_surface_struct) :: segment = segmented_surface_struct()
  real(rp) :: curvature_xy(0:6,0:6) = 0
  logical :: has_curvature = .false.     ! Dependent var. Will be set by Bmad
end type

! Target points are in element coordinates.

type target_point_struct
  real(rp) :: r(3) = 0   ! (x, y, z)
end type

type photon_target_struct
  logical :: deterministic_grid = .false.     ! Use ix/iy_grid instead of random number?  
  integer :: ix_grid = 0, iy_grid = 0         ! Grid pt to go to if deterministic_grid = T.
  integer :: type = off$                      ! or rectangular$, or grid$
  integer :: n_corner = 0
  type (lat_ele_loc_struct) :: ele_loc = lat_ele_loc_struct()
  type (target_point_struct) :: corner(8) = target_point_struct()
  type (target_point_struct) :: center = target_point_struct()
end type

type photon_material_struct
  complex(rp) :: f0_m1 = 0, f0_m2 = 0 ! For multilayer_mirror only.
  complex(rp) :: f_0 = 0
  complex(rp) :: f_h = 0     ! Structure factor for H direction.
  complex(rp) :: f_hbar = 0  ! Structure factor for -H direction.
  complex(rp) :: f_hkl = 0   ! = sqrt(f_h * f_hbar)
  real(rp) :: h_norm(3) = 0  ! Normalized H vector for crystals.
  real(rp) :: l_ref(3) = 0   ! Crystal reference orbit displacement vector in element coords.
end type

! photon_element_struct is an ele_struct component holding photon parameters

type photon_element_struct
  type (photon_surface_struct) :: surface = photon_surface_struct()
  type (photon_target_struct) :: target !! = photon_target_struct() IFORT bug prevents set!
  type (photon_material_struct) :: material = photon_material_struct()
end type

!------------------------------------------------------------------------------
! Beam structures

type bunch_struct
  type (coord_struct), allocatable :: particle(:)
  integer, allocatable :: ix_z(:)  ! bunch%ix_z(1) is index of head particle, etc.
  real(rp) :: charge_tot = 0  ! Total charge in a bunch (Coul).
  real(rp) :: charge_live = 0 ! Charge of live particles (Coul).
  real(rp) :: z_center = 0    ! Longitudinal center of bunch (m). Note: Generally, z_center of 
                              !   bunch #1 is 0 and z_center of the other bunches is negative.
  real(rp) :: t_center = 0    ! Center of bunch creation time relative to head bunch.
  integer :: ix_ele = 0       ! Element this bunch is at.
  integer :: ix_bunch = 0     ! Bunch index. Head bunch = 1, etc.
  integer :: n_live = 0
end type

type beam_struct
  type (bunch_struct), allocatable :: bunch(:)
end type

type ellipse_beam_init_struct
  integer :: part_per_ellipse = 0  ! number of particles per ellipse
  integer :: n_ellipse = 1         ! number of ellipses (>= 1)
  real(rp) :: sigma_cutoff = 0     ! sigma cutoff of the representation
end type

type kv_beam_init_struct
  integer :: part_per_phi(2) = 0    ! number of particles per angle variable.
  integer :: n_I2 = 0               ! number of I2
  real(rp) :: A = 0                 ! A = I1/e
end type

type grid_beam_init_struct
  integer :: n_x = 0         ! Number of columns.
  integer :: n_px = 0        ! Number of rows.
  real(rp) :: x_min = 0      ! Lower x limit.
  real(rp) :: x_max = 0      ! Upper x limit.
  real(rp) :: px_min = 0     ! Lower px limit.
  real(rp) :: px_max = 0     ! Upper px limit,
end type

type beam_init_struct
  character(200) :: file_name = ''           ! For distribution_type(1) = 'FILE'
  character(16) :: distribution_type(3) = '' ! distribution type (in x-px, y-py, and z-pz planes)
                                             ! "ELLIPSE", "KV", "GRID", "FILE", "", or "RAN_GAUSS" 
  type(spin_polar_struct) :: spin = spin_polar_struct()                      ! Initialize the spin
  type (ellipse_beam_init_struct) :: ellipse(3) = ellipse_beam_init_struct() ! Ellipse beam distribution
  type (kv_beam_init_struct) :: KV = kv_beam_init_struct()                   ! KV beam distribution
  type (grid_beam_init_struct) :: grid(3) = grid_beam_init_struct()          ! Grid beam distribution
  !!! The following are for Random distributions
  real(rp) :: center_jitter(6) = 0.0  ! Bunch center rms jitter
  real(rp) :: emit_jitter(2)   = 0.0  ! a and b bunch emittance rms jitter normalized to emittance
  real(rp) :: sig_z_jitter     = 0.0  ! bunch length RMS jitter 
  real(rp) :: sig_e_jitter     = 0.0  ! energy spread RMS jitter 
  integer :: n_particle = 0           ! Number of random particles per bunch.
  logical :: renorm_center = .true.   ! Renormalize centroid?
  logical :: renorm_sigma = .true.    ! Renormalize sigma?
  character(16) :: random_engine = 'pseudo' ! Or 'quasi'. Random number engine to use. 
  character(16) :: random_gauss_converter = 'exact'  
                                            ! Or 'quick'. Uniform to gauss conversion method.
  real(rp) :: random_sigma_cutoff = -1      ! Cut-off in sigmas.
  !!! The following are used  by all distribution types
  real(rp) :: a_norm_emit = 0                ! a-mode normalized emittance (emit * gamma)
  real(rp) :: b_norm_emit = 0                ! b-mode normalized emittance (emit * gamma)
  real(rp) :: a_emit = 0                     ! a-mode emittance
  real(rp) :: b_emit = 0                     ! b-mode emittance
  real(rp) :: dPz_dz = 0                     ! Correlation of Pz with long position.
  real(rp) :: center(6) = 0                  ! Bench center offset relative to reference.
  real(rp) :: dt_bunch = 0                   ! Time between bunches.
  real(rp) :: sig_z = 0                      ! Z sigma in m.
  real(rp) :: sig_e = 0                      ! e_sigma in dE/E.
  real(rp) :: bunch_charge = 1               ! charge (Coul) in a bunch. Default is something non-zero.
  integer :: n_bunch = 1                     ! Number of bunches.
  character(16) :: species = ""              ! "positron", etc. "" => use referece particle.
  logical :: init_spin     = .false.         ! initialize beam spinors
  logical :: full_6D_coupling_calc = .false. ! Use V from 6x6 1-turn mat to match distribution?  
                                             !   Else use 4x4 1-turn mat used.
  logical :: use_lattice_center = .false.    ! Use beam_start[...] in lattice rather than beam_init%center(:)?
  logical :: use_t_coords = .false.          ! If true, the distributions will be taken as in t-coordinates  
  logical :: use_z_as_t   = .false.          ! Only used if  use_t_coords = .true.
                                             !   If true,  z describes the t distribution 
                                             !   If false, z describes the s distribution
end type

! The routines calc_bunch_params and calc_bunch_params_slice calculate bunch parameters.

type bunch_params_struct
  type (twiss_struct) :: x, y, z   ! Projected Twiss parameters
  type (twiss_struct) :: a, b, c   ! Normal mode twiss parameters
  type (coord_struct) :: centroid  ! Lab frame
  type (spin_polar_struct) :: spin ! polarization
  real(rp) sigma(6,6)              ! beam size matrix
  real(rp) rel_max(6)              ! Max orbit relative to centroid
  real(rp) rel_min(6)              ! Min orbit relative to_centroid
  real(rp) s                       ! Longitudinal position.
  real(rp) charge_live             ! Charge of all non-lost particle
  integer n_particle_tot           ! Total number of particles
  integer n_particle_live          ! Number of non-lost particles
  integer n_particle_lost_in_ele   ! Number lost in element (not calculated by Bmad)
end type

!-------------------------------------------------------------------------
! Ele_struct:
! Remember: If this struct is changed you have to:
!     Increase bmad_inc_version by 1.
!     run scripts to regenerate cpp_bmad_interface library.
!     Modify:
!       read_digested_bmad_file
!       write_digested_bmad_file
!       parsing routines to read in modified/new parameters...
!       deallocate_ele_pointers
!       ele_equal_ele
!       type_ele
!       write_bmad_lattice_file
!       pointer_to_attribute
!       pointers_to_attribute


type ele_struct
  character(40) :: name = '<Initialized>'                ! name of element.
  character(40) :: type = ''                             ! type name.
  character(40) :: alias = ''                            ! Another name.
  character(40) :: component_name = ''                   ! Used by overlays, multipass patch, etc.
  character(200), pointer :: descrip => null()           ! Description string.
  type (twiss_struct) :: a = twiss_struct()              ! Twiss parameters at end of element
  type (twiss_struct) :: b = twiss_struct()              ! Twiss parameters at end of element
  type (twiss_struct) :: z = twiss_struct()              ! Twiss parameters at end of element
  type (xy_disp_struct) :: x = xy_disp_struct()          ! Projected dispersions.
  type (xy_disp_struct) :: y = xy_disp_struct()          ! Projected dispersions.
  type (ac_kicker_struct), pointer :: ac_kick => null()  ! ac_kicker element parameters.
  type (bookkeeping_state_struct) :: bookkeeping_state = bookkeeping_state_struct() ! Attribute bookkeeping
  type (branch_struct), pointer :: branch => null()                      ! Pointer to branch containing element.
  type (controller_var_struct), pointer :: control_var(:) => null()      ! group & overlay variables.
  type (cartesian_map_struct), pointer :: cartesian_map(:) => null()     ! Used to define E/M fields
  type (cylindrical_map_struct), pointer :: cylindrical_map(:) => null() ! Used to define E/M fields
  type (ele_struct), pointer :: lord => null()                           ! Pointer to a slice lord.
  type (taylor_field_struct), pointer :: taylor_field(:) => null()       ! Used to define E/M fields.
  type (grid_field_struct), pointer :: grid_field(:) => null()           ! Used to define E/M fields.
  type (fibre), pointer :: ptc_fibre => null()                           ! PTC tracking.
  type (floor_position_struct) :: floor = floor_position_struct(r0_vec, w_unit, 0.0_rp, 0.0_rp, 0.0_rp)
                                                                     ! Reference position in global coords.
  type (ptc_genfield_struct) :: ptc_genfield = ptc_genfield_struct() ! For symp_map$
  type (mode3_struct), pointer :: mode3 => null()                    ! 6D normal mode structure.
  type (photon_element_struct), pointer :: photon => null()
  type (rad_int_ele_cache_struct), pointer :: rad_int_cache => null() 
                                                               ! Radiation integral calc cached values 
  type (space_charge_struct), pointer :: space_charge => null()
  type (taylor_struct) :: taylor(6) = taylor_struct()          ! Phase space Taylor map.
  type (taylor_struct) :: spin_taylor(3,3) = taylor_struct()   ! Spin Taylor map.
  type (wake_struct), pointer :: wake => null()                ! Wakes
  type (wall3d_struct), pointer :: wall3d(:) => null()         ! Chamber or capillary wall
  type (coord_struct) :: map_ref_orb_in = coord_struct()       ! Entrance end transfer map ref orbit
  type (coord_struct) :: map_ref_orb_out = coord_struct()      ! Exit end transfer map ref orbit
  type (coord_struct) :: time_ref_orb_in = coord_struct()      ! Reference orbit at entrance end for ref_time calc.
  type (coord_struct) :: time_ref_orb_out = coord_struct()     ! Reference orbit at exit end for ref_time calc.
  real(rp) :: value(num_ele_attrib$) = 0                       ! attribute values.
  real(rp) :: old_value(num_ele_attrib$) = 0                   ! Used to see if %value(:) array has changed.
  real(rp) :: vec0(6) = 0                                      ! 0th order transport vector.
  real(rp) :: mat6(6,6) = 0                                    ! 1st order transport matrix.
  real(rp) :: c_mat(2,2) = 0                                   ! 2x2 C coupling matrix
  real(rp) :: gamma_c = 1                                      ! gamma associated with C matrix
  real(rp) :: s_start = 0                                      ! longitudinal ref position at entrance_end
  real(rp) :: s = 0                                            ! longitudinal ref position at the exit end.
  real(rp) :: ref_time = 0                                     ! Time ref particle passes exit end.
  real(rp), pointer :: r(:,:,:) => null()                      ! For general use. Not used by Bmad.
  real(rp), pointer :: a_pole(:) => null()                     ! knl for multipole elements.
  real(rp), pointer :: b_pole(:) => null()                     ! tilt for multipole elements.
  real(rp), pointer :: a_pole_elec(:) => null()                ! Electrostatic multipoles.
  real(rp), pointer :: b_pole_elec(:) => null()                ! Electrostatic multipoles.
  integer :: key = 0                              ! key value 
  integer :: sub_key = 0                          ! For wigglers: map_type$, periodic_type$
  integer :: ix_ele = -1                          ! Index in lat%branch(n)%ele(:) array [n = 0 <==> lat%ele(:)].
  integer :: ix_branch = 0                        ! Index in lat%branch(:) array [0 => In lat%ele(:)].
  integer :: lord_status = not_a_lord$            ! Type of lord element this is. overlay_lord$, etc.
  integer :: n_slave = 0                          ! Number of slaves (except field slaves) of this element.
  integer :: n_slave_field = 0                    ! Number of field slaves of this element
  integer :: ix1_slave = 0                        ! Pointer index to this element's slaves.
  integer :: slave_status = free$                 ! Type of slave element this is. multipass_slave$, slice_slave$, etc.
  integer :: n_lord = 0                           ! Number of lords (except field lords).
  integer :: n_lord_field = 0                     ! Number of field lords of this element
  integer :: ic1_lord = 0                         ! Pointer index to this element's lords.
  integer :: ix_pointer = 0                       ! For general use. Not used by Bmad.
  integer :: ixx = 0, iyy = 0                     ! Index for Bmad internal use
  integer :: mat6_calc_method = bmad_standard$    ! taylor$, symp_lie_ptc$, etc.
  integer :: tracking_method = bmad_standard$     ! taylor$, linear$, etc.
  integer :: spin_tracking_method = tracking$     ! symp_lie_ptc$, etc.
  integer :: ptc_integration_type = matrix_kick$  ! drift_kick$, matrix_kick$, or ripken_kick$
  integer :: field_calc = bmad_standard$          ! no_field$, map$, grid$, refer_to_lords$, or custom$
  integer :: aperture_at = exit_end$              ! Aperture location: entrance_end$, ...
  integer :: aperture_type = rectangular$         ! rectangular$, elliptical$, auto_aperture$, ...
  integer :: orientation = 1                 ! -1 -> Element is longitudinally reversed. +1 -> Normal.
  logical :: symplectify = .false.           ! Symplectify mat6 matrices.
  logical :: mode_flip = .false.             ! Have the normal modes traded places?
  logical :: multipoles_on = .true.          ! For turning multipoles on/off
  logical :: scale_multipoles = .true.       ! Are ab_multipoles within other elements (EG: quads, etc.) 
                                             !        scaled by the strength of the element?
  logical :: taylor_map_includes_offsets = .true. ! Taylor map calculated with element misalignments?
  logical :: field_master = .false.          ! Calculate strength from the field value?
  logical :: is_on = .true.                  ! For turning element on/off.
  logical :: logic = .false.                 ! For general use. Not used by Bmad (except during lattice parsing).
  logical :: bmad_logic = .false.            ! For Bmad internal use only.
  logical :: select = .false.                ! For element selection. Used by make_hybrid_ring, etc.
  logical :: csr_calc_on = .true.            ! Coherent synchrotron radiation calculation
  logical :: offset_moves_aperture = .false. ! element offsets affects aperture?
end type

! struct for element to element control

type control_struct
  type (expression_atom_struct), allocatable :: stack(:) ! Evaluation stack
  type (lat_ele_loc_struct) :: slave = lat_ele_loc_struct()
  type (lat_ele_loc_struct) :: lord = lat_ele_loc_struct()
  integer :: ix_attrib = 0       ! Index of attribute controlled
end type

! lat_param_struct should be called branch_param_struct [Present name is an historical artifact.]
! Note that backwards_time_tracking is put in the lat_param_struct rather than begin a global
! for multithreaded applications.

integer, parameter :: incoherent$ = 1, coherent$ = 2
character(16), parameter :: photon_type_name(1:2) = ['Incoherent', 'Coherent  ']

type lat_param_struct
  real(rp) :: n_part = 0                       ! Particles/bunch (for BeamBeam elements).
  real(rp) :: total_length = 0                 ! total_length of branch. Warning: branch may not start at s = 0.
  real(rp) :: unstable_factor = 0              ! If positive: Growth rate/turn if unstable in closed branches or
                                               !   |orbit-aperture|/aperture if particle hits wall. Zero otherwise.
  real(rp) :: t1_with_RF(6,6) = 0              ! Full 1-turn matrix with RF on.
  real(rp) :: t1_no_RF(6,6) = 0                ! Full 1-turn matrix with RF off.
  real(rp) :: spin_tune = 0                    ! Closed orbit spin tune.
  integer :: particle = positron$              ! Reference particle: positron$, electron$, etc.
  integer :: default_tracking_species = ref_particle$  ! Default particle type to use in tracking.
  integer :: geometry = 0                      ! open$ or closed$
  integer :: ixx = 0                           ! Integer for general use
  logical :: stable = .false.                  ! is closed lat stable?
  logical :: live_branch = .true.              ! Should tracking be done on the branch?
  logical :: backwards_time_tracking = .false. ! Internal variable. Do not set.  
  type (bookkeeping_state_struct) :: bookkeeping_state = bookkeeping_state_struct()
                                               ! Overall status for the branch.
  type (beam_init_struct) :: beam_init = beam_init_struct() ! For beam initialization.
end type

! Structure for linking a branch_struct with a collection of ptc layouts

type ptc_layout_pointer_struct
  type (layout), pointer :: ptr => null()
end type

type ptc_branch1_info_struct
  type (ptc_layout_pointer_struct), allocatable :: m_u_layout(:)
  type (layout), pointer :: m_t_layout => null()
end type

!

type mode_info_struct
  real(rp) :: tune   = 0  ! "fractional" tune in radians: 0 < tune < 2pi
  real(rp) :: emit   = 0  ! Emittance.
  real(rp) :: chrom  = 0  ! Chromaticity.
  real(rp) :: sigma  = 0  ! Beam size.
  real(rp) :: sigmap = 0  ! Beam divergence.
end type

type normal_form_struct
  type (taylor_struct) :: M(6)             ! One-turn taylor map: M = A o N o A_inv, N = exp(:h:)
  type (taylor_struct) :: A(6)             ! Map from Floquet -> Lab coordinates
  type (taylor_struct) :: A_inv(6)         ! Map from Lab -> Floquet coordinates
  type (taylor_struct) :: dhdj(6)          ! Nonlinear tune function operating on Floquet coordinates
  type (complex_taylor_struct) :: F(6)     ! Vector field factorization in phasor basis:
  type (complex_taylor_struct) :: L(6)     ! M = A1 o c_inv o L exp(F.grad)Identity o c o A1_inv
                                           ! A1 and L are linear, and c maps to the phasor basis: h+ = x + i p, h- = x - i p
  type (ele_struct), pointer :: ele_origin => null()  ! Element at which the on-turn map was created.
                                           ! See subroutines: normal_form_taylors and normal_form_complex_taylors
end type

!

type branch_struct
  character(40) :: name = ''       ! Set to name of line that defines the branch.
  integer :: ix_branch = -1        ! Index of this branch. 0 => Main branch
  integer :: ix_from_branch = -1   ! -1 => No forking element to beginning of branch.
  integer :: ix_from_ele = -1      ! Forking element which forks to beginning of branch.
  integer, pointer :: n_ele_track => null()
  integer, pointer :: n_ele_max => null()
  type (lat_struct), pointer :: lat => null()
  type (mode_info_struct), pointer :: a => null(), b => null(), z => null() ! Note: Tunes are the fractional part.
  type (ele_struct), pointer :: ele(:) => null()
  type (lat_param_struct), pointer :: param => null()
  type (wall3d_struct), pointer :: wall3d(:) => null()
  type (ptc_branch1_info_struct) ptc
  type (normal_form_struct) normal_form_with_rf, normal_form_no_rf
end type

integer, parameter :: opal$ = 1, impactt$ = 2
character(16) :: pre_tracker_name(0:2) = ['NONE   ', 'OPAL   ', 'IMPACTT']

type pre_tracker_struct
  integer :: who = 0   ! Can be opal$, or impactt$
  integer :: ix_ele_start = 0
  integer :: ix_ele_end = 0
  character(200) :: input_file = ''
end type

! lat_struct
! Remember: If this struct is changed you have to modify:
!     Increase bmad_inc_version by 1.
!     read_digested_bmad_file
!     write_digested_bmad_file
!     transfer_lat_parameters
!     lat_equal_lat
! Rule: When lat2 = lat2, lat2%surface and lat1%surface will point to the same location.

type lat_struct
  character(40) :: use_name = ''                   ! Name of lat given by USE statement
  character(40) :: lattice = ''                    ! Lattice
  character(200) :: input_file_name = ''           ! Name of the lattice input file
  character(80) :: title = ''                      ! General title
  character(60), allocatable :: attribute_alias(:) ! Aliases for custom1$, etc.
  type (mode_info_struct) a, b, z                  ! Tunes, etc. Note: Tunes are the fractional part.
  type (lat_param_struct) param                    ! Parameters
  type (bookkeeping_state_struct) lord_state       ! lord bookkeeping status.
  type (ele_struct) ele_init                       ! For use by any program
  type (ele_struct), pointer ::  ele(:) => null()  ! Array of elements [=> branch(0)].
  type (branch_struct), allocatable :: branch(:)   ! Branch(0:) array
  type (control_struct), allocatable :: control(:) ! Control list
  type (photon_reflect_surface_struct), pointer :: surface(:) => null()
  type (coord_struct) beam_start                  ! Starting coords
  type (pre_tracker_struct) pre_tracker           ! For OPAL/IMPACT-T
  integer :: version = -1                         ! Version number
  integer :: n_ele_track                          ! Number of lat elements to track through.
  integer :: n_ele_max                            ! Index of last valid element in %ele(:) array
  integer :: n_control_max = 0                    ! Last index used in control_array
  integer :: n_ic_max = 0                         ! Last index used in ic_array
  integer :: input_taylor_order = 0               ! As set in the input file
  integer, allocatable :: ic(:)                   ! Index to %control(:)
  integer :: photon_type = incoherent$            ! Or coherent$. For X-ray simulations.
  logical :: absolute_time_tracking = .false.     ! Use abs. time for RF phase? Call autoscale if toggled.
  logical :: ptc_uses_hard_edge_drifts = .false.  ! Does associated ptc layout have hard edge model drifts?
end type

character(2), parameter :: coord_name(6) = ['X ', 'Px', 'Y ', 'Py', 'Z ', 'Pz']

! KEY value definitions
! Note: sbend$ and rbend$ also used for sub_key

integer, parameter :: drift$ = 1, sbend$ = 2, quadrupole$ = 3, group$ = 4
integer, parameter :: sextupole$ = 5, overlay$ = 6, custom$ = 7, taylor$ = 8
integer, parameter :: rfcavity$ = 9
integer, parameter :: elseparator$ = 10, beambeam$ = 11, wiggler$ = 12
integer, parameter :: sol_quad$ = 13, marker$ = 14, kicker$ = 15
integer, parameter :: hybrid$ = 16, octupole$ = 17, rbend$ = 18, multipole$ = 19
integer, parameter :: def_bmad_com$ = 20, def_mad_beam$ = 21, ab_multipole$ = 22, solenoid$ = 23
integer, parameter :: patch$ = 24, lcavity$ = 25, def_parameter$ = 26
integer, parameter :: null_ele$ = 27, beginning_ele$ = 28, line_ele$ = 29
integer, parameter :: match$ = 30, monitor$ = 31, instrument$ = 32
integer, parameter :: hkicker$ = 33, vkicker$ = 34, rcollimator$ = 35
integer, parameter :: ecollimator$ = 36, girder$ = 37, bend_sol_quad$ = 38
integer, parameter :: def_beam_start$ = 39, photon_fork$ = 40
integer, parameter :: fork$ = 41, mirror$ = 42, crystal$ = 43
integer, parameter :: pipe$ = 44, capillary$ = 45, multilayer_mirror$ = 46
integer, parameter :: e_gun$ = 47, em_field$ = 48, floor_shift$ = 49, fiducial$ = 50
integer, parameter :: undulator$ = 51, diffraction_plate$ = 52, photon_init$ = 53
integer, parameter :: sample$ = 54, detector$ = 55, sad_mult$ = 56, mask$ = 57, ac_kicker$ = 58

! "bend_sol_" is used to force the use of at least "bend_sol_q" in defining bend_sol_quad elements

integer, parameter :: n_key$ = 58
character(20), parameter :: key_name(n_key$) = [ &
    'Drift            ', 'Sbend            ', 'Quadrupole       ', 'Group            ', &
    'Sextupole        ', 'Overlay          ', 'Custom           ', 'Taylor           ', &
    'RFcavity         ', 'ELseparator      ', 'BeamBeam         ', 'Wiggler          ', &
    'Sol_Quad         ', 'Marker           ', 'Kicker           ', 'Hybrid           ', &
    'Octupole         ', 'Rbend            ', 'Multipole        ', 'Bend_sol_        ', &
    'Def_Mad_Beam     ', 'AB_multipole     ', 'Solenoid         ', 'Patch            ', &
    'Lcavity          ', 'Def_Parameter    ', 'Null_Ele         ', 'Beginning_Ele    ', &
    'Line_Ele         ', 'Match            ', 'Monitor          ', 'Instrument       ', &
    'Hkicker          ', 'Vkicker          ', 'Rcollimator      ', 'Ecollimator      ', &
    'Girder           ', 'Bend_Sol_Quad    ', 'Def_Beam_Start   ', 'Photon_Fork      ', &
    'Fork             ', 'Mirror           ', 'Crystal          ', 'Pipe             ', &
    'Capillary        ', 'Multilayer_Mirror', 'E_Gun            ', 'EM_Field         ', &
    'Floor_Shift      ', 'Fiducial         ', 'Undulator        ', 'Diffraction_Plate', &
    'Photon_Init      ', 'Sample           ', 'Detector         ', 'Sad_Mult         ', &
    'Mask             ', 'AC_Kicker        ']

! These logical arrays get set in init_attribute_name_array and are used
! to sort elements that have kick or orientation attributes from elements that do not.
! The orientation attributes are: tilt, x/y/z_offset, x/y_pitch, and *_tot versions.
! Note: A solenoid does not formally have a tilt but has everything else.
! Rule: Any element that formally has some but not all orientation attributes is considered
!   internally to have all attributes so any routine can safely work with all the 
!   orientation attributes as a block.

logical has_hkick_attributes(n_key$)
logical has_kick_attributes(n_key$)

! Element attribute name logical definitions

integer, parameter :: n_part$ = 2, taylor_order$ = 3

integer, parameter :: val1$=11, val2$=12, val3$=13, val4$=14, val5$=15, &
          val6$=16, val7$=17, val8$=18, val9$=19, val10$=20, val11$=21, &
          val12$=22

integer, parameter :: beta_a0$ = 2, alpha_a0$ = 3, beta_b0$ = 4
integer, parameter :: alpha_b0$ = 5, beta_a1$ = 6, alpha_a1$ = 7, beta_b1$ = 8
integer, parameter :: alpha_b1$ = 9, dphi_a$ = 10, dphi_b$ = 11
integer, parameter :: eta_x0$ = 12, etap_x0$ = 13, eta_y0$ = 14, etap_y0$ = 15
integer, parameter :: eta_x1$ = 16, etap_x1$ = 17, eta_y1$ = 18, etap_y1$ = 19
integer, parameter :: match_end_input$ = 20
integer, parameter :: match_end$ = 21
integer, parameter :: delta_time$ = 22
integer, parameter :: x0$ = 24, px0$ = 25, y0$ = 26, py0$ = 27, z0$ = 28, pz0$ = 29
integer, parameter :: x1$ = 30, px1$ = 31, y1$ = 32, py1$ = 33, z1$ = 34, pz1$ = 35
integer, parameter :: match_end_orbit_input$ = 36, match_end_orbit$ = 37

integer, parameter :: c_11$ = 40, c_12$ = 41, c_21$ = 42, c_22$ = 43, gamma_c$ = 44

integer, parameter :: x$ = 1, px$ = 2, y$ = 3, py$ = 4, z$ = 5, pz$ = 6
integer, parameter :: t$ = 8
integer, parameter :: field_x$ = 10, field_y$ = 11, phase_x$ = 12, phase_y$ = 13
integer, parameter :: e_photon$ = 9

integer, parameter :: e1$ = 19, e2$ = 20
integer, parameter :: fint$ = 21, fintx$ = 22, hgap$ = 23, hgapx$ = 24, h1$ = 25, h2$ = 26

integer, parameter :: l$ = 1                          ! Assumed unique. Do not assign 1 to another attribute.
integer, parameter :: tilt$ = 2, roll$ = 2  ! Important: tilt$ = roll$
integer, parameter :: ref_tilt$ = 3, rf_frequency$ = 3, direction$ = 3
integer, parameter :: kick$ = 3, x_gain_err$ = 3
integer, parameter :: rf_frequency_err$ = 4, k1$ = 4, harmon$ = 4, h_displace$ = 4, y_gain_err$ = 4
integer, parameter :: critical_angle_factor$ = 4, tilt_corr$ = 4, ref_coordinates$ = 4
integer, parameter :: graze_angle$ = 5, k2$ = 5, b_max$ = 5, v_displace$ = 5, drift_id$ = 5
integer, parameter :: ks$ = 5, flexible$ = 5, crunch$ = 5, ref_orbit_follows$ = 5
integer, parameter :: gradient$ = 6, k3$ = 6, noise$ = 6, new_branch$ = 6
integer, parameter :: g$ = 6, bragg_angle_in$ = 6, symmetry$ = 6, field_scale_factor$ = 6
integer, parameter :: g_err$ = 7, n_pole$ = 7, bbi_const$ = 7, osc_amplitude$ = 7
integer, parameter :: gradient_err$ = 7, critical_angle$ = 7, sad_flag$ = 7
integer, parameter :: bragg_angle_out$ = 7, ix_to_branch$ = 7
integer, parameter :: rho$ = 8, delta_e$ = 8, diffraction_limited$ = 8
integer, parameter :: charge$ = 8, x_gain_calib$ = 8, ix_to_element$ = 8
integer, parameter :: voltage$ = 9, eps_step_scale$ = 9
integer, parameter :: fringe_type$ = 10
integer, parameter :: fringe_at$ = 11, gang$ = 11
integer, parameter :: higher_order_fringe_type$ = 12
integer, parameter :: spin_fringe_on$ = 13
integer, parameter :: fb1$ = 14, sig_x$ = 14, exact_multipoles$ = 14
integer, parameter :: fb2$ = 15, sig_y$ = 15
integer, parameter :: fq1$ = 16, sig_z$ = 16
integer, parameter :: fq2$ = 17, sig_vx$ = 17
integer, parameter :: sig_vy$ = 18, autoscale_amplitude$ = 18
integer, parameter :: sig_e$ = 19, autoscale_phase$ = 19
integer, parameter :: d1_thickness$ = 20, voltage_err$ = 20, default_tracking_species$ = 20, direction_beam_start$ = 20
integer, parameter :: n_slice$ = 20, y_gain_calib$ = 20, bragg_angle$ = 20, E_center$ = 20, constant_ref_energy$ = 20
integer, parameter :: polarity$ = 21, crunch_calib$ = 21, alpha_angle$ = 21, d2_thickness$ = 21
integer, parameter :: e_loss$ = 21, dks_ds$ = 21, gap$ = 21, E_center_relative_to_ref$ = 21, spin_x$ = 21
integer, parameter :: x_offset_calib$ = 22, v1_unitcell$ = 22, psi_angle$ = 22, spatial_distribution$ = 22
integer, parameter :: spin_y$ = 22
integer, parameter :: y_offset_calib$ = 23, v_unitcell$ = 23, v2_unitcell$ = 23, spin_z$ = 23
integer, parameter :: cavity_type$ = 23, beta_a$ = 23, velocity_distribution$ = 23
integer, parameter :: phi0$ = 24, tilt_calib$ = 24, beta_b$ = 24, energy_distribution$ = 24, live_branch$ = 24
integer, parameter :: phi0_err$ = 25, current$ = 25, l_pole$ = 25, particle$ = 25
integer, parameter :: quad_tilt$ = 25, de_eta_meas$ = 25, alpha_a$ = 25, e_field_x$ = 25
integer, parameter :: geometry$ = 26, bend_tilt$ = 26, mode$ = 26, alpha_b$ = 26, e_field_y$ = 26
integer, parameter :: phi0_multipass$ = 26, n_sample$ = 26, origin_ele_ref_pt$ = 26
integer, parameter :: phi0_autoscale$ = 27, dx_origin$ =  27, cmat_11$ = 27, scale_field_to_one$ = 27
integer, parameter :: lattice_type$ = 27, x_quad$ = 27, ds_photon_slice$ = 27
integer, parameter :: phi0_max$ = 28, dy_origin$ = 28, y_quad$ = 28, photon_type$ = 28
integer, parameter :: cmat_12$ = 28
integer, parameter :: floor_set$ = 29, upstream_ele_dir$ = 29, dz_origin$ = 29
integer, parameter :: cmat_21$ = 29, l_sagitta$ = 29
integer, parameter :: dtheta_origin$ = 30, b_param$ = 30, transverse_sigma_cut$ = 30, l_chord$ = 30
integer, parameter :: downstream_ele_dir$ = 30, cmat_22$ = 30, spinor_theta$ = 30
integer, parameter :: l_hard_edge$ = 31, dphi_origin$ = 31, ref_cap_gamma$ = 31, ds_slice$ = 31, spinor_phi$ = 31
integer, parameter :: field_autoscale$ = 32, dpsi_origin$ = 32, t_offset$ = 32, spinor_xi$ = 32
integer, parameter :: angle$ = 33, n_cell$ = 33, x_ray_line_len$ = 33, spinor_polarization$ = 33
integer, parameter :: x_pitch$ = 34
integer, parameter :: y_pitch$ = 35  
integer, parameter :: x_offset$ = 36
integer, parameter :: y_offset$ = 37 
integer, parameter :: z_offset$ = 38 ! Assumed unique. Do not overload further.
integer, parameter :: hkick$ = 39, d_spacing$ = 39, x_offset_mult$ = 39, emittance_a$ = 39
integer, parameter :: vkick$ = 40, y_offset_mult$ = 40, p0c_ref_init$ = 40, emittance_b$ = 40
integer, parameter :: BL_hkick$ = 41, x_pitch_mult$ = 41, e_tot_ref_init$ = 41, emittance_z$ = 41
integer, parameter :: BL_vkick$ = 42, y_pitch_mult$ = 42, darwin_width_sigma$ = 42
integer, parameter :: pendellosung_period_sigma$ = 43, BL_kick$ = 43, B_field$ = 43, E_field$ = 43
integer, parameter :: coupler_phase$ = 44, darwin_width_pi$ = 44, B_field_err$ = 44
integer, parameter :: B1_gradient$ = 45, E1_gradient$ = 45, coupler_angle$ = 45, pendellosung_period_pi$ = 45
integer, parameter :: B2_gradient$ = 46, E2_gradient$ = 46, coupler_strength$ = 46, dbragg_angle_de$ = 46
integer, parameter :: coupler_at$ = 47, E_tot_set$ = 47, ptc_canonical_coords$ = 47
integer, parameter :: B3_gradient$ = 48, E3_gradient$ = 48, ptc_fringe_geometry$ = 48, p0c_set$ = 48
integer, parameter :: Bs_field$ = 49, e_tot_offset$ = 49, ptc_field_geometry$ = 49
integer, parameter :: delta_ref_time$ = 50 ! Assumed unique Do not overload.
integer, parameter :: p0c_start$ = 51      ! Assumed unique. Do not overload.
integer, parameter :: e_tot_start$ = 52    ! Assumed unique. Do not overload.
integer, parameter :: p0c$ = 53            ! Assumed unique. Do not overload.
integer, parameter :: e_tot$ = 54          ! Assumed unique. Do not overload.
integer, parameter :: x_pitch_tot$ = 55, no_end_marker$ = 55
integer, parameter :: y_pitch_tot$ = 56
integer, parameter :: x_offset_tot$ = 57
integer, parameter :: y_offset_tot$ = 58
integer, parameter :: z_offset_tot$ = 59
integer, parameter :: tilt_tot$ = 60, roll_tot$ = 60  ! Important: tilt_tot$ = roll_tot$
integer, parameter :: ref_tilt_tot$ = 61
integer, parameter :: n_ref_pass$ = 62
integer, parameter :: r0_mag$ = 63
integer, parameter :: ref_time_start$ = 64
integer, parameter :: thickness$ = 65, integrator_order$ = 65   ! For Etiennes' PTC: 2, 4, or 6.
integer, parameter :: num_steps$ = 66
integer, parameter :: ds_step$ = 67
integer, parameter :: lord_pad1$ = 68
integer, parameter :: lord_pad2$ = 69, ref_wavelength$ = 69
integer, parameter :: r0_elec$ = 70
integer, parameter :: custom_attribute1$ = 71   ! Assumed unique. Do not overload.
integer, parameter :: custom_attribute2$ = 72   ! Assumed unique. Do not overload.
integer, parameter :: custom_attribute3$ = 73   ! Assumed unique. Do not overload.
integer, parameter :: custom_attribute4$ = 74   ! Assumed unique. Do not overload.
integer, parameter :: custom_attribute5$ = 75, custom_attribute_max$ = 75   ! For general use
integer, parameter :: x1_limit$ = 76   ! Assumed unique. Do not overload.
integer, parameter :: x2_limit$ = 77   ! Assumed unique. Do not overload.
integer, parameter :: y1_limit$ = 78   ! Assumed unique. Do not overload.
integer, parameter :: y2_limit$ = 79   ! Assumed unique. Do not overload.
integer, parameter :: check_sum$ = 80  ! Assumed unique. Do not overload.

!! 81 = 1 + num_ele_attrib$

integer, parameter :: max_aperture_limit$ = 81     ! bmad_com parameters
integer, parameter :: default_ds_step$ = 82
integer, parameter :: significant_length$ = 83
integer, parameter :: rel_tol_tracking$ = 84
integer, parameter :: abs_tol_tracking$ = 85
integer, parameter :: rel_tol_adaptive_tracking$ = 86
integer, parameter :: abs_tol_adaptive_tracking$ = 87
integer, parameter :: init_ds_adaptive_tracking$ = 88
integer, parameter :: min_ds_adaptive_tracking$ = 89
integer, parameter :: fatal_ds_adaptive_tracking$ = 90
integer, parameter :: max_num_runge_kutta_step$ = 91

integer, parameter :: alpha_b_begin$ = 81, use_hard_edge_drifts$ = 81, tt$ = 81, lr_wake_spline$ = 81
integer, parameter :: alias$  = 82, eta_x$ = 82, ptc_max_fringe_order$ = 82
integer, parameter :: eta_y$ = 83, electric_dipole_moment$ = 83, lr_self_wake_on$ = 83, x_ref$ = 83
integer, parameter :: etap_x$ = 84, lr_wake_file$ = 84, px_ref$ = 84
integer, parameter :: etap_y$ = 85, lr_freq_spread$ = 85, y_ref$ = 85
integer, parameter :: lattice$ = 86, phi_a$ = 86, multipoles_on$ = 86, py_ref$ = 86
integer, parameter :: aperture_type$ = 87, eta_z$ = 87
integer, parameter :: taylor_map_includes_offsets$ = 88, cmat_11_begin$ = 88, surface_attrib$ = 88
integer, parameter :: csr_calc_on$ = 89, cmat_12_begin$ = 89, var$ = 89, z_ref$ = 89
integer, parameter :: cmat_21_begin$ = 90, pz_ref$ = 90
integer, parameter :: mat6_calc_method$ = 91, cmat_22_begin$ = 91
integer, parameter :: tracking_method$  = 92, s_long$ = 92
integer, parameter :: ref_time$ = 93, ptc_integration_type$ = 93
integer, parameter :: spin_tracking_method$ = 94, eta_a$ = 94
integer, parameter :: aperture$ = 95, etap_a$ = 95
integer, parameter :: x_limit$ = 96, absolute_time_tracking$ = 96, eta_b$ = 96
integer, parameter :: y_limit$ = 97, etap_b$ = 97
integer, parameter :: offset_moves_aperture$ = 98
integer, parameter :: aperture_limit_on$ = 99

integer, parameter :: ptc_exact_misalign$ = 100, physical_source$ = 100
integer, parameter :: sr_wake_file$ = 100, alpha_a_begin$ = 100
integer, parameter :: term$ = 101, frequencies$ = 101
integer, parameter :: x_position$ = 102, s_spline$ = 102, ptc_exact_model$ = 102
integer, parameter :: symplectify$ = 103, y_position$ = 103, n_slice_spline$ = 103
integer, parameter :: z_position$ = 104, amp_vs_time$ = 104
integer, parameter :: is_on$ = 105, theta_position$ = 105
integer, parameter :: field_calc$ = 106, phi_position$ = 106
integer, parameter :: psi_position$ = 107, wall$ = 107
integer, parameter :: aperture_at$ = 108, beta_a_begin$ = 108
integer, parameter :: ran_seed$ = 109, beta_b_begin$ = 109, origin_ele$ = 109

! Make sure that lat%control(i)%ix_attrib /= field_overlaps$ for all girders, superposition, overlays, and groups.
integer, parameter :: to_line$ = 110, field_overlaps$ = 110 
integer, parameter :: field_master$ = 111, harmon_master$ = 111, to_element$ = 111
integer, parameter :: descrip$ = 112
integer, parameter :: scale_multipoles$ = 113
integer, parameter :: ref_orbit$ = 115
integer, parameter :: phi_b$ = 116, crystal_type$ = 116, material_type$ = 116
integer, parameter :: type$ = 117
integer, parameter :: ref_origin$ = 118
integer, parameter :: ele_origin$ = 119

! superimpose$ through create_jumbo_slave$ assumed unique (or need to modify bmad_parser_mod.f90).

integer, parameter :: superimpose$     = 120   
integer, parameter :: offset$          = 121
integer, parameter :: reference$       = 122
integer, parameter :: cartesian_map$   = 123
integer, parameter :: cylindrical_map$ = 124
integer, parameter :: grid_field$      = 125
integer, parameter :: taylor_field$    = 126
integer, parameter :: create_jumbo_slave$ = 127

integer, parameter :: accordion_edge$  = 128
integer, parameter :: start_edge$  = 129
integer, parameter :: end_edge$  = 130
integer, parameter :: s_position$ = 131

integer, parameter :: a0$  = 140, a21$  = 161
integer, parameter :: b0$  = 162, b21$  = 183

integer, parameter :: k0l$ = 140, k21l$ = 161
integer, parameter :: t0$  = 162, t21$  = 183

integer, parameter :: a0_elec$ = 190, a21_elec$ = 211
integer, parameter :: b0_elec$ = 212, b21_elec$ = 233

integer, parameter :: num_ele_attrib_extended$ = b21_elec$

character(40), parameter :: blank_name$ = ' '

! lattice logical names

integer, parameter :: open$ = 1, closed$ = 2

character(16), parameter :: lattice_type_name(0:2) = ['GARBAGE!        ', 'LINEAR_LATTICE  ', 'CIRCULAR_LATTICE']
character(16), parameter :: geometry_name(0:2) = ['GARBAGE!    ', 'Open        ', 'Closed      ']

! The linac_normal_mode_struct is basically the synchrotron integrals with the
! energy factors thrown in. Useful for linacs.

type anormal_mode_struct
  real(rp) :: emittance = 0       ! Beam emittance
  real(rp) :: synch_int(4:6) = 0  ! Synchrotron integrals
  real(rp) :: j_damp = 0          ! damping partition number
  real(rp) :: alpha_damp = 0      ! damping per turn
  real(rp) :: chrom = 0           ! Chromaticity
  real(rp) :: tune = 0            ! "Fractional" tune in radians
end type

type linac_normal_mode_struct
  real(rp) :: i2_E4 = 0           ! Integral: g^2 * gamma^4
  real(rp) :: i3_E7 = 0           ! Integral: g^3 * gamma^7
  real(rp) :: i5a_E6 = 0          ! Integral: (g^3 * H_a) * gamma^6
  real(rp) :: i5b_E6 = 0          ! Integral: (g^3 * H_b) * gamma^6
  real(rp) :: sig_E1 = 0          ! Energy spread after 1 pass (eV)
  real(rp) :: a_emittance_end = 0 ! a mode emittance at end of linac
  real(rp) :: b_emittance_end = 0 ! b mode emittance at end of linac
end type

type normal_modes_struct
  real(rp) :: synch_int(0:3) = 0  ! Synchrotron integrals I0, I1, I2, and I3
  real(rp) :: sigE_E = 0          ! SigmaE/E
  real(rp) :: sig_z = 0           ! Sigma_Z
  real(rp) :: e_loss = 0          ! Energy loss / turn (eV)
  real(rp) :: rf_voltage = 0      ! Total rfcavity voltage (eV)
  real(rp) :: pz_aperture = 0     ! pz aperture limit
  type (anormal_mode_struct) :: a = anormal_mode_struct()
  type (anormal_mode_struct) :: b = anormal_mode_struct()
  type (anormal_mode_struct) :: z = anormal_mode_struct()
  type (linac_normal_mode_struct) :: lin = linac_normal_mode_struct()
end type

integer, parameter :: bends$ = 201
integer, parameter :: wigglers$ = 202
integer, parameter :: all$ = 203

!---------------------------------------------------------------------------
! Units

integer, parameter :: radians$ = 1, degrees$ = 2, cycles$ = 3, kHz$ = 4
character(8), parameter :: frequency_units_name(4) = ['Radians ', 'Degrees ', 'Cycles  ', 'kHz     ']

! Electric and magnetic fields.

type em_field_struct
  real(rp) :: E(3) = 0        ! electric field
  real(rp) :: B(3) = 0        ! magnetic field
  real(rp) :: dE(3,3) = 0     ! electric field gradient
  real(rp) :: dB(3,3) = 0     ! magnetic field gradient
end type

type em_potential_struct
  real(rp) :: phi = 0         ! Electric scalar potential.
  real(rp) :: phi_B = 0       ! Magnetic scalar potential.
  real(rp) :: A(3) = 0        ! Magnetic vector potential.
end type

! Grid of grid_field information
integer, parameter :: rotationally_symmetric_rz$ = 1, xyz$ = 2
character(30), parameter :: grid_field_geometry_name(0:2) = &
                          [character(30) :: 'GARBAGE!', 'rotationally_symmetric_rz', 'xyz']
integer, parameter :: grid_field_dimension(2) = [2, 3] 


! Structure for saving the track through an element.

type track_map_struct
  real(rp) vec0(6)        ! 0th order part of xfer map from the beginning.
  real(rp) mat6(6,6)      ! 1st order part of xfer map (transfer matrix).
end type

! Valid track%orb(:) points in range 0:track%n_pt

type track_struct
  type (coord_struct), allocatable :: orb(:)      ! An array of track points indexed from 0 (%orb(0:)).
  type (em_field_struct), allocatable:: field(:)  ! An array of em fields indexed from 0 (%field(0:)).
  type (track_map_struct), allocatable :: map(:)  ! An array of linear maps indexed from 0.
  real(rp) :: ds_save = 1d-3                      ! Min distance between points.
  integer :: n_pt = -1                            ! Track upper bound for %orb(0:), etc. arrays.
  integer :: n_bad = 0                            ! Number of bad steps when adaptive tracking is done.
  integer :: n_ok = 0                             ! Number of good steps when adaptive tracking is done.
end type

!------------------------------------------------------------------------------
! Dynamic aperture structures

type aperture_data_struct
  real(rp) x, y     ! (x,y) aperture point
  integer plane     ! plane determining loss
  integer ix_lat    ! ele index lost at
  integer i_turn    ! turn lost at
end type

type aperture_param_struct
  real(rp) :: min_angle = 0
  real(rp) :: max_angle = pi
  integer :: n_angle   = 9
  integer :: n_turn = 100         ! Number of turns a particle must survive
  real(rp) :: x_init = 1e-3_rp    ! Initial x coordinate to start with for theta_xy = 0.
  real(rp) :: y_init = 1e-3_rp    ! Initial y coordinate to start with for theta_xy = pi/2.
  real(rp) :: accuracy = 1e-5_rp  ! Resolution of bracketed aperture (meters.)
end type

type aperture_scan_struct
  type(aperture_data_struct), allocatable :: aperture(:) ! set of apertures at different angles
  type(aperture_param_struct) :: param                   ! parameters used for the scan            
  type (coord_struct) :: ref_orb                         ! Ref orbit around which the scan is made.
end type

!-------------------------------------------------------------------------
! CSR parameters

!+
! Note: Shielding is simulated via the image current due to the 
! top and bottom walls. The side walls are neglected.
!-

type csr_parameter_struct                  ! Common block for csr calc
  real(rp) :: ds_track_step = 0            ! Tracking step size
  real(rp) :: beam_chamber_height = 0      ! Used in shielding calculation.
  real(rp) :: sigma_cutoff = 0.1           ! Cutoff for the lsc calc. If a bin sigma
                                           !  is < cutoff * sigma_ave then ignore.
  integer :: n_bin = 0                     ! Number of bins used
  integer :: particle_bin_span = 2         ! Longitudinal particle length / dz_bin
  integer :: n_shield_images = 0           ! Chamber wall shielding. 0 = no shielding.
  integer :: ix1_ele_csr = -1              ! Start index for csr tracking
  integer :: ix2_ele_csr = -1              ! Stop index for csr tracking
  integer :: sc_min_in_bin = 10            ! Minimum number of particles in a bin for sigmas to be valid.
  logical :: lcsr_component_on = .true.    ! Longitudinal csr component
  logical :: lsc_component_on = .true.     ! Longitudinal space charge component
  logical :: tsc_component_on = .false.    ! Transverse space charge component
  logical :: lsc_kick_transverse_dependence = .false.
  logical :: print_taylor_warning = .true. ! Print warning if Taylor element is present?
  logical :: use_csr_old = .false.         ! Use old CSR tracking? Should only be done for testing.
  logical :: small_angle_approx = .true.   ! Use small angle approximation? ONLY USED WITH OLD CSR.
  logical :: write_csr_wake = .false.      ! Write the CSR wake to csr_wake.dat / old_csr_wake.dat?
end type

type (csr_parameter_struct), save, target :: csr_param

!------------------------------------------------------------------------------
! This is for debugging radiation damping and fluctuations.

type synch_rad_common_struct
  real(rp) :: scale = 1.0               ! used to scale the radiation
  real(rp) :: i2 = 0, i3 = 0            ! radiation integrals
  real(rp) :: i5a = 0, i5b = 0
  logical :: i_calc_on = .false.        ! For calculating i2 and i3    
end type

type (synch_rad_common_struct), save :: synch_rad_com

integer, parameter :: unknown$ = 0, is_logical$ = 1, is_integer$ = 2, is_real$ = 3, is_switch$ = 4, is_string$ = 5
integer, parameter :: is_struct$ = 6 

! For coords_floor_to_curvilinear status argument

integer, parameter :: patch_problem$ = 2, outside$ = 3, cannot_find$ = 4

! extra_parsing_info_struct is used by parsing routines.
! %deterministic:
!   0 = Not, 1 = Ran state on input deterministic, 2 = ran state deterministice
!   will be generated in exactly the same way every time?
! %ran_function_was_called:
!   Only set True when ran function is called with ran_determinisitc = 0.
! %determinisitc_ran_function_was_called:
!   Only set True when ran function is called with ran_determinisitc = 1.

type extra_parsing_info_struct
  type (random_state_struct) :: initial_state       = random_state_struct()
  integer :: deterministic                          = 0
  logical :: ran_function_was_called                = .false.
  logical :: deterministic_ran_function_was_called  = .false.
  logical :: d_orb_set                              = .false.
  logical :: max_aperture_limit_set                 = .false.
  logical :: default_ds_step_set                    = .false.
  logical :: significant_length_set                 = .false.
  logical :: rel_tol_tracking_set                   = .false.
  logical :: abs_tol_tracking_set                   = .false.
  logical :: rel_tol_adaptive_tracking_set          = .false.
  logical :: abs_tol_adaptive_tracking_set          = .false.
  logical :: init_ds_adaptive_tracking_set          = .false.
  logical :: min_ds_adaptive_tracking_set           = .false.
  logical :: fatal_ds_adaptive_tracking_set         = .false.
  logical :: electric_dipole_moment_set             = .false.
  logical :: ptc_cut_factor_set                     = .false.
  logical :: taylor_order_set                       = .false.
  logical :: runge_kutta_order_set                  = .false.
  logical :: default_integ_order_set                = .false.
  logical :: ptc_max_fringe_order_set               = .false.
  logical :: use_hard_edge_drifts_set               = .false.
  logical :: sr_wakes_on_set                        = .false.
  logical :: lr_wakes_on_set                        = .false.
  logical :: mat6_track_symmetric_set               = .false.
  logical :: auto_bookkeeper_set                    = .false.
  logical :: space_charge_on_set                    = .false.
  logical :: coherent_synch_rad_on_set              = .false.
  logical :: spin_tracking_on_set                   = .false.
  logical :: radiation_damping_on_set               = .false.
  logical :: radiation_fluctuations_on_set          = .false.
  logical :: conserve_taylor_maps_set               = .false.
  logical :: absolute_time_tracking_default_set     = .false.
  logical :: convert_to_kinetic_momentum_set        = .false.
  logical :: aperture_limit_on_set                  = .false.
  logical :: debug_set                              = .false.
  logical :: sad_eps_scale_set                      = .false.
  logical :: sad_amp_max_set                        = .false.
  logical :: sad_n_div_max_set                      = .false.
  logical :: max_num_runge_kutta_step_set           = .false.
end type

! ptc_field_geometry for bends

integer, parameter :: sector$ = 1, straight$ = 2, true_rbend$ = 3
character(16), parameter :: ptc_field_geometry_name(0:3) = [ &
                              'Garbage!  ', 'Sector    ', 'Straight  ', 'True_Rbend']

!------------------------------------------------------------------------------
! common stuff

! %max_aperture_limit is used when no limit is specified or when 
!   lat%param%aperture_limit_on = False.
! Remember: Change extra_parsing_info_struct if bmad_common_struct changed.

type bmad_common_struct
  real(rp) :: max_aperture_limit = 1d3            ! Max Aperture.
  real(rp) :: d_orb(6)           = 1d-5           ! Orbit deltas for the mat6 via tracking calc.
  real(rp) :: default_ds_step    = 0.2_rp         ! Integration step size.  
  real(rp) :: significant_length = 1d-10          ! meter 
  real(rp) :: rel_tol_tracking = 1d-8                 ! Closed orbit relative tolerance.
  real(rp) :: abs_tol_tracking = 1d-10                ! Closed orbit absolute tolerance.
  real(rp) :: rel_tol_adaptive_tracking = 1d-8        ! Runge-Kutta tracking relative tolerance.
  real(rp) :: abs_tol_adaptive_tracking = 1d-10       ! Runge-Kutta tracking absolute tolerance.
  real(rp) :: init_ds_adaptive_tracking = 1d-3        ! Initial step size
  real(rp) :: min_ds_adaptive_tracking = 0            ! Min step size to take.
  real(rp) :: fatal_ds_adaptive_tracking = 1d-8       ! If actual step size is below this particle is lost.
  real(rp) :: electric_dipole_moment = 0              ! Particle's EDM. Call set_ptc to transfer value to PTC.
  real(rp) :: ptc_cut_factor = 0.006                  ! Cut factor for PTC tracking
  real(rp) :: sad_eps_scale = 5.0d-3                  ! Used in sad_mult step length calc.
  real(rp) :: sad_amp_max = 5.0d-2                    ! Used in sad_mult step length calc.
  integer :: sad_n_div_max = 1000                     ! Used in sad_mult step length calc.
  integer :: taylor_order = 0                         ! Input Taylor order for maps. 
                                                      !   0 -> default = ptc%taylor_order_saved
                                                      !   ptc_com%taylor_order_ptc gives actual order in use. 
  integer :: runge_kutta_order = 4                    ! Runge Kutta order.
  integer :: default_integ_order = 2                  ! PTC integration order. 
  integer :: ptc_max_fringe_order = 2                 ! PTC max fringe order (2  = > Quadrupole !). 
                                                      !   Must call set_ptc after changing.
  integer :: max_num_runge_kutta_step = 10000         ! Maximum number of RK steps before particle is considered lost.
  logical :: use_hard_edge_drifts = .true.            ! Insert drifts when tracking through cavity?
  logical :: sr_wakes_on = .true.                     ! Short range wakefields?
  logical :: lr_wakes_on = .true.                     ! Long range wakefields
  logical :: mat6_track_symmetric = .true.            ! symmetric offsets
  logical :: auto_bookkeeper = .true.                 ! Automatic bookkeeping?
  logical :: space_charge_on = .false.                ! Space charge switch
  logical :: coherent_synch_rad_on = .false.          ! csr 
  logical :: spin_tracking_on = .false.               ! spin tracking?
  logical :: radiation_damping_on = .false.           ! Damping toggle.
  logical :: radiation_fluctuations_on = .false.      ! Fluctuations toggle.
  logical :: conserve_taylor_maps = .true.            ! Enable bookkeeper to set ele%taylor_map_includes_offsets = F?
  logical :: absolute_time_tracking_default = .false. ! Default for lat%absolute_time_tracking
  logical :: convert_to_kinetic_momentum = .false.    ! Cancel kicks due to finite vector potential when doing symplectic tracking?
                                                      !   Set to True to test symp_lie_bmad against runge_kutta.
  logical :: aperture_limit_on = .true.               ! use apertures in tracking?
  logical :: debug = .false.                          ! Used for code debugging.
end type
  
type (bmad_common_struct), save, target :: bmad_com

! ptc_com common block.
! %taylor_order_saved is what is used if bmad_com%taylor_order is not set ( =  0).
! %taylor_order_saved is initially 3.
! When parsing a lattice file, %taylor_order_saved will be set to the taylor order of the lattice.

type ptc_common_struct
  integer :: taylor_order_ptc = 0       ! What has been set in PTC. 0 -> not yet set
  integer :: taylor_order_saved = 3     ! Default to use.
  logical :: complex_ptc_used = .false. ! Complex PTC code in use? (EG for spin tracking, normal form anal, etc.)
  logical :: use_totalpath = .false.    ! phase space z = time instead of time - ref_time? 
end type

type (ptc_common_struct), save :: ptc_com

real(rp), parameter :: small_rel_change$ = 1d-14

! Some routines need to keep track of where elements are when elements are added or removed from
! the lattice. 

type ele_loc_ele_struct
  type (lat_ele_loc_struct), allocatable :: ele(:)
end type

type ele_loc_branch_struct
  type (ele_loc_ele_struct), allocatable :: branch(:)
end type

type (ele_loc_branch_struct), save, target :: ele_loc_com

! This structure stores the radiation integrals for an individual element except
! lin_norm_emit_a and lin_norm_emit_b are running sums from the beginning of the branch.
! See also rad_int_all_ele_struct.

type rad_int1_struct
  real(rp) :: i0   = 0
  real(rp) :: i1   = 0
  real(rp) :: i2   = 0
  real(rp) :: i3   = 0
  real(rp) :: i4a  = 0
  real(rp) :: i4b  = 0
  real(rp) :: i4z  = 0
  real(rp) :: i5a  = 0
  real(rp) :: i5b  = 0
  real(rp) :: i6b  = 0
  real(rp) :: lin_i2_E4   = 0 
  real(rp) :: lin_i3_E7   = 0
  real(rp) :: lin_i5a_E6  = 0
  real(rp) :: lin_i5b_E6  = 0
  real(rp) :: lin_norm_emit_a = 0 ! Running sum
  real(rp) :: lin_norm_emit_b = 0 ! Running sum
  real(rp) :: n_steps         = 0 ! number of qromb steps needed
end type

! Structure for an array of rad_int1_structs for a single branch

type rad_int_all_ele_struct
  type (rad_int1_struct), allocatable :: ele(:) ! Array is indexed from 0
end type

contains

!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!+
! Function coord_state_name (coord_state) result (state_str)
!
! Routine to return the string representation of a coord%state state.
!
! Input:
!   coord_state -- integer: coord%state value
!
! Output:
!   state_str   -- character(20): String representation.
!-

function coord_state_name (coord_state) result (state_str)

implicit none

integer coord_state
character(20) state_str

!

select case (coord_state)
case (alive$);                 state_str = 'Alive'
case (lost$);                  state_str = 'Lost'
case (not_set$);               state_str = 'Not_Set'
case (lost_neg_x_aperture$);   state_str = 'Hit -X Side'
case (lost_pos_x_aperture$);   state_str = 'Hit +X Side'
case (lost_neg_y_aperture$);   state_str = 'Hit -Y Side'
case (lost_pos_y_aperture$);   state_str = 'Hit +Y Side'
case (lost_pz_aperture$);       state_str = 'Hit Z Side'
case default;                  state_str = 'UNKNOWN!'
end select

end function coord_state_name

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function is_attribute (ix_attrib, which) result (is_attrib)
!
! Routine to determine if an attribute index corresponds to a control variable for overlys/groups.
!
! Input:
!   ix_attrib -- integer: Attribute index.
!   which     -- integer: control_var$, old_control_var$, all_control_var$, multipole$, elec_multipole$
!
! Output:
!   is_attrib -- logical: True if a control variable
!-

function is_attribute (ix_attrib, which) result (is_attrib)

integer ix_attrib, which
logical is_attrib

!

select case (which)
case (control_var$)
  is_attrib = (ix_attrib > var_offset$ .and. ix_attrib < var_offset$+20)
  
case (old_control_var$)
  is_attrib = (ix_attrib > old_control_var_offset$ .and. ix_attrib < old_control_var_offset$+20)

case (all_control_var$)
  is_attrib = ((ix_attrib > var_offset$ .and. ix_attrib < var_offset$+20) .or. &
               (ix_attrib > old_control_var_offset$ .and. ix_attrib < old_control_var_offset$+20))

case (multipole$)
  is_attrib = (ix_attrib >= k0l$ .and. ix_attrib <= t21$)

case (elec_multipole$)
  is_attrib = (ix_attrib >= a0_elec$ .and. ix_attrib <= b21_elec$)

end select

end function is_attribute

end module
