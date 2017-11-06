!+
! See the paper:
!   "An Efficient Formalism for Simulating Coherent Synchrotron Radiation"
!   D. Sagan
!-

module csr_mod

use make_mat6_mod
use beam_utils
use spline_mod
use nr, only: zbrent

! csr_ele_info_struct holds info for a particular lattice element
! The centroid "chord" is the line from the centroid position at the element entrance to
! the centroid position at the element exit end.

type csr_ele_info_struct
  type (ele_struct), pointer :: ele            ! lattice element
  type (coord_struct) orbit0, orbit1           ! centroid orbit at entrance/exit ends
  type (floor_position_struct) floor0, floor1  ! Floor position of centroid at entrance/exit ends
  type (floor_position_struct) e_floor0, e_floor1  ! Floor position of element ref coords at entrance/exit ends
  type (spline_struct) spline                  ! spline for centroid orbit. spline%x = distance along chord
  real(rp) theta_chord                         ! Reference angle of chord in z-x plane
  real(rp) L_chord                             ! Chord Length. Negative if bunch moves backwards in element.
  real(rp) dL_s                                ! L_s(of element) - L_chord
end type

type csr_bunch_slice_struct  ! Structure for a single particle bin.
  real(rp) :: x0 = 0, y0 = 0 ! Transverse center of the particle distrubution
  real(rp) :: z0_edge = 0    ! Left (min z) edge of bin
  real(rp) :: z1_edge = 0    ! Right (max z) edge of bin
  real(rp) :: z_center = 0   ! z at center of bin.
  real(rp) :: sig_x = 0      ! particle's RMS width
  real(rp) :: sig_y = 0      ! particle's RMS width
  real(rp) :: charge = 0     ! charge of the particles
  real(rp) :: dcharge_density_dz = 0      ! Charge density gradient 
  real(rp) :: edge_dcharge_density_dz = 0 ! gradient between this and preceeding bin. [Evaluated at bin edge.]
  real(rp) :: kick_csr = 0                ! CSR kick
  real(rp) :: coef_lsc_plus(0:2,0:2) = 0  ! LSC Kick coefs.
  real(rp) :: coef_lsc_minus(0:2,0:2) = 0 ! LSC Kick coefs.
  real(rp) :: kick_lsc = 0
  real(rp) :: n_particle = 0 ! Number of particles in slice can be a fraction since particles span multiple bins.
end type

! csr_kick1_struct stores the CSR kick, kick integral etc. for a given source and kick positions.
! This structure also holds info on parameters that go into the kick calculation.
! Since an integration step involves one kick position and many source positions,
! the information that is only kick position dependent is held in the csr_struct and
! the csr_struct holds an array of csr_kick1_structs, one for each dz.

type csr_kick1_struct ! Sub-structure for csr calculation cache
  real(rp) I_csr            ! Kick integral.
  real(rp) I_int_csr        ! Integrated Kick integral.
  real(rp) image_kick_csr   ! kick.
  real(rp) L_vec(3)         ! L vector in global coordinates.
  real(rp) L                ! Distance between source and kick points.
  real(rp) dL               ! = epsilon_L = Ls - L
  real(rp) dz_particles     ! Kicked particle - source particle position at constant time.
  real(rp) s_chord_source   ! Source point coordinate along chord.
  real(rp) theta_L          ! Angle of L vector
  real(rp) theta_sl         ! Angle between velocity of particle at source pt and L
  real(rp) theta_lk         ! Angle between L and velocity of kicked particle
  integer ix_ele_source     ! Source element index.
  type (floor_position_struct) floor_s  ! Floor position of source pt
end type

type csr_struct           ! Structurture for binning particle averages
  real(rp) gamma, gamma2        ! Relativistic gamma factor.
  real(rp) rel_mass             ! m_particle / m_electron
  real(rp) beta                 ! Relativistic beta factor.
  real(rp) :: dz_slice = 0      ! Bin width
  real(rp) ds_track_step        ! True step size
  real(rp) s_kick               ! Kick point longitudinal location (element ref coords) from entrance end
  real(rp) s_chord_kick         ! Kick point along beam centroid line
  real(rp) y_source             ! Height of source particle.
  real(rp) kick_factor          ! Coefficient to scale the kick
  real(rp) actual_track_step    ! ds_track_step scalled by Length_centroid_chord / Length_element ratio
  real(rp) x0_bunch, y0_bunch   ! Bunch centroid
  type(floor_position_struct) floor_k   ! Floor coords at kick point
  integer species                       ! Particle type
  integer ix_ele_kick                   ! Same as element being tracked through.
  type (csr_bunch_slice_struct), allocatable :: slice(:)    ! slice(i) refers to the i^th bunch slice.
  type (csr_kick1_struct), allocatable :: kick1(:)          ! kick1(i) referes to the kick between two slices i bins apart.
  type (csr_ele_info_struct), allocatable :: eleinfo(:)     ! Element-by-element information.
  type (ele_struct), pointer :: kick_ele                    ! Element where the kick pt is == ele tracked through.
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine track1_bunch_csr (bunch_start, ele, centroid, bunch_end, err, s_start, s_end)
!
! Routine to track a bunch of particles through an element with csr radiation effects.
!
! Modules needed:
!   use csr_mod
!
! Input:
!   bunch_start  -- Bunch_struct: Starting bunch position.
!   ele          -- Ele_struct: The element to track through. Must be part of a lattice.
!   centroid(0:) -- coord_struct, Approximate beam centroid orbit for the lattice branch.
!                     Hint: Calculate this before beam tracking by tracking a single particle.
!   s_start      -- real(rp), optional: Starting position relative to ele. Default = 0
!   s_end        -- real(rp), optional: Ending position. Default is ele length.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!   err       -- Logical: Set true if there is an error. EG: Too many particles lost.
!-

subroutine track1_bunch_csr (bunch_start, ele, centroid, bunch_end, err, s_start, s_end)

implicit none

type (bunch_struct), target :: bunch_start, bunch_end
type (coord_struct), pointer :: pt
type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (ele_struct), save :: runt
type (ele_struct), pointer :: ele0, s_ele
type (csr_struct), target :: csr
type (csr_ele_info_struct), pointer :: eleinfo
type (coord_struct) :: centroid(0:)
type (floor_position_struct) floor

real(rp), optional :: s_start, s_end
real(rp) s0_step, vec0(6), vec(6), theta_chord, theta0, theta1, L
real(rp) e_tot, f1, x, z, ds_step

integer i, j, n, ie, ns, nb, n_step, n_live, i_step
integer :: iu_wake

character(*), parameter :: r_name = 'track1_bunch_csr'
logical err, auto_bookkeeper, err_flag, parallel0, parallel1

! Init

err = .true.
csr%kick_ele => ele    ! Element where the kick pt is == ele tracked through.
branch => ele%branch

! No CSR for a zero length element.
! And taylor elements get ignored.

if (ele%value(l$) == 0 .or. ele%key == taylor$) then
  ele%csr_calc_on = .false.
  call track1_bunch_hom (bunch_end, ele, branch%param, bunch_end)
  err = .false.
  ! Only do warning if previous element needed
  if (ele%key == taylor$ .and. csr_param%print_taylor_warning) then
    ele0 => pointer_to_next_ele (ele, -1)
    if (ele0%csr_calc_on) call out_io (s_warn$, r_name, &
                        'CSR calc for taylor element not done: ' // ele%name)
  endif
  return
endif

! n_step is the number of steps to take when tracking through the element.
! csr%ds_step is the true step length.

bunch_end = bunch_start

if (csr_param%n_bin <= csr_param%particle_bin_span + 1) then
  call out_io (s_fatal$, r_name, &
            'CSR_PARAM%N_BIN MUST BE GREATER THAN CSR_PARAM%PARTICLE_BIN_SPAN+1!')
  if (global_com%exit_on_error) call err_exit
  return
endif

if (csr_param%ds_track_step == 0) then
  call out_io (s_fatal$, r_name, 'CSR_PARAM%DS_TRACK_STEP NOT SET!')
  if (global_com%exit_on_error) call err_exit
  return
endif

! Calculate beam centroid info at element edges, etc.

n = min(ele%ix_ele+10, branch%n_ele_track)
allocate (csr%eleinfo(0:n))

do i = 0, n
  eleinfo => csr%eleinfo(i)
  eleinfo%ele => branch%ele(i)  ! Pointer to the P' element
  s_ele => eleinfo%ele
  eleinfo%e_floor1 = branch%ele(i)%floor
  eleinfo%e_floor1%r(2) = 0  ! Make sure in horizontal plane

  eleinfo%orbit1 = centroid(i)
  vec = eleinfo%orbit1%vec
  floor%r = [vec(1), vec(3), s_ele%value(l$)]
  eleinfo%floor1 = coords_local_curvilinear_to_floor (floor, s_ele)
  eleinfo%floor1%r(2) = 0  ! Make sure in horizontal plane
  eleinfo%floor1%theta = s_ele%floor%theta + asin(vec(2) / (1+vec(6)))

  if (i == 0) then
    eleinfo%e_floor0 = csr%eleinfo(i)%e_floor1
    eleinfo%floor0   = csr%eleinfo(i)%floor1
    eleinfo%orbit0   = csr%eleinfo(i)%orbit1
  else
    eleinfo%e_floor0 = csr%eleinfo(i-1)%e_floor1
    eleinfo%floor0   = csr%eleinfo(i-1)%floor1
    eleinfo%orbit0   = csr%eleinfo(i-1)%orbit1
  endif

  vec0 = eleinfo%orbit0%vec
  vec = eleinfo%orbit1%vec
  theta_chord = atan2(eleinfo%floor1%r(1)-eleinfo%floor0%r(1), eleinfo%floor1%r(3)-eleinfo%floor0%r(3))
  eleinfo%theta_chord = theta_chord

  ! An element with a negative step length (EG off-center beam in a patch which just has an x-pitch), can
  ! have floor%theta and theta_chord 180 degress off. Hence, restrict theta0 & theta1 to be within [-pi/2,pi/2]
  theta0 = modulo2(eleinfo%floor0%theta - theta_chord, pi/2)
  theta1 = modulo2(eleinfo%floor1%theta - theta_chord, pi/2)
  eleinfo%L_chord = sqrt((eleinfo%floor1%r(1)-eleinfo%floor0%r(1))**2 + (eleinfo%floor1%r(3)-eleinfo%floor0%r(3))**2)
  if (eleinfo%L_chord == 0) then
    eleinfo%spline = spline_struct()
    eleinfo%dL_s = 0
    cycle
  endif

  ! With a negative step length, %L_chord is negative
  parallel0 = (abs(modulo2(eleinfo%floor0%theta - theta_chord, pi)) < pi/2)
  parallel1 = (abs(modulo2(eleinfo%floor1%theta - theta_chord, pi)) < pi/2)
  if (parallel0 .neqv. parallel1) then
    call out_io (s_fatal$, r_name, 'VERY CONFUSED CSR CALCULATION! PLEASE SEAK HELP')
    if (global_com%exit_on_error) call err_exit
    return
  endif
  if (.not. parallel0) eleinfo%L_chord = -eleinfo%L_chord

  call create_a_spline (eleinfo%spline, [0.0_rp, 0.0_rp], [eleinfo%L_chord, 0.0_rp], theta0, theta1)
  eleinfo%dL_s = dspline_len(0.0_rp, eleinfo%L_chord, eleinfo%spline)
enddo

! make sure that ele_len / track_step is an integer.

n_step = max (1, nint(ele%value(l$) / csr_param%ds_track_step))
csr%ds_track_step = ele%value(l$) / n_step
csr%species = bunch_start%particle(1)%species
csr%ix_ele_kick = ele%ix_ele
csr%actual_track_step = csr%ds_track_step * (csr%eleinfo(ele%ix_ele)%L_chord / ele%value(l$))

auto_bookkeeper = bmad_com%auto_bookkeeper ! save state
bmad_com%auto_bookkeeper = .false.   ! make things go faster

!----------------------------------------------------------------------------------------
! Loop over the tracking steps
! runt is the element that is tracked through at each step.

do i_step = 0, n_step

  ! track through the runt

  if (i_step /= 0) then
    call create_uniform_element_slice (ele, branch%param, i_step, n_step, runt, s_start, s_end)
    runt%csr_calc_on = .false.
    call track1_bunch_hom (bunch_end, runt, branch%param, bunch_end)
  endif

  s0_step = i_step * csr%ds_track_step
  if (present(s_start)) s0_step = s0_step + s_start

  ! Cannot do a realistic calculation if there are less particles than bins

  n_live = count(bunch_end%particle%state == alive$)
  if (n_live < csr_param%n_bin) then
    call out_io (s_error$, r_name, 'NUMBER OF LIVE PARTICLES: \i0\ ', &
                          'LESS THAN NUMBER OF BINS FOR CSR CALC.', &
                          'AT ELEMENT: ' // trim(ele%name) // '  [# \i0\] ', &
                          i_array = [n_live, ele%ix_ele ])
    return
  endif

  ! Assume a linear energy gain in a cavity

  f1 = s0_step / ele%value(l$)
  e_tot = f1 * branch%ele(ele%ix_ele-1)%value(e_tot$) + (1 - f1) * ele%value(e_tot$)
  call convert_total_energy_to (e_tot, branch%param%particle, csr%gamma, beta = csr%beta)
  csr%gamma2 = csr%gamma**2
  csr%rel_mass = mass_of(branch%param%particle) / m_electron 

  call csr_bin_particles (bunch_end%particle, csr, err_flag); if (err_flag) return

  csr%s_kick = s0_step
  csr%s_chord_kick = s_ref_to_s_chord (s0_step, csr%eleinfo(ele%ix_ele))
  z = csr%s_chord_kick
  x = spline1(csr%eleinfo(ele%ix_ele)%spline, z)
  theta_chord = csr%eleinfo(ele%ix_ele)%theta_chord
  csr%floor_k%r = [x*cos(theta_chord)+z*sin(theta_chord), 0.0_rp, -x*sin(theta_chord)+z*cos(theta_chord)] + &
                      csr%eleinfo(ele%ix_ele)%floor0%r
  csr%floor_k%theta = theta_chord + spline1(csr%eleinfo(ele%ix_ele)%spline, z, 1)

  ! ns = 0 is the unshielded kick.

  do ns = 0, csr_param%n_shield_images
    ! The factor of -1^ns accounts for the sign of the image currents
    ! Take into account that at the endpoints we are only putting in a half kick.
    ! The factor of two is due to there being image currents both above and below.

    csr%kick_factor = (-1)**ns
    if (i_step == 0 .or. i_step == n_step) csr%kick_factor = csr%kick_factor / 2
    if (ns /= 0) csr%kick_factor = 2* csr%kick_factor

    csr%y_source = ns * csr_param%beam_chamber_height

    call csr_bin_kicks (s0_step, csr, err_flag)
    if (err_flag) return
  enddo

  ! loop over all particles and give them a kick

  do j = 1, size(bunch_end%particle)
    if (bunch_end%particle(j)%state /= alive$) cycle
    call csr_kick_calc (csr, bunch_end%particle(j))
  enddo

  call save_bunch_track (bunch_end, ele, s0_step)

  ! Record to file?

  if (csr_param%write_csr_wake) then
    iu_wake = lunget()
    open (iu_wake, file = 'csr_wake.dat', access = 'append')
    if (i_step == 0) then
      write (iu_wake, '(a)') '!------------------------------------------------------------'
      write (iu_wake, '(a, i6, 2x, a)') '! ', ele%ix_ele, trim(ele%name)
    endif
    write (iu_wake, '(a)') '!#-----------------------------'
    write (iu_wake, '(a, i4, f12.6)') '! ', i_step, ele%s_start + s0_step
    write (iu_wake, '(a)') '!         Z   Charge/Meter    CSR_Kick/m       I_CSR/m      S_Source' 
    ds_step = csr%kick_factor * csr%actual_track_step
    do j = 1, csr_param%n_bin
      ele0 => branch%ele(csr%kick1(j)%ix_ele_source)
      write (iu_wake, '(f12.6, 4es14.6)') csr%slice(j)%z_center, &
                  csr%slice(j)%charge/csr%dz_slice, csr%slice(j)%kick_csr/ds_step, &
                  csr%kick1(j)%I_csr/ds_step, ele0%s_start + csr%kick1(j)%s_chord_source
    enddo
    close (iu_wake)
  endif

enddo

bmad_com%auto_bookkeeper = auto_bookkeeper  ! restore state
err = .false.

end subroutine track1_bunch_csr

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_bin_parcticles (particle, csr)
!
! Routine to bin the particles longitudinally in s. 
!
! To avoid noise in the cacluation, every particle is considered to have a 
! triangular distribution with a base length  given by
!   csr_param%particle_bin_span * csr%dz_slice.
! That is, particles will, in general, overlap multiple bins. 
!
! Input:
!   particle(:)          -- Coord_struct: Array of particles
!   csr_param            -- Csr_parameter_struct: CSR common block (not an argument).
!     %n_bin             -- Number of bins.
!     %particle_bin_span -- Particle length / dz_slice. 
!
! Output:
!   csr         -- Csr_struct: The bin structure.
!     %dz_slice     -- Bin longitudinal length
!     %slice(1:) -- Array of bins.
!-

subroutine csr_bin_particles (particle, csr, err_flag)

implicit none

type this_local_struct   ! Temporary structure 
  real(rp) :: charge     ! how much charge of particle in bin
  real(rp) :: x0, y0     ! particle center
  integer ib             ! bin index
end type

type (coord_struct), target :: particle(:)
type (coord_struct), pointer :: p
type (csr_struct), target :: csr
type (this_local_struct), allocatable :: tloc(:)
type (csr_bunch_slice_struct), pointer :: slice

real(rp) z_center, z_min, z_max, dz_particle, dz, z_maxval, z_minval, c_tot
real(rp) zp_center, zp0, zp1, zb0, zb1, charge, overlap_fraction, f, last_sig_x, last_sig_y

integer i, j, n, ix0, ib, ib2, ic, ib_center, n_bin_eff

logical err_flag

character(*), parameter :: r_name = 'csr_bin_particles'

! Init bins...
! The left edge of csr%slice(1) is at z_min
! The right edge of csr%slice(n_bin) is at z_max
! The first and last bins are empty.

err_flag = .false.

if (.not. csr_param%lcsr_component_on .and. .not. csr_param%lsc_component_on .and. &
    .not. csr_param%tsc_component_on .and. csr_param%n_shield_images == 0) return

n_bin_eff = csr_param%n_bin - 2 - (csr_param%particle_bin_span + 1)
if (n_bin_eff < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF CSR BINS TOO SMALL: \i0\ ', &
              'MUST BE GREATER THAN 3 + PARTICLE_BIN_SPAN.', i_array = [csr_param%n_bin])
  if (global_com%exit_on_error) call err_exit
  particle%state = lost$
  err_flag = .true.
  return
endif

z_maxval = maxval(particle(:)%vec(5), mask = (particle(:)%state == alive$))
z_minval = minval(particle(:)%vec(5), mask = (particle(:)%state == alive$))
dz = z_maxval - z_minval
csr%dz_slice = dz / n_bin_eff
csr%dz_slice = 1.0000001 * csr%dz_slice     ! to prevent round off problems
z_center = (z_maxval + z_minval) / 2
z_min = z_center - csr_param%n_bin * csr%dz_slice / 2
z_max = z_center + csr_param%n_bin * csr%dz_slice / 2
dz_particle = csr_param%particle_bin_span * csr%dz_slice

if (dz == 0) then
  call out_io (s_fatal$, r_name, 'LONGITUDINAL WIDTH OF BEAM IS ZERO!')
  if (global_com%exit_on_error) call err_exit
  err_flag = .true.
  return
endif

! allocate memeory for the bins

if (allocated(csr%slice)) then
  if (size(csr%slice, 1) < csr_param%n_bin) deallocate (csr%slice)
endif

if (.not. allocated(csr%slice)) &
    allocate (csr%slice(csr_param%n_bin), csr%kick1(-csr_param%n_bin:csr_param%n_bin))

csr%slice(:) = csr_bunch_slice_struct()  ! Zero everything

! Fill in some z information

do i = 1, csr_param%n_bin
  csr%slice(i)%z0_edge  = z_min + (i - 1) * csr%dz_slice
  csr%slice(i)%z_center = csr%slice(i)%z0_edge + csr%dz_slice / 2
  csr%slice(i)%z1_edge  = csr%slice(i)%z0_edge + csr%dz_slice
enddo

! Init the tloc structure...
! Each particle is distributed longitudinally in a triangular fashion.
! The tloc records how much charge for a given particle is in a bin.

n = size(particle) * (csr_param%particle_bin_span + 2)
allocate (tloc(n))
tloc%ib = -1

! Compute the particle distribution center in each bin

! The contribution to the charge in a bin from a particle is computed from the overlap
! between the particle and the bin.
 
ic = 0
c_tot = 0    ! Used for debugging sanity check
do i = 1, size(particle)
  p => particle(i)
  if (p%state /= alive$) cycle
  zp_center = p%vec(5) ! center of particle
  zp0 = zp_center - dz_particle / 2       ! particle left edge 
  zp1 = zp_center + dz_particle / 2       ! particle right edge 
  ix0 = nint((zp0 - z_min) / csr%dz_slice)  ! left most bin index
  do j = 0, csr_param%particle_bin_span+1
    ib = j + ix0
    slice => csr%slice(ib)
    zb0 = csr%slice(ib)%z0_edge
    zb1 = csr%slice(ib)%z1_edge   ! edges of the bin
    overlap_fraction = particle_overlap_in_bin (zb0, zb1)
    slice%n_particle = slice%n_particle + overlap_fraction
    charge = overlap_fraction * p%charge
    slice%charge = slice%charge + charge
    slice%x0 = slice%x0 + p%vec(1) * charge
    slice%y0 = slice%y0 + p%vec(3) * charge
    c_tot = c_tot + charge
    ic = ic + 1
    tloc(ic)%charge = charge
    tloc(ic)%x0 = p%vec(1)
    tloc(ic)%y0 = p%vec(3)
    tloc(ic)%ib = ib
  enddo
enddo

do ib = 1, csr_param%n_bin
  if (ib /= 1) csr%slice(ib)%edge_dcharge_density_dz = (csr%slice(ib)%charge - csr%slice(ib-1)%charge) / csr%dz_slice**2
  if (ib == 1) then
    csr%slice(ib)%dcharge_density_dz = (csr%slice(ib+1)%charge - csr%slice(ib)%charge) / csr%dz_slice**2
  elseif (ib == csr_param%n_bin) then
    csr%slice(ib)%dcharge_density_dz = (csr%slice(ib)%charge - csr%slice(ib-1)%charge) / csr%dz_slice**2
  else
    csr%slice(ib)%dcharge_density_dz = (csr%slice(ib+1)%charge - csr%slice(ib-1)%charge) / (2 * csr%dz_slice**2)
  endif

  if (csr%slice(ib)%charge == 0) cycle
  csr%slice(ib)%x0 = csr%slice(ib)%x0 / csr%slice(ib)%charge
  csr%slice(ib)%y0 = csr%slice(ib)%y0 / csr%slice(ib)%charge
enddo

csr%x0_bunch = sum(csr%slice%x0 * csr%slice%charge) / sum(csr%slice%charge)
csr%y0_bunch = sum(csr%slice%y0 * csr%slice%charge) / sum(csr%slice%charge)

! Compute the particle distribution sigmas in each bin.
! Abs(x-x0) is used instead of the usual formula involving (x-x0)^2 to lessen the effect
! of non-Gaussian tails.

do ic = 1, size(tloc)
  if (tloc(ic)%ib < 0) cycle
  slice => csr%slice(tloc(ic)%ib)
  if (slice%n_particle < csr_param%sc_min_in_bin) cycle
  slice%sig_x = slice%sig_x + abs(tloc(ic)%x0 - slice%x0) * tloc(ic)%charge
  slice%sig_y = slice%sig_y + abs(tloc(ic)%y0 - slice%y0) * tloc(ic)%charge
enddo

f = sqrt(pi/2)
do ib = 1, csr_param%n_bin
  slice => csr%slice(ib)
  if (slice%n_particle < csr_param%sc_min_in_bin) cycle
  slice%sig_x = f * slice%sig_x / slice%charge
  slice%sig_y = f * slice%sig_y / slice%charge
enddo

! For bins where there was not enough particles to calculate sigmas, use the sigmas from nearby bins.
! Start at a slice near the center and work outwards.

do ib = 0, csr_param%n_bin/2
  ib_center = csr_param%n_bin/2 + ib
  slice => csr%slice(ib_center)
  if (slice%sig_x /= 0) exit
  ib_center = csr_param%n_bin/2 - ib
  slice => csr%slice(ib_center)
  if (slice%sig_x /= 0) exit
enddo

last_sig_x = slice%sig_x
last_sig_y = slice%sig_y

do ib = ib_center, 1, -1
  slice => csr%slice(ib)
  if (slice%sig_x == 0) then
    ib2 = max(1, ib-1)  ! Bin on other side
    if (csr%slice(ib2)%sig_x == 0) then  ! If this other side bin has no sig_x then ignore
      slice%sig_x = last_sig_x
      slice%sig_y = last_sig_y
    else  ! Else average the bins
      slice%sig_x = (last_sig_x + csr%slice(ib2)%sig_x) / 2
      slice%sig_y = (last_sig_y + csr%slice(ib2)%sig_y) / 2
    endif
  else
    last_sig_x = slice%sig_x
    last_sig_y = slice%sig_y
  endif
enddo

do ib = ib_center, csr_param%n_bin
  slice => csr%slice(ib)
  if (slice%sig_x == 0) then
    ib2 = min(csr_param%n_bin, ib+1) ! Bin on other side
    if (csr%slice(ib2)%sig_x == 0) then  ! If this other side bin has no sig_x then ignore
      slice%sig_x = last_sig_x
      slice%sig_y = last_sig_y
    else  ! Else average the bins
      slice%sig_x = (last_sig_x + csr%slice(ib2)%sig_x) / 2
      slice%sig_y = (last_sig_y + csr%slice(ib2)%sig_y) / 2
    endif
  else
    last_sig_x = slice%sig_x
    last_sig_y = slice%sig_y
  endif
enddo

!---------------------------------------------------------------------------
contains

! computes the contribution to the charge in a bin from a given particle.
! z0_bin, z1_bin are the edge positions of the bin

function particle_overlap_in_bin (z0_bin, z1_bin) result (overlap)

real(rp) z0_bin, z1_bin, overlap, z1, z2

! Integrate over left triangular half of particle distribution

z1 = max(zp0, z0_bin)        ! left integration edge
z2 = min(zp_center, z1_bin)  ! right integration edge
if (z2 > z1) then            ! If left particle half is in bin ...
  overlap = 2 * real(((z2 - zp0)**2 - (z1 - zp0)**2), rp) / dz_particle**2
else
  overlap = 0
endif

! Integrate over right triangular half of particle distribution

z1 = max(zp_center, z0_bin)  ! left integration edge
z2 = min(zp1, z1_bin)        ! right integration edge
if (z2 > z1) then            ! If right particle half is in bin ...
  overlap = overlap + 2 * real(((z1 - zp1)**2 - (z2 - zp1)**2), rp) / dz_particle**2
endif

end function particle_overlap_in_bin

end subroutine csr_bin_particles

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_bin_kicks (ds_kick_pt, csr, err_flag)
!
! Routine to cache intermediate values needed for the csr calculations.
!
! Input:
!   ds_kick_pt   -- real(rp): Distance between the beginning of the element we are
!                    tracking through and the kick point (which is within this element).
!   csr          -- csr_struct: 
!
! Output:
!   csr          -- csr_struct: 
!     %kick1(:)          -- CSR kick calculation bin array. 
!     %slice(:)%kick_csr -- Integrated kick
!   err_flag     -- logical: Set True if there is an error. False otherwise
!-

subroutine csr_bin_kicks (ds_kick_pt, csr, err_flag)

implicit none

type (csr_struct), target :: csr
type (branch_struct), pointer :: branch
type (csr_kick1_struct), pointer :: kick1

real(rp) ds_kick_pt, coef

integer i, n_bin
logical err_flag

character(16) :: r_name = 'csr_bin_kicks'

! The kick point P is fixed.
! Loop over all kick1 bins and compute the kick.

err_flag = .false.

do i = lbound(csr%kick1, 1), ubound(csr%kick1, 1)

  kick1 => csr%kick1(i)
  kick1%dz_particles = i * csr%dz_slice

  if (i == lbound(csr%kick1, 1)) then
    kick1%ix_ele_source = csr%ix_ele_kick  ! Initial guess where source point is
  else
    kick1%ix_ele_source = csr%kick1(i-1)%ix_ele_source
  endif

  ! Calculate what element the source point is in.

  kick1%s_chord_source = s_source_calc(kick1, csr, err_flag)
  if (err_flag) return

  ! calculate csr.
  ! I_csr is only calculated for particles with y = 0 and not for image currents.

  if (csr%y_source == 0) then
    call I_csr (kick1, i, csr)
  else
    call image_charge_kick_calc (kick1, csr)
  endif

enddo

!

coef = csr%actual_track_step * r_e / (csr%rel_mass * e_charge * abs(charge_of(csr%species)) * csr%gamma)
n_bin = csr_param%n_bin

! CSR & Image charge kick

if (csr%y_source == 0) then
  if (csr_param%lcsr_component_on) then
    do i = 1, n_bin
      csr%slice(i)%kick_csr = coef * dot_product(csr%kick1(i:1:-1)%I_int_csr, csr%slice(1:i)%edge_dcharge_density_dz)
    enddo
  endif

else  ! Image charge
  do i = 1, n_bin
    csr%slice(i)%kick_csr = csr%slice(i)%kick_csr + coef * &
                  dot_product(csr%kick1(i-1:i-n_bin:-1)%image_kick_csr, csr%slice(1:n_bin)%charge)
  enddo
endif

! Space charge kick

if (csr_param%lsc_component_on) then
  if (csr%y_source == 0) then    ! image charges do not contribute to LSC.
    call lsc_kick_params_calc (csr)
  endif
endif

end subroutine csr_bin_kicks
  
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function s_source_calc (kick1, csr, err_flag) result (s_source)
!
! Routine to calculate the distance between source and kick points.
!
! Input:
!   kick1         -- csr_kick1_struct:
!   csr           -- csr_struct:
!
! Output:
!   s_source      -- real(rp): source s-position.
!   kick1         -- csr_kick1_struct:
!   err_flag      -- logical: Set True if there is an error. Untouched otherwise.
!-

function s_source_calc (kick1, csr, err_flag) result (s_source)

implicit none

type (csr_kick1_struct), target :: kick1
type (csr_struct), target :: csr
type (csr_ele_info_struct), pointer :: einfo_s, einfo_k
type (floor_position_struct), pointer :: fk, f0, fs
type (ele_struct), pointer :: ele

real(rp) a, b, c, dz, s_source, beta2, L0, Lz, ds_source
real(rp) z0, z1, sz_kick, sz0, Lsz0, ddz0, ddz1

integer i, last_step
logical err_flag

character(*), parameter :: r_name = 's_source_calc'

! Each interation of the loop looks for a possible source point in the lattice element with 
! index kick1%ix_ele_source. If found, return. If not, move on to another element

dz = kick1%dz_particles   ! Target distance.
beta2 = csr%beta**2
last_step = 0
s_source = 0   ! To prevent uninitalized complaints from the compiler.

do

  einfo_s => csr%eleinfo(kick1%ix_ele_source)

  ! If at beginning of lattice assume an infinite drift.
  ! s_source will be negative

  if (kick1%ix_ele_source == 0) then
    fk => csr%floor_k
    fs => kick1%floor_s
    f0 => einfo_s%floor1

    L0 = sqrt((fk%r(1) - f0%r(1))**2 + (fk%r(3) - f0%r(3))**2 + csr%y_source**2)
    ! L_z is the z-component from lat start to the kick point.
    Lz = (fk%r(1) - f0%r(1)) * sin(f0%theta) + (fk%r(3) - f0%r(3)) * cos(f0%theta) 
    ! Lsz0 is Ls from the lat start to the kick point
    Lsz0 = dspline_len(0.0_rp, csr%s_chord_kick, csr%eleinfo(csr%ix_ele_kick)%spline) + csr%s_chord_kick
    do i = 1, csr%ix_ele_kick - 1
      Lsz0 = Lsz0 + csr%eleinfo(i)%dL_s + csr%eleinfo(i)%L_chord
    enddo

    a = 1/csr%gamma2
    b = 2 * (Lsz0 - dz - beta2 * Lz)
    c = (Lsz0 - dz)**2 - beta2 * L0**2
    ds_source = -(-b + sqrt(b**2 - 4 * a * c)) / (2 * a)
    s_source = einfo_s%ele%s + ds_source

    fs%r = [f0%r(1) + ds_source * sin(f0%theta), csr%y_source, f0%r(3) + ds_source * cos(f0%theta)]
    fs%theta = f0%theta
    kick1%L_vec = fk%r - fs%r
    kick1%L = sqrt(dot_product(kick1%L_vec, kick1%L_vec))
    kick1%dL = lsz0 - ds_source - kick1%L  ! Remember s_source is negative
    kick1%theta_L = atan2(kick1%L_vec(1), kick1%L_vec(3))
    kick1%theta_sl = f0%theta - kick1%theta_L
    einfo_k => csr%eleinfo(csr%ix_ele_kick)
    kick1%theta_lk = kick1%theta_L - (spline1(einfo_k%spline, csr%s_chord_kick, 1) + einfo_k%theta_chord)
    return
  endif

  !

  ele => einfo_s%ele

  if (ele%key == match$ .and. (abs(ele%value(x1$) - ele%value(x0$)) > 0 .or. &
                               abs(ele%value(y1$) - ele%value(y0$)) > 0)) then
    call out_io (s_fatal$, r_name, &
        'CSR CALC IS NOT ABLE TO HANDLE A MATCH ELEMENT WITH FINITE ORBIT OFFSETS: ' // ele%name, &
        'WHILE TRACKING THROUGH ELEMENT: ', csr%kick_ele%name)
    if (global_com%exit_on_error) call err_exit
    err_flag = .true.
    return
  endif

  if (ele%key == floor_shift$) then
    call out_io (s_fatal$, r_name, &
        'CSR CALC IS NOT ABLE TO HANDLE A FLOOR_SHIFT ELEMENT: ' // ele%name, &
        'WHILE TRACKING THROUGH ELEMENT: ', csr%kick_ele%name)
    if (global_com%exit_on_error) call err_exit
    err_flag = .true.
    return
  endif

  ! Look at ends of the source element and check if the source point is within the element or not.
  ! Generally dz decreases with increasing s but this may not be true for patch elements.

  ddz0 = ddz_calc_csr(0.0_rp)
  ddz1 = ddz_calc_csr(einfo_s%L_chord)

  if (last_step == -1 .and. ddz1 > 0) then  ! Roundoff error is causing ddz1 to be positive.
    s_source = ddz_calc_csr(einfo_s%L_chord)
    return
  endif

  if (last_step == 1 .and. ddz0 < 0) then  ! Roundoff error is causing ddz0 to be negative.
    s_source = ddz_calc_csr(0.0_rp)
    return
  endif

  if (ddz0 < 0 .and. ddz1 < 0) then
    if (last_step == 1) exit       ! Round off error can cause problems
    last_step = -1
    kick1%ix_ele_source = kick1%ix_ele_source - 1
    cycle
  endif

  if (ddz0 > 0 .and. ddz1 > 0) then
    if (kick1%ix_ele_source == csr%ix_ele_kick) return  ! Source ahead of kick pt => ignore.
    if (last_step == -1) exit      ! Round off error can cause problems
    last_step = 1
    kick1%ix_ele_source = kick1%ix_ele_source + 1
    cycle
  endif

  ! Only possibility left is that root is bracketed.

  s_source = zbrent (ddz_calc_csr, 0.0_rp, einfo_s%L_chord, 1d-8)
  return
    
enddo

call out_io (s_fatal$, r_name, &
    'CSR CALCULATION ERROR. PLEASE REPORT THIS...', &
    'WHILE TRACKING THROUGH ELEMENT: ', csr%kick_ele%name)
if (global_com%exit_on_error) call err_exit
err_flag = .true.

!----------------------------------------------------------------------------
contains

!+
! Function ddz_calc_csr (s_chord_source) result (ddz_this)
!
! Routine to calculate the distance between the source particle and the
! kicked particle at constant time minus the target distance.
!
! Input:
!   s_chord_source  -- real(rp): Chord distance from start of element.
!
! Output:
!   ddz_this        -- real(rp): Distance between source and kick particles: Calculated - Wanted.
!-

function ddz_calc_csr (s_chord_source) result (ddz_this)

implicit none

type (csr_ele_info_struct), pointer :: ce
real(rp), intent(in) :: s_chord_source
real(rp) ddz_this, x, z, c, s, dtheta_L
real(rp) s0, s1, ds, theta_L, dL

integer i

character(*), parameter :: r_name = 'ddz_calc_csr'

! 

x = spline1(einfo_s%spline, s_chord_source)
c = cos(einfo_s%theta_chord)
s = sin(einfo_s%theta_chord)
kick1%floor_s%r = [x*c + s_chord_source*s, csr%y_source, -x*s + s_chord_source*c] + einfo_s%floor0%r    ! Floor coordinates

kick1%L_vec = csr%floor_k%r - kick1%floor_s%r
kick1%L = sqrt(dot_product(kick1%L_vec, kick1%L_vec))
kick1%theta_L = atan2(kick1%L_vec(1), kick1%L_vec(3))

s0 = s_chord_source
s1 = csr%s_chord_kick

if (kick1%ix_ele_source == csr%ix_ele_kick) then
  ds = s1 - s0
  ! dtheta_L = angle of L line in centroid chord ref frame
  dtheta_L = einfo_s%spline%coef(1) + einfo_s%spline%coef(2) * (2*s0 + ds) + einfo_s%spline%coef(3) * (3*s0**2 + 3*s0*ds + ds**2)
  dL = dspline_len(s0, s1, einfo_s%spline, dtheta_L) ! = Ls - L
  ! Ls is negative if the source pt is ahead of the kick pt (ds < 0). But L is always positive. 
  if (ds < 0) dL = dL + 2 * ds    ! Correct for L always being positive.
  kick1%theta_sl = spline1(einfo_s%spline, s0, 1) - dtheta_L
  kick1%theta_lk = dtheta_L - spline1(einfo_s%spline, s1, 1)

else
  ! In an element where the beam centroid takes a backstep, %theta_chord will be anti-parallel to
  ! other angles. This is why pi/2 is used with modulo2 to make sure angles are in the range [-pi/2, pi/2]
  theta_L = kick1%theta_L
  dL = dspline_len(s0, einfo_s%L_chord, einfo_s%spline, modulo2(theta_L-einfo_s%theta_chord, pi/2))
  do i = kick1%ix_ele_source+1, csr%ix_ele_kick-1
    ce => csr%eleinfo(i)
    dL = dL + dspline_len(0.0_rp, ce%L_chord, ce%spline, modulo2(theta_L-ce%theta_chord, pi/2))
  enddo
  ce => csr%eleinfo(csr%ix_ele_kick)
  dL = dL + dspline_len(0.0_rp, s1, ce%spline, modulo2(theta_L-ce%theta_chord, pi/2))
  kick1%theta_sl = modulo2((spline1(einfo_s%spline, s0, 1) + einfo_s%theta_chord) - theta_L, pi/2)
  kick1%theta_lk = modulo2(theta_L - (spline1(ce%spline, s1, 1) + ce%theta_chord), pi/2)
endif

kick1%floor_s%theta = kick1%theta_sl + kick1%theta_L

! The above calc for dL neglected csr%y_source. So must correct for this.

if (csr%y_source /= 0) dL = dL - (kick1%L - sqrt(kick1%L_vec(1)**2 + kick1%L_vec(3)**2))
kick1%dL = dL
ddz_this = kick1%L / (2 * csr%gamma2) + dL
ddz_this = ddz_this - kick1%dz_particles

end function ddz_calc_csr

end function s_source_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine lsc_kick_params_calc (csr)
!
! Routine to cache intermediate values needed for the lsc calculation.
! This routine is not for image currents.
!
! Modules needed:
!   use lsc_mod
!
! Input:
!   csr       -- csr_struct: 
!     %slice(:)   -- bin array of particle averages.
!
! Output:
!   csr     -- csr_struct: Binned particle averages.
!     %slice(:)   -- bin array of particle averages.
!-

subroutine lsc_kick_params_calc (csr)

use da2_mod

implicit none

type (csr_struct), target :: csr
type (csr_bunch_slice_struct), pointer :: slice

real(rp), pointer :: f(:,:)
real(rp) sx, sy, a, b, c, z_slice, factor, sig_x_ave, sig_y_ave, charge_tot, f00
real(rp) radix, sr, z1, z2, rho0, drho_dz, dk0, atz1, atz2, bcd, abcz1, abcz2, b2cz1, b2cz2
real(rp) sx2, sy2, g_z1, g_z2, h_z1, h_z2, alph, bet, ss, f1(0:2,0:2)

integer i, j, sign_of_z_slice

character(*), parameter :: r_name = 'lsc_kick_params_calc'

! If there are too few particles in a bin the sigma calc may give a bad value.
! This can be a problem if the computed sigma is small.
! Therefore  ignore any bins with a very small sigma. 
! To know what is "small" is, compute the average sigma

charge_tot = sum(csr%slice(:)%charge)
if (charge_tot == 0) return

sig_x_ave = dot_product(csr%slice(:)%sig_x, csr%slice(:)%charge) / charge_tot
sig_y_ave = dot_product(csr%slice(:)%sig_y, csr%slice(:)%charge) / charge_tot

if (sig_y_ave == 0 .or. sig_x_ave == 0) return  ! Symptom of not enough particles.
if (.not. csr_param%lsc_component_on) return

factor = csr%kick_factor * csr%actual_track_step * r_e / &
          (csr%rel_mass * e_charge * abs(charge_of(csr%species)) * csr%gamma)

! Compute the kick at the center of each bin
! i = index of slice where kick is computed

do i = 1, csr_param%n_bin

  ! Loop over all slices and calculate kick at slice i due to slice j.

  do j = 1, csr_param%n_bin
    slice => csr%slice(j)
    sx = slice%sig_x
    sy = slice%sig_y
    if (sx < sig_x_ave * csr_param%sigma_cutoff .or. sy < sig_y_ave * csr_param%sigma_cutoff) then
      slice%sig_x = 0  ! Mark for tsc calc.
      slice%sig_y = 0
      cycle
    endif
    sx2 = sx*sx
    sy2 = sy*sy
    a = sx * sy
    b = csr%gamma * (sx2 + sy2) / (sx + sy)
    c = csr%gamma**2

    z_slice = csr%slice(i)%z_center - csr%slice(j)%z_center
    sign_of_z_slice = sign_of(z_slice)

    ! The kick is computed for z_slice positive and then the sign of the kick is corrected.
    radix = -b**2 + 4 * a * c
    sr = sqrt(abs(radix))
    z1 = abs(z_slice) - csr%dz_slice / 2
    z2 = abs(z_slice) + csr%dz_slice / 2
    drho_dz = csr%slice(j)%dcharge_density_dz * sign_of_z_slice
    rho0 = csr%slice(j)%charge / csr%dz_slice - drho_dz * abs(z_slice)

    if (i == j) then ! Self slice kick
      drho_dz = csr%slice(j)%dcharge_density_dz
      rho0 = 0
      z1 = 0
      z2 = csr%dz_slice/2       ! Integrate over 1/2 the slice
      sign_of_z_slice = -2      ! Factor of 2 accounts for 1/2 we did not integrate over.
    endif

    bcd = 2 * c * rho0 - b * drho_dz
    abcz1 = a + b*z1 + c*z1**2
    abcz2 = a + b*z2 + c*z2**2
    b2cz1 = b + 2*c*z1
    b2cz2 = b + 2*c*z2

    if (radix > 0) then
      atz1 = atan (b2cz1/sr) / sr
      atz2 = atan (b2cz2/sr) / sr
    else
      atz1 = log((b2cz1 - sr) / (b2cz1 + sr)) / (2 * sr)
      atz2 = log((b2cz2 - sr) / (b2cz2 + sr)) / (2 * sr)
    endif

    dk0 = factor * ((2 * atz2 * bcd + drho_dz * log(abcz2)) - (2 * atz1 * bcd + drho_dz * log(abcz1))) / (2 * c)
    if (dk0 == 0) cycle

    csr%slice(i)%kick_lsc = csr%slice(i)%kick_lsc + sign_of_z_slice * dk0

    if (csr_param%lsc_kick_transverse_dependence) then
      f00 = 1 / dk0
      g_z1 = a * (b2cz1*bcd + 4*abcz1*atz1*bcd*c - drho_dz*radix) / (2*abcz1*c*radix)
      g_z2 = a * (b2cz2*bcd + 4*abcz2*atz2*bcd*c - drho_dz*radix) / (2*abcz2*c*radix)
      h_z1 = (b*b*b2cz1*bcd - 6*a*b2cz1*bcd*c - 4*abcz1*b2cz1*bcd*c - 24*abcz1**2*atz1*bcd*c**2 + drho_dz*radix**2 - 2*b*b2cz1*bcd*c*z1 - 2*b2cz1*bcd*(c*z1)**2)
      h_z2 = (b*b*b2cz2*bcd - 6*a*b2cz2*bcd*c - 4*abcz2*b2cz2*bcd*c - 24*abcz2**2*atz2*bcd*c**2 + drho_dz*radix**2 - 2*b*b2cz2*bcd*c*z2 - 2*b2cz2*bcd*(c*z2)**2)

      alph = f00**2 * (g_z2 - g_z1)
      bet = f00**3 * (g_z2 - g_z1)**2 + (f00 * a)**2 * (h_z2/abcz2**2 - h_z1/abcz1**2) / (4 * c * radix**2)

      ss = 1.0_rp / sign_of_z_slice
      f1(0,0) = ss * f00
      f1(0,1) = ss * alph / (2 * sy2)
      f1(0,2) = ss * (alph + 2*bet) / (8 * sy2**2)
      f1(1,0) = ss * alph / (2 * sx2)
      f1(1,1) = ss * (alph + 2*bet) / (4 * sx2 * sy2)
      f1(2,0) = ss * (alph + 2*bet) / (8 * sx2**2)

      if (f1(0,0) > 0) then
        f => csr%slice(i)%coef_lsc_plus
      else
        f => csr%slice(i)%coef_lsc_minus
      endif

      if (f(0,0) == 0) then  ! If first time
        f = f1
      else
        f = da2_div(da2_mult(f1, f), f + f1)
      endif
    endif

  enddo

enddo

end subroutine lsc_kick_params_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine I_csr (kick1, i_bin, csr) 
!
! Routine to calculate the CSR kick integral (at y = 0)
!
! Input:
!   kick1      -- csr_kick1_struct: 
!   i_bin      -- integer: Bin index.
!   csr    -- csr_struct:
!
! Output:
!   kick1     -- csr_kick1_struct: 
!     %I_csr     -- real(rp): CSR kick integral.
!     %I_int_csr -- real(rp): Integral of I_csr. Only calculated for i_bin = 1 since it is not needed otherwise.
!-

subroutine I_csr (kick1, i_bin, csr)

implicit none

type (csr_kick1_struct), target :: kick1
type (csr_struct), target :: csr
type (csr_kick1_struct), pointer :: k
type (csr_ele_info_struct), pointer :: eleinfo
type (spline_struct), pointer :: spl

real(rp) g_bend, z, zz, Ls, L, dtheta_L, dL, s_chord_kick
integer i_bin, ix_ele_kick

!

kick1%I_int_csr = 0
kick1%I_csr = 0

! No kick when source particle is ahead of the kicked particle

z = kick1%dz_particles
if (z <= 0) return

!

k => kick1
k%I_csr = -csr%kick_factor * 2 * (k%dL / z + csr%gamma2 * k%theta_sl * k%theta_lk / (1 + csr%gamma2 * k%theta_sl**2)) / k%L

! I_csr Integral 

if (i_bin == 1) then
  ix_ele_kick = csr%ix_ele_kick
  s_chord_kick = csr%s_chord_kick
  do
    if (s_chord_kick /= 0) exit  ! Not at edge of element and element has finite length
    ix_ele_kick = ix_ele_kick - 1
    if (ix_ele_kick == 0) return  ! No kick from before beginning of lattice
    s_chord_kick = csr%eleinfo(ix_ele_kick)%L_chord
  enddo

  spl => csr%eleinfo(ix_ele_kick)%spline
  g_bend = -spline1(spl, s_chord_kick, 2) / sqrt(1 + spline1(spl, s_chord_kick, 1)**2)**3

  if (k%ix_ele_source == ix_ele_kick) then
    Ls = k%L + k%dL
    k%I_int_csr = -csr%kick_factor * ((g_bend * Ls/2)**2 - log(2 * csr%gamma2 * z / Ls) / csr%gamma2)  
  else
    ! Since source pt is in another element, split integral into pieces.
    ! First integrate over element containing the kick point.
    L = s_chord_kick
    eleinfo => csr%eleinfo(ix_ele_kick)
    dtheta_L = eleinfo%spline%coef(1) + eleinfo%spline%coef(2) * L + eleinfo%spline%coef(3) * L**2
    dL = dspline_len(0.0_rp, L, eleinfo%spline, dtheta_L) ! = Ls - L
    Ls = L + dL
    zz = L / (2 * csr%gamma2) + dL
    k%I_int_csr = -csr%kick_factor * ((g_bend * Ls/2)**2 - log(2 * csr%gamma2 * zz / Ls) / csr%gamma2)
    ! Now add on rest of integral
    k%I_int_csr = k%I_int_csr + k%I_csr * (z - zz)
  endif

else
  kick1%I_int_csr = (kick1%I_csr + csr%kick1(i_bin-1)%I_csr) * csr%dz_slice / 2
endif

end subroutine I_csr

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine image_charge_kick_calc (kick1, csr) 
!
! Routine to calculate the image charge kick.
!
! Input:
!   kick1    -- csr_kick1_struct: 
!   csr      -- csr_struct:
!
! Output:
!   kick1             -- csr_kick1_struct: 
!     %image_kick_csr -- real(rp): Image charge kick.
!-

subroutine image_charge_kick_calc (kick1, csr)

implicit none

type (csr_kick1_struct), target :: kick1
type (csr_struct), target :: csr
type (csr_kick1_struct), pointer :: k
type (spline_struct), pointer :: sp

real(rp) N_vec(3), G_vec(3), B_vec(3), Bp_vec(3), NBp_vec(3), NBpG_vec(3), rad_cross_vec(3)
real(rp) z, sin_phi, cos_phi, OneNBp, OneNBp3, radiate, coulomb1, theta, g_bend

!

k => kick1
sp => csr%eleinfo(k%ix_ele_source)%spline

g_bend = -spline1(sp, k%s_chord_source, 2) / sqrt(1 + spline1(sp, k%s_chord_source, 1)**2)**3
theta = k%floor_s%theta

Bp_vec = csr%beta * [sin(theta), 0.0_rp, cos(theta) ]             ! beta vector at source point
G_vec = csr%beta**2 * g_bend * [-cos(theta), 0.0_rp, sin(theta) ] ! Acceleration vector

k%image_kick_csr = 0
z = k%dz_particles

N_vec = (csr%floor_k%r - k%floor_s%r) / k%L

theta = csr%floor_k%theta
B_vec = [sin(theta), 0.0_rp, cos(theta)]        ! Beta vector at kicked point


OneNBp = 1 - sum(N_vec * Bp_vec)
OneNBp3 = OneNBp**3

NBp_vec = N_vec - Bp_vec
NBpG_vec = cross_product(NBp_vec, G_vec)
rad_cross_vec = cross_product(N_vec, NBpG_vec)

radiate  = dot_product (B_vec, rad_cross_vec) / (k%L * OneNBp3)
coulomb1 = dot_product (B_vec, NBp_vec) / (csr%gamma2 * k%L**2 * OneNBp3)
kick1%image_kick_csr = csr%kick_factor * (radiate + coulomb1)

end subroutine image_charge_kick_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine csr_kick_calc (csr, particle)
!
! Routine to calculate the longitudinal coherent synchrotron radiation kick.
!
!   csr  -- csr_struct: 
!   particle -- coord_struct: Particle to kick.
!     %vec(6) -- Initial particle energy.
!
! Output:
!   particle -- Coord_struct: Particle to kick.
!     %vec(6) -- Final particle energy.
!-

subroutine csr_kick_calc (csr, particle)

use da2_mod

implicit none

type (csr_struct), target :: csr
type (coord_struct), target :: particle
type (csr_bunch_slice_struct), pointer :: slice

real(rp) zp, r1, r0, dz, dpz, kx, ky, f0, f, beta0, dpx, dpy, x, y, f_tot
real(rp), pointer :: vec(:)
integer i, i0, i_del

! We use a weighted average between %kick1(j)%I_csr and %kick1(j+1)%I_csr
! so that the integral varies smoothly as a function of particle%vec(5).

if (.not. allocated(csr%slice)) return  ! True if kicks are turned off.

zp = particle%vec(5)
i0 = int((zp - csr%slice(1)%z_center) / csr%dz_slice) + 1
r1 = (zp - csr%slice(i0)%z_center) / csr%dz_slice
r0 = 1 - r1
vec => particle%vec

if (r1 < 0 .or. r1 > 1 .or. i0 < 1 .or. i0 >= csr_param%n_bin) then
  print *, 'CSR INTERNAL ERROR!'
  if (global_com%exit_on_error) call err_exit
endif

vec(6) = vec(6) + r0 * csr%slice(i0)%kick_csr + r1 * csr%slice(i0+1)%kick_csr

! Longitudinal space charge

if (csr_param%lsc_component_on) then

  if (csr_param%lsc_kick_transverse_dependence) then
    x = particle%vec(1)
    y = particle%vec(3)
    dpz = 0

    slice => csr%slice(i0)
    if (slice%coef_lsc_plus(0,0) /= 0)  dpz = dpz + r0 / da2_evaluate(slice%coef_lsc_plus, (x-csr%x0_bunch)**2, (y-csr%y0_bunch)**2)
    if (slice%coef_lsc_minus(0,0) /= 0) dpz = dpz + r0 / da2_evaluate(slice%coef_lsc_minus, (x-csr%x0_bunch)**2, (y-csr%y0_bunch)**2)

    slice => csr%slice(i0+1)
    if (slice%coef_lsc_plus(0,0) /= 0)  dpz = dpz + r1 / da2_evaluate(slice%coef_lsc_plus, (x-csr%x0_bunch)**2, (y-csr%y0_bunch)**2)
    if (slice%coef_lsc_minus(0,0) /= 0) dpz = dpz + r1 / da2_evaluate(slice%coef_lsc_minus, (x-csr%x0_bunch)**2, (y-csr%y0_bunch)**2)

    vec(6) = vec(6) + dpz

  else
    vec(6) = vec(6) + r0 * csr%slice(i0)%kick_lsc + r1 * csr%slice(i0+1)%kick_lsc
  endif

endif

! Must update beta and z due to the energy change

beta0 = particle%beta
call convert_pc_to ((1+vec(6))* particle%p0c, particle%species, beta = particle%beta)
vec(5) = vec(5) * particle%beta / beta0

! Transverse space charge. Like the beam-beam interaction but in this case everything is going
! in the same general direction. In this case, the kicks due to the electric and magnetic fields
! tend to cancel instead of add as in the bbi case.

if (csr_param%tsc_component_on) then
  f0 = csr%kick_factor * csr%actual_track_step * r_e / (twopi * &
           csr%dz_slice * csr%rel_mass * e_charge * abs(charge_of(particle%species)) * csr%gamma**3)

  dpx = 0;  dpy = 0

  slice => csr%slice(i0)
  if (slice%sig_x /= 0 .and. slice%sig_y /= 0) then
    call bbi_kick ((vec(1)-slice%x0)/slice%sig_x, (vec(3)-slice%y0)/slice%sig_y, slice%sig_y/slice%sig_x, kx, ky)
    f = f0 * r0 * slice%charge / (slice%sig_x + slice%sig_y)
    ! The kick is negative of the bbi kick. That is, the kick is outward.
    dpx = -kx * f
    dpy = -ky * f
  endif

  slice => csr%slice(i0+1)
  if (slice%sig_x /= 0 .and. slice%sig_y /= 0) then
    call bbi_kick ((vec(1)-slice%x0)/slice%sig_x, (vec(3)-slice%y0)/slice%sig_y, slice%sig_y/slice%sig_x, kx, ky)
    f = f0 * r1 * slice%charge / (slice%sig_x + slice%sig_y)
    ! The kick is negative of the bbi kick. That is, the kick is outward.
    dpx = dpx - kx * f   
    dpy = dpy - ky * f
  endif

  ! 1/2 of the kick is electric and this leads to an energy change.
  ! The formula for vec(6) assumes that dpx and dpy are small so only the linear terms are kept.

  vec(2) = vec(2) + dpx
  vec(4) = vec(4) + dpy
  vec(6) = vec(6) + (dpx*vec(2) + dpy*vec(4)) / (2 * (1 + vec(6)))

endif

end subroutine csr_kick_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function dspline_len (s_chord0, s_chord1, spline, dtheta_ref) result (dlen)
!
! Routine to calculate the difference in length between the spline curve length and a referece line.
! Referece line is centroid chord (referece system of the spline) rotated by dtheta_ref.
!
! Input:
!   s_chord0    -- real(rp): Start position along centroid chord.
!   s_chord1    -- real(rp): Stop position along central_chord.
!   spline      -- spline_struct: Spline of x-position as a function of s.
!   dtheta_ref  -- real(rp), optional: angle to rotate the reference line from the centroid chord.
!                    Default is 0.
!
! Output:
!   dlen -- real(rp): L_spline - L_chord
!-

function dspline_len (s_chord0, s_chord1, spline, dtheta_ref) result (dlen)

implicit none

type (spline_struct) spline

real(rp) s_chord0, s_chord1, dlen
real(rp), optional :: dtheta_ref
real(rp) c(0:3), s0, ds

! x' = c(1) + 2*c2*s + 3*c3*s^2
! dlen = Integral: x'^2/2 ds

c = spline%coef
s0 = s_chord0
ds = s_chord1 - s_chord0

if (present(dtheta_ref)) then
  c(1) = c(1) - dtheta_ref
endif

dlen = (ds / 2) * ( &
        c(1)**2 + &
        (2*s0 + ds) * 2*c(1)*c(2) + &
        (3*s0*s0 + 3*s0*ds + ds*ds) * (6*c(1)*c(3) + 4*c(2)**2) / 3 + &
        (4*s0**3 + 6*s0*s0*ds + 4*s0*ds*ds + ds**3) * 3*c(2)*c(3) + &
        (5*s0**4 + 10*s0**3*ds + 10*s0*s0*ds*ds + 5*s0*ds**3 + ds**4) * 9*c(3)**2 &
       )

end function dspline_len

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function s_ref_to_s_chord (s_ref, eleinfo) result (s_chord)
!
! Routine to calculate s_chord given s_ref.
!
! Input:
!   s_ref     -- real(rp): s-position along element ref coords.
!   eleinfo     -- csr_ele_info_struct: Element info
!
! Output:
!   s_chord   -- real(rp): s-posijtion along centroid chord.
!-

function s_ref_to_s_chord (s_ref, eleinfo) result (s_chord)

implicit none

type (csr_ele_info_struct), target :: eleinfo
type (ele_struct), pointer :: ele

real(rp) s_ref, s_chord, dtheta, dr(3), x, g, t

!

ele => eleinfo%ele
g = ele%value(g$)

if (ele%key == sbend$ .and. abs(g) > 1d-5) then
  dtheta = eleinfo%floor0%theta - eleinfo%e_floor0%theta
  dr = eleinfo%e_floor0%r - eleinfo%floor0%r
  t = eleinfo%e_floor0%theta + pi/2
  x = dr(1) * cos(t) + dr(3) * sin(t)
  s_chord = abs(ele%value(rho$)) * atan2(s_ref * cos(dtheta), abs(ele%value(rho$) + x + s_ref * sin(dtheta)))

else
  s_chord = s_ref * eleinfo%L_chord /ele%value(l$)
endif

end function s_ref_to_s_chord

end module
