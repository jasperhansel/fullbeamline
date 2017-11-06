!+
! Subroutine twiss_at_start (lat, status, ix_branch)
!
! Subroutine to calculate, for a circular machine, the closed 1-turn 
! solution for the Twiss parameters at the start of the lat.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat         -- lat_struct: Lat
!   ix_branch   -- Integer, option: Branch to use. Default is 0 (main branch).
!
! Output:
!   lat         -- Lat_struct: Lattice with twiss parameters computed.
!     %param%t1_no_RF --  Note: Only the linear part is computed.
!     %ele(0)%a      -- "a" mode Twiss parameters at the start of the lat.
!     %ele(0)%b      -- "b" mode Twiss parameters at the start of the lat.
!     %ele(0)%c_mat  -- Coupling matrix.
!     %a%tune         -- Fractional part of the tune in radians
!     %b%tune         -- Fractional part of the tune in radians
!     %param%stable   -- Set true or false.
!     %param%unstable_factor -- unstable growth rate (= 0 if stable)
!   status      -- Integer, optional: Calculation status:
!                       ok$, in_stop_band$, unstable$, or non_symplectic$
!-

subroutine twiss_at_start (lat, status, ix_branch)

use bookkeeper_mod, except_dummy => twiss_at_start

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch

real(rp) eta_vec(4), t0_4(4,4), mat6(6,6), map0(4), m56
real(rp), allocatable :: on_off_state(:)

integer, optional, intent(in) :: ix_branch
integer, optional, intent(out) ::status
integer i, j, n, iu, n_lines, stat

logical :: debug = .false. 
logical saved_state

character(200), allocatable :: lines(:)

! init one turn. T0 is the transverse part of the matrix

call mat_make_unit (t0_4)       ! form unit matrix
eta_vec = 0
map0 = 0
m56 = 0

branch => lat%branch(integer_option(0, ix_branch))

! Propagate the transfer map around branch. 
! Since the RF is taken to be off we use a trick so we only have to multiply
! 4x4 matrices.

if (debug) then
  iu = lunget()
  open (iu, file = 'twiss_at_start.dat')
endif

!

call set_on_off (rfcavity$, lat, off_and_save$, use_ref_orb = .true., &
                              ix_branch = branch%ix_branch, saved_values = on_off_state)

do n = 1, branch%n_ele_track
  ele => branch%ele(n)
  m56 = m56 + ele%mat6(5,6) + dot_product(ele%mat6(5,1:4), eta_vec)
  eta_vec = matmul (ele%mat6(1:4,1:4), eta_vec) + ele%mat6(1:4,6)
  map0 = matmul (ele%mat6(1:4,1:4), map0) + ele%vec0(1:4)
  t0_4 = matmul (ele%mat6(1:4,1:4), t0_4)
  if (debug) then
    write (iu, *) '!------------------------------------', n
    call type_ele (ele, .false., 0, .false., 0, lines = lines, n_lines = n_lines)
    do i = 1, n_lines
      write (iu, '(a)') trim(lines(i))
    enddo
    write (iu, *) 'Symplectic Check:', mat_symp_error(t0_4)
    call mat_type (t0_4, iu)
    do i = 1, 4
      write (iu, '(es20.12, 2x, es20.12)') eta_vec(i), map0(i)
    enddo
  endif
enddo

call set_on_off (rfcavity$, lat, restore_state$, use_ref_orb = .true., &
                              ix_branch = branch%ix_branch, saved_values = on_off_state)

if (debug) close (iu)

! Put 1-turn matrix into branch%param%t1_no_RF

call mat_make_unit (mat6)
mat6(1:4,1:4) = t0_4

call mat6_dispersion (eta_vec, mat6) ! dispersion to %mat6
mat6(5,6) = m56
branch%param%t1_no_RF = mat6

! Compute twiss parameters

call twiss_from_mat6 (mat6, branch%ele(1)%map_ref_orb_in%vec, branch%ele(0), &
                branch%param%stable, branch%param%unstable_factor, stat, .true.)
if (present(status)) status = stat

lat%a%tune = branch%ele(0)%a%phi
lat%b%tune = branch%ele(0)%b%phi

branch%ele(0)%a%phi = 0
branch%ele(0)%b%phi = 0

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine mat6_dispersion (m_i6, mat6)
!
! Subroutine to set the mat6(5, 1:4) terms given the vector mat6(1:4, 6)
! which is a measure of the dispersion.
!
! Input:
!   m_i6(4)   -- Real(rp): mat6(1:4, 6) components.
!   mat6(6,6) -- Real(rp): Matrix with 4x4 x-y submatrix already made.
!
! Output:
!   mat6(6,6) -- Real(rp): mat6(5, 1:4) components set. 
!-

subroutine mat6_dispersion (m_i6, mat6)

implicit none

real(rp), intent(inout) :: mat6(:,:)
real(rp), intent(in) :: m_i6(:)

real(rp) vec4(4)

!

mat6(1:4, 6) = m_i6(1:4)

vec4(1) = -m_i6(2)
vec4(2) =  m_i6(1)
vec4(3) = -m_i6(4)
vec4(4) =  m_i6(3)

mat6(5,1:4) = matmul (vec4, mat6(1:4,1:4))

end subroutine mat6_dispersion

end subroutine
