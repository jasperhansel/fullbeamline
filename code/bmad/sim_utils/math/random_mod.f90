!+
! Module random_mod
!
! Module for random number generation.
!-

module random_mod

use precision_def
use physical_constants
use output_mod
use sim_utils_interface

!

integer, private, parameter :: kr4b = selected_int_kind(9)
integer(kr4b), private, parameter :: im_nr_ran = 2147483647
integer(i4_b), private, parameter :: sobseq_maxbit = 30, sobseq_maxdim = 6

! common variables for random number generator.

integer, parameter :: pseudo_random$ = 1, quasi_random$ = 2
integer, parameter :: quick_gaussian$ = 3, exact_gaussian$ = 4

type random_state_struct
  integer(kr4b) :: ix = -1, iy = -1
  logical :: number_stored = .false.
  real(rp) :: h_saved = 0
  integer :: engine = pseudo_random$
  ! Params
  integer :: seed = 0
  real(sp) :: am = 0
  integer :: gauss_converter = exact_gaussian$
  real(rp) :: gauss_sigma_cut = -1
  integer(i4_b) :: in_sobseq = 0
  integer(i4_b) :: ix_sobseq(sobseq_maxdim) = 0
  real(rp) :: x_sobseq(sobseq_maxdim) = 0
end type

type (random_state_struct), private, target, save :: ran_state_dflt

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine ran_gauss (harvest, ran_state)
!
! Routine to return a gaussian distributed random number with unit sigma.
! This routine uses the same algorithm as gasdev from Numerical Recipes.
!
! Note: ran_gauss is an overloaded name for:
!     ran_gauss_scalar   ! harvest is a scalar
!     ran_gauss_vector   ! harvest is a 1-D array.
!
! Note: Use ran_seed_put for initialization.
! Note: Use ran_engine to set which random number generator to use.
! Note: Use ran_gauss_converter to set which conversion routine to use.
!
! Module needed:
!   use random_mod
!
! Input:
!   ran_state -- random_state_struct, optional: Internal state.
!                     See the ran_seed_put documentation for more details.
!
! Output:
!   harvest    -- Real(rp): Random number. 
! Or
!   harvest(:) -- Real(rp): Random number array. 
!                  For quasi_random$ numbers, the array size must be less than 6.
!-

interface ran_gauss
  module procedure ran_gauss_scalar
  module procedure ran_gauss_vector
end interface

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine ran_uniform (harvest, ran_state)
!
! Routine to return a random number uniformly distributed in the 
! interval [0, 1]. This routine uses the same algorithm as ran or sobseq
! from Numberical Recipes in Fortran90.
! See ran_engine.
!
! Note: ran_uniform is an overloaded name for:
!     ran_uniform_scalar   ! harvest is a scalar
!     ran_uniform_vector   ! harvest is a 1-D array.
!
! Note: Use ran_seed_put for initialization.
! Note: Use ran_engine to set which random number generator to use.
!
! Modules needed:
!   use random_mod
!
! Input:
!   ran_state -- random_state_struct, optional: Internal state.
!                     See the ran_seed_put documentation for more details.
!
! Output:
!   harvest    -- Real(rp): Random number. 
! Or
!   harvest(:) -- Real(rp): Random number array. 
!                  For quasi_random$ numbers the array size must be less than 6.
!-

interface ran_uniform
  module procedure ran_uniform_scalar
  module procedure ran_uniform_vector
end interface

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_scalar (harvest, ran_state, index_quasi)
!
! Routine to return a gaussian distributed random number with unit sigma.
! See ran_gauss for more details.
!
! Note: The index_quasi argument is used internally for the quasi-random number generator.
!-

subroutine ran_gauss_scalar (harvest, ran_state, index_quasi)

use nr, only: erf_s, erf

implicit none

type (random_state_struct), optional, target :: ran_state
type (random_state_struct), pointer :: r_state

integer, parameter :: sigma_max = 8, n_pts_per_sigma = 25
integer, parameter :: max_g = sigma_max * n_pts_per_sigma

real(rp), intent(out) :: harvest
real(rp) a(2), v1, v2, r, sigma_cut, fac
real(rp), save :: g(0:max_g) = 0

integer, optional :: index_quasi
integer i, ss, ix

! quasi-random must use the quick_gaussian since the exact_gaussian can
! use several uniform random numbers to generate a single Gaussian random number.
! This invalidates the algorithm used to generate a quasi-random Gaussian vector.

! ran_state%g is the normalized error function and maps from the 
! interval [0, 0.5] to [0, infinity].

if (present(ran_state)) then
  r_state => ran_state
else
  r_state => ran_state_dflt
endif

if (r_state%engine == quasi_random$ .or. r_state%gauss_converter == quick_gaussian$) then

  ! Init g

  if (g(1) == 0) then
    fac = 2 * erf_s (sigma_max/sqrt_2)
    do i = 0, max_g-1
      g(i) = erf_s (i / (n_pts_per_sigma * sqrt_2)) / fac
    enddo
    g(max_g) = 0.50000000001_rp
  endif

  !

  sigma_cut = sigma_max
  if (r_state%gauss_sigma_cut > 0) sigma_cut = min(r_state%gauss_sigma_cut, sigma_cut)

  call ran_uniform_scalar (r, r_state, index_quasi)
  if (r > 0.5) then
    r = r - 0.5
    ss = 1
  else
    r = 0.5 - r
    ss = -1
  endif
  call bracket_index(g, 0, max_g, r, ix) 
  harvest = (ix + (r - g(ix)) / (g(ix+1) - g(ix))) * ss / n_pts_per_sigma
  if (harvest >  sigma_cut) harvest =  sigma_cut
  if (harvest < -sigma_cut) harvest = -sigma_cut
  return

endif

! Loop until we get an acceptable number

do 

  ! If we have a stored value then just use it

  if (r_state%number_stored) then
    r_state%number_stored = .false.
    harvest = r_state%h_saved
    if (r_state%gauss_sigma_cut < 0 .or. abs(harvest) < r_state%gauss_sigma_cut) return
  endif

  ! else we generate a number

  do
    call ran_uniform(a, r_state)
    v1 = 2*a(1) - 1
    v2 = 2*a(2) - 1
    r = v1**2 + v2**2
    if (r > 0 .and. r < 1) exit   ! In unit circle
  enddo

  r = sqrt(-2*log(r)/r)
  r_state%h_saved = v2 * r
  r_state%number_stored = .true.

  harvest = v1 * r
  if (r_state%gauss_sigma_cut < 0 .or. abs(harvest) < r_state%gauss_sigma_cut) return

enddo

end subroutine ran_gauss_scalar

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_vector (harvest, ran_state)
!
! Routine to return a gaussian distributed random number with unit sigma.
! See ran_gauss for more details.
!-

subroutine ran_gauss_vector (harvest, ran_state)

implicit none

type (random_state_struct), optional, target :: ran_state

real(rp), intent(out) :: harvest(:)
integer i

!

do i = 1, size(harvest)
  call ran_gauss_scalar (harvest(i), ran_state, i)
enddo

end subroutine ran_gauss_vector

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_engine (set, get, ran_state)
!
! Routine to set what random number generator algorithm is used.
! If this routine is never called then pseudo_random$ is used.
! With sobseq quasi-random numbers the maximum dimension is 6.
!
! Modules needed:
!   use random_mod
! 
! Input:
!   set -- Character(*), optional: Set the random number engine. Possibilities are:
!                'pseudo' -> Uses ran from Numerical Recipies (F90).
!                'quasi'  -> Uses sobseq from Numerical Recipes.
!   get -- Character, optional: Get the current (before any set) random number engine. 
!   ran_state -- random_state_struct, optional: Internal state.
!                     See the ran_seed_put documentation for more details.
!-

subroutine ran_engine (set, get, ran_state)

implicit none

type (random_state_struct), optional, target :: ran_state
type (random_state_struct), pointer :: r_state

character(*), optional :: set, get
character(16) :: r_name = 'ran_engine'

! Set state to use

if (present(ran_state)) then
  r_state => ran_state
else
  r_state => ran_state_dflt
endif

! get

if (present (get)) then
  select case (r_state%engine)
  case (pseudo_random$)
    get = 'pseudo'
  case (quasi_random$)
    get = 'quasi'
  end select
endif

! set

if (present(set)) then
  select case (set)
  case ('pseudo')
    r_state%engine = pseudo_random$
  case ('quasi')
    r_state%engine = quasi_random$
    r_state%number_stored = .false.
  case default
    call out_io (s_error$, r_name, 'BAD RANDOM NUMBER ENGINE NAME: ' // set)
  end select
endif

end subroutine ran_engine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_gauss_converter (set, set_sigma_cut, get, get_sigma_cut, ran_state)
!
! Routine to set what conversion routine is used for converting
! uniformly distributed random numbers to Gaussian distributed random numbers.
!
! If this routine is not called then exact_gaussian$ is used.
!
! exact_gaussian$ is a straight forward converter as explained in Numerical recipes.
!
! quick_gaussian$ is a quick a dirty approximation with a cutoff so that no 
! numbers will be generated beyound what is set for sigma_cut. 
!
! A negative sigma_cut means that the exact_gaussian$ will not be limited
! and the quick_gaussian$ will use a default of 10.0
!
! Note: Because of technical issues, when using the quasi_random$ number generator
! (see the ran_engine routine), the quick_gaussian$ method will automatically be 
! used independent of what was set with this routine.
!
! Modules needed:
!   use random_mod
! 
! Input:
!   set -- Character(*), optional: Set the random number engine. Possibilities are:
!             'exact'
!             'quick'  (Old deprecated: 'limited')
!   set_sigma_cut -- Real(rp), optional: Sigma cutoff. Initially: sigma_cut = -1.
!   ran_state -- random_state_struct, optional: Internal state.
!                     See the ran_seed_put documentation for more details.
!
! Output:
!   get -- Character(*), optional: Get the current (before any set) gaussian converter.
!   get_sigma_cut -- Real(rp), optional: Get the current (before andy set) sigma cutoff.
!-

subroutine ran_gauss_converter (set, set_sigma_cut, get, get_sigma_cut, ran_state)

implicit none

type (random_state_struct), optional, target :: ran_state
type (random_state_struct), pointer :: r_state

real(rp), optional :: set_sigma_cut, get_sigma_cut

character(*), optional :: set, get
character(16) :: r_name = 'ran_gauss_converter'

! Set state to use

if (present(ran_state)) then
  r_state => ran_state
else
  r_state => ran_state_dflt
endif

! get converter

if (present (get)) then
  select case (r_state%gauss_converter)
  case (quick_gaussian$)
    get = 'quick'
  case (exact_gaussian$)
    get = 'exact'
  end select
endif

! get sigma_cut

if (present(get_sigma_cut)) then
  get_sigma_cut = r_state%gauss_sigma_cut
endif

! set converter

if (present(set)) then
  select case (set)
  case ('quick', 'limited')
    r_state%gauss_converter = quick_gaussian$
  case ('exact')
    r_state%gauss_converter = exact_gaussian$
  case default
    call out_io (s_error$, r_name, 'BAD RANDOM NUMBER GAUSS_CONVERTER NAME: ' // set)
  end select
endif

! set sigma_cut

if (present(set_sigma_cut)) then
  r_state%gauss_sigma_cut = set_sigma_cut
endif

end subroutine ran_gauss_converter

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_seed_put (seed, ran_state)
!
! Routine to seed a random number generator. 
!
! Independent random number "generators" are implemented by supplying a ran_state argument.
! The actual ran_state argument for independent generators must be different.
! If the ran_state arg is not present, the default generator is used.
!
! Each generator is independent of the others. That is the sequence of "random" numbers
! generated with a given generator is unaffected by usage of the other generators.
! Multiple generators are useful in cases where you want to maintain the same random number 
! sequence in one part of the program independent of other parts of the program or
! when parallel computing is used.
!
! For example, suppose you were tracking particles whose initial position was determined
! with a random number generator. Additionally suppose that these particles could interact 
! with the residual gas background and that this interaction is modeled using a random
! number generator. If you want the initial particle positions to be independent of whether
! you were simulating the particle-gas interaction or not, you could use different generators 
! for the particle initialization and the particle-gas interaction.
!
! If a program never calls ran_seed_put, or ran_seed_put is called with seed = 0,
! the system clock will be used to generate the seed.
!
! Note: The seed is only used with the pseudo_random$ engine.
!
! Note: Use the subroutine ran_seed_get(seed) to get the seed used.
!
! Modules needed:
!   use random_mod
!
! Intput:
!   seed  -- Integer, optional: Seed number. If seed = 0 then a 
!              seed will be choosen based upon the system clock.
!   ran_state -- random_state_struct, optional: Internal state.
!-

subroutine ran_seed_put (seed, ran_state)

implicit none

type (random_state_struct), pointer :: r_state
type (random_state_struct), optional, target :: ran_state

integer, optional :: seed
integer v(10)

real(rp) dum(2)

! Set state to use

if (present(ran_state)) then
  r_state => ran_state
else
  r_state => ran_state_dflt
endif

! init

r_state%in_sobseq = 0
r_state%ix_sobseq = 0

r_state%am = nearest(1.0,-1.0) / im_nr_ran

if (seed == 0) then
  call date_and_time (values = v)
  r_state%seed = v(2) + 11*v(3) + 111*v(5) + 1111*v(6) + 11111*v(7) + 111111*v(8)
else
  r_state%seed = seed
endif

r_state%iy = ior(ieor(888889999, abs(r_state%seed)), 1)
r_state%ix = ieor(777755555, abs(r_state%seed))

r_state%number_stored = .false.

end subroutine ran_seed_put

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_seed_get (seed, ran_state)
! 
! Routine to return the seed used for the random number generator.
!
! Note: The internal state can be used to put the pseudo-random
! number generator into a known state. See ran_seed_put
!
! Modules needed:
!   use random_mod
!
! Input:
!   ran_state -- random_state_struct, optional: Internal state.
!                     See the ran_seed_put documentation for more details.
!
! Output:
!   seed      -- Integer, optional: Random number seed used.
!   ran_state -- random_state_struct, optional: Internal state.
!-

subroutine ran_seed_get (seed, ran_state)

implicit none

type (random_state_struct), optional, target :: ran_state
type (random_state_struct), pointer :: r_state

integer, optional :: seed

! Set state to use

if (present(ran_state)) then
  r_state => ran_state
else
  r_state => ran_state_dflt
endif

!

if (present(seed)) seed = r_state%seed

end subroutine ran_seed_get 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_default_state (set_state, get_state)
!
! Routine to set or get the state of the default random number generator.
! See the ran_seed_put documentation for more details
!
! Module needed:
!   use random_mod
!
! Input:
!   set_state -- random_state_struct, optional: State to set the default generator to.
!
! Output:
!   get_state -- random_state_struct, optional: Returns the state of the default generator.
!-

subroutine ran_default_state (set_state, get_state)

implicit none

type (random_state_struct), optional :: set_state, get_state

!

if (present(set_state)) ran_state_dflt = set_state
if (present(get_state)) get_state = ran_state_dflt

end subroutine ran_default_state

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_uniform_scalar (harvest, ran_state, index_quasi)
!
! Routine to return a random number uniformly distributed in the 
! interval [0, 1]. 
! See ran_uniform for more details.
!
! Note: The index_quasi argument is used internally for the quasi-random number generator.
!-

subroutine ran_uniform_scalar (harvest, ran_state, index_quasi)

implicit none

type (random_state_struct), pointer :: r_state
type (random_state_struct), optional, target :: ran_state

real(rp), intent(out) :: harvest

integer(kr4b) k, ix_q
integer, optional :: index_quasi


integer(kr4b), parameter :: ia = 16807
integer(kr4b), parameter :: iq = 127773, ir = 2836

character :: r_name = 'ran_uniform_scalar'

! Set state to use

if (present(ran_state)) then
  r_state => ran_state
else
  r_state => ran_state_dflt
endif

! If r_state%iy < 0 then the random number generator has never bee initialized.

if (r_state%iy < 0) call ran_seed_put(r_state%seed)

! quasi-random

if (r_state%engine == quasi_random$) then
  ix_q = integer_option(1, index_quasi)
  if (ix_q == 1) call super_sobseq (r_state%x_sobseq, r_state)
  if (ix_q > sobseq_maxdim) then
    call out_io (s_error$, r_name, 'NUMBER OF DIMENSIONS WANTED IS TOO LARGE!')
    if (global_com%exit_on_error) call err_exit
  endif
  harvest = r_state%x_sobseq(ix_q)
  return
endif

! Pseudo-random
! Marsaglia shift sequence with period 2^32 - 1.

r_state%ix = ieor(r_state%ix, ishft(r_state%ix, 13)) 
r_state%ix = ieor(r_state%ix, ishft(r_state%ix, -17))
r_state%ix = ieor(r_state%ix, ishft(r_state%ix, 5))
k = r_state%iy/iq         ! Park-Miller sequence by Schrage's method,
r_state%iy = ia*(r_state%iy - k*iq) - ir * k       ! period 2^31-2.

if (r_state%iy < 0) r_state%iy = r_state%iy + im_nr_ran

! Combine the two generators with masking to ensure nonzero value.

harvest = r_state%am * ior(iand(im_nr_ran, ieor(r_state%ix, r_state%iy)), 1) 

end subroutine ran_uniform_scalar

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ran_uniform_vector (harvest, ran_state)
!
! Routine to return a vector of random numbers uniformly distributed in the 
! interval [0, 1]. 
! See ran_uniform for more details.
!-

subroutine ran_uniform_vector (harvest, ran_state)

implicit none

type (random_state_struct), optional, target :: ran_state

real(rp), intent(out) :: harvest(:)
integer i

!

do i = 1, size(harvest)
  call ran_uniform_scalar (harvest(i), ran_state, i)
enddo

end subroutine ran_uniform_vector

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine super_sobseq (x, ran_state)
!
! Routine patterened after sobseq in Numerical Recipes.
! Difference is that this version has an argument for the internal state.
!
! Module needed:
!   use random_mod
!
! Input:
!   ran_state -- random_state_struct, optional: Generator state.
!                     See the ran_seed_put documentation for more details.
!
! Output: 
!   x(:)      -- real(dp): Random vector.
!   ran_state -- random_state_struct, optional: Generator state.
!-

subroutine super_sobseq (x, state)

implicit none

type (random_state_struct), optional, target :: state
type (random_state_struct), pointer :: r_state

real(rp), dimension(:), intent(out) :: x
real(dp), parameter :: fac=1.0_dp/2.0_dp**sobseq_maxbit

integer(i4_b) :: im, j
integer(i4_b), parameter :: ip(sobseq_maxdim) = [0,1,1,2,1,4], mdeg(sobseq_maxdim) = [1,2,3,3,4,4] 
! Note: If sobseq_maxbit or sobseq_maxdim are changed, iv needs to be recomputed by
! running the original sobseq routine.
integer(i4_b), parameter :: iv(sobseq_maxdim*sobseq_maxbit) = [ &
											536870912, 536870912, 536870912, 536870912, 536870912, &
											536870912, 805306368, 268435456, 805306368, 805306368, &
											268435456, 268435456, 671088640, 939524096, 939524096, &
											402653184, 402653184, 671088640, 1006632960, 738197504, &
											335544320, 1006632960, 872415232, 603979776, 570425344, &
											436207616, 234881024, 167772160, 838860800, 100663296, &
											855638016, 1023410176, 721420288, 285212672, 150994944, &
											385875968, 713031680, 562036736, 411041792, 713031680, &
											763363328, 1031798784, 1069547520, 331350016, 616562688, &
											566231040, 88080384, 465567744, 538968064, 975175680, &
											920649728, 853540864, 941621248, 497025024, 808452096, &
											756023296, 1062207488, 489684992, 605028352, 198180864, &
											673710080, 431489024, 381157376, 952631296, 706215936, &
											898105344, 1010565120, 1072431104, 258736128, 208928768, &
											1026818048, 804519936, 572653568, 540672000, 771883008, &
											316801024, 531759104, 864944128, 858980352, 271384576, &
											453181440, 758317056, 206110720, 954400768, 715816960, &
											941195264, 545488896, 550076416, 361594880, 238256128, &
											1073725440, 742375424, 817971200, 813154304, 559235072, &
											590921728, 536879104, 438312960, 954261504, 417505280, &
											301998080, 328081408, 805318656, 1024462848, 340963328, &
											1009913856, 419434496, 685969408, 671098880, 565721088, &
											238651392, 172697600, 897587200, 640919552, 1006648320, &
											334244864, 732843008, 297131008, 826291200, 121168896, &
											570434048, 976886272, 417426944, 704744960, 169882112, &
											361637376, 855651072, 757939456, 609285376, 553894656, &
											756025600, 1071843072, 713042560, 429498752, 909831040, &
											847291520, 127413632, 464754048, 1069563840, 1073730496, &
											1068349120, 499194688, 947127616, 486080448, 538976288, &
											536877600, 383778848, 954376224, 663894048, 137247648, &
											808464432, 268451088, 256901168, 204607536, 676930576, &
											875760848, 673720360, 939532728, 783810616, 306915352, &
											1066773016, 775659528, 1010580540, 738202604, 460062740, &
											766893116, 476151092, 856477748, 572662306, 436222522, &
											537001998, 536972810, 229785522, 1000343598, 858993459, &
											1023421741, 805503019, 805552913, 357112905, 214971443]

character(16), parameter :: r_name = 'super_sobseq'

! which state to use?

if (present(state)) then
	r_state => state
else
	r_state => ran_state_dflt
endif

! calc

im = r_state%in_sobseq
do j = 1, sobseq_maxbit
  if (.not. btest(im,0)) exit
  im = im/2
end do

if (j > sobseq_maxbit) then
  call out_io (s_fatal$, r_name, 'SOBSEQ_MAXBIT TOO SMALL')
  if (global_com%exit_on_error) call err_exit
  return
endif

im = (j-1) * sobseq_maxdim
j = min(size(x), sobseq_maxdim)
r_state%ix_sobseq(1:j) = ieor(r_state%ix_sobseq(1:j),iv(1+im:j+im))
x(1:j) = r_state%ix_sobseq(1:j) * fac
r_state%in_sobseq = r_state%in_sobseq + 1

end subroutine super_sobseq

end module

