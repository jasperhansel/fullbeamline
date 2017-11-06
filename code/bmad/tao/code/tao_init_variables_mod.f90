module tao_init_variables_mod

use tao_mod

contains

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!+
! Subroutine tao_init_variables (var_file)
!
! Subroutine to initialize the tao variable structures.
! If var_file is not in the current directory then it 
! will be searched for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   var_file  -- Character(*): Tao variable initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_variables (var_file)

use tao_lattice_calc_mod
use tao_input_struct
use random_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_v1_var_input) v1_var
type (tao_v1_var_struct), pointer :: v1_var_ptr
type (tao_var_input) var(n_var_minn:n_var_maxx)

real(rp) default_weight        ! default merit function weight
real(rp) default_step          ! default "small" step size
real(rp) default_low_lim, default_high_lim, default_key_delta
real(rp), allocatable, save :: default_key_d(:)

integer ios, iu, i, j, j1, j2, k, ix, num
integer n, iostat, n_list
integer ix_min_var, ix_max_var, ix_ele, n_v1, n_v1_var_max

character(*) var_file
character(40) :: r_name = 'tao_init_variables'
character(40) name, universe, default_universe
character(40) default_merit_type, default_attribute
character(40) use_same_lat_eles_as
logical default_key_bound
character(200) file_name
logical, allocatable, save :: default_key_b(:)
character(100) line, search_for_lat_eles

logical err, free, gang
logical searching, limited, do_standard_setup
logical, allocatable, save :: dflt_good_unis(:), good_unis(:)

namelist / tao_var / v1_var, var, default_weight, default_step, default_key_delta, &
                    ix_min_var, ix_max_var, default_universe, default_attribute, &
                    default_low_lim, default_high_lim, default_merit_type, &
                    use_same_lat_eles_as, search_for_lat_eles, default_key_bound

!-----------------------------------------------------------------------
! Init

allocate (s%var(0))
allocate (s%v1_var(0))
s%n_var_used = 0
s%n_v1_var_used = 0

! Call the hook routine and then return if the standard init is not wanted.

call tao_hook_init_var(do_standard_setup) 
if (.not. do_standard_setup) return

! No var file then just setup the key table if needed.

if (var_file == '') then
  call tao_setup_key_table ()
  return
endif

! Standard var init:
! First open the var init file.

call out_io (s_blank$, r_name, '*Init: Opening Variable File: ' // file_name)
call tao_open_file (var_file, iu, file_name, s_fatal$)
if (iu == 0) then
  call out_io (s_fatal$, r_name, 'CANNOT OPEN VARIABLE INIT FILE: ' // var_file)
  call err_exit
endif

! Count how many v1_var definitions there are

allocate (default_key_b(100), default_key_d(100))
n = 0
do
  call set_logical_to_garbage(default_key_bound)
  default_key_delta = 0
  v1_var%name = ''
  read (iu, nml = tao_var, iostat = ios)
  if (ios > 0) then
    call out_io (s_error$, r_name, 'TAO_VAR NAMELIST READ ERROR.')
    rewind (iu)
    do
      read (iu, nml = tao_var)  ! force printing of error message
    enddo
  endif
  if (ios < 0 .and. v1_var%name == '') exit  ! Exit on end-of-file and no namelist read
  n = n + 1
  if (n >= size(default_key_b)) then
    call re_allocate (default_key_b, 2*size(default_key_b))
    call re_allocate (default_key_d, 2*size(default_key_d))
  endif
  call transfer_logical (default_key_bound, default_key_b(n))
  default_key_d(n) = default_key_delta
enddo

call tao_allocate_v1_var (n)
n_list = n

! Now fill in all the information

rewind (iu)

allocate (dflt_good_unis(lbound(s%u,1):ubound(s%u, 1)), good_unis(lbound(s%u,1):ubound(s%u,1)))

n_v1 = 0
do
  n_v1 = n_v1 + 1
  if (n_v1 > n_list) exit

  use_same_lat_eles_as = ''
  search_for_lat_eles  = ''
  v1_var%name        = ''         ! set default
  default_merit_type = 'limit'
  default_weight     = 0     ! set default
  default_step       = 0       ! set default
  default_attribute  = ''
  default_universe   = ''
  default_low_lim    = -1e30
  default_high_lim   = 1e30
  ix_min_var         = 1
  ix_max_var         = 0
  var%ele_name       = ''
  var%merit_type     = ''
  var%weight         = 0         ! set default
  var%step           = 0         ! set default
  var%attribute      = ''
  var%universe       = ''
  var%low_lim        = default_low_lim
  var%high_lim       = default_high_lim
  call set_logical_to_garbage (var%good_user)
  call transfer_logical (default_key_b(n_v1), var(0)%key_bound)
  var%key_bound = var(0)%key_bound
  var%key_delta      = default_key_d(n_v1)

  read (iu, nml = tao_var, iostat = ios)
  if (ios < 0 .and. v1_var%name == '') exit         ! exit on end-of-file
  call out_io (s_blank$, r_name, 'Init: Read tao_var namelist: ' // v1_var%name)

  ! Convert old format to new

  call str_upcase (default_attribute, default_attribute)
  if (var(0)%ele_name(1:7) == 'SEARCH:') then
    call string_trim(var(0)%ele_name(8:), search_for_lat_eles, ix)
  elseif (var(0)%ele_name(1:5) == 'SAME:') then
    call string_trim (var(0)%ele_name(6:), use_same_lat_eles_as, ix)
  endif

  ! Convert to upper case.

  do i = lbound(var, 1), ubound(var, 1)
    call str_upcase (var(i)%attribute, var(i)%attribute)
    call str_upcase (var(i)%ele_name, var(i)%ele_name)
  enddo

  if (v1_var%name == '') then
    call out_io (s_error$, r_name, 'FOUND TAO_VAR NAMELIST WITHOUT V1_VAR%NAME PARAMETER!')
    cycle
  endif

  ! Gang or clone?

  gang = .true.
  if (default_universe(1:5) == 'clone') then
    gang = .false.
    call string_trim (default_universe(6:), default_universe, ix)
  endif
  if (default_universe(1:4) == 'gang') then
    call string_trim (default_universe(5:), default_universe, ix)
  endif

  ! Read universe numbers

  if (default_universe == '*' .or. default_universe == '') then
    dflt_good_unis = .true.
    if (s%com%common_lattice .and. gang) then
      dflt_good_unis = .false.
      dflt_good_unis(ix_common_uni$) = .true.
    endif

  else
    call location_decode (default_universe, dflt_good_unis, 1, num)
    if (num == 0) dflt_good_unis = .true.  ! blank => all
    if (num < 0) then
      call out_io (s_error$, r_name, 'ERROR READING DEFAULT_UNIVERSE FOR: ' // v1_var%name)
      cycle
    endif
  endif

  ! Set up variable lists

  if (gang) then
    call tao_var_stuffit1 (var, v1_var_ptr, v1_var, -1, searching, &
              use_same_lat_eles_as, search_for_lat_eles, ix_min_var, ix_max_var, &
              default_attribute, default_weight, default_step, default_merit_type, &
              default_low_lim, default_high_lim, dflt_good_unis)

    if (.not. searching) then
      do j = lbound(v1_var_ptr%v, 1), ubound(v1_var_ptr%v, 1)

        ! Find which universes
        good_unis = dflt_good_unis

        if (.not. searching) then
          if (var(j)%universe /= '') then
            if (var(j)%universe == '*') then
              good_unis = .true.
            else
              call location_decode (var(j)%universe, good_unis, 1, num)
              if (num < 0) then
                call out_io (s_error$, r_name, 'ERROR READING UNIVERSE FOR: ' // v1_var%name)
                cycle
              endif
            endif
          endif
        endif

        if (count(good_unis) == 0) then
          call out_io (s_error$, r_name, 'ERROR: NO UNIVERSE FOR: ' // v1_var%name)
          call err_exit
        endif

        call tao_var_stuffit2 (good_unis, v1_var_ptr%v(j))
      enddo
    endif

  else   ! If clone...
    do i = lbound(s%u, 1), ubound(s%u, 1)
      if (.not. dflt_good_unis(i)) cycle
      call tao_var_stuffit1 (var, v1_var_ptr, v1_var, i, searching, &
              use_same_lat_eles_as, search_for_lat_eles, ix_min_var, ix_max_var, &
              default_attribute, default_weight, default_step, default_merit_type, &
              default_low_lim, default_high_lim, dflt_good_unis)

      write (v1_var_ptr%name, '(2a, i0)') trim(v1_var_ptr%name), '_u', i
      if (.not. searching) then
        good_unis = .false.
        good_unis(i) = .true.
        do j = lbound(v1_var_ptr%v, 1), ubound(v1_var_ptr%v, 1)
          call tao_var_stuffit2 (good_unis, v1_var_ptr%v(j))
        enddo
      endif
    enddo
  endif

enddo

close (iu)
deallocate (dflt_good_unis, good_unis)
deallocate (default_key_b, default_key_d)

! For a variable pointing to a term in a taylor map, the pointer can get messed up if a later variable points to a taylor term
! that does not exist forcing the reallocation of the map to create the taylor term.
! To get around this just repoint the pointers

call tao_var_repoint()

!

call tao_setup_key_table ()

! Check that all lattice values are within limits.

call tao_limit_calc(limited)
if (limited) then
  call out_io (s_fatal$, r_name, 'VARIABLES OUT OF RANGE ON STARTUP NOT PERMITTED!')
  stop
endif

end subroutine tao_init_variables

!-----------------------------------------------------------------------
!------------------------------------------------------------------------
! stuff common to all universes

subroutine tao_var_stuffit1 (var, v1_var_ptr, v1_var, ix_uni, searching, &
              use_same_lat_eles_as, search_for_lat_eles, ix_min_var, ix_max_var, &
              default_attribute, default_weight, default_step, default_merit_type, &
              default_low_lim, default_high_lim, dflt_good_unis)

use tao_input_struct

implicit none

type (tao_v1_var_array_struct), allocatable, save, target :: v1_array(:)
type (tao_v1_var_struct), pointer :: v1_var_ptr
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: var_ptr
type (ele_pointer_struct), allocatable, save :: eles(:)
type (tao_v1_var_input) v1_var
type (tao_var_input) var(n_var_minn:n_var_maxx)

real(rp) default_weight, default_step, default_low_lim, default_high_lim

character(*) use_same_lat_eles_as, search_for_lat_eles, default_attribute
character(*) default_merit_type
character(20) fmt
character(40) name
character(40), allocatable :: ele_names(:)
character(60) search_string
character(20) :: r_name = 'tao_var_stuffit1'

integer i, iu, ip, j, jj, k, kk, n, nn, n1, n2, ix1, ix2, ix
integer num_ele, ios, ix_uni, ixm, ix2m
integer, allocatable, save :: an_indexx(:)
integer ix_min_var, ix_max_var

logical :: dflt_good_unis(lbound(s%u,1):)
logical searching, grouping, found_one, err, logical_is_garbage

! count number of v1 entries

s%n_v1_var_used = s%n_v1_var_used + 1
nn = s%n_v1_var_used
v1_var_ptr => s%v1_var(nn)
v1_var_ptr%name = v1_var%name

! If reusing a previous element list...

if (use_same_lat_eles_as /= '') then
  call string_trim (use_same_lat_eles_as, name, ix)
  call tao_find_var (err, name, v1_array = v1_array)
  if (err .or. size(v1_array) /= 1) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "USE_SAME_LAT_ELES_AS": ' // name)
    call err_exit
  endif
  v1_ptr => v1_array(1)%v1
  n1 = s%n_var_used + 1
  n2 = s%n_var_used + size(v1_ptr%v)
  call tao_allocate_var_array (n2)

  ix_min_var = lbound(v1_ptr%v, 1)
  ix1 = ix_min_var
  ix2 = ix1 + (n2 - n1)

  s%var(n1:n2)%ele_name    = v1_ptr%v%ele_name
  s%var(n1:n2)%s           = v1_ptr%v%s

  do n = n1, n2
    ix = ix1 + (n - n1)
    ip = 1 + (n - n1) 

    s%var(n)%good_user = v1_ptr%v(ip)%good_user
    if (.not. logical_is_garbage(var(ix)%good_user)) s%var(n)%good_user = var(ix)%good_user

    s%var(n)%key_bound = v1_ptr%v(ip)%key_bound
    if (.not. logical_is_garbage(var(ix)%key_bound)) then
      if (var(ix)%key_bound) s%var(n)%key_bound = var(ix)%key_bound
    endif

    s%var(n)%key_delta = v1_ptr%v(ip)%key_delta
    if (var(ix)%key_delta /= 0) s%var(n)%key_delta = var(ix)%key_delta

    s%var(n)%attrib_name = v1_ptr%v(ip)%attrib_name
    if (default_attribute /= '') s%var(n)%attrib_name = default_attribute
    if (var(ix)%attribute /= '') s%var(n)%attrib_name = var(ix)%attribute

    s%var(n)%weight = v1_ptr%v(ip)%weight
    if (default_weight /= 0) s%var(n)%weight = default_weight
    if (var(ix)%weight /= 0) s%var(n)%weight = var(ix)%weight

    s%var(n)%step = v1_ptr%v(ip)%step
    if (default_step /= 0) s%var(n)%step = default_step
    if (var(ix)%step /= 0) s%var(n)%step = var(ix)%step

    s%var(n)%merit_type = v1_ptr%v(ip)%merit_type
    if (default_merit_type /= '') s%var(n)%merit_type = default_merit_type
    if (var(ix)%merit_type /= '') s%var(n)%merit_type = var(ix)%merit_type

    s%var(n)%low_lim = v1_ptr%v(ip)%low_lim
    if (default_low_lim /= -1e30) s%var(n)%low_lim = default_low_lim
    if (var(ix)%low_lim /= -1e30) s%var(n)%low_lim = var(ix)%low_lim

    s%var(n)%high_lim = v1_ptr%v(ip)%high_lim
    if (default_high_lim /= 1e30) s%var(n)%high_lim = default_high_lim
    if (var(ix)%high_lim /= 1e30) s%var(n)%high_lim = var(ix)%high_lim

    s%var(n)%key_bound = v1_ptr%v(ip)%key_bound

    s%var(n)%ix_v1 = ix_min_var + n - n1
    s%var(n)%v1 => v1_var_ptr
  enddo

  call tao_point_v1_to_var (v1_var_ptr, s%var(n1:n2), ix_min_var)

  return

endif

!------------------------------
! Are we searching for elements?

if (search_for_lat_eles /= '') then

  searching = .true.
  call string_trim (search_for_lat_eles, search_string, ix)
  grouping = .true.  ! Default
  if (search_string(1:12) == '-no_grouping') then
    grouping = .false.
    search_string = search_string(13:)
  endif

  if (index(search_string, '-no_slaves') /= 0) then
    call out_io (s_error$, r_name, &
          'Note: The "-no_slaves" switch is not valid with search_for_lat_eles for *variables*.')
  else
    search_string = '-no_slaves ' // trim(search_string)
  endif

  if (any(var%universe /= '')) then
    call out_io (s_abort$, r_name, &
           "CANNOT SPECIFY INDIVIDUAL UNIVERSES WHEN SEARCHING FOR VARIABLES")
    call err_exit
  endif

  ! find matching elements...
  ! first count how many variables we need

  num_ele = 0
  found_one = .false.
  allocate(ele_names(100), an_indexx(100))

  do iu = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. dflt_good_unis(iu)) cycle
    if (ix_uni > -1 .and. iu /= ix_uni) cycle
    call tao_init_find_elements (s%u(iu), search_string, eles, default_attribute, found_one)

    if (grouping) then
      do kk = 1, size(eles)
        call find_indexx(eles(kk)%ele%name, ele_names, an_indexx, num_ele, ixm, ix2m)
        if (ixm == 0) then
          if (num_ele+1 > size(ele_names)) then
            call re_allocate(ele_names, size(ele_names) + 100)
            call re_allocate(an_indexx, size(an_indexx) + 100)
          endif
          an_indexx(ix2m+1:num_ele+1) = an_indexx(ix2m:num_ele)
          an_indexx(ix2m) = num_ele + 1
          ele_names(num_ele+1) = eles(kk)%ele%name
          num_ele = num_ele + 1
        endif
      enddo
    else
      num_ele = num_ele + size(eles)
    endif
  enddo
  deallocate(ele_names, an_indexx)

  if (.not. found_one) then
    call out_io (s_error$, r_name, &
            'NO ELEMENTS FOUND IN SEARCH FOR: ' // search_for_lat_eles, &
            'WHILE SETTING UP VARIABLE ARRAY FOR: ' // v1_var%name, &
            'THIS V1_VARIABLE WILL NOT BE CREATED!')
    s%n_v1_var_used = s%n_v1_var_used - 1
    return
  endif

  ! Now load the information into the variable structures

  n1 = s%n_var_used + 1
  n2 = s%n_var_used + num_ele
  call tao_allocate_var_array (n2)

  nn = n1 - 1
  do iu = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. dflt_good_unis(iu)) cycle
    if (ix_uni > -1 .and. iu /= ix_uni) cycle
    call tao_init_find_elements (s%u(iu), search_string, eles, default_attribute)
    kk_loop: do kk = 1, size(eles)
      ! If the name matches an existing variable then use that variable
      if (grouping) then
        do j = n1, nn
          if (s%var(j)%ele_name == eles(kk)%ele%name) then
            call tao_pointer_to_var_in_lattice (s%var(j), iu, eles(kk)%ele, err)
            cycle kk_loop
          endif
        enddo
      endif
      ! Here if there is no match...
      ! create a new variable and stuff the info into it.
      nn = nn + 1
      var_ptr => s%var(nn)
      var_ptr%v1 => v1_var_ptr  ! Used for error messages in pointer_to_var
      var_ptr%ix_v1 = ix_min_var + nn - n1
      var_ptr%exists = .true.
      var_ptr%ele_name = eles(kk)%ele%name
      var_ptr%s = eles(kk)%ele%s
      var_ptr%attrib_name = default_attribute
      call tao_pointer_to_var_in_lattice (var_ptr, iu, eles(kk)%ele, err)
    enddo kk_loop
  enddo 

  ix1 = ix_min_var
  ix2 = (n2 - n1) + ix_min_var

!------------------------------
! If not searching or reusing...

else  
  searching = .false.

  ! if ix_max/min_var has not been set then just search for elements that have been named
  ! and use this info to set ix_max/min

  if (ix_max_var < ix_min_var) then
    found_one = .false.
    do i = lbound(var, 1), ubound(var, 1)
      if (var(i)%ele_name /= '') then
        ix_max_var = i
        if (.not. found_one) ix_min_var = i
        found_one = .true.
      endif
    enddo
    if (ix_min_var > 1) ix_min_var = 1 
  endif

  if (ix_max_var < ix_min_var) then
    call out_io (s_abort$, r_name, 'NO ELEMENTS FOR: ' // v1_var%name)
    call err_exit
  endif

  n1 = s%n_var_used + 1
  n2 = s%n_var_used + ix_max_var - ix_min_var + 1
  ix1 = ix_min_var
  ix2 = ix_max_var
 
  call tao_allocate_var_array (n2)

  s%var(n1:n2)%ele_name = var(ix1:ix2)%ele_name

endif

!---------------------------------------------------------------
! Common init stuff

do n = n1, n2
  i = ix1 + n - n1

  s%var(n)%key_bound = .false.
  if (.not. logical_is_garbage(var(i)%key_bound)) then
    s%var(n)%key_bound = var(i)%key_bound
  endif

  if (logical_is_garbage(var(i)%good_user)) then
    s%var(n)%good_user = .true.
  else
    s%var(n)%good_user = var(i)%good_user
  endif

  s%var(n)%ix_v1 = ix_min_var + n - n1
  s%var(n)%v1 => v1_var_ptr

enddo

call tao_point_v1_to_var (v1_var_ptr, s%var(n1:n2), ix_min_var)

s%var(n1:n2)%key_delta = var(ix1:ix2)%key_delta

s%var(n1:n2)%attrib_name = var(ix1:ix2)%attribute
where (s%var(n1:n2)%attrib_name == '') s%var(n1:n2)%attrib_name = default_attribute

s%var(n1:n2)%weight = var(ix1:ix2)%weight
where (s%var(n1:n2)%weight == 0) s%var(n1:n2)%weight = default_weight
 
s%var(n1:n2)%step = var(ix1:ix2)%step
where (s%var(n1:n2)%step == 0) s%var(n1:n2)%step = default_step
 
s%var(n1:n2)%merit_type = var(ix1:ix2)%merit_type
where (s%var(n1:n2)%merit_type == '') s%var(n1:n2)%merit_type = default_merit_type
 
s%var(n1:n2)%low_lim = var(ix1:ix2)%low_lim
where (s%var(n1:n2)%low_lim == -1e30) s%var(n1:n2)%low_lim = default_low_lim
 
s%var(n1:n2)%high_lim = var(ix1:ix2)%high_lim
where (s%var(n1:n2)%high_lim == 1e30) s%var(n1:n2)%high_lim = default_high_lim

end subroutine tao_var_stuffit1 
  
!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine tao_allocate_v1_var (n_v1)

implicit none

integer i, n_v1

!

if (allocated(s%v1_var)) deallocate (s%v1_var)
allocate (s%v1_var(n_v1))

do i = 1, n_v1
  s%v1_var(i)%name = ''  ! blank name means it doesn't (yet) exist
  s%v1_var(i)%ix_v1_var = i
enddo

s%n_v1_var_used = 0
s%n_var_used = 0

end subroutine tao_allocate_v1_var

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine tao_var_stuffit2 (good_unis, var)

implicit none

type (tao_var_struct), target :: var
type (tao_var_slave_struct), pointer :: var_slave
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (all_pointer_struct), allocatable :: a_ptr(:)

integer i, j, n, n1, n2, ie, iu, ib

character(20) :: r_name = 'tao_var_stuffit2'
logical err, good_unis(lbound(s%u, 1):), found

! 

if (var%ele_name == '') then
  var%exists = .false.
  var%key_bound = .false.
  return
endif

found = .false.
do iu = lbound(s%u, 1), ubound(s%u, 1)

  if (.not. good_unis(iu)) cycle
  lat => s%u(iu)%model%lat

  if (var%ele_name == 'BEAM_START') then
    call tao_pointer_to_var_in_lattice2 (var, iu, err)
    if (err) return
    found = .true.
  else
    do ib = 0, ubound(lat%branch, 1)
      do ie = 0, lat%branch(ib)%n_ele_max
        if (var%ele_name /= lat%branch(ib)%ele(ie)%name) cycle
        call tao_pointer_to_var_in_lattice (var, iu, lat%branch(ib)%ele(ie), err)
        if (err) return
        found = .true.
      enddo
    enddo
  endif
enddo

if (.not. found) then
  call out_io (s_error$, r_name, &
            'CANNOT FIND LATTICE ELEMENT WITH NAME: "' // trim(var%ele_name) // '"', &
            'FOR VARIABLE: ' // var%v1%name)
  return
endif

if (size(var%slave) > 0) then
  var%exists = .true.
endif

end subroutine tao_var_stuffit2

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_pointer_to_var_in_lattice (var, ix_uni, ele, err)
! 
! Routine to add a pointer from a given Tao variable
! to the appropriate variable in a lattice.
!
! Input:
!   var       -- Tao_var_struct: Structure has the info of where to point.
!   ix_uni    -- Integer: the universe to use
!   ix_ele    -- Integer: Index of element.
!
! Output:
!   var%slave(ix_slave) -- Tao_var_slave_struct: New component of %slave(:) array is added. 
!     %model_ptr
!     %base_ptr
!     %ix_ele
!     %ix_uni
!   err       -- Logical: Set True if there is an error. False otherwise.
!-

subroutine tao_pointer_to_var_in_lattice (var, ix_uni, ele, err)

implicit none

type (tao_var_struct), target :: var
type (ele_struct) ele
type (ele_struct), pointer :: ele2
type (tao_universe_struct), pointer :: u
type (tao_var_slave_struct), pointer :: var_slave
type (tao_var_slave_struct) :: var_slave_saved(size(var%slave))
type (all_pointer_struct) a_ptr

integer ix, ix_uni, ix_slave
logical :: err
character(30) :: r_name = 'tao_pointer_to_var_in_lattice'

! locate element

err = .true.

u => s%u(ix_uni)
if (s%com%common_lattice) u => s%com%u_working

! allocate space for var%slave.

ix_slave = size(var%slave) + 1
var_slave_saved = var%slave
deallocate (var%slave)  
allocate (var%slave(ix_slave))
var%slave(1:ix_slave-1) = var_slave_saved
var_slave => var%slave(ix_slave)

! locate attribute

ele2 => pointer_to_ele (u%model%lat, ele%ix_ele, ele%ix_branch)
call pointer_to_attribute (ele2, var%attrib_name, .true., a_ptr, err, .false., ix_attrib = var%ix_attrib)
if (err .or. .not. associated(a_ptr%r)) then
  call out_io (s_error$, r_name, &
            'IN VARIBALE: ' // tao_var1_name(var), &
            '  INVALID ATTRIBUTE: ' // var%attrib_name, &
            '  FOR ELEMENT: ' // ele%name)
  var%model_value => var%old_value
  var%base_value  => var%old_value
  var%exists = .false.
  var%key_bound = .false.
  return
endif

var_slave%model_value => a_ptr%r
ele2 => pointer_to_ele (u%base%lat, ele%ix_ele, ele%ix_branch)

call pointer_to_attribute (ele2,  var%attrib_name, .true., a_ptr,  err)

var_slave%base_value => a_ptr%r
var_slave%ix_ele    = ele%ix_ele
var_slave%ix_branch = ele%ix_branch
var_slave%ix_uni    = ix_uni

var%model_value => var%slave(1)%model_value
var%base_value  => var%slave(1)%base_value
var%design_value = var%slave(1)%model_value

! Common pointer

if (associated(u%common)) then
  ele2 => pointer_to_ele (u%common%model%lat, ele%ix_ele, ele%ix_branch)
  call pointer_to_attribute (ele,  var%attrib_name, .true., a_ptr,  err)
  var%common_slave%model_value => a_ptr%r
  ele2 => pointer_to_ele (u%common%base%lat, ele%ix_ele, ele%ix_branch)
  call pointer_to_attribute (ele, var%attrib_name, .true., a_ptr, err)
  var%common_slave%base_value => a_ptr%r
endif

! With unified lattices: model_value and base_value get their own storage
! instead of pointing to var%slave(1). 
! Exception: If variable controls a common parameter

if (s%com%common_lattice) then
  if (var_slave%ix_uni == ix_common_uni$) then
    var%model_value => var%common_slave%model_value
    var%base_value => var%common_slave%base_value
  else
    allocate (var%model_value, var%base_value)
    var%model_value = var_slave%model_value
    var%base_value = var_slave%base_value
  endif
endif

end subroutine tao_pointer_to_var_in_lattice

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_pointer_to_var_in_lattice2 (var, ix_uni, err)
! 
! Routine to add a pointer from a given Tao variable
! to the appropriate variable in a lattice.
!
! Input:
!   var         -- Tao_var_struct: Structure has the info of where to point.
!   ix_uni      -- Integer: the universe to use
!
! Output:
!   var%slave(ix_slave) -- Tao_var_slave_struct: New component of %slave(:) array is added. 
!     %model_ptr
!     %base_ptr
!     %ix_ele
!     %ix_uni
!   err       -- Logical: Set True if there is an error. False otherwise.
!-

subroutine tao_pointer_to_var_in_lattice2 (var, ix_uni, err)

implicit none

type (tao_var_struct), target :: var
type (ele_pointer_struct), allocatable :: eles(:)
type (tao_universe_struct), pointer :: u
type (tao_var_slave_struct), pointer :: var_slave
type (tao_var_slave_struct) :: var_slave_saved(size(var%slave))
type (all_pointer_struct), allocatable :: a_ptr(:), b_ptr(:), cm_ptr(:), cb_ptr(:)

integer ii, ix, ix_uni, n_old
logical :: err
character(30) :: r_name = 'tao_pointer_to_var_in_lattice'

! locate element

err = .true.

u => s%u(ix_uni)
if (s%com%common_lattice) u => s%com%u_working

call pointers_to_attribute (u%model%lat, var%ele_name, var%attrib_name, .true., a_ptr, err, .false., eles, var%ix_attrib)
if (err .or. size(a_ptr) == 0) then
  call out_io (s_error$, r_name, &
            'IN VARIBALE: ' // tao_var1_name(var), &
            '  INVALID ATTRIBUTE: ' // var%attrib_name, &
            '  FOR ELEMENT: ' // var%ele_name)
  var%model_value => var%old_value
  var%base_value  => var%old_value
  var%exists = .false.
  var%key_bound = .false.
  return
endif

call pointers_to_attribute (u%base%lat, var%ele_name, var%attrib_name, .true., b_ptr, err, .false., eles, var%ix_attrib)
if (associated(u%common)) then
  call pointers_to_attribute (u%common%model%lat, var%ele_name, var%attrib_name, .true., cm_ptr, err, .false., eles, var%ix_attrib)
  call pointers_to_attribute (u%common%base%lat,  var%ele_name, var%attrib_name, .true., cb_ptr, err, .false., eles, var%ix_attrib)
endif

! allocate space for var%slave.

n_old = size(var%slave)
var_slave_saved = var%slave
deallocate (var%slave)  
allocate (var%slave(n_old+size(a_ptr)))
var%slave(1:n_old) = var_slave_saved

! locate attribute


do ii = 1, size(a_ptr)
  var_slave => var%slave(n_old+ii)
  var_slave%model_value => a_ptr(ii)%r
  var_slave%base_value  => b_ptr(ii)%r
  var_slave%ix_uni = ix_uni
  if (associated(eles(ii)%ele)) then
    var_slave%ix_ele    = eles(ii)%ele%ix_ele
    var_slave%ix_branch = eles(ii)%ele%ix_branch
  else
    var_slave%ix_ele    = -1
    var_slave%ix_branch = 0
  endif

  var%model_value => var%slave(1)%model_value
  var%base_value  => var%slave(1)%base_value
  var%design_value = var%slave(1)%model_value

  ! Common pointer

  if (associated(u%common)) then
    var%common_slave%model_value => cm_ptr(ii)%r
    var%common_slave%base_value  => cb_ptr(ii)%r
  endif

  ! With unified lattices: model_value and base_value get their own storage
  ! instead of pointing to var%slave(1). 
  ! Exception: If variable controls a common parameter

  if (s%com%common_lattice) then
    if (var_slave%ix_uni == ix_common_uni$) then
      var%model_value => var%common_slave%model_value
      var%base_value => var%common_slave%base_value
    else
      allocate (var%model_value, var%base_value)
      var%model_value = var_slave%model_value
      var%base_value = var_slave%base_value
    endif
  endif

enddo

end subroutine tao_pointer_to_var_in_lattice2

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_allocate_var_array (n_var)
!
! Routine to increase the s%var(:) array size.
!
! Input:
!   n_var -- Integer: Size of s%var(:) wanted.
!-

subroutine tao_allocate_var_array (n_var)

type (tao_var_struct) :: var(size(s%var))
type (tao_v1_var_struct), pointer :: v1

integer i, j1, j2, n0, n_var

! First save the information presently in s%var(:) in the var(:) array.
! s%n_var_used gives the present upper bound to s%var(:).

if (allocated(s%var)) then
  n0 = s%n_var_used
  do i = 1, n0
    allocate(var(i)%slave(size(s%var(i)%slave)))
  enddo
  var(1:n0) = s%var(1:n0)
  deallocate (s%var)
  allocate (s%var(n_var))
  do i = 1, n0
    allocate(s%var(i)%slave(size(var(i)%slave)))
  enddo
  s%var(1:n0) = var(1:n0)
else
  n0 = 0
  allocate (s%var(n_var))
endif

! Since the s%var(:) array has been reallocated, the pointer from s%var(:)%v1 
! to the datums must be reestablished.

j2 = 0
do
  j1 = j2 + 1
  if (j1 > n0) exit
  v1 => s%var(j1)%v1
  do 
    if (j2 == n0) exit
    if (.not. associated(s%var(j2+1)%v1, v1)) exit
    j2 = j2 + 1
  enddo
  call tao_point_v1_to_var (v1, s%var(j1:j2), s%var(j1)%ix_v1)
enddo

! Set the defaults for the newly created slots in the s%var(:) array

do i = n0+1, size(s%var)
  s%var(i)%ix_var    = i
  s%var(i)%good_opt  = .true.
  s%var(i)%exists    = .false.
  s%var(i)%good_var  = .true.
  s%var(i)%good_user = .true.
  s%var(i)%model_value => s%com%dummy_target  ! Just to point to somewhere
  s%var(i)%base_value  => s%com%dummy_target  ! Just to point to somewhere
  s%var(i)%low_lim = -1e30
  s%var(i)%high_lim = 1e30
  allocate(s%var(i)%slave(0))
enddo
  
s%n_var_used = n_var


end subroutine tao_allocate_var_array

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_point_v1_to_var (v1, var, n)
!
! used for arbitrary variable pointer indexing
!
! v1       -- tao_v1_var_struct: Contains the pointer.
! var(n:)  -- tao_var_struct: the variable
! n        -- integer: starting index for the var array.
!-

subroutine tao_point_v1_to_var (v1, var, n)

implicit none

integer n, i, n_var

type (tao_v1_var_struct), target :: v1
type (tao_var_struct), target :: var(n:)

v1%v => var

do i = lbound(var, 1), ubound(var, 1)
  var(i)%v1 => v1
  var(i)%ix_v1 = i
enddo

end subroutine 

end module
