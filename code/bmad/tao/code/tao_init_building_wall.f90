!+
! Subroutine tao_init_building_wall (wall_file)
!
! Subroutine to initialize the building wall.
!
! Input:
!   wall_file  -- character(*): file name containing the wall definition.
!                   If blank then there is no wall.
!-

subroutine tao_init_building_wall (wall_file)

use tao_mod
use tao_input_struct

implicit none

type (tao_building_wall_point_struct) point(100)
type (tao_building_wall_point_struct), pointer :: pt(:)

real(rp) x_mid, z_mid, dx, dz, a, a2
integer i, j, iu, ios, n_wall

character(*) wall_file
character(200) complete_file_name
character(24) :: r_name = 'tao_init_building_wall'
character(16) constraint

namelist / building_wall_section / constraint, point

! Open file

if (wall_file == '') return

call out_io (s_blank$, r_name, '*Init: Opening File: ' // complete_file_name)
call tao_open_file (wall_file, iu, complete_file_name, s_fatal$)
if (iu == 0) then
  call out_io (s_blank$, r_name, 'ERROR OPENING BUILDING WALL FILE. WILL EXIT HERE...')
  call err_exit
endif

! Count the number of walls

n_wall = 0
do
  read (iu, nml = building_wall_section, iostat = ios)
  if (ios > 0) then
    call out_io (s_fatal$, r_name, 'ERROR READING BUILDING_WALL_SECTION NAMELIST')
    do   ! Generate informational message
      read (iu, nml = building_wall_section)
    enddo
  endif
  if (ios < 0) exit
  n_wall = n_wall + 1
enddo

call out_io (s_blank$, r_name, '  Number of building walls: \I2\ ', n_wall)

allocate (s%building_wall%section(n_wall))
if (n_wall == 0) return ! no walls

! Now transfer the information

rewind (iu)
do i = 1, n_wall

  constraint = 'none'
  point%radius = 0
  point%x = real_garbage$
  read (iu, nml = building_wall_section, iostat = ios)

  select case (constraint)
  case ('left_side', 'right_side', 'none')
    s%building_wall%section(i)%constraint = constraint
  case default
    call out_io (s_error$, r_name, 'BAD "CONSTRAINT" FOR BUILDING WALL: ' // constraint)
    call err_exit
  end select

  do j = size(point), 1, -1
    if (point(j)%x == real_garbage$) cycle
    allocate (s%building_wall%section(i)%point(j))
    exit
  enddo

  if (.not. allocated(s%building_wall%section(i)%point)) then
    call out_io (s_fatal$, r_name, 'ERROR READING BUILDING_WALL_SECTION NAMELIST POINT ARRAY', &
                                   'IN BUILDING_WALL_SECTION NAMELIST NUMBER \i0\ ', i_array = [i])
    return
  endif

  if (point(1)%radius /= 0) then
    call out_io (s_fatal$, r_name, 'ERROR IN POINT ARRAY OF BUILDING_WALL_SECTION NAMELIST NUMBER \i0\ ', &
                                   'FIRST POINT HAS NON-ZERO RADIUS', i_array = [i])
    return
  endif

  pt => s%building_wall%section(i)%point

  do j = 1, size(pt)
    
    pt(j)%x = point(j)%x
    pt(j)%z = point(j)%z
    pt(j)%radius = point(j)%radius
    if (pt(j)%radius /= 0) then
      x_mid = (pt(j)%x + pt(j-1)%x) / 2; z_mid = (pt(j)%z + pt(j-1)%z) / 2 
      dx    = (pt(j)%x - pt(j-1)%x) / 2; dz    = (pt(j)%z - pt(j-1)%z) / 2 
      a2 = (pt(j)%radius**2 - dx**2 - dz**2) / (dx**2 + dz**2)
      if (a2 < 0) then
        call out_io (s_fatal$, r_name, 'ERROR IN POINT ARRAY OF BUILDING_WALL_SECTION NAMELIST NUMBER \i0\ ', &
                                       'WALL POINTS TOO FAR APART FOR CIRCLE AT POINT \i0\ ', i_array = [i, j])
        return
      endif
      a = sqrt(a2)
      if (pt(j)%radius < 0) a = -a
      pt(j)%x_center = x_mid - a * dz
      pt(j)%z_center = z_mid + a * dx
    endif
  enddo

enddo

close (iu)

end subroutine tao_init_building_wall
