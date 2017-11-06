module xsif_lat_file_names_mod

character(200), allocatable, private, save :: file_name_list(:)
integer, private, save :: n_file_name_list = -1

contains

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine file_name_list_add (file_name)

  implicit none

  character(*) file_name
  character(200), allocatable :: file_list2(:)
  integer n

! Init if needed

  if (n_file_name_list == -1) then
    allocate (file_name_list(1))
    n_file_name_list = 1
    file_name_list(1) = file_name
    return
  endif

!

  n = n_file_name_list 
  allocate (file_list2(n))
  file_list2 = file_name_list
  deallocate (file_name_list)
  allocate (file_name_list(n+1))
  file_name_list(1:n) = file_list2
  file_name_list(n+1) = file_name 
  n_file_name_list = n + 1

end subroutine

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine file_name_list_show (file_names, n_names)

  implicit none

  character(200), allocatable :: file_names(:)
  integer n, n_names

!

  if (allocated(file_names)) deallocate (file_names)
  n = n_file_name_list 
  allocate (file_names(n))
  file_names = file_name_list
  n_names = n

end subroutine

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine file_name_list_init ()

  implicit none

  if (allocated(file_name_list)) deallocate (file_name_list)
  n_file_name_list = -1

end subroutine

end module
