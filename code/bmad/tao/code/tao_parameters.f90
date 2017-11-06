!+
! module tao_parameters
!
! File of parameters like the size of the u%data array.
!-

module tao_parameters

integer, parameter :: n_data_maxx     = 2000   ! max index of datum per d1_data
integer, parameter :: n_data_minn     = -100   ! min index of datum per d1_data
integer, parameter :: n_var_maxx      = 1000   ! max index of datum per v1_var 
integer, parameter :: n_var_minn      = -100   ! min index of datum per v1_var 
integer, parameter :: n_region_maxx   = 50     ! number of plotting regions.
integer, parameter :: n_curve_maxx    = 20     ! number of curves per graph

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module

