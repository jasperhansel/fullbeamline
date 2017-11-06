module synrad3d_test_mod

use synrad3d_track_mod
use synrad3d_output_mod
use synrad3d_parse_wall

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_reflection_test (param_file, who)
! 
! Routine to proform the reflection test.
!
! Input:
!   param_file    -- character(*): Input parameter file.
!   who           -- character(*): "reflection" or "diffuse_reflection"
!-

subroutine sr3d_reflection_test (param_file, who)

implicit none

type (photon_reflect_surface_struct) surface

real(rp) graze_angle_in, energy
real(rp) p_reflect, rel_p_specular, theta_out, phi_out
real(rp) surface_roughness_rms, roughness_correlation_len

integer n_photons
integer i, ix, ios, random_seed

character(*) param_file, who
character(200) output_file, surface_reflection_file

namelist / reflection_test / graze_angle_in, energy, n_photons, surface_roughness_rms, &
            roughness_correlation_len, surface_reflection_file, output_file, random_seed

! Set defaults

select case (who)
case ('reflection');          output_file = 'test_reflection.dat'
case ('diffuse_reflection');  output_file = 'test_diffuse_reflection.dat'
case default;                 call err_exit
end select

random_seed = 0

! Read parameters

open (1, file = param_file)
read (1, nml = reflection_test, iostat = ios)
if (ios > 0) then
  print *, 'ERROR READING REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
if (ios < 0) then
  print *, 'CANNOT FIND REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
close (1)

!

call ran_seed_put (random_seed)

if (surface_reflection_file == '') then
  call photon_reflection_std_surface_init (surface)
else
  call read_surface_reflection_file (surface_reflection_file, surface)
endif

if (surface_roughness_rms > 0) surface%surface_roughness_rms = surface_roughness_rms
if (roughness_correlation_len > 0) surface%roughness_correlation_len = roughness_correlation_len

call photon_reflectivity (graze_angle_in, energy, surface, p_reflect, rel_p_specular)

!

open (2, file = output_file)

write (2, *) 'Grazing angle in (rad):    ', graze_angle_in
write (2, *) 'Energy (eV):               ', energy
write (2, *) 'surface_roughness_rms:     ', surface_roughness_rms
write (2, *) 'roughness_correlation_len: ', roughness_correlation_len
write (2, *) 'Reflection probability:    ', p_reflect
write (2, *) 'Specular Reflection / Reflection probability ratio:', rel_p_specular
write (2, *) 'surface_reflection_file: "', trim(surface_reflection_file), '"'
write (2, *) 'random_seed:                 "', random_seed

write (2, *) '          #          theta_out                     phi_out'
do i = 1, n_photons
  select case (who)
  case ('reflection')
    call photon_reflection (graze_angle_in, energy, surface, theta_out, phi_out)
  case ('diffuse_reflection')
    call photon_diffuse_scattering (graze_angle_in, energy, surface, theta_out, phi_out)
  end select
  write (2, *) i, pi/2-theta_out, phi_out
enddo

close (2)
print *, 'Output file: ' // trim(output_file)

end subroutine sr3d_reflection_test

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_specular_reflection_test (param_file)
! 
! Routine to proform the reflection test.
!
! Input:
!   param_file    -- Character(*): Input parameter file.
!-

subroutine sr3d_specular_reflection_test (param_file)

implicit none

type (lat_struct), target :: lat
type (sr3d_photon_track_struct) :: photon
type (coord_struct) p
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)
type (random_state_struct) ran_state
type (branch_struct), pointer :: branch

real(rp) vel
integer i, ios, num_ignored, random_seed, n_photon

logical is_inside, err, absorbed

character(*) param_file
character(100) photon_start_input_file, output_file, lattice_file, wall_file

namelist / specular_reflection_test / photon_start_input_file, output_file, lattice_file, wall_file
namelist / start / p, ran_state, random_seed

! Read parameters

open (1, file = param_file)
output_file = ''
read (1, nml = specular_reflection_test, iostat = ios)
if (ios > 0) then
  print *, 'ERROR READING SPECULAR_REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
if (ios < 0) then
  print *, 'CANNOT FIND SPECULAR_REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
close (1)

if (output_file == '') output_file = 'test_specular_reflection.dat'

! Get lattice

if (lattice_file(1:6) == 'xsif::') then
  call xsif_parser(lattice_file(7:), lat)
else
  call bmad_parser (lattice_file, lat)
endif

! Init wall

call sr3d_read_wall_file (wall_file, lat)

! Open photon start input file and count the number of photons

print *, 'Opening photon starting position input file: ', trim(photon_start_input_file)
open (1, file = photon_start_input_file, status = 'old')
open (2, file = output_file)

allocate (wall_hit(0:10))
sr3d_params%specular_reflection_only = .true.
sr3d_params%allow_absorption = .false.
num_ignored = 0
n_photon = 0

branch => lat%branch(0)

do

  read (1, nml = start, iostat = ios)
  if (ios < 0) exit 
  if (ios > 0) then
    print *, 'Error reading photon starting position at photon index:', n_photon
    call err_exit
  endif

  vel = sqrt(p%vec(2)**2 + p%vec(4)**2 + p%vec(6)**2)
  if (abs(vel - 1) > 0.1) then
    print *, 'ERROR: PHOTON VELOCITY NOT PROPERLY NORMALIZED TO 1 FOR PHOTON:', n_photon
    stop
  endif
  p%vec(2:6:2) = p%vec(2:6:2) / vel

  p%p0c = 1000             ! Arbitrary
  p%ix_ele = element_at_s(lat, p%s, .true.)
  photon%start%orb = p
  photon%n_wall_hit = 0

  call sr3d_check_if_photon_init_coords_outside_wall (photon%start, lat, is_inside, num_ignored)

  n_photon = n_photon + 1
  photon%ix_photon_generated = n_photon
  photon%ix_photon = n_photon

  call sr3d_track_photon (photon, lat, wall_hit, err, .true.)
  call sr3d_print_hit_points (2, photon, wall_hit, branch)

enddo

print *, 'Output file: ' // trim(output_file)

end subroutine sr3d_specular_reflection_test 

end module

