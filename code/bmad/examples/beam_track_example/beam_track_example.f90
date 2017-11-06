!+ 
! Program beam_track_example
!
! Example program to track a beam through a lattice
!
! Input (command line)
!   beam_track_params namelist
!		lat_filename lat.bmad
!       beam_init_struct
!-

program beam_track_example

use bmad
use beam_mod

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: ele1, ele2
type (beam_struct), target :: beam
type (beam_init_struct) beam_init
type (coord_struct), pointer :: particle
character(100) in_filename
character(100) lat_filename
character(100) outfile_name
integer :: n_arg, i, j
integer ::  nmlfile = 1
integer ::  outfile = 2
logical :: err

namelist / beam_track_params / lat_filename, beam_init, outfile_name

!------------------------------------------
!Setup

!Get input file from command line						   
n_arg = cesr_iargc()
if (n_arg > 1) then
  print *, 'Usage: beam_track_example <input_file>'
  print *, 'Default: <input_file> = beam_track_example.in'
  stop
endif

in_filename = 'beam_track_example.in' !default
if (n_arg == 1) call cesr_getarg(1, in_filename)

! Defaults of input file
lat_filename = 'lat.bmad'
beam_init%n_particle = 1
beam_init%n_bunch = 1

!read input file
print *, 'Opening: ', trim(in_filename)
open (nmlfile, file = in_filename, status = "old")
read (nmlfile, nml = beam_track_params)
close (nmlfile)

print *, '--------------------------------------'
write (*, '(a, a)') 'lattice:      ', lat_filename
write (*, '(a, i8)') 'n_particle = ', beam_init%n_particle
write (*, '(a, i8)') 'n_bunch    = ', beam_init%n_bunch
print *, '--------------------------------------'

!Parse Lattice
call bmad_parser (lat_filename, lat)

!Initialize beam
ele1 => lat%ele(0) !start element

call init_beam_distribution (ele1, lat%param, beam_init, beam)

!For example, set the macrocharge of all particles
do i = 1, beam_init%n_bunch
	do j = 1, beam_init%n_particle
		beam%bunch(i)%particle(j)%charge = 0.0_rp
	end do 
end do

!Alternatively, call reallocate_beam, which just builds the beam_struct
!call reallocate_beam (beam, beam_init%n_bunch, beam_init%n_particle)


!Track particles to end of lattice
ele2 => lat%ele(lat%n_ele_track)
call track_beam(lat, beam, ele1, ele2, err)

!------------------------------------------
!Write to file

open(outfile, file = outfile_name)

write (outfile, '(a)')  "! End coordinates of all particles"
write (outfile, '(a)')  '! Bunch Particle          X          Px          Y          Py           Z          Pz      State' 
do i = 1, beam_init%n_bunch
  do j = 1, size(beam%bunch(i)%particle)
    particle => beam%bunch(i)%particle(j)
	  write (outfile, '(i7, i9, 6es12.3, 4x, a)') i, j, particle%vec, coord_state_name(particle%state)
  enddo
end do
close(outfile)
print *, "Written: ", trim(outfile_name)
	
	
end program
