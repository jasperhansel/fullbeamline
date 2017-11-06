!+ 
! Program particle_track_example
!
! Example program to track one particle through the first element in a lattice.
!
! Command line syntax:
!   <path-to-bin-dir>/bin/particle_track_example <lattice-file-name>
! Default:
!   <lattice-file-name> = lat.bmad 
!
! Output:
!   track.dat
!-

program particle_track_example

use bmad
use spin_mod

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: ele 
type (coord_struct) particle_beg, particle_end

real(rp) vec0(6)
integer ::  iu = 1

character(100) lat_name
character(100) outfile_name


!------------------------------------------
! Get lattice file name and parse lattice.

lat_name = 'lat.bmad'
if (command_argument_count() > 1) then
  call get_command_argument(1, lat_name)
endif
print *,"Using lattice: ", trim(lat_name)

call bmad_parser (lat_name, lat)

!------------------------------------------
! Track through first element

ele => lat%ele(1) ! First real element in lattice

vec0 = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
call init_coord (particle_beg, vec0, ele, upstream_end$)

particle_beg%spin = polar_to_vec (spin_polar_struct(0.9_rp, 0.1_rp, 0.2_rp, 0.3_rp))
bmad_com%spin_tracking_on = .true.

call track1 (particle_beg, ele, lat%param, particle_end)

!------------------------------------------
! Write to file

outfile_name = 'track.dat'
open(iu, file = outfile_name)

write (iu, '(a)' )  "Start coordinates:"
write (iu, '(a, 6es18.10)') '  Phase Space:', particle_beg%vec
write (iu, '(a, 4f12.6)')   '  Spin:       ', particle_beg%spin
write (iu, *)
write (iu, '(a)' )  "End coordinates:"
write (iu, '(a, 6es18.10)') '  Phase Space:', particle_end%vec
write (iu, '(a, 4f12.6)')   '  Spin:       ', particle_end%spin

close(iu)
print *, "Written: ", trim(outfile_name)

end program
