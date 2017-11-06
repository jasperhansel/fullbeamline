!+
! Function distance_to_aperture_custom (orbit, particle_at, ele, no_aperture_here) result (dist)
!
! Replace this routine with a routine to calculate the percentage distance from the particle to the wall aperture.
! Distances are negative if the particle is inside the wall.
!
! A valid distance_to_aperture_custom routine is needed only if ele%aperture_type is set to custom$.
! General rule: Your code may NOT modify any argument that is not listed as an output agument below."
!
! This routine is called by distance_to_aperture.
!
! Input:
!   orbit         -- coord_struct: Particle position.
!   particle_at   -- integer: 
!   ele           -- ele_struct: Element containing aperture.
!
! Output:
!   no_wall_here  -- logical, optional: True if aperture does not exist at the
!                       longitudinal location of the particle.
!   dist          -- real(rp): Percentage distance to the aperture.
!-

function distance_to_aperture_custom (orbit, particle_at, ele, no_aperture_here) result (dist)

use bmad, dummy => distance_to_aperture_custom

implicit none

type (coord_struct) orbit
type (ele_struct) ele

real(rp) dist, x_particle, y_particle
integer physical_end
logical particle_at, no_aperture_here

character(*), parameter :: r_name = 'distance_to_aperture_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
if (global_com%exit_on_error) call err_exit

dist = real_garbage$

end function distance_to_aperture_custom
