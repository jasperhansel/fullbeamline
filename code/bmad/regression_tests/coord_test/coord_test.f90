!+
! Program coord_test
!
! This program is part of the Bmad regression testing suite.
!-

program coord_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) orbit, orb0
type (branch_struct), pointer :: branch

! 

call bmad_parser ('coord_test.bmad', lat)

open (1, file = 'output.now')

orb0%vec = lat%beam_start%vec
ele => lat%ele(1)
call init_coord (orbit, orb0%vec, ele, upstream_end$, ele%branch%param%particle)
call offset_particle (ele, lat%param, set$, orbit)
call offset_particle (ele, lat%param, unset$, orbit)

write (1, '(a, 6es20.12)') '"orbit1-electron"  ABS 1E-14', orbit%vec
write (1, '(a, 6es20.12)') '"length1-electron" ABS 1E-14', orbit%s, c_light * orbit%t, orbit%path_len


orb0%vec(5) = 0
orb0%vec(6) = sqrt(1 - orb0%vec(2)**2 - orb0%vec(4)**2)

branch => lat%branch(1)
ele => branch%ele(1)
call init_coord (orbit, orb0%vec, ele, upstream_end$, ele%branch%param%particle)
call offset_photon (ele, orbit, set$)
call offset_photon (ele, orbit, unset$)

write (1, '(a, 6es20.12)') '"orbit1-photon"  ABS 1E-14', orbit%vec
write (1, '(a, 6es20.12)') '"length1-photon" ABS 1E-14', orbit%s, c_light * orbit%t, orbit%path_len


ele => branch%ele(2)
orb0%t = 0
call init_coord (orbit, orb0%vec, ele, upstream_end$, ele%branch%param%particle)
call offset_photon (ele, orbit, set$)
call offset_photon (ele, orbit, unset$)

write (1, '(a, 6es20.12)') '"orbit2-photon"  ABS 1E-14', orbit%vec
write (1, '(a, 6es20.12)') '"length2-photon" ABS 1E-14', orbit%s-branch%ele(1)%s, c_light * orbit%t, orbit%path_len

ele => branch%ele(3)
orb0%t = 0
call init_coord (orbit, orb0%vec, ele, upstream_end$, ele%branch%param%particle)
call offset_photon (ele, orbit, set$)
call offset_photon (ele, orbit, unset$)

write (1, '(a, 6es20.12)') '"orbit3-photon"  ABS 1E-14', orbit%vec
write (1, '(a, 6es20.12)') '"length3-photon" ABS 1E-14', orbit%s-branch%ele(1)%s, c_light * orbit%t, orbit%path_len


close(1)

end program
