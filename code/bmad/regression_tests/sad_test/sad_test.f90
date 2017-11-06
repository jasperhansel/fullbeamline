program sad_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: orbit(:)

!

call system_command ('python ../../util_programs/sad_to_bmad/sad_to_bmad.py')

call bmad_parser('sler_1689.bmad', lat)
call twiss_and_track(lat, orbit)

ele => lat%ele(0)

open (1, file = 'output.now')
write (1, '(a, 6es16.8)') '"Orb0"  ABS 3E-12', orbit(0)%vec
write (1, '(a, 2es16.8)') '"Betas"  REL 1E-6', ele%a%beta, ele%b%beta
write (1, '(a, 2es16.8)') '"Alphas" ABS 1E-8', ele%a%alpha, ele%b%alpha



close (1)

end program
