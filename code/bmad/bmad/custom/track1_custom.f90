!+
! Subroutine track1_custom (start_orb, ele, param, end_orb, err_flag, finished, track)
!
! Dummy routine for custom tracking.
! This routine needs to be replaced for a custom calculation.
! If not replaced and this routine is called, this routine will generate an error message.
!
! Also see:
!   track1_preprocess
!   track1_postprocess
!
! If this routine takes into account radiation damping and/or excitation when bmad_com%radiation_damping_on
! and/or bmad_com%radiation_fluctuations_on is True, a custom version of track1_preprocess should be
! constructed to set its radiation_included argument to True.
! If not, the track1 routine will use track1_radiation to include the radiation effects.
! Note: If this routine calles symp_lie_bmad, the symp_lie_bmad routine does take into account radiation effects.
!
! General rule: Your code may NOT modify any argument that is not listed as an output agument below.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- coord_struct: Starting position.
!   ele        -- ele_struct: Element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   end_orb     -- coord_struct: End position.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!   finished    -- logical: When set True, track1 will halt processing and return to its calling routine.
!   track       -- track_struct, optional: Structure holding the track information if the
!                    tracking method does tracking step-by-step.
!-

subroutine track1_custom (start_orb, ele, param, end_orb, err_flag, finished, track)

use bmad, except_dummy => track1_custom

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

logical err_flag, finished

character(*), parameter :: r_name = 'track1_custom'

integer :: fileno, error_indicator

real(rp) :: x, y, bx, by, bz, G, q

! general
finished = .true.

! generate unique file number
fileno = lunget()

! create parameters file
open (fileno, file='__fullbeamline/parameters.txt', status='replace', err=10)

! write custom attributes
if (trim(attribute_name(ele, custom_attribute1$)) /= '!NULL') &
  write (fileno, *, err=10) attribute_name(ele, custom_attribute1$), ' = ', ele%value(custom_attribute1$)
if (trim(attribute_name(ele, custom_attribute2$)) /= '!NULL') &
  write (fileno, *, err=10) attribute_name(ele, custom_attribute2$), ' = ', ele%value(custom_attribute2$)
if (trim(attribute_name(ele, custom_attribute3$)) /= '!NULL') &
  write (fileno, *, err=10) attribute_name(ele, custom_attribute3$), ' = ', ele%value(custom_attribute3$)
if (trim(attribute_name(ele, custom_attribute4$)) /= '!NULL') &
  write (fileno, *, err=10) attribute_name(ele, custom_attribute4$), ' = ', ele%value(custom_attribute4$)
if (trim(attribute_name(ele, custom_attribute5$)) /= '!NULL') &
  write (fileno, *, err=10) attribute_name(ele, custom_attribute5$), ' = ', ele%value(custom_attribute5$)

! close parameters file
close (fileno, err=10)

! call gpt
call system_command('fullbeamline-track-single')

! generate unique file number
fileno = lunget()

! open result file
open (fileno, file='__fullbeamline/result.txt', status='old', err=10)

! read error indicator
read (fileno, *, err=10) error_indicator

! check error indicator
if (error_indicator /= 0) goto 20

! read data
read (fileno, *, err=10) x, y, bx, by, bz, G, q

! close file
close (fileno, err=10, status='delete')

! set particle parameters
end_orb = start_orb
end_orb%beta = sqrt(bx ** 2 + by ** 2 + bz ** 2)
end_orb%vec(1) = x
end_orb%vec(2) = bx * G * m_electron / end_orb%p0c
end_orb%vec(3) = y
end_orb%vec(4) = by * G * m_electron / end_orb%p0c
end_orb%vec(5) = 0.0
end_orb%vec(6) = (bz * G * m_electron / end_orb%p0c) - 1
end_orb%s = ele%s
end_orb%charge = q
end_orb%beta = bz

! return
return

! executed if an io operation fails
10 call out_io (s_fatal$, r_name, 'IO OPERATION FAILED')
call err_exit
return

! executed if tracking failed
20 end_orb = start_orb
end_orb%state = lost$
return

end subroutine
