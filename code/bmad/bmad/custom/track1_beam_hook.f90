!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam_hook (beam_start, lat, ele, beam_end, err, centroid, direction)
!
! Routine that can be customized for tracking a beam through a single element.
!
! Input:
!   beam_start   -- Beam_struct: Starting beam position.
!   lat          -- lat_struct: Lattice containing element to be tracked through.
!   ele          -- Ele_struct: Element to track through.
!   centroid(0:) -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                     Hint: Calculate this before beam tracking by tracking a single particle.
!   direction    -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!
! Output:
!   beam_end    -- beam_struct: Ending beam position.
!   err         -- Logical: Set true if there is an error.
!                    EG: Too many particles lost for a CSR calc.
!   finished    -- logical: When set True, the standard track1_beam code will not be called.
!-

subroutine track1_beam_hook (beam_start, lat, ele, beam_end, err, centroid, direction, finished)

use bmad, dummy => track1_beam_hook

implicit none

type (beam_struct) beam_start
type (beam_struct) :: beam_end
type (lat_struct) :: lat
type (ele_struct) ele
type (coord_struct), optional :: centroid(0:)

integer, optional :: direction
logical err, finished

integer :: fileno, N, i

real(rp) :: t_offset, particles_real, q_total, x, y, z, Bx, By, Bz, G, t

integer :: ios

character(350) :: line

! do nothing unless tracking is through a custom element
if (ele%key /= custom$) then
  finished = .false.
  return
endif

! generate file number
fileno = lunget()

! create parameters file
open (fileno, file='__fullbeamline/parameters.txt', status='new')

! write custom attributes
if (trim(attribute_name(ele, custom_attribute1$)) /= '!NULL') &
  write (fileno, *) attribute_name(ele, custom_attribute1$), ' = ', ele%value(custom_attribute1$)
if (trim(attribute_name(ele, custom_attribute2$)) /= '!NULL') &
  write (fileno, *) attribute_name(ele, custom_attribute2$), ' = ', ele%value(custom_attribute2$)
if (trim(attribute_name(ele, custom_attribute3$)) /= '!NULL') &
  write (fileno, *) attribute_name(ele, custom_attribute3$), ' = ', ele%value(custom_attribute3$)
if (trim(attribute_name(ele, custom_attribute4$)) /= '!NULL') &
  write (fileno, *) attribute_name(ele, custom_attribute4$), ' = ', ele%value(custom_attribute4$)
if (trim(attribute_name(ele, custom_attribute5$)) /= '!NULL') &
  write (fileno, *) attribute_name(ele, custom_attribute5$), ' = ', ele%value(custom_attribute5$)

! close parameters file
close (fileno)

! call gpt
call system_command('fullbeamline-track')

! open statistics.txt
open (fileno, file='__fullbeamline/statistics.txt', status='old')

read (fileno, *) line

read (fileno, *) t_offset, t_offset, particles_real, q_total

close (fileno, status='delete')

N = NINT(particles_real)

call reallocate_beam(beam_end, 1, N)

beam_end%bunch(1)%charge_tot = q_total
beam_end%bunch(1)%charge_live = q_total
beam_end%bunch(1)%z_center = 0.0
beam_end%bunch(1)%t_center = 0.0
beam_end%bunch(1)%ix_ele = ele%ix_ele
beam_end%bunch(1)%ix_bunch = 1
beam_end%bunch(1)%n_live = N


open (fileno, file='__fullbeamline/output.txt', status='old')

do
  read (fileno, *) line
  if (line(1:8) == 'position') exit
enddo

read (fileno, *) line

do i = 1, N
  read (fileno, *, iostat=ios) x, y, z, Bx, By, Bz, G
  beam_end%bunch(1)%particle(i)%vec(1) = x
  beam_end%bunch(1)%particle(i)%vec(2) = Bx * G * m_electron / ele%value(p0c$)
  beam_end%bunch(1)%particle(i)%vec(3) = y
  beam_end%bunch(1)%particle(i)%vec(4) = By * G * m_electron / ele%value(p0c$)
  beam_end%bunch(1)%particle(i)%vec(5) = 0.0
  beam_end%bunch(1)%particle(i)%vec(6) = (Bz * G * m_electron / ele%value(p0c$)) - 1
  beam_end%bunch(1)%particle(i)%s = ele%s
  beam_end%bunch(1)%particle(i)%t = t - t_offset
  beam_end%bunch(1)%particle(i)%charge = q_total / N
  beam_end%bunch(1)%particle(i)%p0c = ele%value(p0c$)
  beam_end%bunch(1)%particle(i)%beta = sqrt(1 - 1 / (G * G))
  beam_end%bunch(1)%particle(i)%ix_ele = ele%ix_ele
  beam_end%bunch(1)%particle(i)%state = alive$
  beam_end%bunch(1)%particle(i)%species = electron$
  beam_end%bunch(1)%particle(i)%location = downstream_end$
enddo

close (fileno, status='delete')


! no need to do regular tracking
finished = .true.

end subroutine
