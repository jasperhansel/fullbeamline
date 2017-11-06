program wake_test

use beam_mod
use wake_mod

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (wake_sr_mode_struct), pointer :: w
type (coord_struct) :: orb0
type (coord_struct) :: p1, p2
type (beam_init_struct) :: beam_init
type (bunch_struct) :: bunch, bunch_init, bunch0, bunch2

real(rp) dz, z0

integer i, nargs, n, ie, i0, i1

logical print_extra, err_flag

character(100) :: lat_file

!------------------------------------------

lat_file = 'wake_test.bmad'

print_extra = .false.
nargs = cesr_iargc()
if (nargs == 1)then
   call cesr_getarg(1, lat_file)
   print *, 'Using ', trim(lat_file)
   print_extra = .true.
elseif (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

call bmad_parser (lat_file, lat)

open (1, file = 'output.now')

! Short range wake test.

ele => lat%ele(1)
write (1, '(a, 2i4, f10.2)') '"SR-Size" ABS 0' , &
                size(ele%wake%sr_long%mode), size(ele%wake%sr_trans%mode), ele%wake%z_sr_max

do i = 1, size(ele%wake%sr_trans%mode)
  w => ele%wake%sr_trans%mode(i)
  write (1, '(a, i0, a, 4es14.6, 2i4)') '"SR-T', i, '" ABS 0', &
                                w%amp, w%damp, w%k, w%phi, w%polarization, w%transverse_dependence
enddo

do i = 1, size(ele%wake%sr_long%mode)
  w => ele%wake%sr_long%mode(i)
  write (1, '(a, i0, a, 4es14.6, 2i4)') '"SR-L', i, '" ABS 0', &
                               w%amp, w%damp, w%k, w%phi, w%polarization, w%transverse_dependence
enddo

!

z0 = 7e-3_rp

p1%vec = 0
p1%vec(1) = 1e-3
p1%vec(3) = -2e-3
p1%vec(5) = z0

call init_coord (p1, p1%vec, lat%ele(0), element_end = upstream_end$)
call init_coord (p2, p1%vec, lat%ele(0), element_end = upstream_end$)

p1%charge = 1
p2%charge = 0

ele => lat%ele(1)

! Make wake from first particle

ele%wake%sr_long%mode%b_sin = 0
ele%wake%sr_long%mode%b_cos = 0
ele%wake%sr_long%mode%a_sin = 0
ele%wake%sr_long%mode%a_cos = 0
ele%wake%sr_long%z_ref = p1%vec(5)

ele%wake%sr_trans%mode%b_sin = 0
ele%wake%sr_trans%mode%b_cos = 0
ele%wake%sr_trans%mode%a_sin = 0
ele%wake%sr_trans%mode%a_cos = 0
ele%wake%sr_trans%z_ref = p1%vec(5)

call sr_long_wake_particle (ele, p1)
call sr_trans_wake_particle (ele, p1)

! Add in a second wake

p1%vec(5) = z0 + lat%beam_start%vec(5)
call sr_long_wake_particle (ele, p1)
call sr_trans_wake_particle (ele, p1)

!

dz = -2e-3

do i = 1, 2
  p2%vec = 0
  p2%vec(5) = i*dz + z0
  call sr_long_wake_particle (ele, p2)
  call sr_trans_wake_particle (ele, p2)
  write (1, '(a, i0, a, 3es20.9)') '"SR', i, '" REL 1E-8' , p2%vec(2:6:2)
enddo

! Beam tracking

beam_init%n_particle = 100
beam_init%random_engine = 'quasi'
beam_init%a_norm_emit = 1e-12
beam_init%b_norm_emit = 1e-12
beam_init%dPz_dz = 0.0
beam_init%bunch_charge = 1 !100.0e-12
beam_init%sig_e = 1e-12
beam_init%sig_z = 5.99585e-3  ! 200 ps * cLight

call init_bunch_distribution (lat%ele(0), lat%param, beam_init, 0, bunch_init)
bunch = bunch_init

call track1_bunch (bunch, lat, ele, bunch, err_flag)

bmad_com%sr_wakes_on = .false.
bunch0 = bunch_init
call track1_bunch (bunch0, lat, ele, bunch0, err_flag)

write (1, '(a, 6es18.9)') '"SR-P20" REL 1E-8' , bunch%particle(20)%vec - bunch0%particle(20)%vec
write (1, '(a, 6es18.9)') '"SR-P40" REL 1E-8' , bunch%particle(40)%vec - bunch0%particle(40)%vec

! Long range wake test

ele => lat%ele(2)

bunch = bunch_init

do n = 1, 3
  call track1_bunch (bunch, lat, ele, bunch, err_flag)
  do i = 1, size(bunch%particle)
    bunch%particle(i)%t = bunch%particle(i)%t + 1e-7
  enddo
enddo

bmad_com%lr_wakes_on = .false.

bunch0 = bunch_init

do n = 1, 3
  call track1_bunch (bunch0, lat, ele, bunch0, err_flag)
  do i = 1, size(bunch0%particle)
    bunch0%particle(i)%t = bunch0%particle(i)%t + 1e-7
  enddo
enddo

write (1, '(a, 6es18.9)') '"LR-P20" REL 1E-8' , bunch%particle(20)%vec - bunch0%particle(20)%vec
write (1, '(a, 6es18.9)') '"LR-P40" REL 1E-8' , bunch%particle(40)%vec - bunch0%particle(40)%vec

! Long range wake with superimposed element.

bmad_com%lr_wakes_on = .true.
ele => lat%ele(3)

bunch2 = bunch_init

do n = 1, 3
  do ie = 3, 5
    call track1_bunch (bunch2, lat, lat%ele(ie), bunch2, err_flag)
  enddo
  do i = 1, size(bunch2%particle)
    bunch2%particle(i)%t = bunch2%particle(i)%t + 1e-7
  enddo
enddo

write (1, '(a, 6es18.9)') '"dB-LR-P20" ABS 1E-19' , bunch2%particle(20)%vec - bunch%particle(20)%vec
write (1, '(a, 6es18.9)') '"dB-LR-P40" ABS 1E-19' , bunch2%particle(40)%vec - bunch%particle(40)%vec

!---------------------------------
! Sort test

bunch%particle(1)%state = lost$
bunch%particle(3)%state = lost$
bunch%particle(10)%state = lost$
call order_particles_in_z (bunch)

do i = 1, size(bunch%particle) - 1
  i0 = bunch%ix_z(i)
  i1 = bunch%ix_z(i+1)

  if (bunch%particle(i1)%state /= alive$) then
    if (i /= bunch%n_live) then
      print *, 'Sort problem1'
      call err_exit
    endif
    exit
  endif

  if (bunch%particle(i0)%vec(5) < bunch%particle(i1)%vec(5)) then
    print *, 'Sort problem2'
    call err_exit
  endif
enddo

end program
