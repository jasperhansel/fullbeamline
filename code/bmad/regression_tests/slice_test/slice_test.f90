program slice_test

use bmad
use transfer_map_mod
use em_field_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele0, ele, zele
type (ele_struct) ele1, ele2a, ele2b, end_ele, ele2
type (coord_struct), allocatable :: ref_orb(:)
type (coord_struct) orb1, orb2a, orb2b, orb2c, orb2d
type (coord_struct) start_orb, end_orb, end_orb2, orbit
type (taylor_struct) t_map(6)
type (em_field_struct) field

real(rp) s0, s1, s2, s_end
real(rp) xmat_c(6,6), vec0_c(6), xmat_d(6,6), vec0_d(6), xmat(6,6), xmat2(6,6)

integer i, j, idum, nargs, is
logical print_extra, err
character(100) lat_file

! Init

open (1, file = 'output.now')

! Test with elements.bmad

lat_file = 'elements.bmad'
print_extra = .false.
nargs = cesr_iargc()
if (nargs == 1) then
   call cesr_getarg(1, lat_file)
   print *, 'Using ', trim(lat_file)
   print_extra = .true.
elseif (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

call bmad_parser (lat_file, lat)

do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (ele%value(l$) == 0) cycle
  call init_coord (start_orb, lat%beam_start, ele, upstream_end$, lat%param%particle)
  call track1 (start_orb, ele, lat%param, end_orb)
  call make_mat6 (ele, lat%param, start_orb)
  xmat = ele%mat6
  s0 = lat%ele(i-1)%s; s1 = ele%s
  end_orb2 = start_orb
  ele2 = lat%ele(i-1)
  ele2%a = lat%ele(0)%a;  ele2%b = lat%ele(0)%b  ! Just to have some non-zero beta
  call mat_make_unit(xmat2)
  do is = 1, 10
    s_end = (s0 * (10 - is) + s1 * is) / 10
    call twiss_and_track_from_s_to_s (lat%branch(0), end_orb2, s_end, end_orb2, ele2, ele2, err)
    xmat2 = matmul(ele2%mat6, xmat2)
    if (err) exit
  enddo

  if (print_extra) then
    print *
    print *, '!--------------------------------------------------'
    print *, ele%name
    print '(a, 6es11.3)', 'Start:', start_orb%vec
    print '(a, 6es11.3)', 'Whole:', end_orb%vec
    print '(a, 6es11.3)', 'Split:', end_orb2%vec
    print '(a, 7es11.3)', 'Diff: ', end_orb2%vec - end_orb%vec, maxval(abs(end_orb2%vec - end_orb%vec))
    print *
    print *, 'Whole:'
    do j = 1, 6
      print '(6f12.6)', ele%mat6(j,:)
    enddo
    print *
    print *, 'Split:'
    do j = 1, 6
      print '(6f12.6)', xmat2(j,:)
    enddo
    print *
    print *, 'Diff:', maxval(abs(ele%mat6-xmat2))
    do j = 1, 6
      print '(6f12.6)', ele%mat6(j,:)-xmat2(j,:)
    enddo
  endif

enddo

if (print_extra) stop

!-------------------------------------------------
!-------------------------------------------------
! Test with slice_test.bmad

call bmad_parser ('slice_test.bmad', lat)

zele => lat%branch(1)%ele(1)
print *, determinant(zele%mat6)

!------------------------------------------------
! Test branch 0

branch => lat%branch(0)
call init_coord (start_orb, ele = branch%ele(1), element_end = upstream_end$)
! Track through split rfcav
call track1 (start_orb, branch%ele(1), branch%param, end_orb)
call track1 (end_orb, branch%ele(3), branch%param, end_orb)

! Track through unsplit rfcav
call track1 (start_orb, branch%ele(4), branch%param, end_orb2)

xmat = matmul(branch%ele(3)%mat6, branch%ele(1)%mat6)

write (1, '(a, 6es18.9)') '"rfcav:dvec" ABS 1e-14', end_orb2%vec - end_orb%vec
write (1, '(a, es18.9)')  '"rfcav:dt"   ABS 1e-8 ', c_light * (end_orb2%t - end_orb%t)
write (1, '(a, 6es18.9)') '"rfcav:dmat" ABS 1e-14', (maxval(abs(xmat(j,:)-branch%ele(4)%mat6(j,:))), j = 1, 6)

! track slice

ele => branch%ele(0)
call init_coord (start_orb, start_orb, branch%ele(5), downstream_end$)
start_orb%s = branch%ele(4)%s + 0.1
start_orb%location = inside$
call twiss_and_track_from_s_to_s (branch, start_orb, branch%ele(6)%s-0.1, end_orb, ele, ele)
write (1, '(a, 2es18.9)') '"rfcav-slice:beta" ABS 1e-8', ele%a%beta, ele%b%beta
write (1, '(a, 2es18.9)') '"rfcav-slice:eta"  ABS 1e-8', ele%a%eta, ele%b%eta

!--

do i = 1, branch%n_ele_max
  if (branch%ele(i)%key /= rfcavity$) cycle
  branch%ele(i)%tracking_method = runge_kutta$
  branch%ele(i)%mat6_calc_method = tracking$
enddo

! Track through split rfcav
call track1 (start_orb, branch%ele(1), branch%param, end_orb)
call track1 (end_orb, branch%ele(3), branch%param, end_orb)

! Track through unsplit rfcav
call track1 (start_orb, branch%ele(4), branch%param, end_orb2)

call lat_make_mat6 (lat, ix_branch = branch%ix_branch)
xmat = matmul(branch%ele(3)%mat6, branch%ele(1)%mat6)

write (1, '(a, 6es18.9)') '"rfcav-rk:dvec" ABS 1e-14', end_orb2%vec - end_orb%vec
write (1, '(a, es18.9)')  '"rfcav-rk:dt"   ABS 1e-8 ', c_light * (end_orb2%t - end_orb%t)
write (1, '(a, 6es18.9)') '"rfcav-rk:dmat" ABS 1e-10', (maxval(abs(xmat(j,:)-branch%ele(4)%mat6(j,:))), j = 1, 6)


!------------------------------------------------
! Test branch 1

branch => lat%branch(1)
call init_coord (start_orb, ele = branch%ele(1), element_end = upstream_end$)

! Track through split bend
call track1 (start_orb, branch%ele(1), branch%param, end_orb)
call track1 (end_orb, branch%ele(3), branch%param, end_orb)

! Track through unsplit bend
call track1 (start_orb, branch%ele(4), branch%param, end_orb2)

xmat = matmul(branch%ele(3)%mat6, branch%ele(1)%mat6)

write (1, *)
write (1, '(a, 6es18.9)') '"bend:dvec" ABS 1e-14', end_orb2%vec - end_orb%vec
write (1, '(a, es18.9)')  '"bend:dt"   ABS 1e-8 ', c_light * (end_orb2%t - end_orb%t)
write (1, '(a, 6es18.9)') '"bend:dmat" ABS 1e-14', (maxval(abs(xmat(j,:)-branch%ele(4)%mat6(j,:))), j = 1, 6)

! track slice

call reallocate_coord (ref_orb, lat, branch%ix_branch)
call twiss_and_track (lat, ref_orb, ix_branch = branch%ix_branch)

end_ele = branch%ele(1)
call init_coord (start_orb, start_orb, branch%ele(1), downstream_end$)
s1 = branch%ele(6)%s-0.1
call twiss_and_track_from_s_to_s (branch, start_orb, s1, end_orb, end_ele, end_ele, err, .true.)
write (1, '(a, 2es18.9)') '"bend-slice:beta"       ABS 1e-8', end_ele%a%beta, end_ele%b%beta
write (1, '(a, 2es18.9)') '"bend-slice:eta"        ABS 1e-8', end_ele%a%eta, end_ele%b%eta
write (1, '(a, 3es18.9)') '"bend-slice:floor"      ABS 1e-10', end_ele%floor%r
write (1, '(a, 3es18.9)') '"bend-slice:floor-ang"  ABS 1e-10', end_ele%floor%theta, end_ele%floor%phi, end_ele%floor%psi

call twiss_and_track_at_s (lat, s1, end_ele, ref_orb, end_orb, branch%ix_branch, err, .false., .true.)
write (1, '(a, 3es18.9)') '"bend-slice:floor2"     ABS 1e-10', end_ele%floor%r
write (1, '(a, 3es18.9)') '"bend-slice:floor-ang2" ABS 1e-10', end_ele%floor%theta, end_ele%floor%phi, end_ele%floor%psi

!------------------------------------------------
! Test branch 2

branch => lat%branch(2)
call reallocate_coord (ref_orb, lat, branch%ix_branch)
ref_orb = lat%beam_start

s1 = 0.5_rp
s2 = 2.5_rp

call track_all (lat, ref_orb, branch%ix_branch)
call lat_make_mat6 (lat, -1, ref_orb, branch%ix_branch)
call twiss_propagate_all (lat, branch%ix_branch)

call twiss_and_track_at_s (lat, s1, ele1, ref_orb, orb1, branch%ix_branch)
call twiss_and_track_at_s (lat, s2, ele2a, ref_orb, orb2a, branch%ix_branch)

orb2c = orb1
call mat6_from_s_to_s (lat, xmat_c, vec0_c, s1, s2, orb2c, branch%ix_branch)

idum = element_at_s (lat, s1, .true., branch%ix_branch, position = orb1)
call twiss_and_track_from_s_to_s (branch, orb1, s2, orb2b, ele1, ele2b)

call transfer_map_from_s_to_s (lat, t_map, s1, s2, ix_branch = branch%ix_branch)
call taylor_to_mat6 (t_map, orb1%vec, vec0_d, xmat_d, orb2d%vec)

write (1, *)
write (1, '(a, es22.12)') '"vec(1)" REL  1E-10', orb2a%vec(1)
write (1, '(a, es22.12)') '"vec(2)" REL  1E-10', orb2a%vec(2)
write (1, '(a, es22.12)') '"vec(3)" REL  1E-10', orb2a%vec(3)
write (1, '(a, es22.12)') '"vec(4)" REL  1E-10', orb2a%vec(4)
write (1, '(a, es22.12)') '"vec(5)" REL  1E-10', orb2a%vec(5)
write (1, '(a, es22.12)') '"vec(6)" REL  1E-10', orb2a%vec(6)
write (1, '(a, es22.12)') '"t"      REL  1E-10', orb2a%t
write (1, '(a, es22.12)') '"s"      REL  1E-10', orb2a%s

write (1, *)
write (1, '(a, es22.12)') '"Db:vec(1)" ABS  1e-14', orb2b%vec(1) - orb2a%vec(1)
write (1, '(a, es22.12)') '"Db:vec(2)" ABS  1e-14', orb2b%vec(2) - orb2a%vec(2)
write (1, '(a, es22.12)') '"Db:vec(3)" ABS  1e-14', orb2b%vec(3) - orb2a%vec(3)
write (1, '(a, es22.12)') '"Db:vec(4)" ABS  1e-14', orb2b%vec(4) - orb2a%vec(4)
write (1, '(a, es22.12)') '"Db:vec(5)" ABS  1e-14', orb2b%vec(5) - orb2a%vec(5)
write (1, '(a, es22.12)') '"Db:vec(6)" ABS  1e-14', orb2b%vec(6) - orb2a%vec(6)
write (1, '(a, es22.12)') '"Db:c*t"    ABS  1e-14', c_light * (orb2b%t - orb2a%t)
write (1, '(a, es22.12)') '"Db:s"      ABS  1e-14', orb2b%s - orb2a%s

write (1, *)
write (1, '(a, es22.12)') '"Dc:vec(1)" ABS  1e-14', orb2c%vec(1) - orb2a%vec(1)
write (1, '(a, es22.12)') '"Dc:vec(2)" ABS  1e-14', orb2c%vec(2) - orb2a%vec(2)
write (1, '(a, es22.12)') '"Dc:vec(3)" ABS  1e-14', orb2c%vec(3) - orb2a%vec(3)
write (1, '(a, es22.12)') '"Dc:vec(4)" ABS  1e-14', orb2c%vec(4) - orb2a%vec(4)
write (1, '(a, es22.12)') '"Dc:vec(5)" ABS  1e-14', orb2c%vec(5) - orb2a%vec(5)
write (1, '(a, es22.12)') '"Dc:vec(6)" ABS  1e-14', orb2c%vec(6) - orb2a%vec(6)
write (1, '(a, es22.12)') '"Dc:c*t"    ABS  1e-14', c_light* (orb2c%t - orb2a%t)
write (1, '(a, es22.12)') '"Dc:s"      ABS  1e-14', orb2c%s - orb2a%s

write (1, *)
write (1, '(a, es22.12)') '"Dd:vec(1)" ABS  1e-14', orb2d%vec(1) - orb2a%vec(1)
write (1, '(a, es22.12)') '"Dd:vec(2)" ABS  1e-14', orb2d%vec(2) - orb2a%vec(2)
write (1, '(a, es22.12)') '"Dd:vec(3)" ABS  1e-14', orb2d%vec(3) - orb2a%vec(3)
write (1, '(a, es22.12)') '"Dd:vec(4)" ABS  1e-14', orb2d%vec(4) - orb2a%vec(4)
write (1, '(a, es22.12)') '"Dd:vec(5)" ABS  1e-14', orb2d%vec(5) - orb2a%vec(5)
write (1, '(a, es22.12)') '"Dd:vec(6)" ABS  1e-14', orb2d%vec(6) - orb2a%vec(6)

write (1, *)
write (1, '(a, es22.12)') '"a%beta " REL  1E-10', ele2a%a%beta
write (1, '(a, es22.12)') '"b%beta " REL  1E-10', ele2b%b%beta
write (1, '(a, es22.12)') '"a%alpha" REL  1E-10', ele2a%a%alpha
write (1, '(a, es22.12)') '"b%alpha" REL  1E-10', ele2b%b%alpha
write (1, '(a, es22.12)') '"a%eta  " REL  1E-10', ele2a%a%eta
write (1, '(a, es22.12)') '"b%eta  " REL  1E-10', ele2b%b%eta

write (1, *)

write (1, '(a, es22.12)') '"D:a%beta"  ABS  1e-14', ele2b%a%beta  - ele2a%a%beta
write (1, '(a, es22.12)') '"D:b%beta"  ABS  1e-14', ele2b%b%beta  - ele2a%b%beta
write (1, '(a, es22.12)') '"D:a%alpha" ABS  1e-14', ele2b%a%alpha - ele2a%a%alpha
write (1, '(a, es22.12)') '"D:b%alpha" ABS  1e-14', ele2b%b%alpha - ele2a%b%alpha
write (1, '(a, es22.12)') '"D:a%eta"   ABS  1e-14', ele2b%a%eta   - ele2a%a%eta
write (1, '(a, es22.12)') '"D:b%eta"   ABS  1e-14', ele2b%b%eta   - ele2a%b%eta

write (1, *)

write (1, '(a, 6f17.12)') '"xmat_c(1,:)" ABS  1e-10', xmat_c(1,:)
write (1, '(a, 6f17.12)') '"xmat_c(2,:)" ABS  1e-10', xmat_c(2,:)
write (1, '(a, 6f17.12)') '"xmat_c(3,:)" ABS  1e-10', xmat_c(3,:)
write (1, '(a, 6f17.12)') '"xmat_c(4,:)" ABS  1e-10', xmat_c(4,:)
write (1, '(a, 6f17.12)') '"xmat_c(5,:)" ABS  1e-10', xmat_c(5,:)
write (1, '(a, 6f17.12)') '"xmat_c(6,:)" ABS  1e-10', xmat_c(6,:)
write (1, '(a, 6f17.12)') '"vec0_c(:)"   ABS  1e-10', vec0_c

write (1, *)

write (1, '(a, 6f17.12)') '"Db:xmat_c(1,:)" ABS  1e-10', xmat_c(1,:) - ele2b%mat6(1,:)
write (1, '(a, 6f17.12)') '"Db:xmat_c(2,:)" ABS  1e-10', xmat_c(2,:) - ele2b%mat6(2,:)
write (1, '(a, 6f17.12)') '"Db:xmat_c(3,:)" ABS  1e-10', xmat_c(3,:) - ele2b%mat6(3,:)
write (1, '(a, 6f17.12)') '"Db:xmat_c(4,:)" ABS  1e-10', xmat_c(4,:) - ele2b%mat6(4,:)
write (1, '(a, 6f17.12)') '"Db:xmat_c(5,:)" ABS  1e-10', xmat_c(5,:) - ele2b%mat6(5,:)
write (1, '(a, 6f17.12)') '"Db:xmat_c(6,:)" ABS  1e-10', xmat_c(6,:) - ele2b%mat6(6,:)
write (1, '(a, 6f17.12)') '"Db:vec0_c(:)"   ABS  1e-10', vec0_c - ele2b%vec0

write (1, *)
write (1, '(a, 6f17.12)') '"Dd:xmat_c(1,:)" ABS  1e-10', xmat_c(1,:) - xmat_d(1,:)
write (1, '(a, 6f17.12)') '"Dd:xmat_c(2,:)" ABS  1e-10', xmat_c(2,:) - xmat_d(2,:)
write (1, '(a, 6f17.12)') '"Dd:xmat_c(3,:)" ABS  1e-10', xmat_c(3,:) - xmat_d(3,:)
write (1, '(a, 6f17.12)') '"Dd:xmat_c(4,:)" ABS  1e-10', xmat_c(4,:) - xmat_d(4,:)
write (1, '(a, 6f17.12)') '"Dd:xmat_c(5,:)" ABS  1e-10', xmat_c(5,:) - xmat_d(5,:)
write (1, '(a, 6f17.12)') '"Dd:xmat_c(6,:)" ABS  1e-10', xmat_c(6,:) - xmat_d(6,:)
write (1, '(a, 6f17.12)') '"Dd:vec0_c(:)"   ABS  1e-10', vec0_c - vec0_d

!------------------------------------------------
! Test branch 3

branch => lat%branch(3)
call init_coord (start_orb, lat%beam_start, branch%ele(0), downstream_end$)

s_end = branch%ele(1)%s
orbit = start_orb
ele1 = branch%ele(0)
do i = 1, 100
  call twiss_and_track_from_s_to_s (branch, orbit, i * s_end / 100, orbit, ele1, ele1)
enddo

call reallocate_coord (ref_orb, lat, branch%ix_branch)
ref_orb(0) = start_orb
call track_all (lat, ref_orb, branch%ix_branch)
call lat_make_mat6(lat, -1, ref_orb, branch%ix_branch)
call twiss_propagate_all (lat, branch%ix_branch)

end_orb = ref_orb(branch%n_ele_track)
end_ele = branch%ele(branch%n_ele_track)

write (1, '(a, 2es22.12)') '"B2:vec(1)" ABS  1e-14', end_orb%vec(1), end_orb%vec(1) - orbit%vec(1)
write (1, '(a, 2es22.12)') '"B2:vec(2)" ABS  1e-14', end_orb%vec(2), end_orb%vec(2) - orbit%vec(2)
write (1, '(a, 2es22.12)') '"B2:vec(3)" ABS  1e-14', end_orb%vec(3), end_orb%vec(3) - orbit%vec(3)
write (1, '(a, 2es22.12)') '"B2:vec(4)" ABS  1e-14', end_orb%vec(4), end_orb%vec(4) - orbit%vec(4)
write (1, '(a, 2es22.12)') '"B2:vec(5)" ABS  2e-14', end_orb%vec(5), end_orb%vec(5) - orbit%vec(5)
write (1, '(a, 2es22.12)') '"B2:vec(6)" ABS  4e-14', end_orb%vec(6), end_orb%vec(6) - orbit%vec(6)
write (1, '(a, 2es22.12)') '"B2:c*t"    ABS  1e-14', c_light * end_orb%t, c_light * (end_orb%t - orbit%t)

write (1, *)
write (1, '(a, 2f22.14)') '"B2:D:a%beta"  ABS  1e-12', end_ele%a%beta,  end_ele%a%beta  - ele1%a%beta
write (1, '(a, 2f22.14)') '"B2:D:b%beta"  ABS  1e-12', end_ele%b%beta,  end_ele%b%beta  - ele1%b%beta
write (1, '(a, 2f22.14)') '"B2:D:a%alpha" ABS  1e-12', end_ele%a%alpha, end_ele%a%alpha - ele1%a%alpha
write (1, '(a, 2f22.14)') '"B2:D:b%alpha" ABS  1e-12', end_ele%b%alpha, end_ele%b%alpha - ele1%b%alpha
write (1, '(a, 2f22.14)') '"B2:D:a%eta"   ABS  1e-12', end_ele%a%eta,   end_ele%a%eta   - ele1%a%eta
write (1, '(a, 2f22.14)') '"B2:D:b%eta"   ABS  1e-12', end_ele%b%eta,   end_ele%b%eta   - ele1%b%eta

!------------------------------------------------
! Test slice_test2.bmad

call bmad_parser ('slice_test2.bmad', lat, .false.)
branch => lat%branch(0)
ele => branch%ele(3)

orb1%vec(1) = 0.8
orb1%vec(2) = 0.2
orb1%vec(3) = -0.2
call em_field_calc (ele, branch%param, 0.1_rp, orb1, .true., field)

write (1, *)
write (1, '(a, 3f14.8)') '"Field%E"  ABS  1e-12', field%e
write (1, '(a, 3f14.8)') '"Field%B"  ABS  1e-12', field%b



!------------------------------------------------

close (1)


end program 
