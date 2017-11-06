!+
! This is a front end for bmad/multiparticle/envelope_mod.f90
!
! The module and this front end currently implement the beam envelope matrix
! based synchrotron radiation and damping calculation contained in
! "From the beam-envelope matrix to synchrotron-radiation integrals" by
! K. Ohmi, K. Hirata, and K. Oide.
!-
program envelope_ibs

use bmad
use transfer_map_mod
use mode3_mod
use envelope_mod

implicit none

character lat_file*200
character in_file*200
character ix_str*5

type(lat_struct) lat
type(ele_struct) temp_ele
type(ele_struct), allocatable :: eles(:)
type(coord_struct), allocatable :: co(:)
type(coord_struct), allocatable :: coos(:)
type(normal_modes_struct) mode

integer i,j
integer nturns, nslices, ns, six
integer status
integer pct_complete, last_display

logical err_flag, include_ibs, tail_cut

real(rp) one_turn_mat(6,6), one_turn_vec(6)
real(rp) alpha(3), emit(3)
real(rp) mat6(6,6), vec0(6)
real(rp) Sigma_ent(6,6), Sigma_exit(6,6)
real(rp) normal(3)
real(rp) current, npart
real(rp) bend_slice_length, slice_length
real(rp) tau(3), tau_max, ey0, Ykick_strength
real(rp) starting_a_emit, starting_b_emit, starting_c_emit
real(rp), allocatable :: M(:,:,:), Bbar(:,:,:), Ybar(:,:,:)

complex(rp) Lambda(6,6), Theta(6,6), Iota_base(6,6), Iota(6,6)

namelist /envelope_tracker/ lat_file, starting_a_emit, starting_b_emit, starting_c_emit, ey0, nturns, include_ibs, current, tail_cut, bend_slice_length, slice_length

!set defaults
ey0 = -1.0 !negative number to disable
bend_slice_length = 100.0d0
slice_length = 100.0d0  !large value results in element length being used
starting_a_emit = 1.0d-9
starting_b_emit = 10.0d-12
starting_c_emit = 0.01 * 1.0d-4
nturns = 40000
include_ibs = .false.
tail_cut = .true.
current = 0.001

call getarg(1, in_file)
open (unit = 20, file = in_file, action='read')
read (20, nml = envelope_tracker)
close (20)

bmad_com%radiation_damping_on = .true.

call bmad_parser(lat_file, lat)
call twiss_and_track(lat,co,status)
call calc_z_tune(lat)

npart = current / e_charge * lat%param%total_length / c_light

!make_smat_from_abc uses mode tunes to label modes
mode%a%tune = lat%a%tune
mode%b%tune = lat%b%tune
mode%z%tune = lat%z%tune
mode%a%emittance = starting_a_emit
mode%b%emittance = starting_b_emit
mode%z%emittance = starting_c_emit
call transfer_matrix_calc(lat, mat6, vec0, ix1=0, one_turn=.true.)
call make_smat_from_abc(mat6, mode, Sigma_ent, err_flag)

write(*,'(a,3es14.5)') "Initial Emittances: ", mode%a%emittance, mode%b%emittance, mode%z%emittance
open(10,file='sigma_start.out')
do i=1,6
  write(10,'(6es14.4)') Sigma_ent(i,:)
enddo
close(10)

! Build element slices and damping matrix D and diffusion matrix B for each slice.
nslices = 0
do i=1,lat%n_ele_track
  if(any(lat%ele(i)%key==(/marker$,monitor$/))) then
    cycle
  else
    if(any(lat%ele(i)%key==(/sbend$,rbend$/))) then
      ns = max(1,ceiling(lat%ele(i)%value(l$) / bend_slice_length))
    else
      ns = max(1,ceiling(lat%ele(i)%value(l$) / slice_length))
    endif
    nslices = nslices + ns
  endif
enddo
write(*,'(a,i8)') "Number of slices: ", nslices
allocate(eles(nslices))
allocate(coos(0:nslices))
allocate(M(6,6,nslices))
allocate(Bbar(6,6,nslices))
allocate(Ybar(6,6,nslices))
six = 0
coos(0)=co(0)
do i=1,lat%n_ele_track
  if(any(lat%ele(i)%key==(/marker$,monitor$/))) then
    cycle
  else
    if(any(lat%ele(i)%key==(/sbend$,rbend$/))) then
      ns = max(1,ceiling(lat%ele(i)%value(l$) / bend_slice_length))
    else
      ns = max(1,ceiling(lat%ele(i)%value(l$) / slice_length))
    endif
    if(ns .gt. 1) then
      do j=1,ns
        six = six + 1
        call create_uniform_element_slice(lat%ele(i),lat%param,j,ns,temp_ele)
        call make_mat6(temp_ele, lat%param, coos(six-1), coos(six))
        eles(six) = temp_ele
        call make_SR_mats(eles(six),coos(six),M(:,:,six),Bbar(:,:,six))
        call make_Ykick_mat(eles(six),Ybar(:,:,six))
      enddo
    else
      six = six + 1
      coos(six) = co(i) !co at element exit
      eles(six) = lat%ele(i)
      call make_SR_mats(eles(six),coos(six),M(:,:,six),Bbar(:,:,six))
      call make_Ykick_mat(eles(six),Ybar(:,:,six))
    endif
  endif
enddo

one_turn_mat = I6
one_turn_vec = 0.0d0
do i=1,size(eles)
  call concat_transfer_mat(eles(i)%mat6,eles(i)%vec0,one_turn_mat,one_turn_vec,one_turn_mat,one_turn_vec)
enddo
open(50,file='one_turn_mat_from_slices.diag')
write(50,*) "mat6:"
do i=1,6
  write(50,'(6es17.8)') one_turn_mat(i,:)
enddo 
write(50,*) "vec0:"
write(50,'(6es17.8)') one_turn_vec(:)
close(50)

call integrated_mats(eles,coos,Lambda,Theta,Iota_base,mode)
Ykick_strength = ey0 * 2.0d0 * real(Lambda(3,3)) / Iota_base(3,3)
Iota = Iota_base * Ykick_strength !for integral calculations
Ybar = Ybar * Ykick_strength !for tracking

if(include_ibs) then
  call envelope_radints_ibs(Lambda,Theta,Iota,eles,alpha,emit,mode,tail_cut,npart)
else
  call envelope_radints(Lambda,Theta,Iota,alpha,emit)
endif
tau = lat%param%total_length / c_light / alpha

open(20,file='emit_from_integrals.out')
do i=6,20,14
  write(i,*) "From envelope-based radiation integrals and IBS calculations: "
  write(i,'(a,3es14.5)') "  Emittance (SR only):        ", emit
  if(include_ibs) write(i,'(a,3es14.5)') "  Emittance (with IBS):       ", mode%a%emittance, mode%b%emittance, mode%z%emittance
  write(i,'(a,3es14.5)') "  Damping decrements (alpha): ", alpha
  write(i,'(a,3es14.5)') "  Damping times:              ", tau
enddo
close(20)

if(nturns .gt. 0) then
  write(*,*) "Starting tracking..."
  tau_max = maxval(tau)
  ! Track element-slice by element-slice for number of turns.
  open(11,file='emit_vs_turn.out')
  write(11,'(a8,3a14)') "# turn", "emit_a", "emit_b", "emit_c"
  last_display = -1.0
  do i=1,nturns
    do j=1,nslices
      if(include_ibs) then
        call transport_with_sr_and_ibs(eles(j),M(:,:,j),Bbar(:,:,j),Ybar(:,:,j),Sigma_ent,Sigma_exit,tail_cut,tau_max,npart)
      else
        call transport_with_sr(eles(j),M(:,:,j),Bbar(:,:,j),Ybar(:,:,j),Sigma_ent,Sigma_exit) 
      endif
      Sigma_ent = Sigma_exit
    enddo
    if( mod(i,10) == 0 ) then
      call get_emit_from_sigma_mat(Sigma_exit, normal, err_flag = err_flag)
      write(11,'(i8,3es14.5)') i, normal(1:3)
    endif
    pct_complete = floor((i*100.0d0)/nturns)
    if( pct_complete .gt. last_display ) then
      write(*,'(i8,a,i8,a)') i, " turns of ", nturns, " complete."
      last_display = pct_complete
    endif
  enddo
  close(11)

  ! Extract normal mode emittances from final sigma matrix.
  call get_emit_from_sigma_mat(Sigma_exit, normal, err_flag = err_flag)

  write(*,*)
  write(*,'(a)') "From tracking beam enevelope: "
  write(*,'(a,3es14.5)') "  Emittance:        ", normal
  write(*,*) "Beam sigma matrix saved as sigma_finish.out"

  open(10,file='sigma_finish.out')
  do i=1,6
    write(10,'(6es14.4)') Sigma_exit(i,:)
  enddo
  close(10)
else
  write(*,*) "Tracking disabled in .in file.  nturns < 0"
endif

deallocate(eles)
deallocate(coos)

end program












