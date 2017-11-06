!+
! Program parse_test
!
! This program is part of the Bmad regression testing suite.
!-

program parse_test

use bmad
use bmad_parser_mod
use write_lat_file_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct) ele2
type (all_pointer_struct) a_ptr
type (coord_struct), allocatable :: orbit(:)
type (coord_struct) orb
type (em_field_struct) field
type (ele_pointer_struct), allocatable :: eles(:)

real(rp) value
integer i, j, inc_version, n_loc
character(200) digested_file
character(1) delim
logical err, delim_found

! 

open (1, file = 'output.now')

!

orb%vec = [0.1_rp, 0.2_rp, 0.3_rp, 0.4_rp, 0.5_rp, 0.6_rp]
orb%species = electron$

call bmad_parser ('em_field.bmad', lat)

do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (ele%key == taylor$) cycle
  if (ele%key == marker$) cycle
  call em_field_calc (ele, lat%param, 0.1_rp, orb, .false., field, .true., rf_time = 1.0_rp)
  write (1, '(a, i0, a, 6es16.8)') '"Field:', i, '" REL 1e-7', field%E, field%B
  do j = 1, 3
    write (1, '(2(a, i0), a, 6es16.8)') '"dField:', i, '-', j, '" REL 1e-7', field%dE(j,:), field%dB(j,:)
  enddo
enddo

call write_bmad_lattice_file ('z.bmad', lat)
call bmad_parser ('z.bmad', lat)
call form_digested_bmad_file_name ('z.bmad', digested_file)
call read_digested_bmad_file (digested_file, lat, inc_version, err)
call write_bmad_lattice_file ('z2.bmad', lat)

do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (ele%key == taylor$) cycle
  if (ele%key == marker$) cycle
  call em_field_calc (ele, lat%param, 0.1_rp, orb, .false., field, .true., rf_time = 1.0_rp)
  write (1, '(a, i0, a, 6es16.8)') '"Field2:', i, '" REL 1e-7', field%E, field%B
  do j = 1, 3
    write (1, '(2(a, i0), a, 6es16.8)') '"dField2:', i, '-', j, '" REL 1e-7', field%dE(j,:), field%dB(j,:)
  enddo
enddo

!!! Test: curved coords...

!

call bmad_parser ('overlap.bmad', lat)

do i = lat%n_ele_track+1, lat%n_ele_max
  write (1, '(a, i0, 3a)')       '"Overlap-Lord', i, '"    STR   "', trim(lat%ele(i)%name), '"'
enddo

lat%branch(0)%ele(1)%key = -1
lat%branch(0)%ele(4)%key = -1
lat%branch(1)%ele(1)%key = -1
lat%branch(2)%ele(1)%key = -1

call remove_eles_from_lat (lat, .true.)
call write_bmad_lattice_file ('overlap_out.bmad', lat)

!

call bmad_parser ('parse_test.bmad', lat)
call write_bmad_lattice_file ('write_parser_test.bmad', lat)
call bmad_parser ('write_parser_test.bmad', lat)

call pointer_to_attribute (lat%ele(1), 'QQQ', .false., a_ptr, err)
write (1, '(a, f8.4)')  '"zzz"                                   ABS 0', a_ptr%r

call set_attribute_alias('CUSTOM_ATTRIBUTE4', 'CALIB', err)
ele2%key = overlay$
call pointer_to_attribute (ele2, 'CALIB', .false., a_ptr, err)
a_ptr%r = 7
write (1, '(a, f8.4)')  '"calib"                                 ABS 0', a_ptr%r

call lat_ele_locator ('gang0', lat, eles, n_loc)
write (1, '(a, 2i4)') '"gang0"  ABS 0 ', n_loc, eles(1)%ele%n_slave

call lat_ele_locator ('gang1', lat, eles, n_loc)
write (1, '(a, 3i4)') '"gang1"  ABS 0 ', n_loc, eles(1)%ele%n_slave, eles(2)%ele%n_slave

!

write (1, '(a, f10.2)') '"bmad_com[max_aperture_limit]"          ABS 0', bmad_com%max_aperture_limit
write (1, '(a, i4)')    '"bmad_com[ptc_max_fringe_order]"        ABS 0', bmad_com%ptc_max_fringe_order
write (1, '(a, l1, a)') '"bmad_com[convert_to_kinetic_momentum]" STR   "', bmad_com%convert_to_kinetic_momentum, '"'
write (1, '(3a)')       '"geometry"                              STR   "', trim(geometry_name(lat%param%geometry)), '"'
write (1, '(3a)')       '"quad-custom_attribute2"                STR   "', trim(attribute_name(quadrupole$, custom_attribute2$)), '"'
write (1, '(3a)')       '"sex-custom_attribute3"                 STR   "', trim(attribute_name(sextupole$, custom_attribute3$)), '"'

bp_com%input_from_file = .false.
bp_com%parse_line = '-2*7)'
call evaluate_value ('ERR', value, lat, delim, delim_found, err, ',)')
write (1, '(a, f10.4)') '"EVAL 1"  ABS 0', value

write (1, '(a, f10.4)') '"1 REL"  ABS 0', lat%ele(1)%value(k1$)
write (1, '(a, f10.4)') '"2 REL"  ABS 0', lat%ele(2)%value(k1$)
write (1, '(a, f10.4)') '"3 REL"  ABS 0', lat%ele(3)%value(k1$)
write (1, '(a, f10.4)') '"4 REL"  ABS 0', lat%ele(4)%value(k1$)
write (1, '(a, f10.4)') '"5 REL"  ABS 0', lat%ele(5)%value(k2$)
write (1, '(4a)')       '"TM1"     STR ', '"', trim(tracking_method_name(lat%ele(1)%tracking_method)), '"'
write (1, '(4a)')       '"TM5"     STR ', '"', trim(tracking_method_name(lat%ele(5)%tracking_method)), '"'
write (1, '(a, i3)')    '"N3"     ABS 0', lat%branch(2)%n_ele_track
write (1, '(4a)')       '"Custom"  STR ', '"', trim(attribute_name(lat%ele(1), custom_attribute1$)), '"'

do i = lbound(mass_of_fundamental, 1), ubound(mass_of_fundamental, 1)
  write (1, '(3a, es20.12, i6)')         '"', trim(species_name(i)), '"  REL 1e-12', mass_of(i), charge_of(i)
enddo

write (1, '(3a)') '"He++"         STR "', trim(species_name(species_id('He++'))), '"'
write (1, '(3a)') '"#12C-5"       STR "', trim(species_name(species_id('#12C-5'))), '"'
write (1, '(3a)') '"CH3@M34.5+2"  STR "', trim(species_name(species_id('CH3@M34.5+2'))), '"'
write (1, '(3a)') '"@M3.45---"    STR "', trim(species_name(species_id('@M3.45---'))), '"'

!

bp_com%always_parse = .true.
call bmad_and_xsif_parser ('DCO4.xsif', lat)
lat%ele(156)%value(hkick$) = 0.00001
lat%ele(156)%value(vkick$) = 0.00002

call twiss_at_start (lat)
call twiss_propagate_all (lat)
call reallocate_coord (orbit, lat)
call closed_orbit_calc (lat, orbit, 4)

write (1, '(a, 2f12.8)')  '"XSIF:Twiss"  REL 1e-8', lat%ele(0)%a%beta, lat%ele(0)%b%beta
write (1, '(a, 2es15.8)') '"XSIF:Orbit"  REL 1e-8', orbit(0)%vec(1), orbit(0)%vec(3)

!

call bmad_parser ('parse_test.bmad', lat, use_line = 'PHOT')

write (1, '(4a)')         '"PHOT-1"    STR ', '"', trim(lat%ele(1)%name), '"'
write (1, '(2a, i0, a)')  '"PHOT-N"    STR ', '"', lat%n_ele_max, '"'
write (1, '(2a, l1, a)')  '"Hard-Edge" STR ', '"', bmad_com%use_hard_edge_drifts, '"'


close(1)

end program
