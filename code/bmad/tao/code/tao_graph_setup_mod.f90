module tao_graph_setup_mod

use tao_mod
use tao_lattice_calc_mod
use tao_command_mod
use tao_data_and_eval_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_setup (plot, graph)

implicit none

type (tao_universe_struct), pointer :: u
type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve

integer i, iu
logical found

!

graph%valid = .true.   ! assume everything OK
graph%why_invalid = ''
graph%text_legend = ''

if (allocated (graph%curve)) then
  do i = 1, size(graph%curve)
    graph%curve(i)%message_text = ''
  enddo
endif

call tao_hook_graph_setup (plot, graph, found)
if (found) return

!

if (allocated (graph%curve)) then
  do i = 1, size(graph%curve)
    curve => graph%curve(i)
    u => tao_pointer_to_universe(tao_curve_ix_uni(curve))
    if (.not. u%is_on) then
      graph%valid = .false.
      write (graph%why_invalid, '(a, i0, a)') 'UNIVERSE ', u%ix_uni, ' IS OFF!'
      return
    endif
  enddo
else
  u => tao_pointer_to_universe(graph%ix_universe)
  if (.not. u%is_on) then
    graph%valid = .false.
    write (graph%why_invalid, '(a, i0, a)') 'UNIVERSE ', u%ix_uni, ' IS OFF!'
    return
  endif
endif

!

select case (graph%type)
case ('phase_space')
  call tao_graph_phase_space_setup (plot, graph)

case ('data', 'lat_layout')
  if (plot%x_axis_type == 'data') then
    call tao_graph_data_slice_setup(plot, graph)
  else
    call tao_graph_data_setup(plot, graph)
  endif

case ('histogram')
  call tao_graph_histogram_setup (plot, graph)

case ('dynamic_aperture')
  call tao_graph_dynamic_aperture_setup (plot, graph)

end select

! Renormalize

if (allocated (graph%curve)) then
  do i = 1, size(graph%curve)
    curve => graph%curve(i)
    if (allocated(curve%x_symb)) then
        curve%x_symb = curve%x_symb * curve%g%x_axis_scale_factor
        curve%y_symb = curve%y_symb * curve%y_axis_scale_factor
    endif
    if (allocated(curve%x_line)) then
      curve%x_line = curve%x_line * curve%g%x_axis_scale_factor
      curve%y_line = curve%y_line * curve%y_axis_scale_factor
    endif
  enddo
endif

call tao_hook_graph_postsetup (plot, graph)

end subroutine tao_graph_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_data_slice_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve

real(rp), pointer :: symb(:)
real(rp) value

integer i, j, k, m, n_symb, ix

character(160) name
character(60) component
character(40) :: r_name = 'tao_graph_data_slice_setup'

logical err

! setup the graph suffix

graph%valid = .false.
if (size(graph%curve) == 0) return

graph%title_suffix = ''

do k = 1, size(graph%curve)
  curve => graph%curve(k)
  if (index(curve%data_type, '#ref') /= 0) graph%title_suffix = &
                  trim(graph%title_suffix) // '[At: ' // trim(curve%ele_ref_name) // ']'
enddo

if (graph%component /= '') graph%title_suffix = trim(graph%title_suffix) // ' [' // trim(graph%component) // ']'

! loop over all curves

curve_loop: do k = 1, size(graph%curve)

  curve => graph%curve(k)
  component = tao_curve_component(curve, graph)

  ! Find data points

  do i = 1, 2  ! x-axis, y-axis
    if (i == 1) name = curve%data_type_x
    if (i == 2) name = curve%data_type
    call tao_data_type_substitute (name, name, curve, graph)
    if (i == 1) call tao_evaluate_expression (name, 0, graph%draw_only_good_user_data_or_vars, scratch%x, scratch%info_x, err, dflt_component = component)
    if (i == 2) call tao_evaluate_expression (name, 0, graph%draw_only_good_user_data_or_vars, scratch%y, scratch%info_y, err, dflt_component = component)
    if (err) then
      call out_io (s_error$, r_name, 'CANNOT FIND DATA ARRAY TO PLOT CURVE: ' // tao_curve_name(curve))   
      return
    endif
  enddo

  ! How many good points?

  if (size(scratch%x) /= size(scratch%y)) then
    call out_io (s_error$, r_name, &
                  'ARRAY SIZES ARE NOT THE SAME FOR BOTH AXES.', &
                  'FOR: ' // tao_curve_name(curve))
    return
  endif

  n_symb = count(scratch%info_x%good .and. scratch%info_y%good)

  call re_allocate (curve%x_symb, n_symb)
  call re_allocate (curve%y_symb, n_symb)
  call re_allocate (curve%ix_symb, n_symb)

  ! Transfer the values

  curve%x_symb = pack (scratch%x, mask = scratch%info_x%good .and. scratch%info_y%good)
  curve%y_symb = pack (scratch%y, mask = scratch%info_x%good .and. scratch%info_y%good)

  ! Calc symbol index

  if (curve%data_index == '') then
    curve%ix_symb = [(i, i = 1, n_symb)]
  else
    call tao_data_type_substitute (curve%data_index, name, curve, graph)
    call tao_evaluate_expression (name, 0, graph%draw_only_good_user_data_or_vars, scratch%x, scratch%info_ix, err, dflt_component = component)
    if (size(scratch%info_x) == size(scratch%info_y)) then
      curve%ix_symb = pack (nint(scratch%x), mask = scratch%info_x%good .and. scratch%info_y%good)
    else
      call out_io (s_error$, r_name, &
          'SIZE OF SYMBOL INDEX ARRAY IS WRONG IN CURVE: ' // tao_curve_name(curve), &
          'CURVE%DATA_INDEX: ' // curve%data_index)
    endif
  endif

  ! The line data just goes through the symbols

  call re_allocate (curve%x_line, n_symb)
  call re_allocate (curve%y_line, n_symb)
  curve%x_line = curve%x_symb
  curve%y_line = curve%y_symb

enddo curve_loop

graph%valid = .true.

end subroutine tao_graph_data_slice_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_data_type_substitute (template, str_out, curve, graph)
!
! Routine substitute the appropriate data type string for instances of "#ref" and 
! "#comp" in template. 
!
! Additionally, if template does not have a "|" character,
! the string "|" + component will be added at the end of str_out.
!
! Input:
!   templace    -- character(*): String template.
!   curve       -- tao_curve_struct: curve%ele_ref_name is substituted for all instances of "#ref".
!   graph       -- tao_graph_struct: 
!
! Output:
!   str_out     -- character(*): String with substitutions.
!-

subroutine tao_data_type_substitute (template, str_out, curve, graph)

implicit none

type (tao_curve_struct) curve
type (tao_graph_struct) graph
character(*) template, str_out
character(60) component
integer ix

!

str_out = template
component = tao_curve_component(curve, graph)

do
  ix = index(str_out, '#ref')
  if (ix == 0) exit
  str_out = trim(str_out(:ix-1)) // trim(curve%ele_ref_name) // trim(str_out(ix+4:))
enddo

do
  ix = index(str_out, '#comp')
  if (ix == 0) exit
  str_out = trim(str_out(:ix-1)) // trim(component) // trim(str_out(ix+5:))
enddo

if (index(str_out, '|') == 0) str_out = trim(str_out) // '|' // component

end subroutine tao_data_type_substitute

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_phase_space_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u
type (ele_struct), pointer :: ele
type (beam_struct), pointer :: beam
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_x, d1_y
type (coord_struct), pointer :: p(:)

real(rp) v_mat(4,4), v_inv_mat(4,4), g_mat(4,4), g_inv_mat(4,4)
real(rp) mat4(4,4), sigma_mat(4,4), theta, theta_xy, rx, ry, phi
real(rp) emit_a, emit_b

integer k, n, m, ib, ix1_ax, ix2_ax, ix3_ax, ix, i
integer n_curve_pts

logical err, same_uni

character(40) name
character(40) :: r_name = 'tao_graph_phase_space_setup'

! Set up the graph suffix

graph%valid = .false.
if (size(graph%curve) == 0) return

n_curve_pts = s%plot_page%n_curve_pts
if (plot%n_curve_pts > 0) n_curve_pts = plot%n_curve_pts

same_uni = .true.
ix = tao_universe_number(tao_curve_ix_uni(graph%curve(1)))
do k = 2, size(graph%curve)
  curve => graph%curve(k)
  if (tao_universe_number(tao_curve_ix_uni(curve)) /= ix) same_uni = .false.
enddo

graph%title_suffix = ''
do k = 1, size(graph%curve)
  curve => graph%curve(k)
  u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
  if (.not. associated(u)) return
  if (curve%ix_ele_ref_track < 0) then
    call out_io (s_error$, r_name, &
                'BAD REFERENCE ELEMENT: ' // curve%ele_ref_name, &
                'CANNOT PLOT PHASE SPACE FOR: ' // tao_curve_name(curve))
    return
  endif
  ele => u%model%lat%ele(curve%ix_ele_ref_track)
  name = curve%ele_ref_name
  if (name == '') name = ele%name
  if (same_uni) then
    write (graph%title_suffix, '(2a, i0, 3a)') trim(graph%title_suffix), &
                                '[', curve%ix_ele_ref, ': ', trim(name), ']'
  else
    write (graph%title_suffix, '(2a, i0, a, i0, 3a)') trim(graph%title_suffix), &
            '[', u%ix_uni, '@', curve%ix_ele_ref, ': ', trim(name), ']'
  endif
enddo

! loop over all curves

do k = 1, size(graph%curve)

  curve => graph%curve(k)
  u => tao_pointer_to_universe (tao_curve_ix_uni(curve))

  ! find phase space axes to plot

  err = .false.
  call tao_phase_space_axis (curve%data_type_x, ix1_ax, err = err); if (err) return
  call tao_phase_space_axis (curve%data_type,   ix2_ax, err = err); if (err) return

  ! fill the curve data arrays

  if (allocated (curve%ix_symb)) deallocate (curve%ix_symb, curve%x_symb, curve%y_symb)
  if (allocated (curve%x_line))  deallocate (curve%x_line, curve%y_line)

  !----------------------------

  if (curve%data_source == 'beam') then
    beam => u%uni_branch(curve%ix_branch)%ele(curve%ix_ele_ref_track)%beam
    ele => u%model%lat%branch(curve%ix_branch)%ele(curve%ix_ele_ref_track)
    if (.not. allocated(beam%bunch)) then
      call out_io (s_abort$, r_name, 'NO ALLOCATED BEAM WITH PHASE_SPACE PLOTTING.')
      if (.not. u%is_on) call out_io (s_blank$, r_name, '   REASON: UNIVERSE IS TURNED OFF!')
      return
    endif

    if (curve%ix_bunch == 0) then
      n = 0
      do ib = 1,  size(beam%bunch)
        n = n + count(beam%bunch(ib)%particle%state == alive$)
      enddo
    else
      n = count(beam%bunch(curve%ix_bunch)%particle%state == alive$)
    endif

    call re_allocate (curve%ix_symb, n)
    call re_allocate (curve%x_symb, n)
    call re_allocate (curve%y_symb, n)
    if (graph%symbol_size_scale > 0) call re_allocate (curve%symb_size, n)
    if (curve%use_z_color) call re_allocate (curve%z_symb, n)

    n = 0
    do ib = 1, size(beam%bunch)
      if (curve%ix_bunch /= 0 .and. curve%ix_bunch /= ib) cycle
      p => beam%bunch(ib)%particle
      m = count(beam%bunch(ib)%particle%state == alive$)
      call tao_phase_space_axis (curve%data_type_x, ix1_ax, p, scratch%axis1, ele)
      call tao_phase_space_axis (curve%data_type,   ix2_ax, p, scratch%axis2, ele)
      curve%x_symb(n+1:n+m) = pack(scratch%axis1, mask = (p%state == alive$))
      curve%y_symb(n+1:n+m) = pack(scratch%axis2, mask = (p%state == alive$))
      if (curve%use_z_color) then
        call tao_phase_space_axis (curve%data_type_z,   ix3_ax, p, scratch%axis3, ele)
        curve%z_symb(n+1:n+m) = pack(scratch%axis3, mask = (p%state == alive$))
      endif
      if (graph%symbol_size_scale > 0) curve%symb_size(n+1:n+m) = pack(graph%symbol_size_scale * &
                           sqrt(p(:)%field(1)**2 + p(:)%field(2)**2), mask = (p%state == alive$))
      curve%ix_symb(n+1:n+m) = pack([(i, i = 1,m)], mask = (p%state == alive$))
      n = n + count(p%state == alive$)
    enddo

  !----------------------------

  elseif (curve%data_source == 'multi_turn_orbit') then
    
    call tao_find_data (err, curve%data_source, d2_array, ix_uni = tao_curve_ix_uni(curve))
    if (err .or. size(d2_array) /= 1) then
      call out_io (s_error$, r_name, &
                'CANNOT FIND DATA ARRAY TO PLOT CURVE: ' // curve%data_type)
      graph%valid = .false.
      return
    endif

    nullify (d1_x, d1_y)
    d2_ptr => d2_array(1)%d2
    do i = 1, size(d2_ptr%d1)
      if (curve%data_type_x == d2_ptr%d1(i)%name) d1_x => d2_ptr%d1(i)
      if (curve%data_type   == d2_ptr%d1(i)%name) d1_y => d2_ptr%d1(i)
    enddo
    if (.not. associated(d1_x)) then
      call out_io (s_error$, r_name, &
              'CANNOT FIND DATA FOR PHASE SPACE COORDINATE: ' // curve%data_type_x, &
              'FOR CURVE: ' // curve%name)
      return
    endif
    if (.not. associated(d1_y)) then
      call out_io (s_error$, r_name, &
              'CANNOT FIND DATA FOR PHASE SPACE COORDINATE: ' // curve%data_type, &
              'FOR CURVE: ' // curve%name)
      return
    endif

    if (lbound(d1_x%d, 1) /= lbound(d1_y%d, 1) .or. &
                                        ubound(d1_x%d, 1) /= ubound(d1_y%d, 1)) then 
      call out_io (s_error$, r_name, &
              'BOUNDS FOR X-AXIS AND Y-AXIS DATA OF PHASE SPACE PLOTTING MISMATCHED.', &
              'FOR CURVE: ' // curve%name)
      return
    endif

    n = size(d1_x%d)
    call re_allocate (curve%ix_symb, n)
    call re_allocate (curve%x_symb, n)
    call re_allocate (curve%y_symb, n)

    do ib = 1, n
      i = ib + lbound(d1_x%d, 1) - 1
      curve%x_symb(ib) = d1_x%d(i)%model_value
      curve%y_symb(ib) = d1_y%d(i)%model_value
    enddo


  elseif (curve%data_source == 'twiss') then

    n = 2 * n_curve_pts
    call re_allocate (curve%x_line, n)
    call re_allocate (curve%y_line, n)

    call make_v_mats (ele, v_mat, v_inv_mat)
    call make_g_mats (ele, g_mat, g_inv_mat)

    mat4 = matmul(v_mat, g_inv_mat)
    emit_a = u%beam%beam_init%a_emit
    if (emit_a == 0) emit_a = 1e-6  ! default value
    emit_b = u%beam%beam_init%b_emit
    if (emit_b == 0) emit_b = 1e-6  ! default value

    sigma_mat =  0
    sigma_mat(1,1) = emit_a
    sigma_mat(2,2) = emit_a
    sigma_mat(3,3) = emit_b
    sigma_mat(4,4) = emit_b
    sigma_mat = matmul (matmul (mat4, sigma_mat), transpose(mat4))

    if (ix1_ax > 4 .or. ix2_ax > 4) then
      call out_io (s_error$, r_name, &
        'Z OR PZ PHASE SPACE PLOTTING NOT YET IMPLEMENTED FOR "twiss" DATA_SOURCE.')
      return
    endif

    rx = sqrt(sigma_mat(ix1_ax, ix1_ax))
    ry = sqrt(sigma_mat(ix2_ax, ix2_ax))
    write (graph%text_legend(1), '(a, es9.2)') 'emit_a:', emit_a
    write (graph%text_legend(2), '(a, es9.2)') 'emit_b:', emit_b

    if(rx == 0 .or. ry == 0) then
      theta_xy = 0
      write (graph%text_legend(3), '(a, f10.4)') 'Theta_tilt (rad):', 0
    else
      theta_xy =  asin(sigma_mat(ix1_ax, ix2_ax) / (rx * ry))
      phi = 0.5 *atan2((rx**2+ry**2) * sin(2*theta_xy), &
                              (rx**2-ry**2) * cos(2*theta_xy)) - theta_xy
      write (graph%text_legend(3), '(a, f10.4)') 'Theta_tilt (rad):', phi
  endif

    n = 2 * n_curve_pts
    call re_allocate (curve%x_line, n)
    call re_allocate (curve%y_line, n)

    do i = 1, n
      theta = (i-1) * twopi / (n-1)
      curve%x_line(i) = rx * cos(theta)
      curve%y_line(i) = ry * sin(theta + theta_xy)
    enddo

  else
    call out_io (s_abort$, r_name, &
        'INVALID CURVE%DATA_SOURCE: ' // curve%data_source, &
        'FOR CURVE: '// curve%name)
    return
  endif

enddo

graph%valid = .true.

end subroutine tao_graph_phase_space_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_dynamic_aperture_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_curve_struct), allocatable :: temp_curve(:)
type (tao_universe_struct), pointer :: u
type (tao_dynamic_aperture_struct), pointer :: da
type (aperture_scan_struct), pointer :: scan
type (coord_struct), allocatable :: orbit(:)
integer :: i, j, k, n_da, n_da_curve, n_pa_curve, n, nc

logical err

character(40) name
character(40) :: r_name = 'tao_graph_dynamic_aperture_setup'

! Valid?

graph%valid = .false.

! Only the graph's universe for now
u => tao_pointer_to_universe (graph%ix_universe)

da => u%dynamic_aperture
n_da = size(da%scan)
nc = size(graph%curve)

! Count dynamic and physical aperture curves 
n_da_curve = 0
n_pa_curve = 0
do i = 1, nc
  curve => graph%curve(i)
  if (is_physical_aperture(curve)) then
    n_pa_curve = n_pa_curve + 1
  else
     n_da_curve = n_da_curve + 1
     ! Assign a default type 
     if (curve%data_type =='') curve%data_type = 'dynamic_aperture'
  endif
enddo 

if (n_da_curve == 0) then
  write (graph%why_invalid, '(a, i0)') 'NO DYNAMIC APERTURES DEFINED FOR UNIVERSE ', u%ix_uni
  return
endif

! Case where only physical aperture curves are defined
if (nc == n_pa_curve) then
  write (graph%why_invalid, '(a, i0)') 'NO DYNAMIC APERTURES CURVES DEFINED FOR UNIVERSE ', u%ix_uni
  return
endif

if (.not. allocated(da%scan(1)%aperture)) then
  write (graph%why_invalid, '(a, i0)') 'DYNAMIC APERTURE NOT CALCULATED FOR UNIVERSE ', u%ix_uni
  return
endif

! If there aren't enough da curves defined for all of the da scan,, 
! automatically create more da curves based on defined curves,
! looping over defined styles
if ( n_da > n_da_curve) then
  allocate(temp_curve(nc))
  call move_alloc(graph%curve, temp_curve)
  allocate(graph%curve(n_da + n_pa_curve))
  graph%curve(1:nc) = temp_curve(1:nc)
  k = nc + 1
  i = 0
  do 
    i = i+1
    if(i>nc) i = 1
    curve => graph%curve(i)
    ! Skip pa curves
    if (is_physical_aperture(curve)) cycle
    graph%curve(k) = temp_curve(i)
    ! Increment name
    write (graph%curve(k)%name, '(a, i0)') 'c', k
    if (k == n_da + n_pa_curve) exit
    k = k+1
  enddo
endif
! new number of curves
nc = n_da_curve + n_pa_curve

! loop over curves
k=0
do i = 1, nc
  
  curve => graph%curve(i) 
  if (is_physical_aperture(curve)) then
    call tao_curve_physical_aperture_setup(curve)
    cycle
  endif
  k = k + 1 
  if (k > size(da%scan)) cycle
  scan => da%scan(k)
  n = size(scan%aperture)
  
  call re_allocate (curve%x_line, n)
  call re_allocate (curve%y_line, n)
  call re_allocate (curve%x_symb, n)
  call re_allocate (curve%y_symb, n)

  write(curve%legend_text, '(a, f6.2, a)') '\gd:', 100*da%pz(k), ' %'

  ! propagate to ix_ele_ref
  if (curve%ix_ele_ref > 0 ) then
    if (curve%ix_ele_ref > u%model%lat%n_ele_track) then
      call out_io (s_error$, r_name, 'IX_ELE_REF out of range for curve: ' // tao_curve_name(curve))  
      return
    endif 
  
    call reallocate_coord(orbit, curve%ix_ele_ref)
    
    
    do j = 1, n
      orbit(0) = scan%ref_orb
      orbit(0)%vec(1) = orbit(0)%vec(1) + scan%aperture(j)%x
      orbit(0)%vec(3) = orbit(0)%vec(3) + scan%aperture(j)%y
      call track_many (u%model%lat, orbit, 0, curve%ix_ele_ref, +1)
    
      curve%x_line(j) = orbit(curve%ix_ele_ref)%vec(1)
      curve%y_line(j) = orbit(curve%ix_ele_ref)%vec(3)      
    enddo
    ! One more track of ref_orb if centered da is requested  
    if (curve%data_type == 'dynamic_aperture_centered') then
      orbit(0) = scan%ref_orb
      call track_many (u%model%lat, orbit, 0, curve%ix_ele_ref, +1)
      curve%x_line(:) = curve%x_line(:) - orbit(curve%ix_ele_ref)%vec(1)
      curve%y_line(:) = curve%y_line(:) - orbit(curve%ix_ele_ref)%vec(3)
    endif
    
  else
  ! use data at ele 0 
    if (curve%data_type == 'dynamic_aperture_centered') then
      curve%x_line(:) = scan%aperture(:)%x
      curve%y_line(:) = scan%aperture(:)%y
    else
      curve%x_line(:) = scan%aperture(:)%x + scan%ref_orb%vec(1)
      curve%y_line(:) = scan%aperture(:)%y + scan%ref_orb%vec(3)
    endif
  endif
  curve%x_symb = curve%x_line
  curve%y_symb = curve%y_line
enddo

graph%valid = .true.

!------------------------------------------------------------
contains

function is_physical_aperture(curve)
logical :: is_physical_aperture
type (tao_curve_struct) :: curve
is_physical_aperture = (curve%data_type == 'physical_aperture')
end function is_physical_aperture

end subroutine tao_graph_dynamic_aperture_setup


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_curve_physical_aperture_setup (curve)

type (tao_curve_struct) :: curve
type (tao_universe_struct), pointer :: u
type (ele_struct), pointer :: ele
real(rp) :: x1, x2, y1, y2, phi
real(rp) :: large = 1 ! m 
integer :: i
integer, parameter :: n = 25 ! points per elliptical quadrant
character(40) :: r_name = 'tao_curve_physical_aperture_setup'

!

u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
! Set to beginning if x_ele_ref out of bounds 
if (curve%ix_ele_ref < 1 .or. curve%ix_ele_ref > u%model%lat%n_ele_track) then
  call out_io (s_warn$, r_name, 'IX_ELE_REF out of bounds for: ' // tao_curve_name(curve) //', setting to 1')
  curve%ix_ele_ref = 1
endif
ele => u%model%lat%ele(curve%ix_ele_ref)


x1 = ele%value(x1_limit$)
x2 = ele%value(x2_limit$)
y1 = ele%value(y1_limit$)
y2 = ele%value(y2_limit$)
if (x1 == 0) x1 = large
if (x2 == 0) x2 = large
if (y1 == 0) y1 = large
if (y2 == 0) y2 = large


select case(ele%aperture_type)
case(elliptical$)
call alloc_curves(4*n)
! draw four quadrants
do i=1, n
  phi = pi/2*(i-1)/(n-1)
  curve%x_line(i) = x1*cos(phi) 
  curve%y_line(i) = y1*sin(phi)
enddo
do i=1, n
  phi = pi/2*(i-1)/(n-1)
  curve%x_line(i+n) = -x2*sin(phi) 
  curve%y_line(i+n) =  y1*cos(phi)
enddo
do i=1, n
  phi = pi/2*(i-1)/(n-1)
  curve%x_line(i+2*n) =  -x2*cos(phi) 
  curve%y_line(i+2*n) =  -y2*sin(phi)
enddo
do i=1, n
  phi = pi/2*(i-1)/(n-1)
  curve%x_line(i+3*n) =   x1*sin(phi) 
  curve%y_line(i+3*n) =  -y2*cos(phi)
enddo

case(rectangular$)
  call alloc_curves(5)
  curve%x_line = [x1, -x2, -x2,  x1, x1]
  curve%y_line = [y1,  y1, -y2, -y2, y1]
case default
  return
end select

curve%x_symb = curve%x_line
curve%y_symb = curve%y_line

!-------------------------------------------------
contains
subroutine alloc_curves(n)
integer :: n
call re_allocate (curve%x_line, n)
call re_allocate (curve%y_line, n)
call re_allocate (curve%x_symb, n)
call re_allocate (curve%y_symb, n)
end subroutine alloc_curves

end subroutine tao_curve_physical_aperture_setup


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_histogram_setup (plot, graph)

use bin_mod

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u
type (ele_struct), pointer :: ele
type (beam_struct), pointer :: beam
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1
type (coord_struct), pointer :: p(:)
type (bin_struct) :: bins

real(rp) v_mat(4,4), v_inv_mat(4,4), g_mat(4,4), g_inv_mat(4,4)
real(rp) mat4(4,4), sigma_mat(4,4), theta, theta_xy, rx, ry, phi
real(rp) emit_a, emit_b, bin_shift
real(rp), allocatable :: data(:), weight(:)

integer k, n, m, ib, ix1_ax, ix, i
integer, allocatable :: number_in_bin(:)

logical err, same_uni

character(40) name
character(40) :: r_name = 'tao_graph_histogram_setup'

! Valid?

graph%valid = .false.
if (size(graph%curve) == 0) return

! Set up the graph suffix

same_uni = .true.
ix = tao_universe_number(tao_curve_ix_uni(graph%curve(1)))
do k = 2, size(graph%curve)
  curve => graph%curve(k)
  if (tao_universe_number(tao_curve_ix_uni(curve)) /= ix) same_uni = .false.
enddo

graph%title_suffix = ''
do k = 1, size(graph%curve)
  curve => graph%curve(k)
  u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
  if (.not. associated(u)) return
  if (curve%ix_ele_ref_track < 0) then
    call out_io (s_error$, r_name, &
                'BAD REFERENCE ELEMENT: ' // curve%ele_ref_name, &
                'CANNOT PLOT HISTOGRAM FOR: ' // tao_curve_name(curve))
    return
  endif
  ele => u%model%lat%ele(curve%ix_ele_ref_track)
  name = curve%ele_ref_name
  if (name == '') name = ele%name
  if (same_uni) then
    write (graph%title_suffix, '(2a, i0, 3a)') trim(graph%title_suffix), &
                                '[', curve%ix_ele_ref, ': ', trim(name), ']'
  else
    write (graph%title_suffix, '(2a, i0, a, i0, 3a)') trim(graph%title_suffix), &
            '[', u%ix_uni, '@', curve%ix_ele_ref, ': ', trim(name), ']'
  endif
enddo

! loop over all curves

do k = 1, size(graph%curve)

  curve => graph%curve(k)
  u => tao_pointer_to_universe (tao_curve_ix_uni(curve))

  ! find phase space axes to plot

  err = .false.
  call tao_phase_space_axis (curve%data_type, ix1_ax, err = err); if (err) return

  ! fill the data array

  if (allocated (curve%ix_symb)) deallocate (curve%ix_symb, curve%x_symb, curve%y_symb)
  if (allocated (curve%x_line))  deallocate (curve%x_line, curve%y_line)

  if (curve%data_source == 'beam') then
    beam => u%uni_branch(curve%ix_branch)%ele(curve%ix_ele_ref_track)%beam
    if (.not. allocated(beam%bunch)) then
      call out_io (s_abort$, r_name, 'NO ALLOCATED BEAM WITH PHASE_SPACE PLOTTING.')
      if (.not. u%is_on) call out_io (s_blank$, r_name, '   REASON: UNIVERSE IS TURNED OFF!')
      return
    endif

    if (curve%ix_bunch == 0) then
      n = 0
      do ib = 1,  size(beam%bunch)
        n = n + count(beam%bunch(ib)%particle%state == alive$)
      enddo
    else
      n = count(beam%bunch(curve%ix_bunch)%particle%state == alive$)
    endif

    allocate (data(n))
    if (curve%hist%weight_by_charge) allocate (weight(n))
    
    n = 0
    do ib = 1, size(beam%bunch)
      p => beam%bunch(ib)%particle
      m = size(p)
      call tao_phase_space_axis (curve%data_type, ix1_ax, p, scratch%axis1, ele)
      data(n+1:n+m) = pack(scratch%axis1, mask = (p%state == alive$))
      if (curve%hist%weight_by_charge) weight(n+1:n+m) = pack(p%charge, mask = (p%state == alive$))
      n = n + count(p%state == alive$)
    enddo

  !----------------------------

  elseif (curve%data_source == 'multi_turn_orbit') then
    
    call tao_find_data (err, curve%data_source, d2_array, ix_uni = tao_curve_ix_uni(curve))
    if (err .or. size(d2_array) /= 1) then
      call out_io (s_error$, r_name, 'CANNOT FIND DATA ARRAY TO PLOT CURVE: ' // curve%data_type)
      graph%valid = .false.
      return
    endif

    nullify (d1)
    d2_ptr => d2_array(1)%d2
    do i = 1, size(d2_ptr%d1)
      if (curve%data_type == d2_ptr%d1(i)%name) d1 => d2_ptr%d1(i)
    enddo
    if (.not. associated(d1)) then
      call out_io (s_error$, r_name, &
              'CANNOT FIND DATA FOR PHASE SPACE COORDINATE: ' // curve%data_type, &
              'FOR CURVE: ' // curve%name)
      return
    endif

    allocate(data(size(d1%d)))
    data = d1%d(i)%model_value

  ! Unrecognized

  else
    call out_io (s_abort$, r_name, &
        'INVALID CURVE%DATA_SOURCE: ' // curve%data_source, &
        'FOR CURVE: '// curve%name)
    return
  endif

  ! Bin the data

   ! Automatically scale data
  curve%hist%minimum = minval(data)
  curve%hist%maximum = maxval(data)

  ! Automatically select the number of bins
  if (curve%hist%number == 0) then
    if (curve%hist%width == 0) then
      curve%hist%number = n_bins_automatic(size((data)))
    else
      ! Set the number according to this width, and set a new maximum
      curve%hist%number = nint((curve%hist%maximum - curve%hist%minimum)/curve%hist%width)
      curve%hist%maximum = curve%hist%number*curve%hist%width + curve%hist%minimum
    endif
  else
    ! Number is set, now set the width
    curve%hist%width = (curve%hist%maximum - curve%hist%minimum)/curve%hist%number
  endif
   
  ! Shift for the center bin
  bin_shift = curve%hist%center - bin_x(bin_index(curve%hist%center, curve%hist%minimum, curve%hist%width), &
                              curve%hist%minimum, curve%hist%width)
  curve%hist%minimum =  curve%hist%minimum + bin_shift         
  curve%hist%maximum =  curve%hist%maximum + bin_shift           
                                
  if (curve%hist%width == 0 ) curve%hist%width = (curve%hist%maximum - curve%hist%minimum)/(curve%hist%number - 1)
                            
  if (curve%hist%weight_by_charge) then
    bins = bin(data, weight = weight, min = curve%hist%minimum, max = curve%hist%maximum, n_bins = curve%hist%number)
  else
    bins = bin(data, min = curve%hist%minimum, max = curve%hist%maximum, n_bins = curve%hist%number)
           
  endif
  ! Set width actually used
 
  curve%hist%width = bins%delta
  
  allocate(curve%x_line(bins%n))
  allocate(curve%y_line(bins%n))
  do i=1, bins%n
    curve%x_line(i) = bin_x(i, bins%min, bins%delta)
    curve%y_line(i) = bins%count(i)
  enddo
  if (curve%hist%density_normalized) curve%y_line = curve%y_line/bins%delta
  
enddo

graph%valid = .true.

end subroutine tao_graph_histogram_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_phase_space_axis (data_type, ix_axis, p, axis, ele, err)

implicit none

type (coord_struct), optional, target :: p(:)
type (coord_struct) :: p1
type (ele_struct), optional :: ele

real(rp), allocatable, optional :: axis(:)

integer ix_axis, i

logical, optional :: err

character(*) data_type
character(16) :: r_name = 'phase_space_axis'

!

if (present(p)) call re_allocate (axis, size(p))

select case (data_type)
case ('x');   ix_axis = 1; if (present(p)) axis = p%vec(1)
case ('px');  ix_axis = 2; if (present(p)) axis = p%vec(2)
case ('y');   ix_axis = 3; if (present(p)) axis = p%vec(3)
case ('py');  ix_axis = 4; if (present(p)) axis = p%vec(4)
case ('z');   ix_axis = 5; if (present(p)) axis = p%vec(5)
case ('pz');  ix_axis = 6; if (present(p)) axis = p%vec(6)
case ('intensity_x'); ix_axis =  7; if (present(p)) axis = p%field(1)**2
case ('intensity_y'); ix_axis =  8; if (present(p)) axis = p%field(2)**2
case ('phase_x');     ix_axis =  9; if (present(p)) axis = p%phase(1)
case ('phase_y');     ix_axis = 10; if (present(p)) axis = p%phase(2)


case ('intensity')
  ix_axis = 11
  if (present(p)) then
    p%charge = p%field(1)**2 + p%field(2)**2
    axis = p%charge
  endif
  
case ('Ja')
  ix_axis = 12
  if (present(p) .and. present(ele)) then
    do i=1, size(p)
      call convert_coords('LAB', p(i), ele, 'ACTION-ANGLE', p1)
      axis(i) = p1%vec(1)
    enddo
  endif

case ('energy')
  ix_axis = 13
  if (present(p)) then
    do i=1, size(p)
      call convert_pc_to((1.0_rp+ p(i)%vec(6))* p(i)%p0c,  p(i)%species, e_tot =  axis(i))
    enddo
  endif
  
case ('t');     ix_axis = 14; if (present(p)) axis = p%t

case default
  call out_io (s_abort$, r_name, 'BAD PHASE_SPACE CURVE DATA_TYPE: ' // data_type)
  if (present(err)) err = .true.
end select

end subroutine tao_phase_space_axis

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_graph_data_setup (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), target :: branch_curve
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u

integer n, ic, ib, n0_line, n0_symb
logical err

!

graph%title_suffix = ''
if (graph%component /= '') graph%title_suffix = '[' // trim(graph%component) // ']'

! Attach x-axis type to title suffix if needed.
! Needed %label is blank and %draw_label = F.
! Note: if %label is blank and %draw_label = T then the x-axis_type is printed elsewhere.
 
if (graph%x%label == '' .and. .not. graph%x%draw_label) then
  if (plot%x_axis_type == "lat" .or. plot%x_axis_type == "var") then
    graph%title_suffix = trim(graph%title_suffix) // ',  X-axis: ' // &
              trim(plot%x_axis_type) // '::' // graph%curve(1)%data_type_x
  else
    graph%title_suffix = trim(graph%title_suffix) // ',  X-axis: ' // plot%x_axis_type
  endif
endif

! Loop over all curves in the graph

graph%valid = .false.

if (allocated(graph%curve)) then
  do ic = 1, size(graph%curve)
    curve => graph%curve(ic)

    ! Floor plan curves use all branches.
    if (graph%type == 'floor_plan') then
      u => tao_pointer_to_universe (tao_curve_ix_uni(graph%curve(ic)))
      if (.not. associated(u)) return
      n0_line = 0
      n0_symb = 0
      !! call deallocate_curve_arrays(curve)
      branch_curve = graph%curve(ic)

      do ib = 0, size(u%model%lat%branch)
        curve%ix_branch = ib
        call tao_curve_data_setup (plot, graph, branch_curve, err)
        if (err) return
        n = n0_line + size(branch_curve%x_line)
        !! call re_allocate (curve%x_line, n);  curve%x_line(n0_line+1:) = branch_curve%x_line
        !! call re_allocate (curve%y_line, n);  curve%y_line(n0_line+1:) = branch_curve%y_line        
      enddo

    else
      call tao_curve_data_setup (plot, graph, graph%curve(ic), err)
      if (err) return
    endif

  enddo
endif

graph%valid = .true.

end subroutine tao_graph_data_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_curve_data_setup (plot, graph, curve, err_flag)

use nrutil, only: swap

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), target :: curve
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: model_lat, base_lat
type (tao_ele_shape_struct), pointer :: ele_shape
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (ele_struct), pointer :: ele, ele1, ele2, slave
type (branch_struct), pointer :: branch

real(rp) f, eps, gs, l_tot, s0, s1, x_max, x_min, val, val0, dx
real(rp), allocatable :: value_arr(:)
real(rp), pointer :: var_ptr

integer ii, k, m, n, n_dat, ie, jj, iv, ic
integer ix, ir, jg, i, j, ix_this, ix_uni, ix1, ix2, n_curve_pts

logical err, err_flag, smooth_curve, found, zero_average_phase, ok
logical straight_line_between_syms, valid, in_graph
logical, allocatable :: good(:)

character(200) data_type, name
character(60) component
character(16) data_source, dflt_index
character(20), parameter :: r_name = 'tao_curve_data_setup'

!

n_curve_pts = s%plot_page%n_curve_pts
if (plot%n_curve_pts > 0) n_curve_pts = plot%n_curve_pts

call re_allocate_eles (scratch%eles, 1, exact = .true.)
err_flag = .true.

if (allocated(curve%x_line)) deallocate (curve%x_line, curve%y_line)
if (allocated(curve%y2_line)) deallocate (curve%y2_line)
if (allocated(curve%ix_line)) deallocate (curve%ix_line)
if (allocated(curve%x_symb)) deallocate (curve%x_symb, curve%y_symb)
if (allocated(curve%z_symb)) deallocate (curve%z_symb)
if (allocated(curve%ix_symb)) deallocate (curve%ix_symb)
if (allocated(curve%symb_size)) deallocate (curve%symb_size)

u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
if (.not. associated(u)) then
  graph%why_invalid = 'NO ASSOCIATED UNIVERSE!'
  return
endif

if (s%com%common_lattice) then
  u%calc%lattice = .true.
  call tao_lattice_calc (ok)
endif

model_lat => u%model%lat
base_lat => u%base%lat
branch => model_lat%branch(curve%ix_branch)

if (curve%ele_ref_name == '') then
  zero_average_phase = .true.
else
  zero_average_phase = .false.
  call tao_locate_elements (curve%ele_ref_name, tao_curve_ix_uni(curve), scratch%eles, err, ignore_blank = .true.)
  if (err) then
    graph%why_invalid = 'CANNOT LOCATE ELEMENT: ' // trim(curve%ele_ref_name)
    return
  endif
  curve%ix_branch  = scratch%eles(1)%ele%ix_branch
  curve%ix_ele_ref = scratch%eles(1)%ele%ix_ele
  call tao_ele_to_ele_track(tao_curve_ix_uni(curve), scratch%eles(1)%ele%ix_branch, &
                                      scratch%eles(1)%ele%ix_ele, curve%ix_ele_ref_track)
endif

component = tao_curve_component(curve, graph)

!----------------------------------------------------------------------------
! Calculate where the symbols are to be drawn on the graph.

data_source = curve%data_source

if (plot%x_axis_type == 'lat' .or. plot%x_axis_type == 'var') data_source = 'plot_x_axis_var'
if (curve%data_type(1:11) == 'expression:' .and. (data_source == 'dat' .or. data_source == 'var')) data_source = 'expression'

select case (data_source)

!----------------------------------------------------------------------------
! Case: expression

case ('expression')

  if (plot%x_axis_type == 'index') then
    n_dat = nint(graph%x%max) - nint(graph%x%min) + 1

    call re_allocate (curve%ix_symb, n_dat)
    call re_allocate (curve%x_symb,  n_dat) ! allocate space for the data
    call re_allocate (curve%y_symb,  n_dat) ! allocate space for the data

    n_dat = 0
    do i = nint(graph%x%min), nint(graph%x%max)
      write (dflt_index, '(i0)') i
      call tao_evaluate_expression  (curve%data_type(12:), 0, graph%draw_only_good_user_data_or_vars, value_arr, &
                   scratch%info, err, .false., scratch%stack, component, curve%data_source, dflt_dat_or_var_index = dflt_index)
      if (err .or. .not. scratch%info(1)%good) cycle
      n_dat = n_dat + 1

      curve%x_symb(n_dat) = i
      curve%y_symb(n_dat) = value_arr(1)
      curve%ix_symb(n_dat) = 0
    enddo

  else
    call tao_evaluate_expression  (curve%data_type(12:), 0, graph%draw_only_good_user_data_or_vars, value_arr, &
                                                 scratch%info, err, .false., scratch%stack, component, curve%data_source)
    n_dat = count(scratch%info%good)
    call re_allocate (curve%ix_symb, n_dat)
    call re_allocate (curve%x_symb,  n_dat) ! allocate space for the data
    call re_allocate (curve%y_symb,  n_dat) ! allocate space for the data

    n_dat = 0
    do i = 1, size(value_arr)
      if (.not. scratch%info(i)%good) cycle
      n_dat = n_dat + 1
      curve%y_symb(n_dat) = value_arr(i)

      dx = (graph%x%max - graph%x%min) / 100.0_rp

      select case (plot%x_axis_type)
      case ('s')
        if (scratch%info(i)%s < graph%x%min - dx .or. scratch%info(i)%s > graph%x%max + dx) cycle
        curve%x_symb(n_dat) = scratch%info(i)%s
      case ('ele_index')
        if (scratch%info(i)%ix_ele < graph%x%min - dx .or. scratch%info(i)%ix_ele > graph%x%max + dx) cycle
        curve%x_symb(n_dat) = scratch%info(i)%ix_ele
      end select
    enddo
  endif

  call re_allocate (curve%ix_symb, n_dat)
  call re_allocate (curve%x_symb,  n_dat) ! allocate space for the data
  call re_allocate (curve%y_symb,  n_dat) ! allocate space for the data

  !

  if (curve%draw_line) then
    call re_allocate (curve%y_line,  n_dat) ! allocate space for the data
    call re_allocate (curve%x_line,  n_dat) ! allocate space for the data
    curve%x_line = curve%x_symb
    curve%y_line = curve%y_symb
  else
    if (allocated (curve%x_line)) deallocate (curve%x_line, curve%y_line)
  endif

!----------------------------------------------------------------------------
! Case: x-axis uses a variable.

case ('plot_x_axis_var')

  call re_allocate (curve%ix_symb, n_curve_pts)
  call re_allocate (curve%x_symb, n_curve_pts)
  call re_allocate (curve%y_symb, n_curve_pts)
  call re_allocate (curve%x_line, n_curve_pts)
  call re_allocate (curve%y_line, n_curve_pts)

  if (plot%x_axis_type == 'lat') then

    call tao_pick_universe (curve%data_type_x, name, scratch%this_u, err, ix_uni)
    if (err .or. count(scratch%this_u) /= 1) then
      graph%why_invalid = 'BAD UNIVERSE CONSTRUCT IN CURVE%DATA_TYPE_X: ' // trim(curve%data_type_x)
      return
    endif

    call upcase_string(name)
    ix1 = index(name, '[')
    ix2 = index(name, ']')
    if (ix1 == 0 .or. ix2 == 0 .or. ix2 /= len_trim(name)) then
      graph%why_invalid = 'BAD VARIABLE CONSTRUCT IN CURVE%DATA_TYPE_X: ' // trim(curve%data_type_x)
      return
    endif

    u => tao_pointer_to_universe(ix_uni)
    call pointers_to_attribute (u%model%lat, name(1:ix1-1), name(ix1+1:ix2-1), .true., scratch%attribs, &
                                                                                 err, eles = scratch%eles)
    if (err .or. size(scratch%attribs) /= 1) then
      graph%why_invalid = 'BAD VARIABLE CONSTRUCT IN CURVE%DATA_TYPE_X: ' // trim(curve%data_type_x)
      return
    endif
    var_ptr => scratch%attribs(1)%r

  else  ! x_axis_type == 'var'
    call tao_find_var (err, curve%data_type_x, v_array = scratch%var_array)
    if (err .or. size(scratch%var_array) /= 1) then
      graph%why_invalid = 'BAD VARIABLE CONSTRUCT IN CURVE%DATA_TYPE_X: ' // trim(curve%data_type_x)
      return
    endif
    var_ptr => scratch%var_array(1)%v%model_value
  endif

  ! Get datum values as a function of the variable

  val0 = var_ptr

  do i = 1, n_curve_pts 
    val = graph%x%min + (graph%x%max - graph%x%min) * (i - 1.0_rp) / (n_curve_pts - 1)
    if (plot%x_axis_type == 'lat')then
      var_ptr = val
      call tao_set_flags_for_changed_attribute (u, name(1:ix1-1), scratch%eles(1)%ele, var_ptr)
      s%u(ix_uni)%calc%lattice = .true.
    else
      call tao_set_var_model_value (scratch%var_array(1)%v, val)
    endif
    call tao_lattice_calc (valid)

    call tao_evaluate_expression (curve%data_type, 0, .false., value_arr, scratch%info, err, &
                          dflt_component = component, dflt_source = curve%data_source)
    if (.not. valid .or. err .or. size(value_arr) /= 1) then
      graph%why_invalid = 'BAD CONSTRUCT IN CURVE%DATA_TYPE: ' // trim(curve%data_type)
      return
    endif

    curve%x_symb(i) = val      
    curve%y_symb(i) = value_arr(1)

    curve%x_line(i) = val      
    curve%y_line(i) = value_arr(1)

  enddo

  ! Reset

  if (plot%x_axis_type == 'lat')then
    var_ptr = val0
    s%u(ix_uni)%calc%lattice = .true.
  else
    call tao_set_var_model_value (scratch%var_array(1)%v, val0)
  endif
  call tao_lattice_calc (valid)

!----------------------------------------------------------------------------
! Case: data_source is a data_array

case ('data')

  ! Calculate values

  call tao_data_type_substitute (curve%data_type, data_type, curve, graph)
  call tao_evaluate_expression  (data_type, 0, graph%draw_only_good_user_data_or_vars, value_arr, scratch%info, err, &
                          stack = scratch%stack, dflt_component = component, dflt_source = 'data')
  if (err) then
    graph%why_invalid = 'CANNOT FIND DATA CORRESPONDING: ' // data_type
    return
  end if

  ! point d1_array to the data to be plotted

  do i = 1, size(scratch%stack)
    if (scratch%stack(i)%type == data_num$) exit
    if (i == size(scratch%stack)) then
      graph%why_invalid = 'CANNOT FIND DATA ARRAY TO PLOT CURVE: ' // curve%data_type
      return
    endif
  enddo

  call tao_find_data (err, scratch%stack(i)%name, d2_array, scratch%d1_array, ix_uni = tao_curve_ix_uni(curve))
  if (err .or. size(scratch%d1_array) /= 1) then
    graph%why_invalid = 'CANNOT FIND VALID DATA ARRAY TO PLOT CURVE: ' // curve%data_type
    return
  endif

  d1_ptr => scratch%d1_array(1)%d1
  if (index(curve%data_type, 'phase') /= 0) then
    if (all(d1_ptr%d(:)%ele_ref_name == '')) then
      zero_average_phase = .true.
    else
      zero_average_phase = .false.
    endif
  endif

  ! Set %good_plot True for all data that is within the x-axis limits.
  ! For a circular lattice "wrap around" at s = 0 may mean 
  !   some data points show up twice.

  d1_ptr => scratch%d1_array(1)%d1
  d1_ptr%d%good_plot = .false.
  if (graph%x%min /= graph%x%max) then
    eps = 1e-4 * (graph%x%max - graph%x%min)
    if (plot%x_axis_type == 'index') then
      where (d1_ptr%d%ix_d1 > graph%x%min-eps .and. &
             d1_ptr%d%ix_d1 < graph%x%max+eps) d1_ptr%d%good_plot = .true.
    elseif (plot%x_axis_type == 'ele_index') then
      where (d1_ptr%d%ix_ele > graph%x%min-eps .and. &
             d1_ptr%d%ix_ele < graph%x%max+eps) d1_ptr%d%good_plot = .true.
    else ! s
      where (d1_ptr%d%s > graph%x%min-eps .and. &
             d1_ptr%d%s < graph%x%max+eps) d1_ptr%d%good_plot = .true.
      if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) then 
        l_tot = branch%param%total_length
        where (d1_ptr%d%s-l_tot > graph%x%min-eps .and. &
               d1_ptr%d%s-l_tot < graph%x%max+eps) d1_ptr%d%good_plot = .true.
        where (d1_ptr%d%s+l_tot > graph%x%min-eps .and. &
               d1_ptr%d%s+l_tot < graph%x%max+eps) d1_ptr%d%good_plot = .true.
      endif

    endif
  endif

  ! make sure %useit_plot up-to-date & count the number of data points

  call tao_data_useit_plot_calc (curve, graph, d1_ptr%d) 
  if (plot%x_axis_type == 's') then
    ! veto non-regular elements when plotting s
    do m = lbound(d1_ptr%d,1), ubound(d1_ptr%d,1)
      if (d1_ptr%d(m)%ix_ele > model_lat%branch(d1_ptr%d(m)%ix_branch)%n_ele_track) then
        d1_ptr%d(m)%useit_plot = .false.
      endif
    enddo
  endif
  n_dat = count (d1_ptr%d%useit_plot)       

  ! resize the curve data arrays

  call re_allocate (curve%ix_symb, n_dat)
  call re_allocate (curve%y_symb, n_dat) ! allocate space for the data
  call re_allocate (curve%x_symb, n_dat) ! allocate space for the data

  ! 

  curve%ix_symb = pack(d1_ptr%d%ix_d1, mask = d1_ptr%d%useit_plot)
  curve%y_symb  = pack(value_arr, mask = d1_ptr%d%useit_plot)

  if (plot%x_axis_type == 'index') then
    curve%x_symb = curve%ix_symb
  elseif (plot%x_axis_type == 'ele_index') then
    curve%x_symb = d1_ptr%d(curve%ix_symb)%ix_ele
  elseif (plot%x_axis_type == 's') then
    curve%x_symb = branch%ele(d1_ptr%d(curve%ix_symb)%ix_ele)%s
    ! If there is a wrap-around then reorder data
    if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) then
      do i = 1, n_dat
        if (curve%x_symb(i) > graph%x%max+eps) curve%x_symb(i) = curve%x_symb(i)-l_tot
        if (curve%x_symb(i) < graph%x%min-eps) curve%x_symb(i) = curve%x_symb(i)+l_tot
      enddo
    endif
    ! Super lords will be out of order so reorder in increasing s.
    do i = 2, n_dat
      do j = i, 2, -1
        if (curve%x_symb(j-1) > curve%x_symb(j)) then
          call swap(curve%x_symb(j-1), curve%x_symb(j))
          call swap(curve%y_symb(j-1), curve%y_symb(j))
          call swap(curve%ix_symb(j-1), curve%ix_symb(j))
        else
          exit
        endif
      enddo
    enddo

  else
    graph%why_invalid = 'UNKNOWN AXIS TYPE!'
    return
  endif


!----------------------------------------------------------------------------
! Case: data_source is a var_array

case ('var')

  call tao_find_var (err, curve%data_type, scratch%v1_array)
  if (err .or. size(scratch%v1_array) /= 1) return
  v1_ptr => scratch%v1_array(1)%v1

  ! find which universe we're viewing
  ix_this = -1
  v_loop: do iv = lbound(v1_ptr%v, 1), ubound(v1_ptr%v,1)
    v_ptr => v1_ptr%v(iv)
    if (.not. v_ptr%exists) cycle
    do jj = 1, size(v_ptr%slave)
      if (v_ptr%slave(jj)%ix_uni .eq. s%com%default_universe) then
        ix_this = jj
        exit v_loop
      endif
    enddo
  enddo v_loop
  if (ix_this .eq. -1) then
    call out_io (s_error$, r_name, &
                   "This variable doesn't point to the currently displayed  universe.")
    return
  endif

  v1_ptr%v%good_plot = .true.
  if (graph%x%min /= graph%x%max) then
    eps = 1e-4 * (graph%x%max - graph%x%min)
    if (plot%x_axis_type == 'index') then
      where (v1_ptr%v%ix_v1 < graph%x%min-eps) v1_ptr%v%good_plot = .false.
      where (v1_ptr%v%ix_v1 > graph%x%max+eps) v1_ptr%v%good_plot = .false.
    elseif (plot%x_axis_type == 'ele_index') then
      do jj = lbound(v1_ptr%v, 1), ubound(v1_ptr%v,1)
        if (v1_ptr%v(jj)%slave(ix_this)%ix_ele < graph%x%min-eps) v1_ptr%v%good_plot = .false.
        if (v1_ptr%v(jj)%slave(ix_this)%ix_ele > graph%x%max+eps) v1_ptr%v%good_plot = .false.
      enddo
    else
      where (v1_ptr%v%s < graph%x%min-eps) v1_ptr%v%good_plot = .false.
      where (v1_ptr%v%s > graph%x%max+eps) v1_ptr%v%good_plot = .false.
    endif
  endif

  call tao_var_useit_plot_calc (graph, v1_ptr%v) ! make sure %useit_plot up-to-date
  n_dat = count(v1_ptr%v%useit_plot)             ! count the number of data points

  call re_allocate (curve%ix_symb, n_dat)
  call re_allocate (curve%y_symb,  n_dat) ! allocate space for the data
  call re_allocate (curve%x_symb,  n_dat) ! allocate space for the data

  curve%ix_symb = pack(v1_ptr%v%ix_v1, mask = v1_ptr%v%useit_plot)

  graph%x%label = plot%x_axis_type

  if (plot%x_axis_type == 'index') then
    curve%x_symb = curve%ix_symb
  elseif (plot%x_axis_type == 'ele_index') then
    do jj = lbound(curve%ix_symb,1), ubound(curve%ix_symb,1)
      curve%x_symb(jj) = v1_ptr%v(curve%ix_symb(jj))%slave(ix_this)%ix_ele
    enddo
  elseif (plot%x_axis_type == 's') then
    do jj = lbound(curve%ix_symb,1), ubound(curve%ix_symb,1)
      ele => branch%ele(v1_ptr%v(curve%ix_symb(jj))%slave(ix_this)%ix_ele)
      if (ele%lord_status == multipass_lord$) ele => pointer_to_slave(ele, 1)
      curve%x_symb(jj) = ele%s
    enddo
  endif

  ! calculate the y-axis data point values.

  data_type = trim(curve%data_type) // '|' // trim(component)
  call tao_evaluate_expression (data_type, 0, graph%draw_only_good_user_data_or_vars, value_arr, scratch%info, err)
  if (err) then
    call out_io (s_error$, r_name, 'BAD PLOT COMPONENT: ' // data_type)
    return
  end if

  curve%y_symb = pack(value_arr, mask = v1_ptr%v%useit_plot)

!----------------------------------------------------------------------------
! Case: data_source is from lattice, or beam

case ('lat', 'beam')

  ! Find how many symbol points there are...
  ! Here 'index' and 'ele_index' mean the same thing.

  select case (plot%x_axis_type)
  case ('index', 'ele_index')
    x_min = 1
    x_max = branch%n_ele_track
    if (graph%x%min /= graph%x%max) then
      x_min = max(x_min, graph%x%min)
      x_max = min(x_max, graph%x%max)
    endif 
    n_dat = max(0, nint(x_max+1-x_min))
    call re_allocate_eles (scratch%eles, n_dat, exact = .true.)
    do i = 1, n_dat
      scratch%eles(i)%ele => pointer_to_ele (model_lat, nint(i+x_min-1), curve%ix_branch)
    enddo

  ! x_axis_type == 's':

  case ('s')
    ! Symbols are to be put at the ends of displayed elements in the lat_layout
    eps = 1e-4 * (graph%x%max - graph%x%min)       ! a small number
    branch%ele%logic = .false.                     ! Mark if ele is in the graph

    ! Mark all eles in branch if they match a shape.
    do i = 0, branch%n_ele_track
      ele => branch%ele(i)
      ele_shape => tao_pointer_to_ele_shape (u%ix_uni, ele, s%plot_page%lat_layout%ele_shape)
      if (.not. associated(ele_shape)) cycle
      if (.not. ele_shape%draw) cycle
      call find_element_ends (ele, ele1, ele2)
      ele1%logic = .true.
      ele2%logic = .true.
    enddo

    ! Mark slaves of lord elements that match a shape.
    do i = model_lat%n_ele_track+1, model_lat%n_ele_max
      ele => model_lat%ele(i)
      ele_shape => tao_pointer_to_ele_shape (u%ix_uni, ele, s%plot_page%lat_layout%ele_shape)
      if (.not. associated(ele_shape)) cycle
      if (.not. ele_shape%draw) cycle
      if (ele%lord_status == multipass_lord$) then
        do j = 1, ele%n_slave
          slave => pointer_to_slave (ele, j)
          call find_element_ends (slave, ele1, ele2)
          ele1%logic = .true.
          ele2%logic = .true.
        enddo
      else
        call find_element_ends (ele, ele1, ele2)
        ele1%logic = .true.
        ele2%logic = .true.
      endif
    enddo

    ! Now unmark all elements in the branch that are not within the graph boundries.
    do i = 0, branch%n_ele_track
      ele => branch%ele(i)
      if (.not. ele%logic) cycle
      if (graph%x%min == graph%x%max) cycle
      s0 = ele%s_start
      s1 = ele%s
      in_graph = (s0 >= graph%x%min-eps) .and. (s1 <= graph%x%max+eps)
      l_tot = branch%param%total_length
      if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) in_graph = in_graph .or. &
                      (s0-l_tot >= graph%x%min-eps) .and. (s1-l_tot <= graph%x%max+eps)
      ele%logic = ele%logic .and. in_graph                                 
    enddo

    ! Allocate eles(:) array and set the eles(:)%ele pointers to point to the marked elements.
    n_dat = count (branch%ele(:)%logic)
    call re_allocate_eles (scratch%eles, n_dat, exact = .true.)

    n = 0
    do i = 0, branch%n_ele_max
      ele => branch%ele(i)
      if (.not. ele%logic) cycle
      n = n + 1
      scratch%eles(n)%ele => ele 
    enddo      

  ! Error for x_axis_type unrecognized.

  case default
    call out_io (s_error$, r_name, 'BAD PLOT%X_AXIS_TYPE: ' // plot%x_axis_type)
    return
  end select

if (curve%draw_symbols) then
  call tao_curve_datum_calc (scratch%eles, plot, curve, 'SYMBOL', valid)
  if (.not. valid) return
endif

!----------------------------------------------------------------------------
! Case: Bad data_source

case default
  call out_io (s_error$, r_name, 'UNKNOWN DATA_SOURCE: ' // curve%data_source)
  return
end select

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Now calculate the points for drawing the curve through the symbols...

! If the x-axis is by index or ele_index then these points are the same as the symbol points.
! That is, for x-axis = 'index' or 'ele_index' the line is piece-wise linear between the symbols.

select case (plot%x_axis_type)
case ('index', 'ele_index')
  call re_allocate (curve%y_line, size(curve%x_symb)) ! allocate space for the data
  call re_allocate (curve%x_line, size(curve%y_symb)) ! allocate space for the data
  curve%x_line = curve%x_symb
  curve%y_line = curve%y_symb

! If the axis is by s-value then, if possible, the line is a "smooth" curve with n_curve_pts points.

case ('s')

  ! beam data_source is not interpolated.
  smooth_curve = (curve%data_source == 'lat' .and. curve%smooth_line_calc .and. .not. s%global%disable_smooth_line_calc)
  if (index(curve%data_type, 'emit.') /= 0) smooth_curve = .false.

  if (index(component, 'meas') /= 0 .or. index(component, 'ref') /= 0 .or. &
      curve%data_source == 'data' .or. curve%data_source == 'var') then
    straight_line_between_syms = .true.
    smooth_curve = .false.
  else
    straight_line_between_syms = .false.
  endif

  if (smooth_curve) then

    ! allocate data space

    call re_allocate (curve%y_line, n_curve_pts) 
    call re_allocate (curve%x_line, n_curve_pts) 
    call re_allocate (good, n_curve_pts) 
    curve%y_line = 0
    good = .true.

    call tao_split_component(component, scratch%comp, err)
    if (err) then
      graph%valid = .false.
      graph%why_invalid = 'Bad Graph or Curve Component Expression'
      return
    endif
    if (component == '') then
      graph%valid = .false.
      graph%why_invalid = 'No Graph or Curve Component.'
      return
    endif

    do m = 1, size(scratch%comp)
      select case (scratch%comp(m)%name)
      case ('') 
        cycle
      case ('model')
        call tao_calc_data_at_s (u%model, curve, scratch%comp(m)%sign, good)
      case ('base')  
        call tao_calc_data_at_s (u%base, curve, scratch%comp(m)%sign, good)
      case ('design')  
        call tao_calc_data_at_s (u%design, curve, scratch%comp(m)%sign, good)
      case default
        call out_io (s_error$, r_name, &
                     'BAD PLOT COMPONENT WITH "S" X-AXIS: ' // scratch%comp(m)%name)
        return
      end select
    enddo

    !! if (all(.not. good)) exit
    n_dat = count(good)
    curve%x_line(1:n_dat) = pack(curve%x_line, mask = good)
    curve%y_line(1:n_dat) = pack(curve%y_line, mask = good)
    call re_allocate (curve%y_line, n_dat) ! allocate space for the data
    call re_allocate (curve%x_line, n_dat) ! allocate space for the data

  ! For non-smooth curves: Draw straight lines through the symbols if
  ! the data uses "ref" or "meas" values. Else evaluate at the element ends.

  else if (straight_line_between_syms) then
    if (allocated (curve%x_symb)) then
      n_dat = size(curve%x_symb)
      call re_allocate (curve%y_line, n_dat) 
      call re_allocate (curve%x_line, n_dat) 
      curve%x_line = curve%x_symb 
      curve%y_line = curve%y_symb 
    else
      call tao_curve_datum_calc (scratch%eles, plot, curve, 'LINE', valid)
    endif

  ! Evaluate at element ends

  else

    eps = 1e-4 * (graph%x%max - graph%x%min)             ! a small number
    l_tot = branch%param%total_length
    branch%ele%logic = .false.
    do i = 0, branch%n_ele_track
      ele => branch%ele(i)
      if (graph%x%min == graph%x%max) cycle
      s0 = ele%s_start
      s1 = ele%s
      ele%logic = (s0 >= graph%x%min-eps) .and. (s1 <= graph%x%max+eps)
      if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) then
        ele%logic = ele%logic .or. ((s0-l_tot >= graph%x%min-eps) .and. (s1-l_tot <= graph%x%max+eps))
      endif
    enddo
    n_dat = count (branch%ele(:)%logic)
    call re_allocate_eles (scratch%eles, n_dat, exact = .true.)
    i = 0
    do j = 0, ubound(branch%ele, 1)
      if (.not. branch%ele(j)%logic) cycle
      i = i + 1
      scratch%eles(i)%ele => branch%ele(j)
    enddo

    ! If there is a wrap-around then reorder the data

    ix1 = 0
    if (branch%param%geometry == closed$ .and. graph%allow_wrap_around) then
      do i = 1, n_dat
        if (ix1 == 0 .and. branch%ele(scratch%eles(i)%ele%ix_ele)%s - l_tot > graph%x%min) ix1 = i
        if (branch%ele(scratch%eles(i)%ele%ix_ele)%s < graph%x%max+eps) ix2 = i
      enddo
      if (ix1 /= 0) then
        call re_allocate_eles(scratch%eles, n_dat + ix2 + 1 - ix1, .true., .true.)
        scratch%eles = [scratch%eles(ix1:n_dat), scratch%eles(1:ix2)]
      endif
    endif

    if (curve%draw_line) then
      call tao_curve_datum_calc (scratch%eles, plot, curve, 'LINE', graph%valid)
      if (.not. graph%valid) return
    endif

    do i = 1, size(curve%x_line)
      curve%x_line(i) = branch%ele(scratch%eles(i)%ele%ix_ele)%s
      if (ix1 /= 0 .and. i <= n_dat - ix1 + 1) curve%x_line(i) = curve%x_line(i) - l_tot
    enddo

  endif

end select

!----------------------------------------------------------------------------
! Note: Since there is an arbitrary overall phase, the phase data 
! gets renormalized so that the average value is zero.

if (index(curve%data_type, 'phase') /= 0 .and. n_dat /= 0 .and. zero_average_phase) then
  if (allocated(curve%y_symb)) then
    f = sum(curve%y_symb) / n_dat
    curve%y_symb = curve%y_symb - f
    if (allocated(curve%y_line)) curve%y_line = curve%y_line - f 
  elseif (allocated(curve%y_line)) then
    f = sum(curve%y_line) / n_dat
    curve%y_line = curve%y_line - f
  endif
endif 

err_flag = .false.

end subroutine tao_curve_data_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_calc_data_at_s (tao_lat, curve, comp_sign, good)

use transfer_map_mod

implicit none

type (tao_lattice_struct), target :: tao_lat
type (tao_curve_struct) curve
type (tao_lattice_branch_struct), pointer :: tao_branch
type (bunch_params_struct), pointer :: bunch_params0, bunch_params1
type (bunch_params_struct) :: bunch_params
type (coord_struct), pointer :: orb(:), orb_ref
type (coord_struct) orbit_end, orbit_last, orbit
type (lat_struct), pointer :: lat
type (ele_struct) ele, ele_dum
type (ele_struct), pointer :: ele_ref
type (taylor_struct) t_map(6)
type (branch_struct), pointer :: branch
type (all_pointer_struct) a_ptr
type (em_field_struct) field
type (spin_polar_struct) polar

real(rp) x1, x2, cbar(2,2), s_last, s_now, value, mat6(6,6), vec0(6)
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), one_pz, gamma, len_tot
real(rp) comp_sign, vec3(3), r_bunch, amp, phase, ds

integer i, ii, ix, j, k, expnt(6), ix_ele, ix_ref, ix_branch, idum, n_ele_track
integer cache_status   
integer, parameter :: cache_off$ = 0, loading_cache$ = 1, using_cache$ = 2

character(40) data_type, name
character(40) data_type_select, data_source
character(*), parameter :: r_name = 'tao_calc_data_at_s'
logical err, good(:), first_time, radiation_fluctuations_on

! Some init

data_type = curve%data_type

ix_branch = curve%ix_branch
lat => tao_lat%lat
tao_branch => tao_lat%tao_branch(ix_branch)
orb => tao_branch%orbit
branch => lat%branch(ix_branch)
first_time = .true.
n_ele_track = branch%n_ele_track

ix_ref = curve%ix_ele_ref_track
if (ix_ref < 0) ix_ref = 0

if (lat%param%geometry == closed$ .and. .not. lat%param%stable) then
  curve%g%why_invalid = 'Unstable Lattice'
  good = .false.
  return
endif

if (curve%data_source == 'lat') then
  select case (data_type(1:5))
  case ('sigma', 'emitt', 'norm_')
    curve%g%why_invalid = 'curve%data_source = "lat" is not compatable with data_type: ' // data_type
    call out_io (s_warn$, r_name, curve%g%why_invalid)
    call out_io (s_blank$, r_name, "Will not perform any plot smoothing")
    good = .false.
    return
  end select 
endif

! x1 and x2 are the longitudinal end points of the plot

radiation_fluctuations_on = bmad_com%radiation_fluctuations_on
bmad_com%radiation_fluctuations_on = .false.

x1 = branch%ele(0)%s
x2 = branch%ele(n_ele_track)%s
len_tot = x2 - x1
if (curve%g%x%min /= curve%g%x%max) then
  if (branch%param%geometry == closed$) then
    x1 = min(branch%ele(n_ele_track)%s, max(curve%g%x%min, x1-len_tot))
    x2 = min(x2, max(curve%g%x%max, branch%ele(0)%s-len_tot))
  else
    x1 = min(branch%ele(n_ele_track)%s, max(curve%g%x%min, x1))
    x2 = min(x2, max(curve%g%x%max, branch%ele(0)%s))
  endif
endif
ele_ref => branch%ele(ix_ref)
orb_ref => orb(ix_ref)
s_last = ele_ref%s
orbit_last = orbit

ix = index(data_type, '.')
if (ix == 0) then
  data_type_select = data_type
else
  data_type_select = data_type(:ix)
endif

! Only cache plot data if the number of points is equal to s%plot_page%n_curve_pts

if (curve%data_source == 'lat') then
  if (tao_branch%plot_cache_valid .and. tao_branch%cache_x_min == x1 .and. &
            tao_branch%cache_x_max == x2 .and. tao_branch%cache_n_pts == size(curve%x_line)) then
    cache_status = using_cache$

  else
    cache_status = loading_cache$
    if (allocated(tao_branch%plot_cache)) then
      if (size(tao_branch%plot_cache) /= s%plot_page%n_curve_pts) deallocate(tao_branch%plot_cache)
    endif
    if (.not. allocated(tao_branch%plot_cache)) allocate (tao_branch%plot_cache(s%plot_page%n_curve_pts))
    tao_branch%cache_x_min = x1
    tao_branch%cache_x_max = x2 
    tao_branch%cache_n_pts = size(curve%x_line)
    tao_branch%plot_cache_valid = .true.
  endif
endif

!

do ii = 1, size(curve%x_line)

  ! Good(ii) may be false if this is not first time tao_calc_data_at_s is called from tao_curve_data_setup.
  ! For example, tao_calc_data_at_s is called twice when plotting "meas - design".

  if (.not. good(ii)) then
    first_time = .true.
    cycle
  endif

  s_now = x1 + (ii-1) * (x2-x1) / (size(curve%x_line)-1)
  if (s_now > branch%ele(n_ele_track)%s) s_now = branch%ele(n_ele_track)%s
  curve%x_line(ii) = s_now
  value = 0

  ! Check if in a hybrid or taylor element within which interpolation cannot be done.

  ix_ele = element_at_s (lat, s_now, .true., ix_branch, err)
  if (branch%ele(ix_ele)%key == hybrid$ .or. branch%ele(ix_ele)%key == taylor$ .or. err) then
    good(ii) = .false.
    first_time = .true.
    cycle
  endif

  !-----------------------------

  select case (curve%data_source)
  case ('beam')
    if (.not. allocated(tao_branch%bunch_params)) then
      call out_io (s_fatal$, r_name, 'BUNCH_PARAMS NOT ALLOCATED.')
      return
    endif
 
    call bracket_index (tao_branch%bunch_params(:)%s, 0, n_ele_track, s_now, ix)
    bunch_params0 => tao_branch%bunch_params(ix)
    bunch_params1 => tao_branch%bunch_params(min(ix,n_ele_track))

    if (bunch_params0%s == bunch_params1%s) then
      r_bunch = 0
    else
      r_bunch = (s_now - bunch_params0%s) / (bunch_params1%s - bunch_params0%s)
    endif

    ! Cannot return if no particles since bunch may be injected not at the start of the branch.
    if (bunch_params0%n_particle_live == 0) then
      good(ii) = .true.
      cycle
    endif

    ix_ele = element_at_s (lat, s_now, .true.)
    ele = branch%ele(ix_ele)
    orbit%vec = (1-r_bunch) * bunch_params0%centroid%vec + r_bunch * bunch_params1%centroid%vec

    bunch_params%sigma = (1-r_bunch) * bunch_params0%sigma + (1-r_bunch) * bunch_params1%sigma

    bunch_params%a%emit = (1-r_bunch) * bunch_params0%a%emit + (1-r_bunch) * bunch_params1%a%emit
    bunch_params%b%emit = (1-r_bunch) * bunch_params0%b%emit + (1-r_bunch) * bunch_params1%b%emit

    bunch_params%x%emit = (1-r_bunch) * bunch_params0%x%emit + (1-r_bunch) * bunch_params1%x%emit
    bunch_params%y%emit = (1-r_bunch) * bunch_params0%y%emit + (1-r_bunch) * bunch_params1%y%emit

    bunch_params%a%norm_emit = (1-r_bunch) * bunch_params0%a%norm_emit + (1-r_bunch) * bunch_params1%a%norm_emit
    bunch_params%b%norm_emit = (1-r_bunch) * bunch_params0%b%norm_emit + (1-r_bunch) * bunch_params1%b%norm_emit
    bunch_params%z%norm_emit = (1-r_bunch) * bunch_params0%z%norm_emit + (1-r_bunch) * bunch_params1%z%norm_emit

  case ('lat')
    if (cache_status == using_cache$) then
      ele = tao_branch%plot_cache(ii)%ele
      orbit = tao_branch%plot_cache(ii)%orbit
      mat6 = tao_branch%plot_cache(ii)%ele%mat6
      vec0 = tao_branch%plot_cache(ii)%ele%vec0

    else
      if (first_time) then
        call twiss_and_track_at_s (lat, s_now, ele, orb, orbit, ix_branch, err)
        orbit_end = orbit
        first_time = .false.
      else
        call twiss_and_track_from_s_to_s (branch, orbit, s_now, orbit_end, ele, ele, err)
        orbit = orbit_end
      endif

      if (err) then
        good(ii:) = .false.
        bmad_com%radiation_fluctuations_on = radiation_fluctuations_on
        return
      endif

      if (orbit_end%state /= alive$) then
        write (curve%message_text, '(f10.3)') s_now
        curve%message_text = trim(curve%data_type) // ': Particle lost at s = ' // &
                             trim(adjustl(curve%message_text))
        bmad_com%radiation_fluctuations_on = radiation_fluctuations_on
        return
      endif

      if (data_type == 'momentum_compaction' .or. data_type == 'r56_compaction') then
        if (first_time) then
          call mat6_from_s_to_s (lat, mat6, vec0, ele_ref%s, s_now, orb_ref, ix_branch, err_flag = err)
          first_time = .false.
        else
          mat6 = matmul(ele%mat6, mat6)
          vec0 = matmul(ele%mat6, vec0) + ele%vec0
        endif
      endif

      if (cache_status == loading_cache$) then
        tao_branch%plot_cache(ii)%ele = ele
        tao_branch%plot_cache(ii)%orbit = orbit
        tao_branch%plot_cache(ii)%ele%mat6  = mat6
        tao_branch%plot_cache(ii)%ele%vec0 = vec0
      endif
    endif

  case default
    call out_io (s_fatal$, r_name, 'I DO NOT KNOW HOW TO HANDLE THIS curve%data_source: ' // curve%data_source)
    return
  end select

  !-------------------------------

  select case (data_type_select)
  case ('alpha.')
    select case (data_type)
    case ('alpha.a')
      value = ele%a%alpha
    case ('alpha.b')
      value = ele%b%alpha
    case default
      goto 9000  ! Error message & Return
    end select

  case ('apparent_emit.', 'norm_apparent_emit.')
    select case (data_type)
    case ('apparent_emit.x', 'norm_apparent_emit.x')
      if (curve%data_source == 'beam') then
        value = tao_beam_emit_calc (x_plane$, apparent_emit$, ele, bunch_params)
      else
        value = tao_lat_emit_calc (x_plane$, apparent_emit$, ele, tao_branch%modes)
      endif
      if (data_type_select(1:4) == 'norm') value = value * ele%value(E_tot$) / mass_of(branch%param%particle)
    case ('apparent_emit.y', 'norm_apparent_emit.y')
      if (curve%data_source == 'beam') then
        value = tao_beam_emit_calc (y_plane$, apparent_emit$, ele, bunch_params)
      else
        value = tao_lat_emit_calc (y_plane$, apparent_emit$, ele, tao_branch%modes)
      endif
      if (data_type_select(1:4) == 'norm') value = value * ele%value(E_tot$) / mass_of(branch%param%particle)
    case default
      goto 9000  ! Error message & Return
    end select

  case ('b_curl.')
    call em_field_derivatives (ele, branch%param, orbit%s-ele%s_start, orbit, .false., field)
    select case (data_type)
    case ('b_curl.x')
      value = field%dB(2,3) - field%dB(3,2)
    case ('b_curl.y')
      value = field%dB(3,1) - field%dB(1,3)
    case ('b_curl.z')
      value = field%dB(1,2) - field%dB(2,1)
    case default
      goto 9000  ! Error message & Return
    end select

  case ('b_div')
    call em_field_derivatives (ele, branch%param, orbit%s-ele%s_start, orbit, .false., field)
    value = field%dB(1,1) + field%dB(2,2) + field%dB(3,3)

  case ('b_field.')
    call em_field_calc (ele, branch%param, orbit%s-ele%s_start, orbit, .false., field)
    select case (data_type)
    case ('b_field.x')
      value = field%b(1)
    case ('b_field.y')
      value = field%b(2)
    case ('b_field.z')
      value = field%b(3)
    case default
      goto 9000  ! Error message & Return
    end select

  case ('beta.')
    select case (data_type)
    case ('beta.a')
      value = ele%a%beta
    case ('beta.b')
      value = ele%b%beta
    case default
      goto 9000  ! Error message & Return
    end select

  case ('cbar.')
    select case (data_type)
    case ('cbar.11')
      call c_to_cbar (ele, cbar)
      value = cbar(1,1)
    case ('cbar.12')
      call c_to_cbar (ele, cbar)
      value = cbar(1,2)
    case ('cbar.21')
      call c_to_cbar (ele, cbar)
      value = cbar(2,1)
    case ('cbar.22')
      call c_to_cbar (ele, cbar)
      value = cbar(2,2)
    case default
      goto 9000  ! Error message & Return
    end select

  case ('coupling.')
    select case (data_type)
    case ('coupling.11b')
      call c_to_cbar (ele, cbar)
      value = cbar(1,1) * sqrt(ele%a%beta/ele%b%beta) / ele%gamma_c
    case ('coupling.12a')
      call c_to_cbar (ele, cbar)
      value = cbar(1,2) * sqrt(ele%b%beta/ele%a%beta) / ele%gamma_c
    case ('coupling.12b')
      call c_to_cbar (ele, cbar)
      value = cbar(1,2) * sqrt(ele%a%beta/ele%b%beta) / ele%gamma_c
    case ('coupling.22a')
      call c_to_cbar (ele, cbar)
      value = cbar(2,2)* sqrt(ele%b%beta/ele%a%beta) / ele%gamma_c
    case default
      goto 9000  ! Error message & Return
    end select

  case ('e_curl.')
    call em_field_derivatives (ele, branch%param, orbit%s-ele%s_start, orbit, .false., field)
    select case (data_type)
    case ('e_curl.x')
      value = field%dE(2,3) - field%dE(3,2)
    case ('e_curl.y')
      value = field%dE(3,1) - field%dE(1,3)
    case ('e_curl.z')
      value = field%dE(1,2) - field%dE(2,1)
    case default
      goto 9000  ! Error message & Return
    end select

  case ('e_div')
    call em_field_derivatives (ele, branch%param, orbit%s-ele%s_start, orbit, .false., field)
    value = field%dE(1,1) + field%dE(2,2) + field%dE(3,3)

  case ('e_field.')
    call em_field_calc (ele, branch%param, orbit%s-ele%s_start, orbit, .false., field)
    select case (data_type)
    case ('e_field.x')
      value = field%e(1)
    case ('e_field.y')
      value = field%e(2)
    case ('e_field.z')
      value = field%e(3)
    case default
      goto 9000  ! Error message & Return
    end select

  case ('element_attrib.')
    name = upcase(curve%data_source(16:))
    ele_dum%key = overlay$  ! so entire attribute name table will be searched
    i = attribute_index(ele_dum, name)
    if (i < 1) then
      good = .false.
      return  ! Bad attribute name
    endif
    call pointer_to_attribute (ele_ref, name, .false., a_ptr, err, .false.)
    if (associated (a_ptr%r)) value = a_ptr%r

  case ('e_tot_ref')
    value = ele%value(E_tot$)

  case ('emit.')
    select case (data_type)
    case ('emit.a')
      value = bunch_params%a%emit
    case ('emit.b')
      value = bunch_params%b%emit
    case ('emit.x', 'norm_emit.x')
      if (curve%data_source == 'beam') then
        value = bunch_params%x%emit
      else
        value = tao_lat_emit_calc (x_plane$, projected_emit$, ele, tao_branch%modes)
      endif
      if (data_type_select(1:4) == 'norm') value = value * ele%value(E_tot$) / mass_of(branch%param%particle)
    case ('emit.y', 'norm_emit.y')
      if (curve%data_source == 'beam') then
        value = bunch_params%y%emit
      else
        value = tao_lat_emit_calc (y_plane$, projected_emit$, ele, tao_branch%modes)
      endif
      if (data_type_select(1:4) == 'norm') value = value * ele%value(E_tot$) / mass_of(branch%param%particle)
    case default
      goto 9000  ! Error message & Return
    end select

  case ('eta.')
    select case (data_type)
    case ('eta.x')
      value = ele%x%eta
    case ('eta.y')
      value = ele%y%eta
    case ('eta.z')
      value = ele%z%eta
    case ('eta.a')
      value = ele%a%eta
    case ('eta.b')
      value = ele%b%eta
    case default
      goto 9000  ! Error message & Return
    end select

  case ('etap.')
    select case (data_type)
    case ('etap.x')
      value = ele%x%etap
    case ('etap.y')
      value = ele%y%etap
    case ('etap.a')
      value = ele%a%etap
    case ('etap.b')
      value = ele%b%etap
    case default
      goto 9000  ! Error message & Return
    end select

  case ('ref_time')
    value = ele%ref_time

  case ('floor.')
    select case (data_type)
    case ('floor.x')
      value = ele%floor%r(1)
    case ('floor.y')
      value = ele%floor%r(2)
    case ('floor.z')
      value = ele%floor%r(3)
    case ('floor.theta')
      value = ele%floor%theta
    case ('floor.phi')
      value = ele%floor%phi
    case ('floor.psi')
      value = ele%floor%psi
    case default
      goto 9000  ! Error message & Return
    end select

  case ('momentum')
    value = orbit%p0c * (1 + orbit%vec(6)) 

  case ('momentum_compaction')
    call make_v_mats (ele_ref, v_mat, v_inv_mat)
    eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
    eta_vec = matmul (v_mat, eta_vec)
    one_pz = 1 + orb_ref%vec(6)
    eta_vec(2) = eta_vec(2) * one_pz + orb_ref%vec(2) / one_pz
    eta_vec(4) = eta_vec(4) * one_pz + orb_ref%vec(4) / one_pz
    ds = ele%s - branch%ele(0)%s
    if (ds == 0) then
      value = 0
    else
      value = -(sum(mat6(5,1:4) * eta_vec) + mat6(5,6)) / ds
    endif

  case ('norm_emit.')
    select case (data_type)
    case ('norm_emit.a')
      value = bunch_params%a%norm_emit
    case ('norm_emit.b')
      value = bunch_params%b%norm_emit
    case ('norm_emit.z')
      value = bunch_params%z%norm_emit
    case default
      goto 9000  ! Error message & Return
    end select

  case ('orbit.')
    select case (data_type)
    case ('orbit.e_tot')
      if (orbit%beta == 0) then
        value = mass_of(branch%param%particle)
      else
        value = orbit%p0c * (1 + orbit%vec(6)) / orbit%beta
      endif
    case ('orbit.x')
      value = orbit%vec(1)
    case ('orbit.y')
      value = orbit%vec(3)
    case ('orbit.z')
      value = orbit%vec(5)
    case ('orbit.px')
      value = orbit%vec(2)
    case ('orbit.py')
      value = orbit%vec(4)
    case ('orbit.pz')
      value = orbit%vec(6)
    case ('orbit.amp_a')
      call orbit_amplitude_calc (ele, orbit, amp_a = value)
    case ('orbit.amp_b')
      call orbit_amplitude_calc (ele, orbit, amp_b = value)
    case ('orbit.norm_amp_a')
      call orbit_amplitude_calc (ele, orbit, amp_na = value)
    case ('orbit.norm_amp_b')
      call orbit_amplitude_calc (ele, orbit, amp_nb = value)
    case default
      goto 9000  ! Error message & Return
    end select

  case ('ping_a.')
    select case (data_type)
    case ('ping_a.amp_x')
      value = ele%gamma_c * sqrt(ele%a%beta)
    case ('ping_a.phase_x')
      value = ele%a%phi
    case ('ping_a.sin_y')
      call c_to_cbar (ele, cbar)
      amp = sqrt(ele%b%beta * (cbar(1,2)**2 + cbar(2,2)**2))
      phase = ele%a%phi + atan2(cbar(1,2), - cbar(2,2))
      value = amp * sin(phase)
    case ('ping_a.cos_y')
      call c_to_cbar (ele, cbar)
      amp = sqrt(ele%b%beta * (cbar(1,2)**2 + cbar(2,2)**2))
      phase = ele%a%phi + atan2(cbar(1,2), - cbar(2,2))
      value = amp * cos(phase)
    end select

  case ('ping_b.')
    select case (data_type)
    case ('ping_b.amp_y')
      value = ele%gamma_c * sqrt(ele%b%beta)
    case ('ping_b.phase_y')
      value = ele%b%phi
    case ('ping_b.sin_x')
      call c_to_cbar (ele, cbar)
      amp = sqrt(ele%a%beta * (cbar(1,2)**2 + cbar(1,1)**2))
      phase = ele%b%phi + atan2(cbar(1,2), cbar(1,1))
      value = amp * sin(phase)
    case ('ping_b.cos_x')
      call c_to_cbar (ele, cbar)
      amp = sqrt(ele%a%beta * (cbar(1,2)**2 + cbar(1,1)**2))
      phase = ele%b%phi + atan2(cbar(1,2), cbar(1,1))
      value = amp * cos(phase)
    end select

  case ('phase.')
    select case (data_type)
    case ('phase.a')
      value = ele%a%phi
    case ('phase.b')
      value = ele%b%phi
    case default
      goto 9000  ! Error message & Return
    end select

  case ('r.')
    if (ii == 1) call mat_make_unit (mat6)
    if (s_now < s_last) cycle
    i = tao_read_this_index (data_type, 3); if (i == 0) return
    j = tao_read_this_index (data_type, 4); if (j == 0) return
    call mat6_from_s_to_s (lat, mat6, vec0, s_last, s_now, orbit_last, ix_branch, unit_start = .false.)
    value = mat6(i, j)

  case ('r56_compaction')
    call make_v_mats (ele_ref, v_mat, v_inv_mat)
    eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
    eta_vec = matmul (v_mat, eta_vec)
    one_pz = 1 + orb_ref%vec(6)
    eta_vec(2) = eta_vec(2) * one_pz + orb_ref%vec(2) / one_pz
    eta_vec(4) = eta_vec(4) * one_pz + orb_ref%vec(4) / one_pz
    value = sum(mat6(5,1:4) * eta_vec) + mat6(5,6)

  case ('sigma.')
    select case (data_type)
    case ('sigma.x')
      value = sqrt(bunch_params%sigma(1,1))
    case ('sigma.px')
      value = sqrt(bunch_params%sigma(2,2))
    case ('sigma.y')
      value = sqrt(bunch_params%sigma(3,3))
    case ('sigma.py')
      value = sqrt(bunch_params%sigma(4,4))
    case ('sigma.z')
      value = sqrt(bunch_params%sigma(5,5))
    case ('sigma.pz')
      value = sqrt(bunch_params%sigma(6,6))
    case default
      goto 9000  ! Error message & Return
    end select

  case ('spin.')
    select case (data_type)
    case ('spin.theta')
      polar = vec_to_polar(orbit%spin)
      value = polar%theta
    case ('spin.phi')
      polar = vec_to_polar(orbit%spin)
      value = polar%phi
    case ('spin.amp')
      polar = vec_to_polar(orbit%spin)
      value = polar%polarization
    case ('spin.x')
      value = orbit%spin(1)
    case ('spin.y')
      value = orbit%spin(2)
    case ('spin.z')
      value = orbit%spin(3)
    case default
      goto 9000  ! Error message & Return
    end select

  case ('time')
    value = orbit%t

  case ('t.')
    if (ii == 1) then
      call twiss_and_track_at_s (lat, s_last, orb = orb, orb_at_s = orbit, ix_branch = ix_branch)
      call taylor_make_unit (t_map, orbit%vec)
    endif
    if (s_now < s_last) cycle
    i = tao_read_this_index (data_type, 3); if (i == 0) goto 9000
    j = tao_read_this_index (data_type, 4); if (j == 0) goto 9000
    k = tao_read_this_index (data_type, 5); if (k == 0) goto 9000
    call transfer_map_from_s_to_s (lat, t_map, s_last, s_now, ix_branch = ix_branch, unit_start = .false.)
    value = taylor_coef (t_map(i), taylor_expn([j, k]))

  case ('tt.')
    if (ii == 1) then
      call twiss_and_track_at_s (lat, s_last, orb = orb, orb_at_s = orbit, ix_branch = ix_branch)
      call taylor_make_unit (t_map, orbit%vec)
    endif
    if (s_now < s_last) cycle
    expnt = 0
    i = tao_read_this_index (data_type, 4); if (i == 0) goto 9000
    do j = 5, 15
      if (data_type(j:j) == '') exit
      k = tao_read_this_index (data_type, j); if (k == 0) goto 9000
      expnt(k) = expnt(k) + 1
    enddo
    call transfer_map_from_s_to_s (lat, t_map, s_last, s_now, ix_branch = ix_branch, unit_start = .false.)
    value = taylor_coef (t_map(i), expnt)
  
  case default
    goto 9000  ! Error message & Return
  end select

  curve%y_line(ii) = curve%y_line(ii) + comp_sign * value
  s_last = s_now

enddo

bmad_com%radiation_fluctuations_on = radiation_fluctuations_on

return

! Error message

9000 continue
call out_io (s_warn$, r_name, &
                    'For the smooth curve calculation: I do not know about this data_type: ' // data_type)
call out_io (s_blank$, r_name, "Will not perform any smoothing.")
good = .false.
bmad_com%radiation_fluctuations_on = radiation_fluctuations_on


end subroutine tao_calc_data_at_s

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_data_useit_plot_calc (curve, graph, data)
!
! Subroutine to set the data for plotting.
!
! Input:
!   graph   -- tao_graph_struct
!   curve   -- tao_curve_struct
!
! Output:
!   data     -- Tao_data_struct:
!     %useit_plot -- True if good for plotting.
!                  = %exists & %good_plot (w/o measured & reference data)
!                  = %exists & %good_plot & %good_user & %good_meas (w/ meas data)
!                  = %exists & %good_plot & %good_user & %good_ref (w/ reference data)
!                  = %exists & %good_plot & %good_user & %good_meas & %good_ref 
!                                                        (w/ measured & reference data)
!-

subroutine tao_data_useit_plot_calc (curve, graph, data)

implicit none

type (tao_curve_struct) curve
type (tao_graph_struct) graph
type (tao_data_struct) data(:)
character(60) component

!

component = tao_curve_component(curve, graph)
data%useit_plot = data%exists .and. data%good_plot .and. &
                  (data%good_user .or. .not. graph%draw_only_good_user_data_or_vars)

if (index(component, 'meas') /= 0)   data%useit_plot = data%useit_plot .and. data%good_meas
if (index(component, 'ref') /= 0)    data%useit_plot = data%useit_plot .and. data%good_ref
if (index(component, 'model') /= 0)  data%useit_plot = data%useit_plot .and. data%good_model

end subroutine tao_data_useit_plot_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_curve_datum_calc (eles, plot, curve, who, valid)
!
! Routine to calculate datum values. 
! The values are calculated at the end of each eles(:)%ele element.
!
! Input:
!   eles(:)   -- ele_pointer_struct: Array of elements.
!   plot      -- Tao_plot_struct:
!   curve     -- Tao_curve_struct:
!   who       -- Character(*): Where to put the data. 
!                  Either: "SYMBOL" or "LINE".
!
! Output:
!   curve -- Tao_curve_struct: Structure holding the datum values
!   valid -- Logical: Set True is OK. False otherwise.
!-

subroutine tao_curve_datum_calc (eles, plot, curve, who, valid)

implicit none

type (tao_plot_struct) plot
type (tao_curve_struct) curve
type (tao_universe_struct), pointer :: u
type (tao_data_struct) datum
type (taylor_struct) t_map(6)
type (ele_pointer_struct), allocatable :: eles(:)

real(rp) y_val

logical valid, err
logical, allocatable :: good(:)

character(*) who
character(20) :: r_name = 'tao_curve_datum_calc'
character(80) why_invalid

integer i, j, m, ie, n_dat

! calculate the y-axis data point values.

u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
n_dat = size(eles)

call re_allocate (good, n_dat)   ! allocate space for the data
call re_allocate (scratch%y_value, n_dat) ! allocate space for the data

scratch%y_value = 0
good = .true.

datum%exists         = .true.
datum%ix_ele_ref     = curve%ix_ele_ref_track
datum%ix_ele_start   = -1
datum%ele_start_name = ''
datum%merit_type     = 'target'
datum%data_type      = curve%data_type
datum%ele_ref_name   = curve%ele_ref_name
datum%data_source    = curve%data_source
datum%ix_branch      = curve%ix_branch

call tao_split_component (tao_curve_component(curve, curve%g), scratch%comp, err)
if (err) then
  valid = .false.
  curve%g%why_invalid = 'Bad Graph / Curve Component Expression.'
  return
endif
if (tao_curve_component(curve, curve%g) == '') then
  valid = .false.
  curve%g%why_invalid = 'No Graph / Curve Components.'
  return
endif

do m = 1, size(scratch%comp)

  do ie = 1, n_dat

    datum%ix_ele = eles(ie)%ele%ix_ele
    datum%ix_branch = eles(ie)%ele%ix_branch
    datum%exists = .true.

    select case (scratch%comp(m)%name)
    case ('') 
      cycle
    case ('model')   
      call tao_evaluate_a_datum (datum, u, u%model, y_val, valid, why_invalid)
    case ('base')  
      call tao_evaluate_a_datum (datum, u, u%base, y_val, valid, why_invalid)
    case ('design')  
      call tao_evaluate_a_datum (datum, u, u%design, y_val, valid, why_invalid)
    case ('ref', 'meas')
      call out_io (s_error$, r_name, &
              'PLOT COMPONENT WHICH IS: ' // scratch%comp(m)%name, &
              '    FOR DATA_TYPE: ' // curve%data_type, &
              '    NOT ALLOWED SINCE DATA_SOURCE IS SET TO: ' // curve%data_source)
      return
    case default
      call out_io (s_error$, r_name, &
              'BAD PLOT COMPONENT: ' // scratch%comp(m)%name, &
              '    FOR DATA_TYPE: ' // curve%data_type)
      return
    end select
    scratch%y_value(ie) = scratch%y_value(ie) + scratch%comp(m)%sign * y_val
    if (.not. valid) good(ie) = .false.
    if (datum%data_type(1:3) == 'tt.' .or. datum%data_type(1:2) == 't.') then
      if (datum%ix_ele < datum%ix_ele_ref) datum%ix_ele_ref = datum%ix_ele
    endif

  enddo
enddo

if (n_dat > 0 .and. all(.not. good)) then
  valid = .false.
  curve%g%why_invalid = why_invalid
  return
endif

n_dat = count(good)

if (who == 'SYMBOL') then
  call re_allocate (curve%x_symb, n_dat) ! allocate space for the data
  call re_allocate (curve%y_symb, n_dat) ! allocate space for the data
  call re_allocate (curve%ix_symb, n_dat)
  j = 0
  do i = 1, size(eles)
    if (.not. good(i)) cycle
    j = j + 1
    if (plot%x_axis_type == 's') then
      curve%x_symb(j)  = eles(i)%ele%s
    else  ! 'index' or 'ele_index'
      curve%x_symb(j)  = eles(i)%ele%ix_ele
    endif
    curve%ix_symb(j) = eles(i)%ele%ix_ele
    curve%y_symb(j)  = scratch%y_value(i)
  enddo

else
  call re_allocate (curve%x_line, n_dat) ! allocate space for the data
  call re_allocate (curve%y_line, n_dat) ! allocate space for the data
  j = 0
  do i = 1, size(eles)
    if (.not. good(i)) cycle
    j = j + 1
    curve%x_line(j)  = eles(i)%ele%ix_ele
    curve%y_line(j)  = scratch%y_value(i)
    eles(j) = eles(i)
  enddo
  call re_allocate_eles (eles, n_dat, .true., .true.)
endif

valid = .true.

end subroutine

end module
