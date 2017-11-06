module tao_scale_mod

use tao_mod
use quick_plot

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_scale_cmd (where, y_min_in, y_max_in, axis, gang, turn_autoscale_off)
!
! Routine to scale a plot. If y_min = y_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   where    -- Character(*): Region to scale. Eg: "top:x"
!   y_min_in -- Real(rp): Plot y-axis min value.
!   y_max_in -- Real(rp): Plot y-axis max value.
!   axis     -- Character(*), optional: 'y', 'y2', or '' (both). Default = ''.
!   gang     -- Character(*), optional: 'gang', 'nogang', ''. Default = ''.
!   turn_autoscale_off
!            -- Logical, optional: If True (default) then turn off plot%autoscale_y logical
!               for all plots that are scaled.
!-

subroutine tao_scale_cmd (where, y_min_in, y_max_in, axis, gang, turn_autoscale_off)

implicit none

type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)

real(rp) y_min_in, y_max_in, y_min, y_max

integer i, j, ix, places, p1, p2

character(*) where
character(*), optional :: axis, gang
character(8) this_axis, this_gang
character(20) :: r_name = 'tao_scale_cmd'

logical, optional :: turn_autoscale_off
logical err

! Error check

if (present(axis)) then
  if (axis /= 'y' .and. axis /= 'y2' .and. axis /= '') then
    call out_io (s_error$, r_name, 'BAD AXIS NAME: ' // axis)
    return
  endif
endif

! Use local vars for y_min and y_max in case the actual args are something 
! like graph%y%min, etc.

y_min = y_min_in
y_max = y_max_in

! If the where argument is blank or '*', then scale all plots.

if (len_trim(where) == 0 .or. where == '*' .or. where == 'all') then
  do j = 1, size(s%plot_page%region)
    if (.not. s%plot_page%region(j)%visible) cycle
    call tao_scale_plot (s%plot_page%region(j)%plot, y_min, y_max, axis, gang)
  enddo
  return
endif

! Locate the plot by the region name given by the where argument.
! If no graph is specified then we scale all the graphs of the plot.

call tao_find_plots (err, where, 'REGION', plot, graph)
if (err) return

if (allocated(graph)) then                ! If all the graphs of a plot...
  do j = 1, size(graph)
    call tao_scale_graph (graph(j)%g, y_min, y_max, axis)
    if (logic_option(.true., turn_autoscale_off)) graph(j)%g%p%autoscale_y = .false.
  enddo
else                          ! else just the one graph...
  do i = 1, size(plot)
    call tao_scale_plot (plot(i)%p, y_min, y_max, axis, gang)
    if (logic_option(.true., turn_autoscale_off)) plot(i)%p%autoscale_y = .false.
  enddo
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_scale_plot (plot, y_min_in, y_max_in, axis, gang)

type (tao_plot_struct), target :: plot
type (tao_graph_struct), pointer :: graph
type (qp_axis_struct) axis_save

real(rp) y_min_in, y_max_in, y_min, y_max
real(rp) y_range(2), y2_range(2)

character(*), optional :: axis, gang
character(16) this_axis
character(16) :: r_name = 'tao_scale_plot'
integer i, p1, p2

logical do_gang

! Use local vars in case the actual args are something like graph%y%min, etc.

y_min = y_min_in
y_max = y_max_in

! If we scale a whole plot with auto scale then at the end all graphs
! are adjusted to have the same scale such that all the data fits on
! all the graphs.

if (.not. allocated (plot%graph)) return

y_range = [1d30, -1d30]
y2_range = [1d30, -1d30]

do i = 1, size(plot%graph)
  call tao_scale_graph (plot%graph(i), y_min, y_max, axis, y_range, y2_range)
enddo

! if auto scale was done...

this_axis = string_option ('', axis)

do_gang = plot%autoscale_gang_y
if (present(gang)) then
  if (gang == 'gang') then
    do_gang = .true.
  elseif (gang == 'nogang') then
    do_gang = .false.
  elseif (gang /= '') then
    call out_io (s_error$, r_name, 'BAD GANG SWITCH: ' // gang)
    call err_exit
  endif
endif

if (y_min == y_max .and. do_gang) then

  if (this_axis == '' .or. this_axis == 'y') then
    do i = 1, size(plot%graph)
      graph => plot%graph(i)
      if (graph%y%major_div_nominal > 0) then
        p1 = nint(0.7 * graph%y%major_div_nominal)  
        p2 = nint(1.3 * graph%y%major_div_nominal)
        call qp_calc_and_set_axis ('Y', y_range(1), y_range(2), p1, p2, 'GENERAL', graph%y%type)
        call qp_get_axis_attrib ('Y', graph%y%min, graph%y%max, graph%y%major_div, graph%y%places)
      else
        call qp_calc_axis_scale (y_range(1), y_range(2), graph%y)
      endif
    enddo
  endif

  if (this_axis == '' .or. this_axis == 'y2') then
    do i = 1, size(plot%graph)
      graph => plot%graph(i)
      if (graph%y2_mirrors_y) then
        axis_save = graph%y2
        graph%y2 = graph%y
        graph%y2%label        = axis_save%label
        graph%y2%draw_label   = axis_save%draw_label
        graph%y2%draw_numbers = axis_save%draw_numbers
      else
        call qp_calc_axis_scale (y2_range(1), y2_range(2), graph%y2)
        if (graph%y2%major_div_nominal > 0) then
          p1 = nint(0.7 * graph%y2%major_div_nominal)  
          p2 = nint(1.3 * graph%y2%major_div_nominal)
          call qp_calc_and_set_axis ('Y', y2_range(1), y2_range(2), p1, p2, 'GENERAL', graph%y2%type)
          call qp_get_axis_attrib ('Y', graph%y2%min, graph%y2%max, graph%y2%major_div, graph%y2%places)
        else
          call qp_calc_axis_scale (y2_range(1), y2_range(2), graph%y2)
        endif
      endif
    enddo
  endif

endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_scale_graph (graph, y_min, y_max, axis, y_range, y2_range)

type (tao_graph_struct), target :: graph
type (lat_struct), pointer :: lat
type (tao_ele_shape_struct), pointer :: shape
type (tao_curve_struct), pointer :: curve
type (floor_position_struct) floor, end
type (qp_axis_struct) axis_save

real(rp), optional :: y_range(2), y2_range(2)
real(rp) y_min, y_max, this_min, this_max, this_min2, this_max2, del

integer i, j, k, ix, ib, p1, p2
character(*), optional :: axis
character(4) this_axis
logical found_data, found_data2

! If specific min/max values are given then life is easy.

this_axis = string_option ('', axis)

if (y_min /= y_max) then

  if (this_axis == '' .or. this_axis == 'y' .or. graph%y2_mirrors_y) then
    graph%y%min = y_min
    graph%y%max = y_max
    if (graph%y%major_div_nominal> 0) then
      p1 = nint(0.7 * graph%y%major_div_nominal)  
      p2 = nint(1.3 * graph%y%major_div_nominal)  
    else
      p1 = graph%y%major_div
      p2 = p1
    endif
    graph%y%min = y_min
    graph%y%max = y_max
    call qp_calc_axis_divisions (y_min, y_max, p1, p2, graph%y%major_div)
    call qp_calc_axis_places (graph%y)
  endif

  ! y2

  if (graph%y2_mirrors_y) then
    axis_save = graph%y2
    graph%y2 = graph%y
    graph%y2%label        = axis_save%label
    graph%y2%draw_label   = axis_save%draw_label
    graph%y2%draw_numbers = axis_save%draw_numbers

  elseif (this_axis == '' .or. this_axis == 'y2') then
    graph%y2%min = y_min
    graph%y2%max = y_max
    if (graph%y2%major_div_nominal> 0) then
      p1 = nint(0.7 * graph%y2%major_div_nominal)  
      p2 = nint(1.3 * graph%y2%major_div_nominal)  
    else
      p1 = graph%y2%major_div
      p2 = p1
    endif
    graph%y2%min = y_min
    graph%y2%max = y_max
    call qp_calc_axis_divisions (y_min, y_max, p1, p2, graph%y2%major_div)
    call qp_calc_axis_places (graph%y2)
  endif

  return

endif

! Since y_min = y_max then autoscale: That is we need to find the 
! min/max so all the data points are within bounds.

! For a floor plan 

if (graph%type == 'floor_plan') then
  ix = tao_universe_number(graph%ix_universe)
  lat => s%u(ix)%model%lat
  this_min = 1e30
  this_max = -1e30

  do ib = 0, ubound(lat%branch, 1)
    do i = 0, lat%branch(ib)%n_ele_track
      call floor_to_screen_coords (graph, lat%branch(ib)%ele(i)%floor, end)
      if (end%r(1) > graph%x%max .or. end%r(1) < graph%x%min) cycle
      this_min = min(this_min, end%r(2))
      this_max = max(this_max, end%r(2))
    enddo
  enddo

  if (allocated(s%building_wall%section)) then
    do i = 1, size(s%plot_page%floor_plan%ele_shape)
      shape => s%plot_page%floor_plan%ele_shape(i)
      if (shape%ele_id /= 'wall::building') cycle
      if (.not. shape%draw) cycle
      do j = 1, size(s%building_wall%section)
        do k = 1, size(s%building_wall%section(j)%point)
          floor%r(1) = s%building_wall%section(j)%point(k)%x
          floor%r(2) = 0
          floor%r(3) = s%building_wall%section(j)%point(k)%z
          floor%theta = 0
          call floor_to_screen_coords (graph, floor, end)
          if (end%r(1) > graph%x%max .or. end%r(1) < graph%x%min) cycle
          this_min = min(this_min, end%r(2))
          this_max = max(this_max, end%r(2))
        enddo
      enddo
    enddo
  endif

! For a lat_layout: Default is [-100, 100]

elseif (graph%type == 'lat_layout') then
    this_min = -1
    this_max = 1

! Not a floor plan

else
  if (.not. allocated (graph%curve)) return
  if (.not. graph%valid) return  ! Do not scale until there is valid data

  this_min =  1e30
  this_max = -1e30
  this_min2 =  1e30
  this_max2 = -1e30
  found_data = .false.
  found_data2 = .false.

  do i = 1, size(graph%curve)
    curve => graph%curve(i)

    if (allocated(curve%y_symb) .and. curve%draw_symbols) then
      if (size(curve%y_symb) > 0) then
        if (curve%use_y2) then
          this_min2 = min(this_min2, minval(curve%y_symb))
          this_max2 = max(this_max2, maxval(curve%y_symb))
          found_data2 = .true.
        else
          this_min = min(this_min, minval(curve%y_symb))
          this_max = max(this_max, maxval(curve%y_symb))
          found_data = .true.
        endif
      endif
    endif

    if (allocated(curve%y_line) .and. curve%draw_line) then
      if (size(curve%y_line) > 0) then
        if (curve%use_y2) then
          this_min2 = min(this_min2, minval(curve%y_line))
          this_max2 = max(this_max2, maxval(curve%y_line))
          found_data2 = .true.
        else
          this_min = min(this_min, minval(curve%y_line))
          this_max = max(this_max, maxval(curve%y_line))
          found_data = .true.
        endif
      endif
    endif

  enddo

  if (.not. found_data) then
    this_max = 10
    this_min = -10
  endif

  if (this_max >  1d252) this_max =  1d252
  if (this_min < -1d252) this_min = -1d252

  if (.not. found_data2) then
    this_max2 = 10
    this_min2 = -10
  endif

  if (this_max2 >  1d252) this_max2 =  1d252
  if (this_min2 < -1d252) this_min2 = -1d252

endif

if (graph%y%major_div_nominal > 0) then
  p1 = nint(0.7 * graph%y%major_div_nominal)  
  p2 = nint(1.3 * graph%y%major_div_nominal)
else
  p1 = graph%y%major_div
  p2 = p1
endif

del = this_max - this_min
this_min = this_min - graph%scale_margin%y1 * del
this_max = this_max + graph%scale_margin%y2 * del

if (axis == '' .or. axis == 'y' .or. graph%y2_mirrors_y) then
  call qp_calc_and_set_axis ('Y', this_min, this_max, p1, p2, 'GENERAL', graph%y%type)
  if (present(y_range)) y_range = [min(y_range(1), this_min), max(y_range(2), this_max)]
  call qp_get_axis_attrib ('Y', graph%y%min, graph%y%max, graph%y%major_div, graph%y%places)
endif

! y2

if (graph%y2_mirrors_y) then
  axis_save = graph%y2
  graph%y2 = graph%y
  graph%y2%label        = axis_save%label
  graph%y2%draw_label   = axis_save%draw_label
  graph%y2%draw_numbers = axis_save%draw_numbers

elseif (axis == '' .or. axis == 'y2') then
  call qp_calc_and_set_axis ('Y2', this_min, this_max, p1, p2, 'GENERAL', graph%y2%type)
  if (present(y2_range)) y2_range = [min(y2_range(1), this_min), max(y2_range(2), this_max)]
  call qp_get_axis_attrib ('Y2', graph%y2%min, graph%y2%max, graph%y2%major_div, graph%y2%places)
endif

end subroutine

end module
