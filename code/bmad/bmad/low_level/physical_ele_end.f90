!+
! Function physical_ele_end (track_end, track_direction, ele_orientation, return_stream_end) result (physical_end)
!
! Rotine to determine which physical end of an element a particle is at given 
! the position in terms of upstream/downstream and the element's orientation
!
! Input:
!   track_end         -- Integer: first_track_edge$, second_track_edge$, surface$, or in_between$
!   track_direction   -- Integer: +1 or -1
!   ele_orientation   -- Integer: Either 1 = Normal or -1 = element reversed.
!   return_stream_end -- logical, optional: If True return the stream end instead of the physical end.
!                          Default is False.
!
! Output:
!   physical_end     -- Integer: Return_stream_end ->  Possibilities
!                                False             ->  entrance_end$, exit_end$, surface$, or in_between$
!                                True              ->  upstream_end$, downstream_end$ 
!-

function physical_ele_end (track_end, track_direction, ele_orientation, return_stream_end) result (physical_end)

use bmad_routine_interface, dummy => physical_ele_end

implicit none

integer track_end, track_direction, ele_orientation, physical_end
logical, optional :: return_stream_end
character(*), parameter :: r_name  = 'physical_ele_end'

!

if (track_end == surface$ .or. track_end == in_between$) then
  physical_end = track_end
  return
endif

!

if (logic_option(.false., return_stream_end)) then
  select case (ele_orientation * track_direction)
  case (1) 
    select case (track_end)
    case (first_track_edge$);   physical_end = upstream_end$
    case (second_track_edge$);  physical_end = downstream_end$
    end select

  case (-1)
    select case (track_end)
    case (first_track_edge$);   physical_end = downstream_end$
    case (second_track_edge$);  physical_end = upstream_end$
    end select
  end select

  return
endif

!

select case (ele_orientation * track_direction)

case (1) 
  select case (track_end)
  case (first_track_edge$);   physical_end = entrance_end$
  case (second_track_edge$); physical_end = exit_end$
  case default;
    call out_io (s_fatal$, r_name, 'BAD TRACK_END: \i0\ ', i_array = [track_end])
    if (global_com%exit_on_error) call err_exit
  end select

case (-1) 
  select case (track_end)
  case (first_track_edge$);   physical_end = exit_end$
  case (second_track_edge$); physical_end = entrance_end$
  case default;
    call out_io (s_fatal$, r_name, 'BAD TRACK_END: \i0\ ', i_array = [track_end])
    if (global_com%exit_on_error) call err_exit
  end select

case default
  call out_io (s_fatal$, r_name, 'BAD ELEMENT ORIENTATION: \2i4\ ', i_array = [ele_orientation, track_direction])
  if (global_com%exit_on_error) call err_exit

end select

end function physical_ele_end 
