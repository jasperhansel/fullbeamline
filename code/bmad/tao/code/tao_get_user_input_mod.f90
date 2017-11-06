module tao_get_user_input_mod

use tao_mod
use input_mod

character(5), parameter, private :: sub_str(9) = ['[[1]]', '[[2]]', '[[3]]', '[[4]]', '[[5]]', &
                            '[[6]]', '[[7]]', '[[8]]', '[[9]]']

private tao_alias_translate

contains

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!+
! Subroutine tao_get_user_input (cmd_out, prompt_str, wait_flag, cmd_in, will_need_input)
!
! Subroutine to get input from the terminal.
!
! Input:
!   prompt_str -- Character(*), optional: Primpt string to print at terminal. If not
!                   present then s%global%prompt_string will be used.
!   wait_flag  -- logical, optional: Used for single mode: Wait state for get_a_char call.
!   cmd_in     -- Character(*), optional: Command. Used, for example, when using Python.
!
! Output:
!   cmd_out    -- Character(*): Command from the user.
!   will_need_input
!              -- logical, optional :: This argument is set to True when the local command buffer is 
!                   empty and more input will be needed when this routine is next called.
!-

subroutine tao_get_user_input (cmd_out, prompt_str, wait_flag, cmd_in, will_need_input)

implicit none

type do_loop_struct
  character(20) :: name ! do loop index name
  integer index, start, end, step ! for do loops
  integer n_line_start, n_line_end ! lines in each nested loop
  integer value
end type

type (do_loop_struct), allocatable, save :: loop(:)

integer i, j, ix, ix1, ix2
integer, save :: in_loop = 0 ! in loop nest level
integer ios, n_level

character(*) :: cmd_out
character(*), optional :: prompt_str, cmd_in
character(80) prompt_string, color_prompt_string
character(40) name
character(200), save :: saved_line
character(40) :: r_name = 'tao_get_user_input'

logical, optional :: wait_flag, will_need_input
logical err, wait, flush, boldit
logical, save :: init_needed = .true.

! Init single char input

if (init_needed) then
#ifdef CESR_WINCVF
#else
  call init_tty_char
#endif
  init_needed = .false.
endif

! Init

if (present(cmd_in)) cmd_out = cmd_in

prompt_string = s%global%prompt_string
if (present(prompt_str)) prompt_string = prompt_str
color_prompt_string = prompt_string

select case (s%global%prompt_color)
case ('BLACK');   call add_color (prompt_string, color_prompt_string, black_color)
case ('RED');     call add_color (prompt_string, color_prompt_string, red_color)
case ('GREEN');   call add_color (prompt_string, color_prompt_string, green_color)
case ('YELLOW');  call add_color (prompt_string, color_prompt_string, yellow_color)
case ('BLUE');    call add_color (prompt_string, color_prompt_string, blue_color)
case ('MAGENTA'); call add_color (prompt_string, color_prompt_string, magenta_color)
case ('CYAN');    call add_color (prompt_string, color_prompt_string, cyan_color)
case ('GRAY');    call add_color (prompt_string, color_prompt_string, gray_color)
end select

! If single character input wanted then...

if (s%com%single_mode) then
  if (s%global%wait_for_CR_in_single_mode) then
    if (s%com%single_mode_buffer == '') then
      do
        read '(a)', s%com%single_mode_buffer
        if (s%com%single_mode_buffer /= '') exit
      enddo
    endif
    cmd_out(1:1) = s%com%single_mode_buffer(1:1)
    s%com%single_mode_buffer = s%com%single_mode_buffer(2:)
  else
    wait = logic_option(.true., wait_flag)
    call get_a_char (cmd_out(1:1), wait, [' '])  ! ignore blanks
    s%com%cmd_from_cmd_file = .false.
  endif
  return
endif

s%com%single_mode_buffer = '' ! Reset buffer when not in single mode

! check if we still have something from a line with multiple commands

if (s%com%multi_commands_here) then
  call string_trim (saved_line, saved_line, ix)
  if (ix == 0) then
    s%com%multi_commands_here = .false.
  else
    cmd_out = saved_line
  endif
endif

! If recalling a command from the cmd history stack...

if (s%com%use_cmd_here) then
  cmd_out = s%com%cmd
  call out_io (s_blank$, r_name, '  ' // cmd_out)
  cmd_out = tao_alias_translate (cmd_out, err)
  s%com%use_cmd_here = .false.
  return
endif

! If a command file is open then read a line from the file.

n_level = s%com%cmd_file_level
if (n_level == 0) s%global%quiet = .false.  ! Disable if not running from a command file

if (n_level /= 0 .and. .not. s%com%cmd_file(n_level)%paused) then

  call output_direct (do_print = s%global%command_file_print_on)

  if (.not. s%com%multi_commands_here) then
    do
      read (s%com%cmd_file(n_level)%ix_unit, '(a)', end = 8000, iostat = ios) cmd_out
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'CANNOT READ LINE FROM FILE: ' // s%com%cmd_file(n_level)%name)
        goto 8000
      endif

      s%com%cmd_file(n_level)%n_line = s%com%cmd_file(n_level)%n_line + 1
      s%com%cmd_from_cmd_file = .true.
      call string_trim (cmd_out, cmd_out, ix)
      if (ix /= 0) exit
    enddo

    ! nothing more to do if an alias definition

    if (cmd_out(1:5) == 'alias') then
      if (.not. s%global%quiet) call out_io (s_blank$, r_name, '', trim(color_prompt_string) // ': ' // trim(cmd_out))
      return
    endif

    ! replace argument variables

    do i = 1, 9
      do j = 1, 10
        ix = index (cmd_out, sub_str(i))
        if (ix == 0) exit
        cmd_out = cmd_out(1:ix-1) // trim(s%com%cmd_file(n_level)%cmd_arg(i)) // cmd_out(ix+5:)
      enddo
    enddo

    loop1: do
      ix1 = index(cmd_out, '[[')
      ix2 = index(cmd_out, ']]')
      if (ix1 == 0) exit
      if (.not. allocated(loop) .or. ix2 < ix1) then
        call out_io (s_error$, r_name, 'MALFORMED LINE IN COMMAND FILE: ' // cmd_out)
        call tao_abort_command_file()
        cmd_out = ''
        return
      endif
      name = cmd_out(ix1+2:ix2-1)

      do i = 1, size(loop)
        if (name == loop(i)%name) then
          write (cmd_out, '(a, i0, a)') cmd_out(1:ix1-1), loop(i)%value, cmd_out(ix2+2:) 
          cycle loop1
        endif
      enddo

      call out_io (s_error$, r_name, 'CANNOT MATCH NAME IN [[...]] CONSTRUCT: ' // cmd_out)
      call tao_abort_command_file()
      cmd_out = ''
      return
      
    enddo loop1

    if (.not. s%global%quiet) call out_io (s_blank$, r_name, '', trim(color_prompt_string) // ': ' // trim(cmd_out))
    
    ! Check if in a do loop
    call do_loop(cmd_out)
    
  endif

  cmd_out = tao_alias_translate (cmd_out, err)
  call check_for_multi_commands

  if (s%com%multi_commands_here) then
    call out_io (s_blank$, r_name, '', trim(color_prompt_string) // ': ' // trim(cmd_out))
  endif

  return

  !-----------------------------------------

  8000 continue
  close (s%com%cmd_file(n_level)%ix_unit)
  s%com%cmd_file(n_level)%ix_unit = 0 
  s%com%cmd_file_level = n_level - 1 ! signal that the file has been closed
  cmd_out = ''
  ! If still lower nested command file to complete then return
  if (s%com%cmd_file_level /= 0) then
    if (s%com%cmd_file(n_level-1)%paused) then
      call out_io (s_info$, r_name, 'To continue the paused command file type "continue".')
    else
      return 
    endif
  endif
  call output_direct (do_print = s%com%print_to_terminal)
endif

! Here if no command file is being used.

if (.not. present(cmd_in)) then
  if (.not. s%com%multi_commands_here) then
    cmd_out = ''
    s%com%cmd_from_cmd_file = .false.
    boldit = (s%global%prompt_color /= '' .and. s%global%prompt_color /= 'DEFAULT')
    call out_io (s_blank$, r_name, '')
    call read_a_line (trim(prompt_string) // '> ', cmd_out, prompt_color = s%global%prompt_color, prompt_bold = boldit)
    if (cmd_out == achar(24)) cmd_out = 'exit'   ! Cntl-D pressed
  endif
endif

if (present (will_need_input)) then
  will_need_input = .false.
  if (.not. s%com%multi_commands_here .and. (s%com%cmd_file_level == 0 .or. s%com%cmd_file(n_level)%paused)) then
     will_need_input = .true.
  endif
endif

cmd_out = tao_alias_translate (cmd_out, err)
call check_for_multi_commands

if (s%com%multi_commands_here) then
  call out_io (s_blank$, r_name, '', trim(color_prompt_string) // ': ' // trim(cmd_out))
endif

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
contains

subroutine check_for_multi_commands

integer ix
character(1) quote

!

if (cmd_out(1:5) == 'alias') return

quote = ''   ! Not in quote string
do ix = 1, len(cmd_out)

  if (quote == '') then
    select case (cmd_out(ix:ix))
    case ('!')
      exit

    case (';')
      s%com%multi_commands_here = .true.
      saved_line = cmd_out(ix+1:)
      cmd_out = cmd_out(:ix-1)
      return

    case ("'", '"')
     quote = cmd_out(ix:ix)
    end select

  else ! quote /= ''
    if (cmd_out(ix:ix) == quote .and. cmd_out(ix-1:ix-1) /= '\') quote = ''           ! '
  endif

enddo

saved_line = ' '

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains

subroutine do_loop (cmd_out)

use tao_command_mod

integer ix

character(*) cmd_out
character(8) :: r_name = "do_loop"
character(20) cmd_word(9)

logical err

! Check if a "do" statement

call tao_cmd_split (cmd_out, 9, cmd_word, .false., err, '=,')
if (cmd_word(1) /= 'do' .and. cmd_word(1) /= 'enddo') return

!

if (cmd_word(1) == 'do') then

  if (cmd_word(3) /= '=' .or. .not. is_integer(cmd_word(4)) .or. &
      cmd_word(5) /= ',' .or. .not. is_integer(cmd_word(6)) .or. &
     (cmd_word(7) /= '' .and. ( &
      cmd_word(7) /= ',' .or. .not. is_integer(cmd_word(8)) .or. cmd_word(9) /= ''))) then
    call out_io (s_error$, r_name, 'MALFORMED DO STATEMENT.')
    call tao_abort_command_file()
    return
  endif

  call set_loop_level (in_loop + 1)
  loop(in_loop)%name = cmd_word(2)
  read (cmd_word(4), *) loop(in_loop)%start
  read (cmd_word(6), *) loop(in_loop)%end
  loop(in_loop)%step = 1
  if (cmd_word(7) /= '') read (cmd_word(8), *) loop(in_loop)%step

  loop(in_loop)%n_line_start = s%com%cmd_file(n_level)%n_line
  loop(in_loop)%value = loop(in_loop)%start
  call out_io (s_blank$, r_name, 'Loop: ' // trim(loop(in_loop)%name) // ' = \i0\ ', &
                                                  i_array = (/ loop(in_loop)%value /) )
  cmd_out = '' ! So tao_command will not try to parse this line.

! Check if hit 'enddo'.

elseif (cmd_word(1) == 'enddo') then
  cmd_out = ''  ! so tao_command will not try to process this.
  if (in_loop == 0) then
    call out_io (s_error$, r_name, 'ENDDO FOUND WITHOUT CORRESPODING DO STATEMENT')
    call tao_abort_command_file()
    return
  endif

  loop(in_loop)%value = loop(in_loop)%value + loop(in_loop)%step
  if ((loop(in_loop)%value <= loop(in_loop)%end .and. loop(in_loop)%step > 0) .or. &
      (loop(in_loop)%value >= loop(in_loop)%end .and. loop(in_loop)%step < 0)) then
    ! rewind
    do i = s%com%cmd_file(n_level)%n_line, loop(in_loop)%n_line_start+1, -1
      backspace (s%com%cmd_file(s%com%cmd_file_level)%ix_unit)
      s%com%cmd_file(n_level)%n_line = s%com%cmd_file(n_level)%n_line - 1
    enddo
    call out_io (s_blank$, r_name, 'Loop: ' // trim(loop(in_loop)%name) // ' = \i0\ ', &
                                                  i_array = (/ loop(in_loop)%value /) )
  else
    ! looped correct number of times, now exit loop
    call set_loop_level (in_loop-1)
  endif
  ! read next line
endif
    
end subroutine do_loop

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains

subroutine set_loop_level (level)

implicit none

type (do_loop_struct), allocatable, save :: temp(:)
integer level

! If going down then only need to set in_loop

if (level < in_loop) then
  in_loop = level
  return
endif

!

in_loop = level

if (allocated(loop)) then
  allocate (temp(level-1))
  temp = loop(1:level-1)
  deallocate (loop)
endif

allocate (loop(level))
if (allocated(temp)) then
  loop(1:size(temp)) = temp
  deallocate(temp)
endif

end subroutine set_loop_level

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains

subroutine add_color (str, c_str, color)
character(*) str, c_str, color
c_str = trim(color) // trim(str) // trim(reset_color)
end subroutine add_color

end subroutine tao_get_user_input

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

recursive function tao_alias_translate (cmd_in, err, depth) result(cmd_out)

character(*) cmd_in
character(:), allocatable :: cmd_out, tail, alias_cmd

integer, optional :: depth
integer ic, i, j, ix, this_depth
logical err
character(*), parameter :: r_name = 'tao_alias_translate'

!

if (present(depth)) then
  this_depth = depth + 1
  if (this_depth > 100) then
    call out_io (s_error$, r_name, 'INFINITE RECURSION LOOP TRANSLATING ALIASES.')
    cmd_out = ''
    err = .true.
    return
  endif
else
  this_depth = 1
  err = .false.
endif

! Look multiple commands

ix = index(cmd_in, ';')
if (ix /= 0) then
  cmd_out = tao_alias_translate(cmd_in(1:ix-1), err, this_depth)
  tail = cmd_in(ix+1:)

  do 
    ix = index(cmd_out, ';')
    if (ix == 0) exit
    cmd_out = trim(cmd_out) // ';' // tao_alias_translate(tail(1:ix-1), err, depth)
    tail = tail(ix+1:)
  enddo

  cmd_out = trim(cmd_out) // ';' // tao_alias_translate(tail, err, depth)
  return
endif

! Here if there is just a single command to translate

tail = cmd_in
call string_trim (tail, tail, ic)

do i = 1, s%com%n_alias

  if (tail(1:ic) /= s%com%alias(i)%name) cycle

  ! We have a match...
  ! Now get the actual arguments and replace dummy args with actual args.

  alias_cmd = trim(s%com%alias(i)%expanded_str)

  do j = 1, 9
    call string_trim (tail(ic+1:), tail, ic)
    do
      ix = index (alias_cmd, sub_str(j))
      if (ix == 0) exit
      alias_cmd = alias_cmd(1:ix-1) // trim(tail(1:ic)) // alias_cmd(ix+5:)
    enddo
  enddo

  ! The translation may need to be translated.

  cmd_out = tao_alias_translate (alias_cmd, err, this_depth) ! Translation is an alias?
  return
enddo

! No translation needed

cmd_out = trim(cmd_in)

end function tao_alias_translate

end module

