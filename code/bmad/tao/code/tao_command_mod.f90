module tao_command_mod

use output_mod
use tao_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_cmd_history_record (cmd)
!
! Subroutine to record a cmd in the command history stack
!-

subroutine tao_cmd_history_record (cmd)

implicit none

character(*) cmd

!

s%com%ix_history = s%com%ix_history + 1
if (s%com%ix_history > size(s%history)) s%com%ix_history = 1
s%com%n_history = s%com%n_history + 1
s%history(s%com%ix_history)%ix = s%com%n_history
if (s%com%cmd_from_cmd_file) then
  s%history(s%com%ix_history)%cmd = '  ! ' // trim(cmd)
else
  s%history(s%com%ix_history)%cmd = trim(cmd)
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_re_exectue (string, err)
!
! Subroutine to execute a previous command.
!-

subroutine tao_re_execute (string, err)

implicit none

integer ios, ix1, ix, ix_rec
character(100) line
character(*) string
character(*), parameter :: r_name = 'tao_history_cmd'
logical err

!

err = .true.

if (is_integer(string)) then
  call string_trim(string, line, ix)
  if (line(ix+1:) /= '') then
    call out_io (s_error$, r_name, 'EXTRA STUFF AFTER INTEGER INDEX.')
    return
  endif

  read (string, *, iostat = ios) ix_rec
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'ERROR READING HISTORY NUMBER')
    return
  endif

  if (ix_rec > 0) then
    if (ix_rec > s%com%n_history .or. ix_rec < s%com%n_history - (size(s%history) - 1)) then
      call out_io (s_error$, r_name, 'INVALID INDEX FOR THE HISTORY LIST.')
      return
    endif
    ix = ix_rec + s%com%ix_history - s%com%n_history
  else
    if (-ix_rec > size(s%history) - 1 .or. -ix_rec > s%com%n_history - 1) then 
      call out_io (s_error$, r_name, 'INVALID INDEX FOR THE HISTORY LIST.')
      return
    endif
    ix = s%com%ix_history + ix_rec
  endif

  if (ix < 1) ix = ix + size(s%history)

!

else

  ix = s%com%ix_history
  do

    if (index(s%history(ix)%cmd, trim(string)) == 1) exit

    ix = ix - 1
    if (ix < 1) ix = ix + size(s%history)

    if (ix == s%com%ix_history .or. s%history(ix)%ix == 0) then
      call out_io (s_error$, r_name, 'COMMAND NOT FOUND IN THE HISTORY LIST.')
      return
    endif

  enddo

endif

! put the command in the common area so it can be used next.

call string_trim(s%history(ix)%cmd, s%com%cmd, ix)
if (s%com%cmd(1:1) == '!') s%com%cmd = s%com%cmd(2:)
s%com%use_cmd_here = .true.

err = .false.

end subroutine tao_re_execute

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_cmd_split (cmd_line, n_word, cmd_word, no_extra_words, err, separator)
!
! This routine splits the command line into words.
!
! Input: 
!   cmd_line       -- Character(*): The command line.
!   n_word         -- Integer: Maximum number of words to split command line into.
!   no_extra_words -- Logical: are extra words allowed at the end?
!                        If False then cmd_word(n_word) will contain everything after 
!                        the n_word-1 word.
!   separator      -- Character(*), optional: a list of characters that,
!                        besides a blank space, signify a word boundary. 
!
! Output:
!   cmd_word(n_word) -- Character(*): The individual words.
!   err              -- Logical: error in splitting words
!
! For example: 
!   separator = '-+' 
!   cmd_line = 'model-design'
! Restult:
!   cmd_word(1) = 'model'
!   cmd_word(2) = '-'
!   cmd_word(3) = 'design'
!
! Anything between single or double quotes is treated as a single word.
! Whitespace or a separator inside of "{}", "()", or "[]" is ignored.
!-

subroutine tao_cmd_split (cmd_line, n_word, cmd_word, no_extra_words, err, separator)

integer i, ix, nw, n_word, ix_b1, ix_b2, ix_b3

character(*) cmd_line
character(*), optional :: separator
character(*) cmd_word(:)
character(16) :: r_name = 'tao_cmd_split'
character(300) line
character(1), parameter :: tab = char(9)

logical err
logical no_extra_words

!

err = .false.
line = cmd_line
cmd_word(:) = ''
nw = 0

nw_loop: do 
  call string_trim (line, line, ix)

  if (ix == 0) exit

  ! If extra words allowed, everything left goes into cmd_word(n_word)
  if (nw == n_word - 1 .and. .not. no_extra_words) then
    nw=nw+1; cmd_word(nw) = trim(line)
    return
  endif

  if (line(1:1) == '"') then
    ix = index(line(2:), '"')
    if (ix == 0) ix = len(line)
    nw=nw+1; cmd_word(nw) = line(2:ix)
    line = line(ix+1:)
    cycle
  elseif (line(1:1) == "'") then
    ix = index(line(2:), "'")
    if (ix == 0) ix = len(line)
    nw=nw+1; cmd_word(nw) = line(2:ix)
    line = line(ix+1:)
    cycle
  endif

  ix_b1 = 0; ix_b2 = 0; ix_b3 = 0
  do i = 1, len(line)

    if (line(i:i) == '{') ix_b1 = ix_b1 + 1
    if (line(i:i) == '}') ix_b1 = ix_b1 - 1
    if (line(i:i) == '(') ix_b2 = ix_b2 + 1
    if (line(i:i) == ')') ix_b2 = ix_b2 - 1
    if (line(i:i) == '[') ix_b3 = ix_b3 + 1
    if (line(i:i) == ']') ix_b3 = ix_b3 - 1

    if (ix_b1 /= 0 .or. ix_b2 /= 0 .or. ix_b3 /= 0) cycle

    if (present(separator)) then
      if (index(separator, line(i:i)) /= 0) then
        if (i /= 1) then
          nw=nw+1; cmd_word(nw) = line(1:i-1)
          line = line(i:)
        endif
        nw=nw+1; cmd_word(nw) = line(1:1)
        line = line(2:)
        cycle nw_loop
      endif
    endif

    if (line(i:i) == ' ' .or. line(i:i) == tab) then
      nw=nw+1; cmd_word(nw) = line(1:i-1)
      line = line(i+1:)
      cycle nw_loop
    endif

  enddo

  if (ix_b1 /= 0 .or. ix_b2 /= 0 .or. ix_b3 /= 0) then
    call out_io (s_error$, r_name, 'MISMATCHED "{...}", "(...)", OR "[...]".')
    err = .true.
    return
  endif

  call out_io (s_fatal$, r_name, 'INTERNAL ERROR!')
  call err_exit

enddo nw_loop

!  

if (no_extra_words .and. nw > n_word) then
  call out_io (s_error$, r_name, 'EXTRA STUFF ON COMMAND LINE: ' // line)
  err = .true.
endif

end subroutine tao_cmd_split

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine tao_next_switch (line, switch_list, return_next_word, switch, err, ix_word)
!
! Subroutine look at the next word on the command line and match this word to a list of "switches"
! given by the switch_list argument.
! 
! Switch abbreviations are permitted.
!
! If return_next_word = True then, when a non-switch word is encountered, the switch argument 
! will be set to that word and that word will be removed from the line argument.
!
! If return_next_word = False then, when a non-switch word is encountered, the switch argument 
! will be set to '' and the non-switch word will be left on the line argument.
!
! Input:
!   line              -- character(*): Command line
!   switch_list(*)    -- character(*): List of valid switches. 
!   return_next_word  -- logical: See above.
!
! Output:
!   line            -- character(*): Line with switch word removed. If the first word
!                       does not look like a switch then nothing is removed.
!   switch          -- character(*): Switch found. This is the full name
!                       even if what was on the command line was an abbreviation.
!                       See above for more details.
!   err             -- logical: Set True if the next word begins with '-' but there is no match
!                       to anything in switch_list. An error message will be printed.
!   ix_word         -- integer: Character length of first word left on line.
!-

subroutine tao_next_switch (line, switch_list, return_next_word, switch, err, ix_word)

character(*) line, switch, switch_list(:)
character(20) :: r_name = 'tao_next_switch'
logical err

integer ix, n, ix_word
logical return_next_word

!

err = .false.
switch = ''

call string_trim(line, line, ix_word)
if (ix_word == 0) return
if (line(1:1) /= '-') then
  if (return_next_word) then
    switch = line(1:ix_word)
    call string_trim(line(ix_word+1:), line, ix_word)
  endif
  return
endif

!

call match_word (line(:ix_word), switch_list, n, .true., matched_name=switch)
if (n < 1) then
  err = .true.
  if (n == 0) then 
    call out_io (s_error$, r_name, 'UNKNOWN SWITCH: ' // line(:ix_word))
  else
    call out_io (s_error$, r_name, 'AMBIGUOUS SWITCH: ' // line(:ix_word))
  endif
  return
endif

call string_trim(line(ix_word+1:), line, ix_word)

end subroutine

end module

