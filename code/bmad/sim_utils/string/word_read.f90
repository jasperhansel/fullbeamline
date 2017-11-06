!+
! Subroutine word_read (in_str, delim_list, word, ix_word, delim, delim_found, out_str, ignore_interior)
!
! Subroutine to extract the first word and its length from a string
! Also: Subroutine returns the delimiter between the words
! Note:
!     1) Leading blanks will be trimmed.
!     2) If a blank is in delim_list then it is considered an
!      "inferior" delimiter: Leading blanks are not counted as a delimiter
!      and trailing blanks are not counted as a delimiter if there is
!      any other non-blank delimiter before the next non-delim character or
!      all the trailing characters are blank
!
! If the ignore_interior argument is True, "interior" delimitors will be ignored.
! An "interior" character is a character that is enclosed in brackets: "(...)", "{...}", or "[...]".
! For example: 
!     call word_read "(;)a;b", ";", word, ix_word, delim, delim_found, out_str, .true.)
! Would result in:
!     word = "(;)a"
!     out_str = "b"
! Note: If there is no matching end bracket then all characters after the opening bracket will be considered 
! interior characters. For example, with "(abc]def", all characters after the "(" are interior characters.
!
! Input:
!   in_str          -- character(*): String constant to be parsed.
!   delim_list      -- character(*): String constant containing delimiters to be used
!                        by the program.
!   ignore_interior -- logical, optional: See above. Default is False.
!
! Output:
!   word        -- Character(*): First word with leading blanks trimmed
!   ix_word     -- Integer: Index in WORD of last character in the first word.
!                    Set to 0 if word is blank.
!   delim       -- Character(1): Delimiter found. Set to ' ' if no delimiter found
!   delim_found -- Logical: Set to true if delimiter found. False otherwise.
!   out_str     -- Character(*):  Rest of string after the delimiter.
!
! 1) If in_str has no non-blank characters then ix_word = 0 and delim_found = False.
!    
! 2) Example:
!     call word_read ("  to be : or not", ":", word, ix_word, delim, delim_found, out_str)
!
! Output:
!     word = 'to be    '
!     ix_word:  5
!     delim = ':'
!     delim_found = .true.
!     out_str = ' or not'
!
! 3) Example:
!     call word_read (',,,', ', ', word, ix_word, delim, delim_found, out_str)
!
! Output:
!     word = ' '
!     ix_word:  0
!     delim = ','
!     delim_found = .true.
!     out_str = ',,    '
!-

subroutine word_read (in_str, delim_list, word, ix_word, delim, delim_found, out_str, ignore_interior)

use utilities_mod

implicit none

character(*) in_str, out_str, word, delim_list, delim
character(1) tab, ch
parameter (tab = char(9))

integer i, j, ix_word, n_len, ix1, ix2
integer ix_b1, ix_b2, ix_b3

logical, optional :: ignore_interior
logical blank_delim_in_list, non_blank_delim_in_list
logical delim_found, non_blank_found, exterior

! Init

non_blank_delim_in_list = .false.
blank_delim_in_list = .false.
do i = 1, len(delim_list)
  if (delim_list(i:i) == ' ' .or. delim_list(i:i) == tab) then
    blank_delim_in_list = .true.
  else
    non_blank_delim_in_list = .true.
  endif
enddo

n_len = len(in_str)
non_blank_found = .false.
word    = ' '          ! default
delim = ' '            ! default if no delim found
ix_word = 0
delim_found = .true.   ! assume this for now
ix1 = 0
ix2 = 0
ix_b1 = 0
ix_b2 = 0
ix_b3 = 0
exterior = .true.

! loop over all characters

do i = 1, n_len

  ch = in_str(i:i)

  ! if a blank character...

  if (ch == tab .or. ch == ' ') then

    if (blank_delim_in_list .and. non_blank_found .and. exterior) then
      ix_word = ix2 - ix1 + 1
      word = in_str(ix1:ix2)
      goto 1000
    endif

  ! else if (non-blank) character is a delimiter

  elseif (index(delim_list, ch) /= 0 .and. exterior) then

    if (non_blank_found) then
      word = in_str(ix1:ix2)
      ix_word  = ix2 - ix1 + 1
    endif
    delim = ch
    if (i /= n_len) then
      out_str = in_str  ! Needed in case both actual args are the same.
      out_str = out_str(i+1:)
    else
      out_str = ' '
    endif
    return

  ! else if this is the first non_blank character found then start counting

  elseif (.not. non_blank_found) then
    ix1 = i     ! index for first non-blank
    ix2 = i
    non_blank_found = .true.

  ! else we are in the middle of a word so update end pointer

  else
    ix2 = i
  endif

  ! Ignore interior?

  if (logic_option(.false., ignore_interior)) then
    select case (ch)
    case ('(');  ix_b1 = ix_b1 + 1
    case (')');  ix_b1 = ix_b1 - 1
    case ('[');  ix_b2 = ix_b2 + 1
    case (']');  ix_b2 = ix_b2 - 1
    case ('{');  ix_b3 = ix_b3 + 1
    case ('}');  ix_b3 = ix_b3 - 1
    end select
    exterior = (ix_b1 == 0 .and. ix_b2 == 0 .and. ix_b3 == 0)
  endif

enddo

! here if no delim found

if (ix1 == 0) then
  word = in_str
  ix_word = 0
else
  word = in_str(ix1:ix2)
  ix_word = ix2 - ix1 + 1
endif
delim_found = .false.
out_str = ' '
return

! We have a blank delim. See if there is a "true" delimiter.

1000  continue

do j = i+1, n_len
  if (in_str(j:j) /= ' ' .and. in_str(j:j) /= tab) then
    if (index(delim_list, in_str(j:j)) /= 0) then
      delim = in_str(j:j)
      if (j == n_len) then
        out_str = ' '
      else
        out_str = in_str  ! Needed in case both actual args are the same.
        out_str = out_str(j+1:)
      endif
    else
      out_str = in_str  ! Needed in case both actual args are the same.
      out_str = out_str(j:)
    endif

    return
  endif
enddo

! we are here only if rest of string is blank

delim_found = .false.
out_str = ' '

end subroutine
