! Original PGPLOT code had pixmap declared as an integer and then used the C routines grgmem 
! and grfmem to malloc and free memory using the pixmap value as an *address*.
! This meant that pixmap had to be a 32-bit integer for 32-bit machines and a 64-bit integer
! for 64-bit machines.
! Note: The following routines have been converted to use pgplot_mod:
!   gidriv
!   wddriv
! The following routines have not yet been converted:
!   xedriv
!   ppdriv
!
! David Sagan 9/20/2011.

      module pgplot_mod
        BYTE, allocatable :: PIXMAP(:,:)
      end module

