!****************************************************************************
!*GREXEC -- PGPLOT device handler dispatch routine
! Code converted using TO_F90 by Alan Miller
! Date: 2009-01-28  Time: 14:34:28
!+  Author: Clive Page 2009 Jan 28
      SUBROUTINE grexec(idev,ifunc,rbuf,nbuf,chr,lchr)
      use gidriv_mod
      INTEGER, INTENT(IN OUT)                  :: idev
      INTEGER, INTENT(IN OUT)                  :: ifunc
      REAL, INTENT(OUT)                        :: rbuf(*)
      INTEGER, INTENT(OUT)                     :: nbuf
      CHARACTER (LEN=*), INTENT(IN OUT)        :: chr
      INTEGER, INTENT(IN OUT)                  :: lchr
!---
      INTEGER, PARAMETER :: ndev=13
      CHARACTER (LEN=10) :: msg
!
      SELECT CASE ( idev )
      CASE(0)
          rbuf(1) = ndev
          nbuf = 1
      CASE ( 1)
          CALL nudriv(ifunc,rbuf,nbuf,chr,lchr)
      CASE ( 2)
          CALL psdriv(ifunc,rbuf,nbuf,chr,lchr,1)
      CASE ( 3)
          CALL psdriv(ifunc,rbuf,nbuf,chr,lchr,2)
      CASE ( 4)
          CALL psdriv(ifunc,rbuf,nbuf,chr,lchr,3)
      CASE ( 5)
          CALL psdriv(ifunc,rbuf,nbuf,chr,lchr,4)
      CASE (6)
          CALL gidriv(ifunc, rbuf, nbuf, chr, lchr, 1)
      CASE (7)
          CALL gidriv(ifunc, rbuf, nbuf, chr, lchr, 2)
      CASE DEFAULT
          WRITE (msg,'(I10)') idev
          CALL grwarn('Unknown device code in GREXEC: '//msg)
      END SELECT
      
      END SUBROUTINE grexec

