      SUBROUTINE PRLINE(IFILE)
C
C     member of MAD INPUT PARSER
C
C---- PRINT A LINE OF DASHES                                             
C----------------------------------------------------------------------- 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      SAVE
C----------------------------------------------------------------------- 
      integer*4 ifile   ! dcs
C----------------------------------------------------------------------- 
      WRITE (IFILE,910)                                                  
      RETURN                                                             
C----------------------------------------------------------------------- 
  910 FORMAT(' ',11('------------'))                                     
C----------------------------------------------------------------------- 
      END                                                                