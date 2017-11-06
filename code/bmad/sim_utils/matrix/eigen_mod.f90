module eigen_mod

use precision_def
use output_mod

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine mat_eigen (mat, eigen_val, eigen_vec, error, print_err)
!
! Routine for determining the eigen vectors and eigen values of a matrix.
!
! Modules needed:
!   use eigen_mod
!
! Input:
!   mat(n,n)  -- Real(rp): Matrix. n must be even.
!   print_err -- Logical, optional: If present and False then suppress 
!                  "no eigen-system found" messages.
!
! Output:
!   eigen_val(n)   -- Complex(rp): Eigenvalue.
!   eigen_vec(n,n) -- Complex(rp): Eigenvector.
!   error          -- Logical: Set True on error. False otherwise.
!-

subroutine mat_eigen (mat, eigen_val, eigen_vec, error, print_err)

implicit none

real(rp) mat(:,:)
real(rp) :: val(size(mat, 1)), vec(size(mat, 1), size(mat, 1))
integer :: iv(size(mat, 1))

complex(rp) eigen_val(:), eigen_vec(:,:)

integer i, n, ier

logical, optional :: print_err
logical error

!

n = ubound(mat, 1)
error = .false.

call eigensys (mat, val, vec, iv, n, ier, print_err)
if (ier /= 0) then
  error = .true.
  return
endif

!

call ordersys (val, vec, iv, n)

do i = 2, n, 2
  if (iv(i-1) == 0) then
    eigen_val(i-1) = val(i-1)
    eigen_val(i)   = val(i)
    eigen_vec(i-1, :) = vec(:, i-1)
    eigen_vec(i, :)   = vec(:, i)

  elseif (iv(i-1) == 1) then
    eigen_val(i-1) = cmplx(val(i-1), val(i))
    eigen_val(i)   = cmplx(val(i-1), -val(i))
    eigen_vec(i-1, :) = cmplx(vec(:, i-1), vec(:, i))
    eigen_vec(i, :)   = cmplx(vec(:, i-1), -vec(:, i))

  else
    call out_io (s_fatal$, 'mat_eigen', 'BAD IV FROM EIGENSYS')
    if (global_com%exit_on_error) call err_exit
  endif
enddo


end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      SUBROUTINE EIGENSYS(DN,DV,DM,IV,N,IER, print_err)
!     finds the N dimensional eigenvalue vector DV, the N*N dimensional
!     eigenvectormatrix DM, and the indication vector IV for
!     real/complex of a N*N dimensional matrix DN.
!
!     DN (input) : N*N dimensional matrix of real*8
!
!     DV (input) : N dimensional matrix of real*8
!                  contains the eigenvalues sorted by real and imaginary
!                  parts. Which part DV(I) actually is, is described in
!                  IV(I).
!
!     DM (input) : N*N dimensional matrix of real*8
!                  contains the eigenvectors sorted by real and imaginary
!                  parts. Which part DM(:,I) actually is, is described in
!                  IV(I).
!
!     IV (output): N dimensional vector of integers describing the type
!                  of eigenvalues and eigenvectors.  The types are not
!                  ordered in any way.
!                  IV(I)=0 indicates that DV(I) is a real eigenvalue and
!                          DM(:,I) is the corresponding real eigenvector
!                  IV(I)=1 indicates that DV(I)+I*DV(I+1) and DV(I)-I*DV(I+1)
!                          are complex eigenvalues and DM(:,I)+I*DM(:,I+1),
!                          DM(:,I)-I*DM(:,I+1) are the corresponding complex
!                          eigenvectors. The imaginary part DV(I+1) is always
!                          chosen to be positive. IV(I)=1 is always followed
!                          by IV(I+1)=2.
!                  IV(I)=2 Indicates the imaginary part of eigenvalues and
!                          eigenvectors.  It is always preceedid by IV(I-1)=1.
!
!     N  (input) : N specifies the dimensions of the matrixes DN, DM and the
!                  vectores DV, IV.  It has to be smaller than the local
!                  parameter NMX.
!
!     IER(output): Error flag, is zero when no error occured and 1, if an
!                  error occured.

!=====local stuff

      integer, parameter :: NMX=8
      integer n
      integer i, j, ii, ij

      real(rp) DN(N,N), DV(N), DM(N,N), DNS(NMX,NMX), Z(NMX,NMX), &
          ORT(NMX), VR(NMX), VI(NMX)
      integer IV(N), ier, ierr
      logical, optional :: print_err

!%%      write(*,*)'eigensys'

      IF (N.GT.NMX) then
        CALL out_io(s_fatal$, 'EIGENSYS', 'dimension too high, increase NMX')
        if (global_com%exit_on_error) call err_exit
      endif

      DNS(1:n,1:n) = DN(1:n,1:n)

!=====Produce an upper Hessenberg form of DNS by
!=====orthogonal similarity transformations and store it in DNS.
!=====The transformations are stored in the rest of DNS and in ORT.
      CALL ETY (NMX,N,1,N,DNS,ORT)

!=====Accumulate the orthogonal transformations in Z.
      CALL ETYT(NMX,N,1,N,DNS,ORT,Z)

!=====Find the real eigenvalues VR(I) and the corresponding eigenvectors
!=====Z(:,I) for real eigenvalues.  For complex eigenvalues find
!=====the eigenvector Z(:,I)+I*Z(:,I+1) corresponding to the eigenvalue
!=====VR(I)+I*VI(I) with VI(I)>0. VR(I+1)+I*VI(I+1) is then the conjugate.
      CALL ETY2(NMX,N,1,N,DNS,VR,VI,Z,IERR)

      IF(IERR.NE.0) THEN
         if (logic_option(.true., print_err)) &
              call out_io(s_error$, 'EIGENSYS', 'no eigensystem found ; skipping')
         IER = 1
         RETURN
      ELSE
         IER = 0
      ENDIF

      II = 1
      IJ = 2

      DO 20 I=1,N

         IF(VI(II).NE.0.D+0) THEN
!=====Eigenvalue is complex:
            IV(II) = 1
            IV(IJ) = 2
            DV(II) = VR(II)
            IF(VI(II).LT.0.D+0) CALL out_io(s_fatal$,  &
                'EIGENSYS', 'eigenvectors missinterpreted')
            DV(IJ) = VI(II)

            DO 120 J=1,N
               DM(J,II) = Z(J,II)
               DM(J,IJ) = Z(J,IJ)
 120        CONTINUE

            II = II+2
            IJ = II+1
         ELSE
!=====Eigenvalue is real:
            IV(II) = 0
            DV(II) = VR(II)

            DO 220 J=1,N
               DM(J,II) = Z(J,II)
 220        CONTINUE

            II = II+1
            IJ = II+1
         ENDIF

         IF(II.GT.N) GOTO 21
 20   CONTINUE
 21   CONTINUE

      RETURN
      end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      subroutine ETY(NM,N,LOW,IGH,A,ORT)
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTHES,
!     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
!     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
!     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
!     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
!          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
!          SET LOW=1, IGH=N,
!
!        A CONTAINS THE INPUT MATRIX.
!
!     ON OUTPUT-
!
!        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT
!          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION
!          IS STORED IN THE REMAINING TRIANGLE UNDER THE
!          HESSENBERG MATRIX,
!
!        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
!          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
!
!     FORTRAN ROUTINE BY B. S. GARBOW
!
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      REAL*8 A(NM,NM),ORT(IGH)
      REAL*8 F,G,H,SCALE

!

      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GOTO 200
!
      DO 180 M = KP1, LA
         H = 0.0
         ORT(M) = 0.0
         SCALE = 0.0
!     SCALE COLUMN (ALGOL TOL THEN NOT NEEDED)
         DO 90 I = M, IGH
   90    SCALE = SCALE + abs(A(I,M-1))
!
         IF (SCALE .EQ. 0.0) GOTO 180
         MP = M + IGH
!     FOR I=IGH STEP -1 UNTIL M DO --
         DO 100 II = M, IGH
            I = MP - II
            ORT(I) = A(I,M-1) / SCALE
            H = H + ORT(I) * ORT(I)
  100    CONTINUE
!
         G = -DSIGN(SQRT(H),ORT(M))
         H = H - ORT(M) * G
         ORT(M) = ORT(M) - G
!     FORM (I-(U*UT)/H) * A
         DO 130 J = M, N
            F = 0.0
            DO 110 II = M, IGH
               I = MP - II
               F = F + ORT(I) * A(I,J)
  110       CONTINUE
!
            F = F / H
!
            DO 120 I = M, IGH
  120       A(I,J) = A(I,J) - F * ORT(I)
!
  130    CONTINUE
!     FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
         DO 160 I = 1, IGH
            F = 0.0
            DO 140 JJ = M, IGH
               J = MP - JJ
               F = F + ORT(J) * A(I,J)
  140       CONTINUE
!
            F = F / H
!
            DO 150 J = M, IGH
  150       A(I,J) = A(I,J) - F * ORT(J)
!
  160    CONTINUE
!
         ORT(M) = SCALE * ORT(M)
         A(M,M-1) = SCALE * G
  180 CONTINUE
!
  200 RETURN
      end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      subroutine ETYT(NM,N,LOW,IGH,A,ORT,Z)
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTRANS,
!     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!     THIS SUBROUTINE ACCUMULATES THE ORTHOGONAL SIMILARITY
!     TRANSFORMATIONS USED IN THE REDUCTION OF A REAL GENERAL
!     MATRIX TO UPPER HESSENBERG FORM BY  ETY.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
!          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
!          SET LOW=1, IGH=N,
!
!        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
!          FORMATIONS USED IN THE REDUCTION BY  ORTHES
!          IN ITS STRICT LOWER TRIANGLE,
!
!          ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANS-
!          FORMATIONS USED IN THE REDUCTION BY  ETY.
!          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
!
!     ON OUTPUT-
!
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTION BY  ETY,
!
!        ORT HAS BEEN ALTERED.
!
!     FORTRAN ROUTINE BY B. S. GARBOW.
!
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL*8 A(NM,IGH),ORT(IGH),Z(NM,N)
!check ??? should this not be A(NM,NM),ORT(NM),Z(NM,NM) ???
      REAL*8 G

!%%      write(*,*)'etyt'
!
!     INITIALIZE Z TO IDENTITY MATRIX
      DO 80 I = 1, N
!
         DO 60 J = 1, N
   60    Z(I,J) = 0.0
!
         Z(I,I) = 1.0
   80 CONTINUE
!
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GOTO 200
      DO 140 MM = 1, KL
         MP = IGH - MM
         IF (A(MP,MP-1) .EQ. 0.0) GOTO 140
         MP1 = MP + 1
!
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
!
         DO 130 J = MP, IGH
            G = 0.0
!
            DO 110 I = MP, IGH
  110       G = G + ORT(I) * Z(I,J)
!     DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
!     DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW
            G = (G / ORT(MP)) / A(MP,MP-1)
!
            DO 120 I = MP, IGH
  120       Z(I,J) = Z(I,J) + G * ORT(I)
!
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
      end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      subroutine ETY2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
!     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
!     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
!     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
!     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
!     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
!          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
!          SET LOW=1, IGH=N,
!
!        H CONTAINS THE UPPER HESSENBERG MATRIX,
!
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
!          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
!          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
!          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
!          IDENTITY MATRIX.
!
!     ON OUTPUT-
!
!        H HAS BEEN DESTROYED,
!
!        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
!          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
!          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
!          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
!          FOR INDICES IERR+1,...,N,
!
!        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
!          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
!          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
!          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
!          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
!          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
!          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND,
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER n_max_iter ITERATIONS.
!
!     ARITHMETIC IS REAL*8. COMPLEX DIVISION
!     IS SIMULATED BY ROUTIN ETDIV.
!
!     FORTRAN ROUTINE BY B. S. GARBOW.
!
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN, &
             IGH,ITS,LOW,MP2,ENM2,IERR, n_max_iter
      REAL*8 H(NM,N),WR(N),WI(N),Z(NM,N)
!check ??? should this not be H(NM,NM),WR(NM),WI(NM),Z(NM,NM) ???
      REAL*8 P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP
      LOGICAL NOTLAS
      REAL*8 Z3R,Z3I

!     COMPLEX CMPLX
!
!     MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!     THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
      MACHEP = 10.0**(-precision(p))
!     MACHEP = R1MACH(4)
!
      IERR = 0
      n_max_iter = 50  ! Originally: 30
      NORM = 0.0
      K = 1
!     STORE ROOTS ISOLATED BY BALANC AND COMPUTE MATRIX NORM
      DO 50 I = 1, N
!
         DO 40 J = K, N
   40    NORM = NORM + abs(H(I,J))
!
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GOTO 50
         WR(I) = H(I,I)
         WI(I) = 0.0
   50 CONTINUE
!
      EN = IGH
      T = 0.0
!     SEARCH FOR NEXT EIGENVALUES
   60 IF (EN .LT. LOW) GOTO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
!     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GOTO 100
         S = abs(H(L-1,L-1)) + abs(H(L,L))
         IF (S .EQ. 0.0) S = NORM
         IF (abs(H(L,L-1)) .LE. MACHEP * S) GOTO 100
   80 CONTINUE
!     FORM SHIFT
  100 X = H(EN,EN)
      IF (L .EQ. EN) GOTO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GOTO 280
      IF (ITS .EQ. n_max_iter) GOTO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GOTO 130
!     FORM EXCEPTIONAL SHIFT
      T = T + X
!
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
!
      S = abs(H(EN,NA)) + abs(H(NA,ENM2))
      X = 0.75 * S
      Y = X
      W = -0.4375 * S * S
  130 ITS = ITS + 1
!     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS.
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = abs(P) + abs(Q) + abs(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GOTO 150
         IF (abs(H(M,M-1)) * (abs(Q) + abs(R)) .LE. MACHEP * abs(P) &
               * (abs(H(M-1,M-1)) + abs(ZZ) + abs(H(M+1,M+1)))) GOTO 150
  140 CONTINUE
!
  150 MP2 = M + 2
!
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0
         IF (I .EQ. MP2) GOTO 160
         H(I,I-3) = 0.0
  160 CONTINUE
!     DOUBLE QR STEP INVOLVING ROWS L TO EN AND COLUMNS M TO EN
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GOTO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0
         IF (NOTLAS) R = H(K+2,K-1)
         X = abs(P) + abs(Q) + abs(R)
         IF (X .EQ. 0.0) GOTO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = DSIGN(sqrt(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GOTO 180
         H(K,K-1) = -S * X
         GOTO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
!     ROW MODIFICATION
         DO 210 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            IF (.NOT. NOTLAS) GOTO 200
            P = P + R * H(K+2,J)
            H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
            H(K,J) = H(K,J) - P * X
  210    CONTINUE
!
         J = MIN0(EN,K+3)
!     COLUMN MODIFICATION
         DO 230 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            IF (.NOT. NOTLAS) GOTO 220
            P = P + ZZ * H(I,K+2)
            H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
            H(I,K) = H(I,K) - P
  230    CONTINUE
!     ACCUMULATE TRANSFORMATIONS
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            IF (.NOT. NOTLAS) GOTO 240
            P = P + ZZ * Z(I,K+2)
            Z(I,K+2) = Z(I,K+2) - P * R
  240       Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K) = Z(I,K) - P
  250    CONTINUE
!
  260 CONTINUE
!
      GOTO 70
!     ONE ROOT FOUND
  270 H(EN,EN) = X + T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0
      EN = NA
      GOTO 60
!     TWO ROOTS FOUND
  280 P = (Y - X) / 2.0
      Q = P * P + W
      ZZ = sqrt(abs(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0) GOTO 320
!     REAL PAIR
      ZZ = P + DSIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0
      WI(EN) = 0.0
      X = H(EN,NA)
      S = abs(X) + abs(ZZ)
      P = X / S
      Q = ZZ / S
      R = sqrt(P*P+Q*Q)
      P = P / R
      Q = Q / R
!     ROW MODIFICATION
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
!     COLUMN MODIFICATION
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
!     ACCUMULATE TRANSFORMATIONS
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
!
      GOTO 330
!     COMPLEX PAIR
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GOTO 60
!     ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND VECTORS OF UPPER TRIANGULAR FORM
  340 IF (NORM .EQ. 0.0) GOTO 1001
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
!     REAL VECTOR
  600    M = EN
         H(EN,EN) = 1.0
         IF (NA .EQ. 0) GOTO 800
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = H(I,EN)
            IF (M .GT. NA) GOTO 620
!
            DO 610 J = M, NA
  610       R = R + H(I,J) * H(J,EN)
!
  620       IF (WI(I) .GE. 0.0) GOTO 630
            ZZ = W
            S = R
            GOTO 700
  630       M = I
            IF (WI(I) .NE. 0.0) GOTO 640
            T = W
            IF (W .EQ. 0.0) T = MACHEP * NORM
            H(I,EN) = -R / T
            GOTO 700
!     SOLVE REAL EQUATIONS
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
            H(I,EN) = T
            IF (abs(X) .LE. abs(ZZ)) GOTO 650
            H(I+1,EN) = (-R - W * T) / X
            GOTO 700
  650       H(I+1,EN) = (-S - Y * T) / ZZ
  700    CONTINUE
!     END REAL VECTOR
         GOTO 800
!     COMPLEX VECTOR
  710    M = NA
!     LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT EIGENVECTOR MATRIX IS
!     TRIANGULAR
         IF (abs(H(EN,NA)) .LE. abs(H(NA,EN))) GOTO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
         GOTO 730
! 720    Z3 = CMPLX(0.0,-H(NA,EN)) / CMPLX(H(NA,NA)-P,Q)
!        H(NA,NA) = REAL(Z3)
!        H(NA,EN) = AIMAG(Z3)
  720    CALL ETDIV(Z3R,Z3I,0.D0,-H(NA,EN),H(NA,NA)-P,Q)
         H(NA,NA) = Z3R
         H(NA,EN) = Z3I
  730    H(EN,NA) = 0.0
         H(EN,EN) = 1.0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GOTO 800
         DO 790 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0
            SA = H(I,EN)
!
            DO 760 J = M, NA
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
!
            IF (WI(I) .GE. 0.0) GOTO 770
            ZZ = W
            R = RA
            S = SA
            GOTO 790
  770       M = I
            IF (WI(I) .NE. 0.0) GOTO 780
!           Z3 = CMPLX(-RA,-SA) / CMPLX(W,Q)
!           H(I,NA) = REAL(Z3)
!           H(I,EN) = AIMAG(Z3)
            CALL ETDIV(Z3R,Z3I,-RA,-SA,W,Q)
            H(I,NA) = Z3R
            H(I,EN) = Z3I
            GOTO 790
!     SOLVE COMPLEX EQUATIONS
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0 * Q
            IF (VR .EQ. 0.0 .AND. VI .EQ. 0.0) VR = MACHEP * NORM &
                 * (abs(W) + abs(Q) + abs(X) + abs(Y) + abs(ZZ))
!           Z3 = CMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA) / CMPLX(VR,VI)
!           H(I,NA) = REAL(Z3)
!           H(I,EN) = AIMAG(Z3)
            CALL ETDIV(Z3R,Z3I,X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI)
            H(I,NA) = Z3R
            H(I,EN) = Z3I
            IF (abs(X) .LE. abs(ZZ) + abs(Q)) GOTO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GOTO 790
! 785       Z3 = CMPLX(-R-Y*H(I,NA),-S-Y*H(I,EN)) / CMPLX(ZZ,Q)
!           H(I+1,NA) = REAL(Z3)
!           H(I+1,EN) = AIMAG(Z3)
  785       CALL ETDIV(Z3R,Z3I,-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q)
            H(I+1,NA) = Z3R
            H(I+1,EN) = Z3I
  790    CONTINUE
!     END COMPLEX VECTOR
  800 CONTINUE
!     END BACK SUBSTITUTION
!     VECTORS OF ISOLATED ROOTS
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GOTO 840
!
         DO 820 J = I, N
  820    Z(I,J) = H(I,J)
!
  840 CONTINUE
!     MULTIPLY BY TRANSFORMATION MATRIX TO GIVE VECTORS OF ORIGINAL FULL MATRIX

      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)
         DO 880 I = LOW, IGH
            ZZ = 0.0
            DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)
            Z(I,J) = ZZ
  880 CONTINUE
!
      GOTO 1001
!     SET ERROR -- NO CONVERGENCE TO AN EIGENVALUE AFTER n_max_iter ITERATIONS
 1000 IERR = EN

 1001 RETURN
      end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      subroutine ETDIV(A,B,C,D,E,F)
!     **************************************
!
!     HIGH-PRECISION DIVISION
!
      REAL*8 A,B,C,D,E,F
      REAL*8 S,T
      REAL*8 CC,DD,EE,FF
      REAL*8 TEMP
      INTEGER FLIP
      FLIP = 0
      CC = C
      DD = D
      EE = E
      FF = F
      IF( abs(F).GE.abs(E) ) THEN
        EE = F
        FF = E
        CC = D
        DD = C
        FLIP = 1
      ENDIF
      S = 1.D0/EE
      T = 1.D0/(EE+ FF*(FF*S))
      IF ( abs(FF) .GE. abs(S) ) THEN
        TEMP = FF
        FF = S
        S = TEMP
      ENDIF
      IF( abs(DD) .GE. abs(S) ) THEN
        A = T*(CC + S*(DD*FF))
      ELSE IF ( abs(DD) .GE. abs(FF) ) THEN
        A = T*(CC + DD*(S*FF))
      ELSE
        A = T*(CC + FF*(S*DD))
      ENDIF
      IF ( abs(CC) .GE. abs(S)) THEN
        B = T*(DD - S*(CC*FF))
      ELSE IF ( abs(CC) .GE. abs(FF)) THEN
        B = T*(DD - CC*(S*FF))
      ELSE
        B = T*(DD - FF*(S*CC))
      ENDIF
      IF (FLIP.NE.0 ) THEN
        B = -B
      ENDIF

      RETURN
      end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      subroutine ORDERSYS(DV,DM,IV,N)
!     orders the eigenvalue vector DV, the eigenvector matrix DM, and the
!     indication vector IV of an N dimensional system.
!
!     DV (in/out): N dimensional matrix of real*8
!                  contains the eigenvalues sorted by real and imaginary
!                  parts. Which part DV(I) actually is, is described in
!                  IV(I).
!
!     DM (in/out): N*N dimensional matrix of real*8
!                  contains the eigenvectors sorted by real and imaginary
!                  parts. Which part DM(:,I) actually is, is described in
!                  IV(I).
!
!     IV (in/out): N dimensional vector of integers describing the type
!                  of eigenvalues and eigenvectors.  The types are not
!                  ordered in any way.
!                  IV(I)=0 indicates that DV(I) is a real eigenvalue and
!                          DM(:,I) is the corresponding real eigenvector
!                  IV(I)=1 indicates that DV(I)+I*DV(I+1) and DV(I)-I*DV(I+1)
!                          are complex eigenvalues and DM(:,I)+I*DM(:,I+1),
!                          DM(:,I)-I*DM(:,I+1) are the corresponding complex
!                          eigenvectors. The imaginary part DV(I+1) is always
!                          chosen to be positive. IV(I)=1 is always followed
!                          by IV(I+1)=2.
!                  IV(I)=2 Indicates the imaginary part of eigenvalues and
!                          eigenvectors.  It is always preceedid by IV(I-1)=1.
!
!     N  (input) : N specifies the dimensions of matrix DM and the
!                  vectores DV, IV.  It has to be smaller than the local
!                  parameter NMX.

!=====job parameters
      real(rp), parameter :: TINY=1.D-99

!=====local stuff
      integer, parameter ::NMX=8,NMX2=4
      integer n
      integer n2, i, j, k, jm1, im, id, jm2, jd1, jd2 
      real(rp) dum

      real(rp) DV(N), DVN(NMX), DM(N,N), DN(NMX,NMX), &
          ORDER(NMX2,NMX)
      integer IV(N), IVN(NMX), JDM(NMX), IDM(NMX2)

!%%      write(*,*)'ordersys'

      IF(N.GT.NMX) CALL &
          out_io(s_fatal$, 'ORDERSYS', 'increase NMX')
      IF(MOD(N,2).NE.0) CALL &
          out_io(s_fatal$, 'ORDERSYS', 'eigensystem not even dimensional')

!=====get matrix ORDER which specifies the size of two succesive components
!=====for each eigenvector :
      N2 = N/2

      DO 10 J=1,N
         DO 110 I=1,N2
            ORDER(I,J) = DM(2*I-1,J)**2+DM(2*I,J)**2
 110     CONTINUE
 10   CONTINUE

      DO 20 I=1,N2
         IDM(I) = I
 20   CONTINUE

      DO 30 J=1,N
         JDM(J) = J
 30   CONTINUE

      DO 40 K=1,N2
         DO 140 J=2*K-1,N
            DUM = 0

            DO 1140 I=K,N2
               DUM = DUM + ORDER(IDM(I),JDM(J))
 1140       CONTINUE

            IF(DUM.LT.TINY) CALL out_io(s_fatal$, 'ORDERSYS', 'column of zeros')

            DO 2140 I=K,N2
               ORDER(IDM(I),JDM(J)) = ORDER(IDM(I),JDM(J))/DUM
 2140       CONTINUE

 140     CONTINUE

!========find maximum in row IM and column JM1 of ORDER :
         DUM = -1.D+38

         DO 240 I=K,N2
            DO 1240 J=2*K-1,N
               IF(ORDER(IDM(I),JDM(J)).GT.DUM) THEN
                  DUM = ORDER(IDM(I),JDM(J))
                  JM1 = J
                  IM  = I
               ENDIF
 1240       CONTINUE
 240     CONTINUE

         ID = IDM(IM)

         DO 340 I=IM,K+1,-1
            IDM(I) = IDM(I-1)
 340     CONTINUE

         IDM(K) = ID

         IF(IV(JDM(JM1)).EQ.2) THEN
            JM2 = JM1
            JM1 = JM1-1
         ELSEIF(IV(JDM(JM1)).EQ.1) THEN
            JM2 = JM1+1
         ELSEIF(IV(JDM(JM1)).EQ.0) THEN
!========find conjugated real column :
            JM2 = 0
            DUM = -1.D+38

            DO 440 J=K*2-1,N
               IF(IV(JDM(J)).EQ.0.AND.J.NE.JM1.AND. &
                   ORDER(IM,JDM(J)).GT.DUM) THEN
                  DUM = ORDER(IM,JDM(J))
                  JM2 = J
               ENDIF
 440        CONTINUE

            IF(JM2.EQ.0) &
                CALL out_io(s_fatal$, 'ORDERSYS', 'odd number of real columns')
!===========place conjugate real columns together :
            J   = MAX(JM1,JM2)
            JM1 = MIN(JM1,JM2)
            JM2 = J

            DO 540 J=JM2,JM1+2,-1
               JDM(J) = JDM(J-1)
 540        CONTINUE

            JM2 = JM1 + 1

            IF(abs(DV(JDM(JM1))).LT.abs(DV(JDM(JM2)))) THEN
               JD1 = JDM(JM2)
               JDM(JM2) = JDM(JM1)
               JDM(JM1) = JD1
            ENDIF
         ENDIF

!========place pair of columns to the front :
         JD1 = JDM(JM1)
         JD2 = JDM(JM2)

         DO 640 J=JM2,2*K+1,-1
            JDM(J) = JDM(J-2)
 640     CONTINUE

         JDM(2*K-1) = JD1
         JDM(2*K)   = JD2
 40   CONTINUE

!=====reorder DV, DM, and IV  ::
!=====start with local copies :
      DO 60 K=1,N2

         DO 160 I=1,N
            DN(I,2*IDM(K)-1) = DM(I,JDM(2*K-1))
            DN(I,2*IDM(K)  ) = DM(I,JDM(2*K  ))
 160     CONTINUE

         IVN(2*IDM(K)-1) = IV(JDM(2*K-1))
         IVN(2*IDM(K)  ) = IV(JDM(2*K  ))
         DVN(2*IDM(K)-1) = DV(JDM(2*K-1))
         DVN(2*IDM(K)  ) = DV(JDM(2*K  ))
 60   CONTINUE

!=====write back originals :
      DO 80 J=1,N

         DO 180 I=1,N
            DM(I,J) = DN(I,J)
 180     CONTINUE

         IV(J) = IVN(J)
         DV(J) = DVN(J)
 80   CONTINUE

      RETURN
      end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      subroutine NORM66(DV,DM,IV,N)
!     normalizes the N columns of the N*N dimensional matrix DM
!     in order to let it be symplectic, if possible, and stores it in DV.
!     should be called after ORDERSYS.

!=====job parameters
      real(rp), parameter :: TINY=1.D-99

!=====local stuff
      integer n
      integer j, k
      real(rp) DV(N), DM(N,N), dum
      integer IV(N)
!%%      write(*,*)'norm66'

      IF(MOD(N,2).NE.0) CALL &
          out_io(s_fatal$, 'NORM66', 'eigensystem not even dimensional')

      DO 10 J=1,N,2
         DUM = 0.D+0

         DO 110 K=1,N,2
            DUM = DUM + DM(K,J)*DM(K+1,J+1) - DM(K+1,J)*DM(K,J+1)
 110     CONTINUE

         IF(DUM.GE.0.D+0) THEN
            DUM = SQRT(DUM)
            IF(DUM.LT.TINY) CALL out_io(s_fatal$, 'NORM66', 'non symplectic')
            DM(1:n,J+1) = DM(1:n,J+1) / DUM
         ELSE
!=====take complex conjugate eigenvectore and value
            DUM = SQRT(-DUM)
            IF(DUM.LT.TINY) CALL out_io(s_fatal$, 'NORM66', 'non symplectic')
            DM(1:n,J+1) = -DM(1:n,J+1) / DUM
            IF(IV(J+1).EQ.2) DV(J+1) = -DV(J+1)
         ENDIF

         DM(1:n,J) = DM(1:n,J) / DUM
 10   CONTINUE

      RETURN
      end subroutine

end module
