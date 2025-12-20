!** DPCHFD
PURE SUBROUTINE DPCHFD(N,X,F,D,Incfd,Skip,Ne,Xe,Fe,De,Ierr)
  !> Evaluate a piecewise cubic Hermite function and its first
  !  derivative at an array of points.  May be used by itself
  !  for Hermite interpolation, or as an evaluator for DPCHIM or DPCHIC.
  !  If only function values are required, use DPCHFE instead.
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Category:**  E3, H1
  !***
  ! **Type:**      DOUBLE PRECISION (PCHFD-S, DPCHFD-D)
  !***
  ! **Keywords:**  CUBIC HERMITE DIFFERENTIATION, CUBIC HERMITE EVALUATION,
  !             HERMITE INTERPOLATION, PCHIP, PIECEWISE CUBIC EVALUATION
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             P.O. Box 808  (L-316)
  !             Livermore, CA  94550
  !             FTS 532-4275, (510) 422-4275
  !***
  ! **Description:**
  !
  !          DPCHFD:  Piecewise Cubic Hermite Function and Derivative evaluator
  !
  !     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
  !     gether with its first derivative, at the points  XE(J), J=1(1)NE.
  !
  !***
  ! **References:**  Fritsch, F.N. and Carlson, R.E., "Monotone Piecewise
  !                 Cubic Interpolation", SIAM J. Numer. Anal. 17(2), 1980
  !***
  ! **Revision history:**
  !   811020  DATE WRITTEN (Fritsch, LLNL)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251220  Eliminated GOTOs using Jacob Williams pattern. (ZH)
  !           Ref: github.com/jacobwilliams/PCHIP
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER, INTENT(IN) :: N, Incfd, Ne
  INTEGER, INTENT(OUT) :: Ierr
  REAL(DP), INTENT(IN) :: X(N), F(Incfd,N), D(Incfd,N), Xe(Ne)
  REAL(DP), INTENT(OUT) :: Fe(Ne), De(Ne)
  LOGICAL, INTENT(INOUT) :: Skip
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER :: i, ierc, ir, j, jfirst, next(2), nj
  LOGICAL :: located, error
  !
  !  VALIDITY-CHECK ARGUMENTS.
  !
  IF( .NOT. Skip ) THEN
    IF( N < 2 ) THEN
      ERROR STOP 'DPCHFD : NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN
    ELSEIF( Incfd < 1 ) THEN
      ERROR STOP 'DPCHFD : INCREMENT LESS THAN ONE'
      RETURN
    ELSE
      DO i = 2, N
        IF( X(i) <= X(i-1) ) EXIT
      END DO
      IF( i <= N ) THEN
        ERROR STOP 'DPCHFD : X-ARRAY NOT STRICTLY INCREASING'
        RETURN
      END IF
    END IF
  END IF
  !
  IF( Ne < 1 ) THEN
    ERROR STOP 'DPCHFD : NUMBER OF EVALUATION POINTS LESS THAN ONE'
    RETURN
  END IF
  !
  Ierr = 0
  Skip = .TRUE.
  !
  !  LOOP OVER INTERVALS.  (Interval index is IL = IR-1)
  !                        (Interval is X(IL) <= X < X(IR))
  !
  jfirst = 1
  ir = 2
  error = .FALSE.
  !
  main_loop: do
    !
    !  Exit when all evaluation points processed.
    !
    if (jfirst > Ne) return
    !
    !  Locate all points in current interval.
    !
    located = .FALSE.
    do j = jfirst, Ne
      if (Xe(j) >= X(ir)) then
        located = .TRUE.
        exit
      end if
    end do
    !
    !  Determine number of points in this interval.
    !
    if (located) then
      if (ir == N) j = Ne + 1  ! Last interval - extrapolation points
    else
      j = Ne + 1  ! All remaining points in current interval
    end if
    !
    nj = j - jfirst
    !
    !  Evaluate cubic if there are points in this interval.
    !
    if (nj /= 0) then
      !
      call DCHFDV(X(ir-1), X(ir), F(1,ir-1), F(1,ir), D(1,ir-1), D(1,ir), &
                  nj, Xe(jfirst), Fe(jfirst), De(jfirst), next, ierc)
      !
      if (ierc < 0) then
        error = .TRUE.
        exit main_loop
      end if
      !
      !  Handle points to the right of X(IR).
      !
      if (next(2) /= 0) then
        if (ir < N) then
          error = .TRUE.
          exit main_loop
        end if
        Ierr = Ierr + next(2)  ! Extrapolation points
      end if
      !
      !  Handle points to the left of X(IR-1).
      !
      if (next(1) /= 0) then
        if (ir > 2) then
          !  XE is not ordered relative to X, must adjust interval.
          !  First, locate first point to left of X(IR-1).
          located = .FALSE.
          do i = jfirst, j - 1
            if (Xe(i) < X(ir-1)) then
              located = .TRUE.
              exit
            end if
          end do
          !
          if (.NOT. located) then
            !  NOTE: Cannot get here unless error in DCHFDV.
            error = .TRUE.
            exit main_loop
          end if
          !
          !  Reset J (this will be the new JFIRST).
          j = i
          !
          !  Find how far to back up in the X-array.
          do i = 1, ir - 1
            if (Xe(j) < X(i)) exit
          end do
          !  At this point: XE(J) < X(1) or X(I-1) <= XE(J) < X(I)
          !  Reset IR (will be incremented before cycling).
          ir = MAX(1, i - 1)
        else
          !  IR == 2: these are extrapolation points at left.
          Ierr = Ierr + next(1)
        end if
      end if
      !
      jfirst = j
    end if
    !
    !  Move to next interval.
    !
    ir = ir + 1
    if (ir > N) return
    !
  end do main_loop
  !
  !  Error exit.
  !
  if (error) then
    ERROR STOP 'DPCHFD : ERROR RETURN FROM DCHFDV -- FATAL'
  end if
  !
END SUBROUTINE DPCHFD
