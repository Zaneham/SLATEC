!** PCHFE
PURE SUBROUTINE PCHFE(N,X,F,D,Incfd,Skip,Ne,Xe,Fe,Ierr)
  !> Evaluate a piecewise cubic Hermite function at an array of points.
  !  May be used by itself for Hermite interpolation, or as an evaluator for PCHIM or PCHIC.
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Category:**  E3
  !***
  ! **Type:**      SINGLE PRECISION (PCHFE-S, PCHFE-S)
  !***
  ! **Keywords:**  CUBIC HERMITE EVALUATION, HERMITE INTERPOLATION, PCHIP,
  !             PIECEWISE CUBIC EVALUATION
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             P.O. Box 808  (L-316)
  !             Livermore, CA  94550
  !             FTS 532-4275, (510) 422-4275
  !***
  ! **Description:**
  !
  !          PCHFE:  Piecewise Cubic Hermite Function Evaluator
  !
  !     Evaluates the cubic Hermite function defined by  N, X, F, D  at
  !     the points  XE(J), J=1(1)NE.
  !
  !     To provide compatibility with PCHIM and PCHIC, includes an
  !     increment between successive values of the F- and D-arrays.
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER, INTENT(IN) :: N, Incfd, Ne
  INTEGER, INTENT(OUT) :: Ierr
  REAL(SP), INTENT(IN) :: X(N), F(Incfd,N), D(Incfd,N), Xe(Ne)
  REAL(SP), INTENT(OUT) :: Fe(Ne)
  LOGICAL, INTENT(INOUT) :: Skip
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER :: i, ierc, ir, j, jfirst, next(2), nj
  LOGICAL :: found_left_point
  !
  !  VALIDITY-CHECK ARGUMENTS.
  !
  !* FIRST EXECUTABLE STATEMENT  PCHFE
  IF( .NOT. Skip ) THEN
    IF( N < 2 ) THEN
      ERROR STOP 'PCHFE : NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN
    ELSEIF( Incfd < 1 ) THEN
      ERROR STOP 'PCHFE : INCREMENT LESS THAN ONE'
      RETURN
    ELSE
      DO i = 2, N
        IF( X(i) <= X(i-1) ) EXIT
      END DO
      IF( i <= N ) THEN
        ERROR STOP 'PCHFE : X-ARRAY NOT STRICTLY INCREASING'
        RETURN
      END IF
    END IF
  END IF
  !
  !  FUNCTION DEFINITION IS OK, GO ON.
  !
  IF( Ne < 1 ) THEN
    ERROR STOP 'PCHFE : NUMBER OF EVALUATION POINTS LESS THAN ONE'
    RETURN
  END IF

  Ierr = 0
  Skip = .TRUE.
  !
  !  LOOP OVER INTERVALS.
  !
  jfirst = 1
  ir = 2
  !
  ir_loop: DO WHILE( ir <= N )
    !
    IF( jfirst > Ne ) EXIT ir_loop
    !
    !  LOCATE ALL POINTS IN INTERVAL.
    !
    DO j = jfirst, Ne
      IF( Xe(j) >= X(ir) ) EXIT
    END DO
    !
    IF( j > Ne ) THEN
      j = Ne + 1
    ELSEIF( ir == N ) THEN
      j = Ne + 1
    END IF
    !
    nj = j - jfirst
    !
    IF( nj /= 0 ) THEN
      !
      CALL CHFEV(X(ir-1),X(ir),F(1,ir-1),F(1,ir),D(1,ir-1),D(1,ir),nj,&
        Xe(jfirst),Fe(jfirst),next,ierc)
      !
      IF( ierc < 0 ) THEN
        ERROR STOP 'PCHFE : ERROR RETURN FROM CHFEV -- FATAL'
        RETURN
      END IF
      !
      IF( next(2) /= 0 ) THEN
        IF( ir < N ) THEN
          ERROR STOP 'PCHFE : ERROR RETURN FROM CHFEV -- FATAL'
          RETURN
        END IF
        Ierr = Ierr + next(2)
      END IF
      !
      IF( next(1) /= 0 ) THEN
        IF( ir > 2 ) THEN
          found_left_point = .FALSE.
          DO i = jfirst, j - 1
            IF( Xe(i) < X(ir-1) ) THEN
              found_left_point = .TRUE.
              EXIT
            END IF
          END DO
          IF( .NOT. found_left_point ) THEN
            ERROR STOP 'PCHFE : ERROR RETURN FROM CHFEV -- FATAL'
            RETURN
          END IF
          j = i
          DO i = 1, ir - 1
            IF( Xe(j) < X(i) ) EXIT
          END DO
          ir = MAX(1, i-1)
        ELSE
          Ierr = Ierr + next(1)
        END IF
      END IF
      !
      jfirst = j
    END IF
    !
    ir = ir + 1
  END DO ir_loop
  !
  RETURN
END SUBROUTINE PCHFE
