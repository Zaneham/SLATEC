!** DSORT
PURE SUBROUTINE DSORT(Dx,Dy,N,Kflag)
  !> Sort an array and optionally make the same interchanges in an auxiliary array.
  !  The array may be sorted in increasing or decreasing order.
  !  A slightly modified QUICKSORT algorithm is used.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N6A2B
  !***
  ! **Type:**      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
  !***
  ! **Keywords:**  SINGLETON QUICKSORT, SORT, SORTING
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !   DSORT sorts array DX and optionally makes the same interchanges in
  !   array DY.  The array DX may be sorted in increasing order or
  !   decreasing order.  A slightly modified quicksort algorithm is used.
  !
  !   Description of Parameters
  !      DX - array of values to be sorted   (usually abscissas)
  !      DY - array to be (optionally) carried along
  !      N  - number of values in array DX to be sorted
  !      KFLAG - control parameter
  !            =  2  means sort DX in increasing order and carry DY along.
  !            =  1  means sort DX in increasing order (ignoring DY)
  !            = -1  means sort DX in decreasing order (ignoring DY)
  !            = -2  means sort DX in decreasing order and carry DY along.
  !
  !***
  ! **References:**  R. C. Singleton, Algorithm 347, An efficient algorithm
  !                 for sorting with minimal storage, Communications of
  !                 the ACM, 12, 3 (1969), pp. 185-187.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   761101  DATE WRITTEN
  !   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
  !   2024xx  Eliminated GOTOs using state machine pattern. (ZH)

  INTEGER, INTENT(IN) :: Kflag, N
  REAL(DP), INTENT(INOUT) :: Dx(N), Dy(N)
  REAL(DP) :: r, t, tt, tty, ty
  INTEGER :: i, ij, j, k, kk, l, m, nn
  INTEGER :: il(21), iu(21)
  INTEGER, PARAMETER :: ST_MAIN = 1, ST_PIVOT = 2, ST_POP = 3, ST_DECIDE = 4
  INTEGER :: istate
  INTRINSIC ABS, INT

  nn = N
  IF( nn<1 ) ERROR STOP 'DSORT : The number of values to be sorted is not positive.'
  kk = ABS(Kflag)
  IF( kk/=1 .AND. kk/=2 ) THEN
    ERROR STOP 'DSORT : The sort control parameter, K, is not 2, 1, -1, or -2.'
    RETURN
  END IF

  IF( Kflag<=-1 ) THEN
    DO i = 1, nn
      Dx(i) = -Dx(i)
    END DO
  END IF

  m = 1; i = 1; j = nn; r = 0.375_DP

  IF( kk==2 ) THEN
    istate = ST_MAIN
    quicksort_xy: DO
      SELECT CASE (istate)
      CASE (ST_MAIN)
        IF( i==j ) THEN; istate = ST_POP; CYCLE quicksort_xy; END IF
        IF( r<=0.5898437_DP ) THEN; r = r + 3.90625E-2_DP
        ELSE; r = r - 0.21875_DP; END IF
        istate = ST_PIVOT
      CASE (ST_PIVOT)
        k = i; ij = i + INT((j-i)*r); t = Dx(ij); ty = Dy(ij)
        IF( Dx(i)>t ) THEN
          Dx(ij) = Dx(i); Dx(i) = t; t = Dx(ij)
          Dy(ij) = Dy(i); Dy(i) = ty; ty = Dy(ij)
        END IF
        l = j
        IF( Dx(j)<t ) THEN
          Dx(ij) = Dx(j); Dx(j) = t; t = Dx(ij)
          Dy(ij) = Dy(j); Dy(j) = ty; ty = Dy(ij)
          IF( Dx(i)>t ) THEN
            Dx(ij) = Dx(i); Dx(i) = t; t = Dx(ij)
            Dy(ij) = Dy(i); Dy(i) = ty; ty = Dy(ij)
          END IF
        END IF
        partition_xy: DO
          l = l - 1
          IF( Dx(l)<=t ) THEN
            inner_xy: DO
              k = k + 1
              IF( Dx(k)>=t ) THEN
                IF( k<=l ) THEN
                  tt = Dx(l); Dx(l) = Dx(k); Dx(k) = tt
                  tty = Dy(l); Dy(l) = Dy(k); Dy(k) = tty
                  EXIT inner_xy
                END IF
                IF( l-i>j-k ) THEN
                  il(m) = i; iu(m) = l; i = k; m = m + 1
                ELSE
                  il(m) = k; iu(m) = j; j = l; m = m + 1
                END IF
                istate = ST_DECIDE; EXIT partition_xy
              END IF
            END DO inner_xy
          END IF
        END DO partition_xy
        IF( istate==ST_PIVOT ) istate = ST_DECIDE
      CASE (ST_POP)
        m = m - 1
        IF( m==0 ) EXIT quicksort_xy
        i = il(m); j = iu(m); istate = ST_DECIDE
      CASE (ST_DECIDE)
        IF( j-i>=1 ) THEN; istate = ST_PIVOT
        ELSEIF( i==1 ) THEN; istate = ST_MAIN
        ELSE
          i = i - 1
          insertion_xy: DO
            i = i + 1
            IF( i==j ) THEN; istate = ST_POP; EXIT insertion_xy; END IF
            t = Dx(i+1); ty = Dy(i+1)
            IF( Dx(i)>t ) THEN
              k = i
              insert_inner_xy: DO
                Dx(k+1) = Dx(k); Dy(k+1) = Dy(k); k = k - 1
                IF( t>=Dx(k) ) THEN; Dx(k+1) = t; Dy(k+1) = ty; EXIT insert_inner_xy; END IF
              END DO insert_inner_xy
            END IF
          END DO insertion_xy
        END IF
      END SELECT
    END DO quicksort_xy
  ELSE
    istate = ST_MAIN
    quicksort_x: DO
      SELECT CASE (istate)
      CASE (ST_MAIN)
        IF( i==j ) THEN; istate = ST_POP; CYCLE quicksort_x; END IF
        IF( r<=0.5898437_DP ) THEN; r = r + 3.90625E-2_DP
        ELSE; r = r - 0.21875_DP; END IF
        istate = ST_PIVOT
      CASE (ST_PIVOT)
        k = i; ij = i + INT((j-i)*r); t = Dx(ij)
        IF( Dx(i)>t ) THEN; Dx(ij) = Dx(i); Dx(i) = t; t = Dx(ij); END IF
        l = j
        IF( Dx(j)<t ) THEN
          Dx(ij) = Dx(j); Dx(j) = t; t = Dx(ij)
          IF( Dx(i)>t ) THEN; Dx(ij) = Dx(i); Dx(i) = t; t = Dx(ij); END IF
        END IF
        partition_x: DO
          l = l - 1
          IF( Dx(l)<=t ) THEN
            inner_x: DO
              k = k + 1
              IF( Dx(k)>=t ) THEN
                IF( k<=l ) THEN; tt = Dx(l); Dx(l) = Dx(k); Dx(k) = tt; EXIT inner_x; END IF
                IF( l-i>j-k ) THEN
                  il(m) = i; iu(m) = l; i = k; m = m + 1
                ELSE
                  il(m) = k; iu(m) = j; j = l; m = m + 1
                END IF
                istate = ST_DECIDE; EXIT partition_x
              END IF
            END DO inner_x
          END IF
        END DO partition_x
        IF( istate==ST_PIVOT ) istate = ST_DECIDE
      CASE (ST_POP)
        m = m - 1
        IF( m==0 ) EXIT quicksort_x
        i = il(m); j = iu(m); istate = ST_DECIDE
      CASE (ST_DECIDE)
        IF( j-i>=1 ) THEN; istate = ST_PIVOT
        ELSEIF( i==1 ) THEN; istate = ST_MAIN
        ELSE
          i = i - 1
          insertion_x: DO
            i = i + 1
            IF( i==j ) THEN; istate = ST_POP; EXIT insertion_x; END IF
            t = Dx(i+1)
            IF( Dx(i)>t ) THEN
              k = i
              insert_inner_x: DO
                Dx(k+1) = Dx(k); k = k - 1
                IF( t>=Dx(k) ) THEN; Dx(k+1) = t; EXIT insert_inner_x; END IF
              END DO insert_inner_x
            END IF
          END DO insertion_x
        END IF
      END SELECT
    END DO quicksort_x
  END IF

  IF( Kflag<=-1 ) THEN
    DO i = 1, nn
      Dx(i) = -Dx(i)
    END DO
  END IF

END SUBROUTINE DSORT
