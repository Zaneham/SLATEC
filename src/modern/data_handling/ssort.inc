!** SSORT
PURE SUBROUTINE SSORT(X,Y,N,Kflag)
  !> Sort an array and optionally make the same interchanges in an auxiliary array.
  !  The array may be sorted in increasing or decreasing order.
  !  A slightly modified QUICKSORT algorithm is used.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N6A2B
  !***
  ! **Type:**      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
  !***
  ! **Keywords:**  SINGLETON QUICKSORT, SORT, SORTING
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !   SSORT sorts array X and optionally makes the same interchanges in
  !   array Y.  The array X may be sorted in increasing order or
  !   decreasing order.  A slightly modified quicksort algorithm is used.
  !
  !   Description of Parameters
  !      X - array of values to be sorted   (usually abscissas)
  !      Y - array to be (optionally) carried along
  !      N - number of values in array X to be sorted
  !      KFLAG - control parameter
  !            =  2  means sort X in increasing order and carry Y along.
  !            =  1  means sort X in increasing order (ignoring Y)
  !            = -1  means sort X in decreasing order (ignoring Y)
  !            = -2  means sort X in decreasing order and carry Y along.
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
  REAL(SP), INTENT(INOUT) :: X(N), Y(N)
  REAL(SP) :: r, t, tt, tty, ty
  INTEGER :: i, ij, j, k, kk, l, m, nn
  INTEGER :: il(21), iu(21)
  INTEGER, PARAMETER :: ST_MAIN = 1, ST_PIVOT = 2, ST_POP = 3, ST_DECIDE = 4
  INTEGER :: istate
  INTRINSIC ABS, INT

  nn = N
  IF( nn<1 ) ERROR STOP 'SSORT : The number of values to be sorted is not positive.'
  kk = ABS(Kflag)
  IF( kk/=1 .AND. kk/=2 ) ERROR STOP 'SSORT : The sort control parameter, K, is not 2, 1, -1, or -2.'

  IF( Kflag<=-1 ) THEN
    DO i = 1, nn
      X(i) = -X(i)
    END DO
  END IF

  m = 1; i = 1; j = nn; r = 0.375_SP

  IF( kk==2 ) THEN
    istate = ST_MAIN
    quicksort_xy: DO
      SELECT CASE (istate)
      CASE (ST_MAIN)
        IF( i==j ) THEN; istate = ST_POP; CYCLE quicksort_xy; END IF
        IF( r<=0.5898437_SP ) THEN; r = r + 3.90625E-2_SP
        ELSE; r = r - 0.21875_SP; END IF
        istate = ST_PIVOT
      CASE (ST_PIVOT)
        k = i; ij = i + INT((j-i)*r); t = X(ij); ty = Y(ij)
        IF( X(i)>t ) THEN
          X(ij) = X(i); X(i) = t; t = X(ij)
          Y(ij) = Y(i); Y(i) = ty; ty = Y(ij)
        END IF
        l = j
        IF( X(j)<t ) THEN
          X(ij) = X(j); X(j) = t; t = X(ij)
          Y(ij) = Y(j); Y(j) = ty; ty = Y(ij)
          IF( X(i)>t ) THEN
            X(ij) = X(i); X(i) = t; t = X(ij)
            Y(ij) = Y(i); Y(i) = ty; ty = Y(ij)
          END IF
        END IF
        partition_xy: DO
          l = l - 1
          IF( X(l)<=t ) THEN
            inner_xy: DO
              k = k + 1
              IF( X(k)>=t ) THEN
                IF( k<=l ) THEN
                  tt = X(l); X(l) = X(k); X(k) = tt
                  tty = Y(l); Y(l) = Y(k); Y(k) = tty
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
            t = X(i+1); ty = Y(i+1)
            IF( X(i)>t ) THEN
              k = i
              insert_inner_xy: DO
                X(k+1) = X(k); Y(k+1) = Y(k); k = k - 1
                IF( t>=X(k) ) THEN; X(k+1) = t; Y(k+1) = ty; EXIT insert_inner_xy; END IF
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
        IF( r<=0.5898437_SP ) THEN; r = r + 3.90625E-2_SP
        ELSE; r = r - 0.21875_SP; END IF
        istate = ST_PIVOT
      CASE (ST_PIVOT)
        k = i; ij = i + INT((j-i)*r); t = X(ij)
        IF( X(i)>t ) THEN; X(ij) = X(i); X(i) = t; t = X(ij); END IF
        l = j
        IF( X(j)<t ) THEN
          X(ij) = X(j); X(j) = t; t = X(ij)
          IF( X(i)>t ) THEN; X(ij) = X(i); X(i) = t; t = X(ij); END IF
        END IF
        partition_x: DO
          l = l - 1
          IF( X(l)<=t ) THEN
            inner_x: DO
              k = k + 1
              IF( X(k)>=t ) THEN
                IF( k<=l ) THEN; tt = X(l); X(l) = X(k); X(k) = tt; EXIT inner_x; END IF
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
            t = X(i+1)
            IF( X(i)>t ) THEN
              k = i
              insert_inner_x: DO
                X(k+1) = X(k); k = k - 1
                IF( t>=X(k) ) THEN; X(k+1) = t; EXIT insert_inner_x; END IF
              END DO insert_inner_x
            END IF
          END DO insertion_x
        END IF
      END SELECT
    END DO quicksort_x
  END IF

  IF( Kflag<=-1 ) THEN
    DO i = 1, nn
      X(i) = -X(i)
    END DO
  END IF

END SUBROUTINE SSORT
