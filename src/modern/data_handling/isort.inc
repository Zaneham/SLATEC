!** ISORT
PURE SUBROUTINE ISORT(Ix,Iy,N,Kflag)
  !> Sort an array and optionally make the same interchanges in an auxiliary array.
  !  The array may be sorted in increasing or decreasing order.
  !  A slightly modified QUICKSORT algorithm is used.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N6A2A
  !***
  ! **Type:**      INTEGER (SSORT-S, DSORT-D, ISORT-I)
  !***
  ! **Keywords:**  SINGLETON QUICKSORT, SORT, SORTING
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !           Kahaner, D. K., (NBS)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !   ISORT sorts array IX and optionally makes the same interchanges in
  !   array IY.  The array IX may be sorted in increasing order or
  !   decreasing order.  A slightly modified quicksort algorithm is used.
  !
  !   Description of Parameters
  !      IX - integer array of values to be sorted
  !      IY - integer array to be (optionally) carried along
  !      N  - number of values in integer array IX to be sorted
  !      KFLAG - control parameter
  !            =  2  means sort IX in increasing order and carry IY along.
  !            =  1  means sort IX in increasing order (ignoring IY)
  !            = -1  means sort IX in decreasing order (ignoring IY)
  !            = -2  means sort IX in decreasing order and carry IY along.
  !
  !***
  ! **References:**  R. C. Singleton, Algorithm 347, An efficient algorithm
  !                 for sorting with minimal storage, Communications of
  !                 the ACM, 12, 3 (1969), pp. 185-187.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   761118  DATE WRITTEN
  !   810801  Modified by David K. Kahaner.
  !   2024xx  Eliminated GOTOs using state machine pattern. (ZH)

  INTEGER, INTENT(IN) :: Kflag, N
  INTEGER, INTENT(INOUT) :: Ix(N), Iy(N)
  REAL(SP) :: r
  INTEGER :: i, ij, j, k, kk, l, m, nn, t, tt, tty, ty
  INTEGER :: il(21), iu(21)
  INTEGER, PARAMETER :: ST_MAIN = 1, ST_PIVOT = 2, ST_POP = 3, ST_DECIDE = 4
  INTEGER :: istate
  INTRINSIC ABS, INT

  nn = N
  IF( nn<1 ) ERROR STOP 'ISORT : The number of values to be sorted is not positive.'
  kk = ABS(Kflag)
  IF( kk/=1 .AND. kk/=2 ) ERROR STOP 'ISORT : The sort control parameter, K, is not 2, 1, -1, or -2.'

  IF( Kflag<=-1 ) THEN
    DO i = 1, nn
      Ix(i) = -Ix(i)
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
        k = i; ij = i + INT((j-i)*r); t = Ix(ij); ty = Iy(ij)
        IF( Ix(i)>t ) THEN
          Ix(ij) = Ix(i); Ix(i) = t; t = Ix(ij)
          Iy(ij) = Iy(i); Iy(i) = ty; ty = Iy(ij)
        END IF
        l = j
        IF( Ix(j)<t ) THEN
          Ix(ij) = Ix(j); Ix(j) = t; t = Ix(ij)
          Iy(ij) = Iy(j); Iy(j) = ty; ty = Iy(ij)
          IF( Ix(i)>t ) THEN
            Ix(ij) = Ix(i); Ix(i) = t; t = Ix(ij)
            Iy(ij) = Iy(i); Iy(i) = ty; ty = Iy(ij)
          END IF
        END IF
        partition_xy: DO
          l = l - 1
          IF( Ix(l)<=t ) THEN
            inner_xy: DO
              k = k + 1
              IF( Ix(k)>=t ) THEN
                IF( k<=l ) THEN
                  tt = Ix(l); Ix(l) = Ix(k); Ix(k) = tt
                  tty = Iy(l); Iy(l) = Iy(k); Iy(k) = tty
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
            t = Ix(i+1); ty = Iy(i+1)
            IF( Ix(i)>t ) THEN
              k = i
              insert_inner_xy: DO
                Ix(k+1) = Ix(k); Iy(k+1) = Iy(k); k = k - 1
                IF( t>=Ix(k) ) THEN; Ix(k+1) = t; Iy(k+1) = ty; EXIT insert_inner_xy; END IF
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
        k = i; ij = i + INT((j-i)*r); t = Ix(ij)
        IF( Ix(i)>t ) THEN; Ix(ij) = Ix(i); Ix(i) = t; t = Ix(ij); END IF
        l = j
        IF( Ix(j)<t ) THEN
          Ix(ij) = Ix(j); Ix(j) = t; t = Ix(ij)
          IF( Ix(i)>t ) THEN; Ix(ij) = Ix(i); Ix(i) = t; t = Ix(ij); END IF
        END IF
        partition_x: DO
          l = l - 1
          IF( Ix(l)<=t ) THEN
            inner_x: DO
              k = k + 1
              IF( Ix(k)>=t ) THEN
                IF( k<=l ) THEN; tt = Ix(l); Ix(l) = Ix(k); Ix(k) = tt; EXIT inner_x; END IF
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
            t = Ix(i+1)
            IF( Ix(i)>t ) THEN
              k = i
              insert_inner_x: DO
                Ix(k+1) = Ix(k); k = k - 1
                IF( t>=Ix(k) ) THEN; Ix(k+1) = t; EXIT insert_inner_x; END IF
              END DO insert_inner_x
            END IF
          END DO insertion_x
        END IF
      END SELECT
    END DO quicksort_x
  END IF

  IF( Kflag<=-1 ) THEN
    DO i = 1, nn
      Ix(i) = -Ix(i)
    END DO
  END IF

END SUBROUTINE ISORT
