!** PROCP
PURE SUBROUTINE PROCP(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Y,M,A,B,C,D,U,W)
  !> Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      COMPLEX (PRODP-C, PROCP-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  ! PROCP applies a sequence of matrix operations to the vector X and
  ! stores the result in Y (periodic boundary conditions).
  !
  ! BD,BM1,BM2 are arrays containing roots of certain B polynomials.
  ! ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
  ! AA         Array containing scalar multipliers of the vector X.
  ! NA         is the length of the array AA.
  ! X,Y        The matrix operations are applied to X and the result is Y.
  ! A,B,C      are arrays which contain the tridiagonal matrix.
  ! M          is the order of the matrix.
  ! D,U,W      are working arrays.
  ! IS         determines whether or not a change in sign is made.
  !
  !***
  ! **See also:**  CBLKTR
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   202412  Eliminated GOTOs using nested DO loops with CYCLE/EXIT.  (modernization)

  INTEGER, INTENT(IN) :: M, Na, Nd, Nm1, Nm2
  REAL(SP), INTENT(IN) :: Aa(Na), Bd(Nd), Bm1(Nm1), Bm2(Nm2)
  COMPLEX(SP), INTENT(IN) :: X(M), A(M), B(M), C(M)
  COMPLEX(SP), INTENT(OUT) :: Y(M), D(M), W(M), U(M)
  !
  INTEGER :: ia, ibr, id, j, k, m1, m2, mm, mm2
  REAL(SP) :: rt
  COMPLEX(SP) :: den, ym, v, bh, am
  !* FIRST EXECUTABLE STATEMENT  PROCP
  DO j = 1, M
    Y(j) = X(j)
    W(j) = Y(j)
  END DO
  mm = M - 1
  mm2 = M - 2
  id = Nd
  ibr = 0
  m1 = Nm1
  m2 = Nm2
  ia = Na
  ! Main iteration loop (replaces labels 100-500)
  main_loop: DO
  inner_loop: DO
    IF( ia>0 ) THEN
      rt = Aa(ia)
      IF( Nd==0 ) rt = -rt
      ia = ia - 1
      DO j = 1, M
        Y(j) = rt*W(j)
      END DO
    END IF
    IF( id<=0 ) RETURN
    rt = Bd(id)
    id = id - 1
    IF( id==0 ) ibr = 1
    !
    ! BEGIN SOLUTION TO SYSTEM
    !
    bh = B(M) - rt
    ym = Y(M)
    den = B(1) - rt
    D(1) = C(1)/den
    U(1) = A(1)/den
    W(1) = Y(1)/den
    v = C(M)
    IF( mm2>=2 ) THEN
      DO j = 2, mm2
        den = B(j) - rt - A(j)*D(j-1)
        D(j) = C(j)/den
        U(j) = -A(j)*U(j-1)/den
        W(j) = (Y(j)-A(j)*W(j-1))/den
        bh = bh - v*U(j-1)
        ym = ym - v*W(j-1)
        v = -v*D(j-1)
      END DO
    END IF
    den = B(M-1) - rt - A(M-1)*D(M-2)
    D(M-1) = (C(M-1)-A(M-1)*U(M-2))/den
    W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/den
    am = A(M) - v*D(M-2)
    bh = bh - v*U(M-2)
    ym = ym - v*W(M-2)
    den = bh - am*D(M-1)
    IF( ABS(den)/=0 ) THEN
      W(M) = (ym-am*W(M-1))/den
    ELSE
      W(M) = (1._SP,0._SP)
    END IF
    W(M-1) = W(M-1) - D(M-1)*W(M)
    DO j = 2, mm
      k = M - j
      W(k) = W(k) - D(k)*W(k+1) - U(k)*W(M)
    END DO
    IF( Na<=0 ) THEN
      IF( m1<=0 ) THEN
        IF( m2>0 ) THEN
          ! Was GOTO 400: subtract Bm2
          IF( ibr<=0 ) THEN
            IF( ABS(Bm2(m2)-Bd(id))<ABS(Bm2(m2)-rt) ) THEN
              Y(1:M) = W(1:M)
              ibr = 1
              CYCLE main_loop
            END IF
          END IF
          rt = rt - Bm2(m2)
          m2 = m2 - 1
          Y(1:M) = Y(1:M) + rt*W(1:M)
          CYCLE main_loop
        END IF
        EXIT inner_loop
      ELSE
        IF( m2<=0 .OR. ABS(Bm1(m1))>ABS(Bm2(m2)) ) THEN
          ! Was GOTO 300: subtract Bm1
          IF( ibr<=0 ) THEN
            IF( ABS(Bm1(m1)-Bd(id))<ABS(Bm1(m1)-rt) ) THEN
              Y(1:M) = W(1:M)
              ibr = 1
              CYCLE main_loop
            END IF
          END IF
          rt = rt - Bm1(m1)
          m1 = m1 - 1
          Y(1:M) = Y(1:M) + rt*W(1:M)
          CYCLE main_loop
        ELSE
          ! Was GOTO 400: subtract Bm2
          IF( ibr<=0 ) THEN
            IF( ABS(Bm2(m2)-Bd(id))<ABS(Bm2(m2)-rt) ) THEN
              Y(1:M) = W(1:M)
              ibr = 1
              CYCLE main_loop
            END IF
          END IF
          rt = rt - Bm2(m2)
          m2 = m2 - 1
          Y(1:M) = Y(1:M) + rt*W(1:M)
          CYCLE main_loop
        END IF
      END IF
    END IF
  END DO inner_loop
  EXIT main_loop
  END DO main_loop
  !
  RETURN
END SUBROUTINE PROCP