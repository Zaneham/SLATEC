!** DSOSEQ
PURE SUBROUTINE DSOSEQ(FNC,N,S,Rtolx,Atolx,Tolf,Iflag,Mxit,Ncjs,Nsrrc,Nsri,&
    Fmax,C,Nc,B,P,Temp,X,Y,Fac,Is)
  !> Subsidiary to DSOS
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SOSEQS-S, DSOSEQ-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     DSOSEQ solves a system of N simultaneous nonlinear equations.
  !     See the comments in the interfacing routine DSOS for a more
  !     detailed description of some of the items in the calling list.
  !
  !***
  ! **See also:**  DSOS
  !***
  ! **Routines called:**  D1MACH, DSOSSL, I1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   2024xx  Eliminated GOTOs using structured control flow.  (ZH)
  USE ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT
  USE service, ONLY : eps_dp, tiny_dp

  INTERFACE
    REAL(DP) PURE FUNCTION FNC(X,K)
      IMPORT DP
      INTEGER, INTENT(IN) :: K
      REAL(DP), INTENT(IN) :: X(:)
    END FUNCTION FNC
  END INTERFACE
  INTEGER, INTENT(IN) :: N, Nc, Ncjs, Nsri, Nsrrc
  INTEGER, INTENT(INOUT) :: Mxit, Is(N)
  INTEGER, INTENT(OUT) :: Iflag
  REAL(DP), INTENT(IN) :: Atolx, Rtolx, Tolf
  REAL(DP), INTENT(INOUT) :: C(Nc), B(N), Fac(N), P(N), S(N), Temp(N), X(N), Y(N)
  REAL(DP), INTENT(OUT) :: Fmax

  INTEGER :: ic, icr, isj, isv, it, item, itry, j, jk, js, k, kd, kj, kk, km1, &
    kn, ksv, l, loun, ls, m, mit, np1
  REAL(DP) :: csv, f, fact, fdif, fmin, fmxs, fn1, fn2, fp, h, hx, pmax, re, &
    sruro, test, uro, xnorm, yj, yn1, yn2, yn3, ynorm, yns, zero
  LOGICAL :: singular_error, do_exit, converged, retry_perturb

  !* FIRST EXECUTABLE STATEMENT  DSOSEQ
  uro = eps_dp
  loun = OUTPUT_UNIT
  zero = tiny_dp
  re = MAX(Rtolx,uro)
  sruro = SQRT(uro)

  Iflag = 0
  np1 = N + 1
  icr = 0
  ic = 0
  itry = Ncjs
  yn1 = 0._DP
  yn2 = 0._DP
  yn3 = 0._DP
  yns = 0._DP
  mit = 0
  fn1 = 0._DP
  fn2 = 0._DP
  fmxs = 0._DP
  singular_error = .FALSE.
  do_exit = .FALSE.

  DO k = 1, N
    Is(k) = k
    X(k) = S(k)
    Temp(k) = X(k)
  END DO

  m = 0
  fmin = 0._DP
  main_loop: DO WHILE( m < Mxit .AND. .NOT. do_exit .AND. .NOT. singular_error )
    m = m + 1

    DO k = 1, N
      Fac(k) = sruro
    END DO

    retry_perturb = .TRUE.
    perturb_loop: DO WHILE( retry_perturb .AND. .NOT. singular_error )
      retry_perturb = .FALSE.
      kn = 1
      Fmax = 0._DP

      triangular_loop: DO k = 1, N
        km1 = k - 1

        IF( km1/=0 ) THEN
          CALL DSOSSL(k,N,km1,Y,C,B,kn)
          DO j = 1, km1
            js = Is(j)
            X(js) = Temp(js) + Y(j)
          END DO
        END IF

        f = FNC(X,k)
        Fmax = MAX(Fmax,ABS(f))

        IF( itry>=Ncjs ) THEN
          it = 0

          DO j = k, N
            item = Is(j)
            hx = X(item)
            h = Fac(item)*hx
            IF( ABS(h)<=zero ) h = Fac(item)
            X(item) = hx + h
            IF( km1/=0 ) THEN
              Y(j) = h
              CALL DSOSSL(k,N,j,Y,C,B,kn)
              DO l = 1, km1
                ls = Is(l)
                X(ls) = Temp(ls) + Y(l)
              END DO
            END IF
            fp = FNC(X,k)
            X(item) = hx
            fdif = fp - f
            IF( ABS(fdif)<=uro*ABS(f) ) THEN
              fdif = 0._DP
              it = it + 1
            END IF
            P(j) = fdif/h
          END DO

          IF( it>(N-k) ) THEN
            DO j = k, N
              isj = Is(j)
              fact = 100._DP*Fac(isj)
              IF( fact>1.0E10_DP ) THEN
                singular_error = .TRUE.
                EXIT triangular_loop
              END IF
              Fac(isj) = fact
            END DO
            retry_perturb = .TRUE.
            EXIT triangular_loop
          ELSEIF( k/=N ) THEN
            pmax = 0._DP
            DO j = k, N
              test = ABS(P(j))
              IF( test>pmax ) THEN
                pmax = test
                isv = j
              END IF
            END DO
            IF( pmax==0._DP ) THEN
              singular_error = .TRUE.
              EXIT triangular_loop
            END IF

            pmax = P(isv)
            kk = kn
            DO j = k, N
              IF( j/=isv ) C(kk) = -P(j)/pmax
              kk = kk + 1
            END DO
            P(k) = pmax

            IF( isv/=k ) THEN
              ksv = Is(k)
              Is(k) = Is(isv)
              Is(isv) = ksv
              kd = isv - k
              kj = k
              DO j = 1, k
                csv = C(kj)
                jk = kj + kd
                C(kj) = C(jk)
                C(jk) = csv
                kj = kj + N - j
              END DO
            END IF
          END IF
        END IF

        kn = kn + np1 - k
        B(k) = -f/P(k)
      END DO triangular_loop
    END DO perturb_loop

    IF( singular_error ) EXIT main_loop

    kn = kn - 1
    Y(N) = B(N)
    IF( N>1 ) CALL DSOSSL(N,N,N,Y,C,B,kn)
    xnorm = 0._DP
    ynorm = 0._DP
    DO j = 1, N
      yj = Y(j)
      ynorm = MAX(ynorm,ABS(yj))
      js = Is(j)
      X(js) = Temp(js) + yj
      xnorm = MAX(xnorm,ABS(X(js)))
    END DO

    converged = .TRUE.
    DO j = 1, N
      js = Is(j)
      IF( ABS(Y(j))>re*ABS(X(js))+Atolx ) THEN
        converged = .FALSE.
        EXIT
      END IF
    END DO
    IF( converged .AND. Fmax<=fmxs ) Iflag = 1

    IF( Fmax<=Tolf ) Iflag = Iflag + 2
    IF( Iflag>0 ) EXIT main_loop

    IF( m>1 ) THEN
      IF( Fmax<fmin ) THEN
        mit = m + 1
        yn1 = ynorm
        yn2 = yns
        fn1 = fmxs
        fmin = Fmax
        DO j = 1, N
          S(j) = X(j)
        END DO
        ic = 0
      END IF

      IF( ynorm<=sruro*xnorm ) THEN
        IF( Fmax>=0.2_DP*fmxs .AND. Fmax<=5._DP*fmxs ) THEN
          IF( ynorm>=0.2_DP*yns .AND. ynorm<=5._DP*yns ) THEN
            icr = icr + 1
            IF( icr>=Nsrrc ) THEN
              Iflag = 4
              Fmax = fmin
              do_exit = .TRUE.
              CYCLE main_loop
            ELSE
              ic = 0
            END IF
          END IF
        END IF
      END IF
      IF( .NOT. do_exit ) icr = 0

      IF( ynorm>2._DP*yns .OR. Fmax>2._DP*fmxs ) THEN
        ic = ic + 1
        IF( ic>=Nsri ) THEN
          Iflag = 7
          EXIT main_loop
        END IF
      ELSE
        ic = 0
      END IF
    ELSE
      fmin = Fmax
    END IF

    itry = itry - 1
    IF( itry==0 ) THEN
      itry = Ncjs
    ELSEIF( 20._DP*ynorm>xnorm ) THEN
      itry = Ncjs
    ELSEIF( ynorm>2._DP*yns ) THEN
      itry = Ncjs
    ELSEIF( Fmax>=2._DP*fmxs ) THEN
      itry = Ncjs
    END IF

    DO j = 1, N
      Temp(j) = X(j)
    END DO
    IF( m==mit ) THEN
      fn2 = Fmax
      yn3 = ynorm
    END IF
    fmxs = Fmax
    yns = ynorm
  END DO main_loop

  IF( singular_error ) THEN
    Iflag = 8
    DO j = 1, N
      S(j) = Temp(j)
    END DO
  ELSEIF( .NOT. do_exit .AND. Iflag==0 ) THEN
    Iflag = 5
    IF( yn1>10._DP*yn2 .OR. yn3>10._DP*yn1 ) Iflag = 6
    IF( fn1>5._DP*fmin .OR. fn2>5._DP*fmin ) Iflag = 6
    IF( Fmax>5._DP*fmin ) Iflag = 6
    DO j = 1, N
      S(j) = X(j)
    END DO
  ELSEIF( .NOT. do_exit ) THEN
    DO j = 1, N
      S(j) = X(j)
    END DO
  END IF

  Mxit = m

END SUBROUTINE DSOSEQ
