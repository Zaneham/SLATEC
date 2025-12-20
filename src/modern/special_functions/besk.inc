!** BESK
PURE SUBROUTINE BESK(X,Fnu,Kode,N,Y,Nz)
  !> Implement forward recursion on the three term recursion
  !  relation for a sequence of non-negative order Bessel
  !  functions K_{FNU+I-1}(X), or scaled Bessel functions
  !  EXP(X)*K_{FNU+I-1}(X), I=1,...,N for REAL, positive
  !  X and non-negative orders FNU.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C10B3
  !***
  ! **Type:**      SINGLE PRECISION (BESK-S, DBESK-D)
  !***
  ! **Keywords:**  K BESSEL FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         BESK implements forward recursion on the three term
  !         recursion relation for a sequence of non-negative order Bessel
  !         functions K_{FNU+I-1}(X), or scaled Bessel functions
  !         EXP(X)*K_{FNU+I-1}(X), I=1,..,N for real X > 0.0 and
  !         non-negative orders FNU.  If FNU < NULIM, orders FNU and
  !         FNU+1 are obtained from BSKNU to start the recursion.  If
  !         FNU >= NULIM, the uniform asymptotic expansion is used for
  !         orders FNU and FNU+1 to start the recursion.  NULIM is 35 or
  !         70 depending on whether N=1 or N >= 2.
  !
  !     Description of Arguments
  !
  !         Input
  !           X      - X > 0.0
  !           FNU    - order of the initial K function, FNU >= 0.0
  !           KODE   - a parameter to indicate the scaling option
  !                    KODE=1 returns Y(I)=       K_{FNU+I-1}(X), I=1,...,N
  !                    KODE=2 returns Y(I)=EXP(X)*K_{FNU+I-1}(X), I=1,...,N
  !           N      - number of members in the sequence, N >= 1
  !
  !         Output
  !           Y      - a vector whose first N components contain values
  !                    for the sequence Y(I)=K_{FNU+I-1}(X), I=1,...,N or
  !                    Y(I)=EXP(X)*K_{FNU+I-1}(X), I=1,...,N depending on KODE
  !           NZ     - number of components of Y set to zero due to
  !                    underflow with KODE=1
  !
  !     Error Conditions
  !         Improper input arguments - a fatal error
  !         Overflow - a fatal error
  !         Underflow with KODE=1 - a non-fatal error (NZ /= 0)
  !
  !***
  ! **References:**  F. W. J. Olver, Tables of Bessel Functions of Moderate
  !                 or Large Orders, NPL Mathematical Tables 6, Her
  !                 Majesty's Stationery Office, London, 1962.
  !               N. M. Temme, On the numerical evaluation of the modified
  !                 Bessel function of the third kind, Journal of
  !                 Computational Physics 19, (1975), pp. 324-337.
  !***
  ! **Routines called:**  ASYIK, BESK0, BESK1, BSK0E, BSK1E, BSKNU, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790201  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   2024xx  Eliminated GOTOs using structured control flow.  (ZH)
  USE service, ONLY : log10_radix_sp, tiny_sp, min_exp_sp
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(SP), INTENT(IN) :: Fnu, X
  REAL(SP), INTENT(OUT) :: Y(N)
  INTEGER :: i, j, k, mz, nb, nd, nn, nud
  REAL(SP) :: cn, dnu, elim, etx, flgik, fn, fnn, gln, gnu, rtz, s, s1, s2, t, tm, &
    trx, w(2), xlim, zn
  INTEGER, PARAMETER :: nulim(2) = [ 35, 70 ]
  LOGICAL :: do_underflow_update, do_overflow_test, do_bsknu, do_asymp
  LOGICAL :: skip_to_y1, skip_to_forward_recur, done
  !* FIRST EXECUTABLE STATEMENT  BESK
  nn = -min_exp_sp
  elim = 2.303_SP*(nn*log10_radix_sp-3._SP)
  xlim = tiny_sp*1.E+3_SP
  IF( Kode<1 .OR. Kode>2 ) THEN
    ERROR STOP 'BESK : SCALING OPTION, KODE, NOT 1 OR 2'
    RETURN
  ELSEIF( Fnu<0._SP ) THEN
    ERROR STOP 'BESK : ORDER, FNU, < 0'
    RETURN
  ELSEIF( X<=0._SP ) THEN
    ERROR STOP 'BESK : X <= 0 '
    RETURN
  ELSEIF( X<xlim ) THEN
    ERROR STOP 'BESK : OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL'
    RETURN
  ELSEIF( N<1 ) THEN
    ERROR STOP 'BESK : N < 1'
    RETURN
  END IF

  etx = Kode - 1
  nd = N
  Nz = 0
  nud = INT(Fnu)
  dnu = Fnu - nud
  gnu = Fnu
  nn = MIN(2,nd)
  fn = Fnu + N - 1
  fnn = fn

  ! Initialize control flags
  do_underflow_update = .FALSE.
  do_overflow_test = .FALSE.
  do_bsknu = .FALSE.
  do_asymp = .FALSE.
  done = .FALSE.

  IF( fn<2._SP ) THEN
    ! UNDERFLOW TEST FOR KODE=1
    IF( Kode==2 ) THEN
      do_overflow_test = .TRUE.
    ELSEIF( X<=elim ) THEN
      do_overflow_test = .TRUE.
    ELSE
      do_underflow_update = .TRUE.
    END IF
  ELSE
    ! OVERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
    zn = X/fn
    IF( zn==0._SP ) THEN
      ERROR STOP 'BESK : OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL'
      RETURN
    END IF
    rtz = SQRT(1._SP+zn*zn)
    gln = LOG((1._SP+rtz)/zn)
    t = rtz*(1._SP-etx) + etx/(zn+rtz)
    cn = -fn*(t-gln)
    IF( cn>elim ) THEN
      ERROR STOP 'BESK : OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL'
      RETURN
    ELSEIF( nud<nulim(nn) ) THEN
      IF( Kode==2 ) THEN
        do_bsknu = .TRUE.
      ELSEIF( X<=elim ) THEN
        do_bsknu = .TRUE.
      ELSE
        do_underflow_update = .TRUE.
      END IF
    ELSEIF( nn==1 ) THEN
      ! Skip underflow test, go directly to check cn < -elim
      do_asymp = .TRUE.
    ELSE
      ! Need underflow test for first order (label 100)
      fn = gnu
      zn = X/fn
      rtz = SQRT(1._SP+zn*zn)
      gln = LOG((1._SP+rtz)/zn)
      t = rtz*(1._SP-etx) + etx/(zn+rtz)
      cn = -fn*(t-gln)
      do_asymp = .TRUE.
    END IF
  END IF

  ! UNDERFLOW UPDATE LOOP (was label 700)
  IF( do_underflow_update ) THEN
    underflow_loop: DO
      nud = nud + 1
      nd = nd - 1
      IF( nd==0 ) EXIT underflow_loop
      nn = MIN(2,nd)
      gnu = gnu + 1._SP
      IF( fnn>=2._SP ) THEN
        IF( nud>=nulim(nn) ) THEN
          ! Compute underflow test (was label 100)
          fn = gnu
          zn = X/fn
          rtz = SQRT(1._SP+zn*zn)
          gln = LOG((1._SP+rtz)/zn)
          t = rtz*(1._SP-etx) + etx/(zn+rtz)
          cn = -fn*(t-gln)
          IF( cn>=(-elim) ) THEN
            do_asymp = .TRUE.
            EXIT underflow_loop
          END IF
        END IF
      END IF
    END DO underflow_loop
    done = (nd==0) .AND. .NOT. do_asymp
  END IF

  ! ASYMPTOTIC EXPANSION SECTION (was labels 200, 500)
  IF( do_asymp .AND. .NOT. done ) THEN
    IF( cn<(-elim) ) THEN
      ! Need underflow update
      underflow_loop2: DO
        nud = nud + 1
        nd = nd - 1
        IF( nd==0 ) EXIT underflow_loop2
        nn = MIN(2,nd)
        gnu = gnu + 1._SP
        IF( fnn>=2._SP ) THEN
          IF( nud>=nulim(nn) ) THEN
            fn = gnu
            zn = X/fn
            rtz = SQRT(1._SP+zn*zn)
            gln = LOG((1._SP+rtz)/zn)
            t = rtz*(1._SP-etx) + etx/(zn+rtz)
            cn = -fn*(t-gln)
            IF( cn>=(-elim) ) EXIT underflow_loop2
          END IF
        END IF
      END DO underflow_loop2
      done = (nd==0)
    END IF

    IF( .NOT. done ) THEN
      flgik = -1._SP
      CALL ASYIK(X,gnu,Kode,flgik,rtz,cn,nn,Y)
      IF( nn/=1 ) THEN
        trx = 2._SP/X
        tm = (gnu+gnu+2._SP)/X
        ! Forward recurrence (was label 500)
        IF( nd/=2 ) THEN
          DO i = 3, nd
            Y(i) = tm*Y(i-1) + Y(i-2)
            tm = tm + trx
          END DO
        END IF
      END IF
      done = .TRUE.
    END IF
  END IF

  ! OVERFLOW TEST SECTION (was label 600)
  IF( do_overflow_test .AND. .NOT. done ) THEN
    IF( fn>1._SP ) THEN
      IF( -fn*(LOG(X)-0.693_SP)>elim ) THEN
        ERROR STOP 'BESK : OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL'
        RETURN
      END IF
    END IF
    IF( dnu==0._SP ) THEN
      j = nud
      IF( j/=1 ) THEN
        j = j + 1
        IF( Kode==2 ) THEN
          Y(j) = BSK0E(X)
        ELSE
          Y(j) = BESK0(X)
        END IF
        IF( nd==1 ) THEN
          done = .TRUE.
        ELSE
          j = j + 1
        END IF
      END IF
      IF( .NOT. done ) THEN
        IF( Kode==2 ) THEN
          Y(j) = BSK1E(X)
        ELSE
          Y(j) = BESK1(X)
        END IF
        done = .TRUE.
      END IF
    ELSE
      CALL BSKNU(X,Fnu,Kode,nd,Y,mz)
      done = .TRUE.
    END IF
  END IF

  ! BSKNU SECTION (was label 300)
  IF( do_bsknu .AND. .NOT. done ) THEN
    skip_to_y1 = .FALSE.
    skip_to_forward_recur = .FALSE.

    IF( dnu/=0._SP ) THEN
      nb = 2
      IF( nud==0 .AND. nd==1 ) nb = 1
      CALL BSKNU(X,dnu,Kode,nb,w,Nz)
      s1 = w(1)
      IF( nb==1 ) THEN
        skip_to_y1 = .TRUE.
      ELSE
        s2 = w(2)
      END IF
    ELSE
      IF( Kode==2 ) THEN
        s1 = BSK0E(X)
      ELSE
        s1 = BESK0(X)
      END IF
      IF( nud==0 .AND. nd==1 ) THEN
        skip_to_y1 = .TRUE.
      ELSE
        IF( Kode==2 ) THEN
          s2 = BSK1E(X)
        ELSE
          s2 = BESK1(X)
        END IF
      END IF
    END IF

    IF( .NOT. skip_to_y1 ) THEN
      trx = 2._SP/X
      tm = (dnu+dnu+2._SP)/X
      ! Forward recur from DNU to FNU+1 to get Y(1) and Y(2)
      IF( nd==1 ) nud = nud - 1
      IF( nud>0 ) THEN
        DO i = 1, nud
          s = s2
          s2 = tm*s2 + s1
          s1 = s
          tm = tm + trx
        END DO
        IF( nd==1 ) s1 = s2
      ELSEIF( nd<=1 ) THEN
        s1 = s2
      END IF
    END IF

    ! Y(1) = s1 (was label 400)
    Y(1) = s1
    IF( nd/=1 ) THEN
      Y(2) = s2
      ! Forward recurrence (was label 500)
      IF( nd/=2 ) THEN
        DO i = 3, nd
          Y(i) = tm*Y(i-1) + Y(i-2)
          tm = tm + trx
        END DO
      END IF
    END IF
    done = .TRUE.
  END IF

  ! FINAL RETURN SECTION (was label 800)
  Nz = N - nd
  IF( Nz/=0 ) THEN
    IF( nd/=0 ) THEN
      DO i = 1, nd
        j = N - i + 1
        k = nd - i + 1
        Y(j) = Y(k)
      END DO
    END IF
    DO i = 1, Nz
      Y(i) = 0._SP
    END DO
  END IF

END SUBROUTINE BESK
