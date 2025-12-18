!** BESJ
SUBROUTINE BESJ(X,Alpha,N,Y,Nz)
  !> Compute an N member sequence of J Bessel functions
  !  J_{ALPHA+K-1}(X), K=1,...,N for non-negative ALPHA and X.
  !***
  ! **Library:**   SLATEC
  ! **Category:**  C10A3
  ! **Type:**      SINGLE PRECISION (BESJ-S, DBESJ-D)
  ! **Keywords:**  J BESSEL FUNCTION, SPECIAL FUNCTIONS
  ! **Author:**  Amos, D. E., (SNLA)
  !           Daniel, S. L., (SNLA)
  !           Weston, M. K., (SNLA)
  !* REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   251219  Eliminated 23 GOTOs per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.9 (SELECT CASE), S11.2.1 (DO)
  !           Pattern: State machine with istate + internal subroutines.
  USE service, ONLY : eps_2_sp, log10_radix_sp, tiny_sp, digits_sp, min_exp_sp
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(OUT) :: Nz
  REAL(SP), INTENT(IN) :: Alpha, X
  REAL(SP), INTENT(OUT) :: Y(N)
  INTEGER :: i, ialp, idalp, iflw, in, is, i1, i2, k, kk, km, kt, nn, ns
  REAL(SP) :: ak, akm, ans, ap, arg, coef, dalpha, dfn, dtm, earg, elim1, &
    etx, fidal, flgjy, fn, fnf, fni, fnp1, fnu, gln, rden, relb, &
    rtx, rzden, s, sa, sb, sxo2, s1, s2, t, ta, tau, tb, temp(3), &
    tfn, tm, tol, tolln, trx, tx, t1, t2, wk(7), xo2, xo2l, slim, rtol
  REAL(SP), PARAMETER :: rtwo = 1.34839972492648E+00_SP, pdf = 7.85398163397448E-01_SP, &
    rttp = 7.97884560802865E-01_SP, pidt = 1.57079632679490_SP
  REAL(SP), PARAMETER :: pp(4) = [ 8.72909153935547E+00_SP, 2.65693932265030E-01_SP, &
    1.24578576865586E-01_SP, 7.70133747430388E-04_SP ]
  INTEGER, PARAMETER :: inlim = 150
  REAL(SP), PARAMETER :: fnulim(2) = [ 100._SP, 60._SP ]
  INTEGER :: istate
  INTEGER, PARAMETER :: ST_DONE = 0, ST_ASYMP_NU = 1, ST_SERIES = 2
  INTEGER, PARAMETER :: ST_SERIES_LOOP = 3, ST_SERIES_UFL = 4, ST_BACKWARD = 5
  INTEGER, PARAMETER :: ST_ASYMP_X = 6, ST_NORM_RECUR = 7
  !* FIRST EXECUTABLE STATEMENT  BESJ
  Nz = 0
  kt = 1
  ns = 0
  ta = eps_2_sp
  tol = MAX(ta,1.E-15_SP)
  i1 = digits_sp + 1
  i2 = min_exp_sp
  tb = log10_radix_sp
  elim1 = -2.303_SP*(i2*tb+3._SP)
  rtol = 1._SP/tol
  slim = tiny_sp*rtol*1.E+3_SP
  tolln = 2.303_SP*tb*i1
  tolln = MIN(tolln,34.5388_SP)
  IF( N<1 ) THEN
    ERROR STOP 'BESJ : N < 1'
  ELSEIF( N==1 ) THEN
    kt = 2
  END IF
  nn = N
  IF( X<0 ) THEN
    ERROR STOP 'BESJ : X < 0'
  ELSEIF( Alpha<0._SP ) THEN
    ERROR STOP 'BESJ : ORDER, ALPHA, < 0'
  ELSEIF( X==0 ) THEN
    IF( Alpha==0 ) THEN
      Y(1) = 1._SP
      IF( N==1 ) RETURN
      i1 = 2
    ELSE
      i1 = 1
    END IF
    DO i = i1, N
      Y(i) = 0._SP
    END DO
    RETURN
  ELSE
    ialp = INT(Alpha)
    fni = ialp + N - 1
    fnf = Alpha - ialp
    dfn = fni + fnf
    fnu = dfn
    xo2 = X*0.5_SP
    sxo2 = xo2*xo2
    istate = ST_DONE
    ! Decision tree for region selection
    IF( sxo2<=(fnu+1._SP) ) THEN
      fn = fnu
      fnp1 = fn + 1._SP
      xo2l = LOG(xo2)
      is = kt
      IF( X<=0.50_SP ) THEN
        istate = ST_SERIES
      ELSE
        ns = 0
        fni = fni + ns
        dfn = fni + fnf
        fn = dfn
        fnp1 = fn + 1._SP
        is = kt
        IF( N-1+ns>0 ) is = 3
        istate = ST_SERIES
      END IF
    ELSE
      ta = MAX(20._SP,fnu)
      IF( X>ta ) THEN
        rtx = SQRT(X)
        tau = rtwo*rtx
        ta = tau + fnulim(kt)
        IF( fnu<=ta ) THEN
          in = INT(Alpha-tau+2._SP)
          IF( in<=0 ) THEN
            idalp = ialp
            in = 0
          ELSE
            idalp = ialp - in - 1
            kt = 1
          END IF
          is = kt
          fidal = idalp
          dalpha = fidal + fnf
          arg = X - pidt*dalpha - pdf
          sa = SIN(arg)
          sb = COS(arg)
          coef = rttp/rtx
          etx = 8._SP*X
          istate = ST_ASYMP_X
        ELSE
          fn = fnu
          is = kt
          istate = ST_ASYMP_NU
        END IF
      ELSEIF( X>12._SP ) THEN
        ans = MAX(36._SP-fnu,0._SP)
        ns = INT(ans)
        fni = fni + ns
        dfn = fni + fnf
        fn = dfn
        is = kt
        IF( N-1+ns>0 ) is = 3
        istate = ST_ASYMP_NU
      ELSE
        xo2l = LOG(xo2)
        ns = INT(sxo2-fnu) + 1
        fni = fni + ns
        dfn = fni + fnf
        fn = dfn
        fnp1 = fn + 1._SP
        is = kt
        IF( N-1+ns>0 ) is = 3
        istate = ST_SERIES
      END IF
    END IF
    ! Main state machine loop
    state_machine: DO WHILE( istate /= ST_DONE )
      SELECT CASE( istate )
      CASE( ST_ASYMP_NU )
        asymp_nu_loop: DO
          i1 = ABS(3-is)
          i1 = MAX(i1,1)
          flgjy = 1._SP
          CALL ASYJY(JAIRY,X,fn,flgjy,i1,temp(is),wk,iflw)
          IF( iflw/=0 ) THEN
            Y(nn) = 0._SP
            nn = nn - 1
            fni = fni - 1._SP
            dfn = fni + fnf
            fn = dfn
            IF( nn<1 ) THEN
              Nz = N - nn
              istate = ST_DONE
              EXIT asymp_nu_loop
            END IF
            IF( nn==1 ) THEN
              kt = 2
              is = 2
            END IF
            CYCLE asymp_nu_loop
          END IF
          SELECT CASE (is)
            CASE (1)
              is = 2
              fni = fni - 1._SP
              dfn = fni + fnf
              fn = dfn
              IF( i1==2 ) THEN
                istate = ST_BACKWARD
                EXIT asymp_nu_loop
              END IF
              CYCLE asymp_nu_loop
            CASE (2)
              istate = ST_BACKWARD
              EXIT asymp_nu_loop
            CASE (3)
              gln = wk(3) + wk(2)
              IF( wk(6)>30._SP ) THEN
                ta = 0.5_SP*tolln/wk(4)
                ta = ((0.0493827160_SP*ta-0.1111111111_SP)*ta+0.6666666667_SP)*ta*wk(6)
                IF( wk(1)<0.10_SP ) THEN
                  tb = (1.259921049_SP+(0.1679894730_SP+0.0887944358_SP*wk(1))*wk(1))/wk(7)
                ELSE
                  tb = gln/wk(5)
                END IF
              ELSE
                rden = (pp(4)*wk(6)+pp(3))*wk(6) + 1._SP
                rzden = pp(1) + pp(2)*wk(6)
                ta = rzden/rden
                IF( wk(1)<0.10_SP ) THEN
                  tb = (1.259921049_SP+(0.1679894730_SP+0.0887944358_SP*wk(1))*wk(1))/wk(7)
                ELSE
                  tb = gln/wk(5)
                END IF
              END IF
              in = INT(ta/tb+1.5_SP)
              IF( in<=inlim ) THEN
                istate = ST_NORM_RECUR
                EXIT asymp_nu_loop
              END IF
              temp(1) = temp(3)
              kt = 1
              is = 2
              fni = fni - 1._SP
              dfn = fni + fnf
              fn = dfn
              IF( i1==2 ) THEN
                istate = ST_BACKWARD
                EXIT asymp_nu_loop
              END IF
              CYCLE asymp_nu_loop
            CASE DEFAULT
              temp(1) = temp(3)
              kt = 1
              is = 2
              fni = fni - 1._SP
              dfn = fni + fnf
              fn = dfn
              IF( i1==2 ) THEN
                istate = ST_BACKWARD
                EXIT asymp_nu_loop
              END IF
              CYCLE asymp_nu_loop
          END SELECT
        END DO asymp_nu_loop
      CASE( ST_SERIES )
        gln = LOG_GAMMA(fnp1)
        arg = fn*xo2l - gln
        IF( arg<(-elim1) ) THEN
          istate = ST_SERIES_UFL
        ELSE
          earg = EXP(arg)
          istate = ST_SERIES_LOOP
        END IF
      CASE( ST_SERIES_LOOP )
        s = 1._SP
        IF( X>=tol ) THEN
          ak = 3._SP
          t2 = 1._SP
          t = 1._SP
          s1 = fn
          DO k = 1, 17
            s2 = t2 + s1
            t = -t*sxo2/s2
            s = s + t
            IF( ABS(t)<tol ) EXIT
            t2 = t2 + ak
            ak = ak + 2._SP
            s1 = s1 + fn
          END DO
        END IF
        temp(is) = s*earg
        SELECT CASE (is)
          CASE (2)
            istate = ST_BACKWARD
          CASE (3)
            akm = MAX(3._SP-fn,0._SP)
            km = INT(akm)
            tfn = fn + km
            ta = (gln+tfn-0.9189385332_SP-0.0833333333_SP/tfn)/(tfn+0.5_SP)
            ta = xo2l - ta
            tb = -(1._SP-1.5_SP/tfn)/tfn
            akm = tolln/(-ta+SQRT(ta*ta-tolln*tb)) + 1.5_SP
            in = km + INT(akm)
            istate = ST_NORM_RECUR
          CASE DEFAULT
            earg = earg*fn/xo2
            fni = fni - 1._SP
            dfn = fni + fnf
            fn = dfn
            is = 2
        END SELECT
      CASE( ST_SERIES_UFL )
        series_ufl_loop: DO
          Y(nn) = 0._SP
          nn = nn - 1
          fnp1 = fn
          fni = fni - 1._SP
          dfn = fni + fnf
          fn = dfn
          IF( nn<1 ) THEN
            Nz = N - nn
            istate = ST_DONE
            EXIT series_ufl_loop
          END IF
          IF( nn==1 ) THEN
            kt = 2
            is = 2
          END IF
          IF( sxo2>fnp1 ) THEN
            istate = ST_ASYMP_NU
            EXIT series_ufl_loop
          END IF
          arg = arg - xo2l + LOG(fnp1)
          IF( arg>=(-elim1) ) THEN
            fnp1 = fn + 1._SP
            gln = LOG_GAMMA(fnp1)
            earg = EXP(arg)
            istate = ST_SERIES_LOOP
            EXIT series_ufl_loop
          END IF
        END DO series_ufl_loop
      CASE( ST_BACKWARD )
        CALL backward_recursion()
        istate = ST_DONE
      CASE( ST_ASYMP_X )
        asymp_x_loop: DO
          dtm = fidal + fidal
          dtm = dtm*dtm
          tm = 0._SP
          IF( fidal/=0._SP .OR. ABS(fnf)>=tol ) tm = 4._SP*fnf*(fidal+fidal+fnf)
          trx = dtm - 1._SP
          t2 = (trx+tm)/etx
          s2 = t2
          relb = tol*ABS(t2)
          t1 = etx
          s1 = 1._SP
          fn = 1._SP
          ak = 8._SP
          DO k = 1, 13
            t1 = t1 + etx
            fn = fn + ak
            trx = dtm - fn
            ap = trx + tm
            t2 = -t2*ap/t1
            s1 = s1 + t2
            t1 = t1 + etx
            ak = ak + 8._SP
            fn = fn + ak
            trx = dtm - fn
            ap = trx + tm
            t2 = t2*ap/t1
            s2 = s2 + t2
            IF( ABS(t2)<=relb ) EXIT
            ak = ak + 8._SP
          END DO
          temp(is) = coef*(s1*sb-s2*sa)
          IF( is==2 ) THEN
            CALL forward_recursion()
            istate = ST_DONE
            EXIT asymp_x_loop
          END IF
          fidal = fidal + 1._SP
          dalpha = fidal + fnf
          is = 2
          tb = sa
          sa = -sb
          sb = tb
        END DO asymp_x_loop
      CASE( ST_NORM_RECUR )
        CALL normalization_recurrence()
        istate = ST_DONE
      END SELECT
    END DO state_machine
  END IF

CONTAINS

  SUBROUTINE backward_recursion()
    IF( ns==0 ) THEN
      Nz = N - nn
      IF( kt==2 ) THEN
        Y(1) = temp(2)
        RETURN
      END IF
      Y(nn) = temp(1)
      Y(nn-1) = temp(2)
      IF( nn==2 ) RETURN
    END IF
    trx = 2._SP/X
    dtm = fni
    tm = (dtm+fnf)*trx
    ak = 1._SP
    ta = temp(1)
    tb = temp(2)
    IF( ABS(ta)<=slim ) THEN
      ta = ta*rtol
      tb = tb*rtol
      ak = tol
    END IF
    kk = 2
    in = ns - 1
    IF( in==0 ) THEN
      CALL final_backward_recur()
      RETURN
    END IF
    IF( ns/=0 ) THEN
      CALL backward_recur_unindexed()
      CALL final_backward_recur()
      RETURN
    END IF
    k = nn - 2
    DO i = 3, nn
      s = tb
      tb = tm*tb - ta
      ta = s
      Y(k) = tb*ak
      dtm = dtm - 1._SP
      tm = (dtm+fnf)*trx
      k = k - 1
    END DO
  END SUBROUTINE backward_recursion

  SUBROUTINE forward_recursion()
    IF( kt==2 ) THEN
      Y(1) = temp(2)
      RETURN
    END IF
    s1 = temp(1)
    s2 = temp(2)
    tx = 2._SP/X
    tm = dalpha*tx
    IF( in/=0 ) THEN
      DO i = 1, in
        s = s2
        s2 = tm*s2 - s1
        tm = tm + tx
        s1 = s
      END DO
      IF( nn==1 ) THEN
        Y(1) = s2
        RETURN
      END IF
      s = s2
      s2 = tm*s2 - s1
      tm = tm + tx
      s1 = s
    END IF
    Y(1) = s1
    Y(2) = s2
    IF( nn==2 ) RETURN
    DO i = 3, nn
      Y(i) = tm*Y(i-1) - Y(i-2)
      tm = tm + tx
    END DO
  END SUBROUTINE forward_recursion

  SUBROUTINE backward_recur_unindexed()
    norm_loop: DO
      DO i = 1, in
        s = tb
        tb = tm*tb - ta
        ta = s
        dtm = dtm - 1._SP
        tm = (dtm+fnf)*trx
      END DO
      IF( kk/=1 ) EXIT norm_loop
      s = temp(3)
      sa = ta/tb
      ta = s
      tb = s
      IF( ABS(s)<=slim ) THEN
        ta = ta*rtol
        tb = tb*rtol
        ak = tol
      END IF
      ta = ta*sa
      kk = 2
      in = ns
      IF( ns==0 ) EXIT norm_loop
    END DO norm_loop
  END SUBROUTINE backward_recur_unindexed

  SUBROUTINE final_backward_recur()
    Y(nn) = tb*ak
    Nz = N - nn
    IF( nn==1 ) RETURN
    k = nn - 1
    s = tb
    tb = tm*tb - ta
    ta = s
    Y(k) = tb*ak
    IF( nn==2 ) RETURN
    dtm = dtm - 1._SP
    tm = (dtm+fnf)*trx
    k = nn - 2
    DO i = 3, nn
      s = tb
      tb = tm*tb - ta
      ta = s
      Y(k) = tb*ak
      dtm = dtm - 1._SP
      tm = (dtm+fnf)*trx
      k = k - 1
    END DO
  END SUBROUTINE final_backward_recur

  SUBROUTINE normalization_recurrence()
    dtm = fni + in
    trx = 2._SP/X
    tm = (dtm+fnf)*trx
    ta = 0._SP
    tb = tol
    kk = 1
    ak = 1._SP
    CALL backward_recur_unindexed()
    CALL final_backward_recur()
  END SUBROUTINE normalization_recurrence

END SUBROUTINE BESJ
