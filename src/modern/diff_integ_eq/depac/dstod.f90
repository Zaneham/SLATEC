!** DSTOD
SUBROUTINE DSTOD(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Acor,Wm,Iwm,DF,DJAC)
  !> Subsidiary to DDEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (STOD-S, DSTOD-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !   DSTOD integrates a system of first order odes over one step in the
  !   integrator package DDEBDF.
  ! ----------------------------------------------------------------------
  ! DSTOD  performs one step of the integration of an initial value
  ! problem for a system of ordinary differential equations.
  ! Note.. DSTOD  is independent of the value of the iteration method
  ! indicator MITER, when this is /= 0, and hence is independent
  ! of the type of chord method used, or the Jacobian structure.
  ! Communication with DSTOD  is done with the following variables..
  !
  ! Y      = An array of length >= N used as the Y argument in
  !          all calls to DF and DJAC.
  ! NEQ    = Integer array containing problem size in NEQ(1), and
  !          passed as the NEQ argument in all calls to DF and DJAC.
  ! YH     = An NYH by LMAX array containing the dependent variables
  !          and their approximate scaled derivatives, where
  !          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate
  !          J-th derivative of Y(I), scaled by H**J/FACTORIAL(J)
  !          (J = 0,1,...,NQ).  On entry for the first step, the first
  !          two columns of YH must be set from the initial values.
  ! NYH    = A constant integer >= N, the first dimension of YH.
  ! YH1    = A one-dimensional array occupying the same space as YH.
  ! EWT    = An array of N elements with which the estimated local
  !          errors in YH are compared.
  ! SAVF   = An array of working storage, of length N.
  ! ACOR   = A work array of length N, used for the accumulated
  !          corrections.  On a successful return, ACOR(I) contains
  !          the estimated one-step local error in Y(I).
  ! WM,IWM = DOUBLE PRECISION and INTEGER work arrays associated with
  !          matrix operations in chord iteration (MITER /= 0).
  ! DPJAC   = Name of routine to evaluate and preprocess Jacobian matrix
  !          if a chord method is being used.
  ! DSLVS   = Name of routine to solve linear system in chord iteration.
  ! H      = The step size to be attempted on the next step.
  !          H is altered by the error control algorithm during the
  !          problem.  H can be either positive or negative, but its
  !          sign must remain constant throughout the problem.
  ! HMIN   = The minimum absolute value of the step size H to be used.
  ! HMXI   = Inverse of the maximum absolute value of H to be used.
  !          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
  !          HMIN and HMXI may be changed at any time, but will not
  !          take effect until the next change of H is considered.
  ! TN     = The independent variable. TN is updated on each step taken.
  ! JSTART = An integer used for input only, with the following
  !          values and meanings..
  !               0  Perform the first step.
  !           >0  Take a new step continuing from the last.
  !              -1  Take the next step with a new value of H, MAXORD,
  !                    N, METH, MITER, and/or matrix parameters.
  !              -2  Take the next step with a new value of H,
  !                    but with other inputs unchanged.
  !          On return, JSTART is set to 1 to facilitate continuation.
  ! KFLAG  = a completion code with the following meanings..
  !               0  The step was successful.
  !              -1  The requested error could not be achieved.
  !              -2  Corrector convergence could not be achieved.
  !          A return with KFLAG = -1 or -2 means either
  !          ABS(H) = HMIN or 10 consecutive failures occurred.
  !          On a return with KFLAG negative, the values of TN and
  !          the YH array are as of the beginning of the last
  !          step, and H is the last step size attempted.
  ! MAXORD = The maximum order of integration method to be allowed.
  ! METH/MITER = The method flags.  See description in driver.
  ! N      = The number of first-order differential equations.
  ! ----------------------------------------------------------------------
  !
  !***
  ! **See also:**  DDEBDF
  !***
  ! **Routines called:**  DCFOD, DPJAC, DSLVS, DVNRMS
  !***
  ! COMMON BLOCKS    DDEBD1

  !* REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920422  Changed DIMENSION statement.  (WRB)
  USE DDEBD1, ONLY : conit_com, crate_com, el_com, elco_com, hold_com, rc_com, &
    rmax_com, tesco_com, el0_com, h_com, hmin_com, hmxi_com, hu_com, tn_com, &
    ksteps_com, ialth_com, ipup_com, lmax_com, meo_com, nqnyh_com, nstepj_com, &
    ier_com, jstart_com, kflag_com, l_com, meth_com, miter_com, maxord_com, n_com, &
    nq_com, nst_com, nfe_com, nqu_com
  !
  INTERFACE
    SUBROUTINE DF(X,U,Uprime)
      IMPORT DP
      REAL(DP), INTENT(IN) :: X
      REAL(DP), INTENT(IN) :: U(:)
      REAL(DP), INTENT(OUT) :: Uprime(:)
    END SUBROUTINE DF
    PURE SUBROUTINE DJAC(X,U,Pd,Nrowpd)
      IMPORT DP
      INTEGER, INTENT(IN) :: Nrowpd
      REAL(DP), INTENT(IN) :: X
      REAL(DP), INTENT(IN) :: U(:)
      REAL(DP), INTENT(OUT) :: Pd(:,:)
    END SUBROUTINE DJAC
  END INTERFACE
  INTEGER, INTENT(IN) :: Neq, Nyh
  INTEGER, INTENT(INOUT) :: Iwm(:)
  REAL(DP), INTENT(IN) :: Ewt(n_com)
  REAL(DP), INTENT(INOUT) :: Yh(Nyh,maxord_com+1), Yh1(Nyh*maxord_com+Nyh), Wm(:)
  REAL(DP), INTENT(OUT) :: Y(n_com), Savf(n_com), Acor(n_com)
  !
  INTEGER :: i, i1, iredo, j, jb, m, ncf, newq
  REAL(DP) :: dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup, r, rh, rhdn, &
    rhsm, rhup, told
  LOGICAL :: do_rescale, do_prediction, do_corrector, do_exit
  LOGICAL :: update_coefficients, change_order
  INTEGER :: istate  ! State for coefficient update: 1=normal, 2=order change, 3=first step
  !
  !* FIRST EXECUTABLE STATEMENT  DSTOD
  kflag_com = 0
  told = tn_com
  ncf = 0
  delp = 0._DP
  do_exit = .FALSE.

  ! Initialisation based on jstart
  IF( jstart_com > 0 ) THEN
    ! Continue from last step - go directly to prediction
    do_rescale = .FALSE.
    do_prediction = .TRUE.
  ELSEIF( jstart_com == -1 ) THEN
    ! New H, MAXORD, N, METH, MITER, and/or matrix parameters
    ipup_com = miter_com
    lmax_com = maxord_com + 1
    IF( ialth_com == 1 ) ialth_com = 2

    IF( meth_com /= meo_com ) THEN
      CALL DCFOD(meth_com, elco_com, tesco_com)
      meo_com = meth_com
      IF( nq_com <= maxord_com ) THEN
        ialth_com = l_com
        istate = 1
        update_coefficients = .TRUE.
        do_rescale = .FALSE.
        do_prediction = .FALSE.
      ELSE
        ! Reduce order to maxord
        nq_com = maxord_com
        l_com = lmax_com
        el_com(1:l_com) = elco_com(1:l_com, nq_com)
        nqnyh_com = nq_com * Nyh
        rc_com = rc_com * el_com(1) / el0_com
        el0_com = el_com(1)
        conit_com = 0.5_DP / (nq_com + 2)
        ddn = DVNRMS(n_com, Savf, Ewt) / tesco_com(1, l_com)
        exdn = 1._DP / l_com
        rhdn = 1._DP / (1.3_DP * ddn**exdn + 0.0000013_DP)
        rh = MIN(rhdn, 1._DP)
        iredo = 3
        IF( h_com == hold_com ) THEN
          rh = MAX(rh, hmin_com / ABS(h_com))
        ELSE
          rh = MIN(rh, ABS(h_com / hold_com))
          h_com = hold_com
        END IF
        update_coefficients = .FALSE.
        do_rescale = .TRUE.
        do_prediction = .FALSE.
      END IF
    ELSEIF( nq_com <= maxord_com ) THEN
      ! Check if H changed
      IF( h_com == hold_com ) THEN
        do_rescale = .FALSE.
        do_prediction = .TRUE.
      ELSE
        rh = h_com / hold_com
        h_com = hold_com
        iredo = 3
        do_rescale = .TRUE.
        do_prediction = .FALSE.
      END IF
      update_coefficients = .FALSE.
    ELSE
      ! Reduce order to maxord
      nq_com = maxord_com
      l_com = lmax_com
      el_com(1:l_com) = elco_com(1:l_com, nq_com)
      nqnyh_com = nq_com * Nyh
      rc_com = rc_com * el_com(1) / el0_com
      el0_com = el_com(1)
      conit_com = 0.5_DP / (nq_com + 2)
      ddn = DVNRMS(n_com, Savf, Ewt) / tesco_com(1, l_com)
      exdn = 1._DP / l_com
      rhdn = 1._DP / (1.3_DP * ddn**exdn + 0.0000013_DP)
      rh = MIN(rhdn, 1._DP)
      iredo = 3
      IF( h_com == hold_com ) THEN
        rh = MAX(rh, hmin_com / ABS(h_com))
      ELSE
        rh = MIN(rh, ABS(h_com / hold_com))
        h_com = hold_com
      END IF
      update_coefficients = .FALSE.
      do_rescale = .TRUE.
      do_prediction = .FALSE.
    END IF
  ELSEIF( jstart_com == -2 ) THEN
    ! New H only
    IF( h_com == hold_com ) THEN
      do_rescale = .FALSE.
      do_prediction = .TRUE.
    ELSE
      rh = h_com / hold_com
      h_com = hold_com
      iredo = 3
      do_rescale = .TRUE.
      do_prediction = .FALSE.
    END IF
    update_coefficients = .FALSE.
  ELSE
    ! First call (jstart == 0)
    lmax_com = maxord_com + 1
    nq_com = 1
    l_com = 2
    ialth_com = 2
    rmax_com = 10000._DP
    rc_com = 0._DP
    el0_com = 1._DP
    crate_com = 0.7_DP
    delp = 0._DP
    hold_com = h_com
    meo_com = meth_com
    nstepj_com = 0
    CALL DCFOD(meth_com, elco_com, tesco_com)
    istate = 3
    update_coefficients = .TRUE.
    do_rescale = .FALSE.
    do_prediction = .FALSE.
  END IF

  ! Main integration loop
  main_loop: DO WHILE( .NOT. do_exit )

    ! Update coefficients if needed
    IF( update_coefficients ) THEN
      el_com(1:l_com) = elco_com(1:l_com, nq_com)
      nqnyh_com = nq_com * Nyh
      rc_com = rc_com * el_com(1) / el0_com
      el0_com = el_com(1)
      conit_com = 0.5_DP / (nq_com + 2)

      SELECT CASE( istate )
        CASE( 2 )
          ! Order change - need rescale
          rh = MAX(rh, hmin_com / ABS(h_com))
          do_rescale = .TRUE.
        CASE( 3 )
          ! First step or retry after failure - go to prediction
          do_prediction = .TRUE.
        CASE DEFAULT
          ! Normal continuation
          do_prediction = .TRUE.
      END SELECT
      update_coefficients = .FALSE.
    END IF

    ! Rescale YH array if needed
    IF( do_rescale ) THEN
      rh = MIN(rh, rmax_com)
      rh = rh / MAX(1._DP, ABS(h_com) * hmxi_com * rh)
      r = 1._DP
      DO j = 2, l_com
        r = r * rh
        DO i = 1, n_com
          Yh(i,j) = Yh(i,j) * r
        END DO
      END DO
      h_com = h_com * rh
      rc_com = rc_com * rh
      ialth_com = l_com

      IF( iredo == 0 ) THEN
        rmax_com = 10._DP
        r = 1._DP / tesco_com(2, nqu_com)
        Acor(1:n_com) = Acor(1:n_com) * r
        do_exit = .TRUE.
        CYCLE main_loop
      END IF
      do_rescale = .FALSE.
      do_prediction = .TRUE.
    END IF

    ! Prediction step
    IF( do_prediction ) THEN
      IF( ABS(rc_com - 1._DP) > 0.3_DP ) ipup_com = miter_com
      IF( nst_com >= nstepj_com + 20 ) ipup_com = miter_com
      tn_com = tn_com + h_com
      i1 = nqnyh_com + 1
      DO jb = 1, nq_com
        i1 = i1 - Nyh
        DO i = i1, nqnyh_com
          Yh1(i) = Yh1(i) + Yh1(i+Nyh)
        END DO
      END DO
      ksteps_com = ksteps_com + 1
      do_prediction = .FALSE.
      do_corrector = .TRUE.
    ELSE
      do_corrector = .FALSE.
    END IF

    ! Corrector iteration
    corrector_outer: DO WHILE( do_corrector )
      m = 0
      Y(1:n_com) = Yh(1:n_com, 1)
      CALL DF(tn_com, Y, Savf)
      nfe_com = nfe_com + 1

      IF( ipup_com > 0 ) THEN
        ipup_com = 0
        rc_com = 1._DP
        nstepj_com = nst_com
        crate_com = 0.7_DP
        CALL DPJAC(Neq, Y, Yh, Nyh, Ewt, Acor, Savf, Wm, Iwm, DF, DJAC)
        IF( ier_com /= 0 ) THEN
          ! Jacobian evaluation failed
          tn_com = told
          ncf = ncf + 1
          rmax_com = 2._DP
          i1 = nqnyh_com + 1
          DO jb = 1, nq_com
            i1 = i1 - Nyh
            DO i = i1, nqnyh_com
              Yh1(i) = Yh1(i) - Yh1(i+Nyh)
            END DO
          END DO
          IF( ABS(h_com) <= hmin_com * 1.00001_DP ) THEN
            kflag_com = -2
            do_exit = .TRUE.
            EXIT corrector_outer
          ELSEIF( ncf /= 10 ) THEN
            rh = 0.25_DP
            ipup_com = miter_com
            iredo = 1
            rh = MAX(rh, hmin_com / ABS(h_com))
            do_rescale = .TRUE.
            EXIT corrector_outer
          ELSE
            kflag_com = -2
            do_exit = .TRUE.
            EXIT corrector_outer
          END IF
        END IF
      END IF

      Acor(1:n_com) = 0._DP

      ! Inner corrector loop
      corrector_loop: DO
        IF( miter_com /= 0 ) THEN
          ! Chord method
          DO i = 1, n_com
            Y(i) = h_com * Savf(i) - (Yh(i,2) + Acor(i))
          END DO
          CALL DSLVS(Wm, Iwm, Y)
          IF( ier_com /= 0 ) THEN
            ! Linear solve failed - check if Jacobian update needed
            IF( ipup_com /= 0 ) THEN
              ipup_com = miter_com
              CYCLE corrector_outer
            END IF
            ! Retract YH and try with smaller H
            tn_com = told
            ncf = ncf + 1
            rmax_com = 2._DP
            i1 = nqnyh_com + 1
            DO jb = 1, nq_com
              i1 = i1 - Nyh
              DO i = i1, nqnyh_com
                Yh1(i) = Yh1(i) - Yh1(i+Nyh)
              END DO
            END DO
            IF( ABS(h_com) <= hmin_com * 1.00001_DP ) THEN
              kflag_com = -2
              do_exit = .TRUE.
              EXIT corrector_outer
            ELSEIF( ncf /= 10 ) THEN
              rh = 0.25_DP
              ipup_com = miter_com
              iredo = 1
              rh = MAX(rh, hmin_com / ABS(h_com))
              do_rescale = .TRUE.
              EXIT corrector_outer
            ELSE
              kflag_com = -2
              do_exit = .TRUE.
              EXIT corrector_outer
            END IF
          END IF
          del = DVNRMS(n_com, Y, Ewt)
          DO i = 1, n_com
            Acor(i) = Acor(i) + Y(i)
            Y(i) = Yh(i,1) + el_com(1) * Acor(i)
          END DO
        ELSE
          ! Functional iteration
          DO i = 1, n_com
            Savf(i) = h_com * Savf(i) - Yh(i,2)
            Y(i) = Savf(i) - Acor(i)
          END DO
          del = DVNRMS(n_com, Y, Ewt)
          DO i = 1, n_com
            Y(i) = Yh(i,1) + el_com(1) * Savf(i)
            Acor(i) = Savf(i)
          END DO
        END IF

        ! Test for convergence
        IF( m /= 0 ) crate_com = MAX(0.2_DP * crate_com, del / delp)
        dcon = del * MIN(1._DP, 1.5_DP * crate_com) / (tesco_com(2, nq_com) * conit_com)

        IF( dcon <= 1._DP ) THEN
          ! Corrector converged
          IF( miter_com /= 0 ) ipup_com = -1
          IF( m == 0 ) dsm = del / tesco_com(2, nq_com)
          IF( m > 0 ) dsm = DVNRMS(n_com, Acor, Ewt) / tesco_com(2, nq_com)

          IF( dsm > 1._DP ) THEN
            ! Error test failed
            kflag_com = kflag_com - 1
            tn_com = told
            i1 = nqnyh_com + 1
            DO jb = 1, nq_com
              i1 = i1 - Nyh
              DO i = i1, nqnyh_com
                Yh1(i) = Yh1(i) - Yh1(i+Nyh)
              END DO
            END DO
            rmax_com = 2._DP

            IF( ABS(h_com) <= hmin_com * 1.00001_DP ) THEN
              kflag_com = -1
              do_exit = .TRUE.
              EXIT corrector_outer
            ELSEIF( kflag_com > -3 ) THEN
              iredo = 2
              rhup = 0._DP
              ! Compute step size factors and possibly change order
              CALL compute_rh_factors()
              EXIT corrector_outer
            ELSEIF( kflag_com /= -10 ) THEN
              ! 3 or more failures - reduce to order 1
              rh = 0.1_DP
              rh = MAX(hmin_com / ABS(h_com), rh)
              h_com = h_com * rh
              Y(1:n_com) = Yh(1:n_com, 1)
              CALL DF(tn_com, Y, Savf)
              nfe_com = nfe_com + 1
              DO i = 1, n_com
                Yh(i,2) = h_com * Savf(i)
              END DO
              ipup_com = miter_com
              ialth_com = 5
              IF( nq_com == 1 ) THEN
                do_prediction = .TRUE.
                EXIT corrector_outer
              END IF
              nq_com = 1
              l_com = 2
              istate = 3
              update_coefficients = .TRUE.
              EXIT corrector_outer
            ELSE
              kflag_com = -1
              do_exit = .TRUE.
              EXIT corrector_outer
            END IF
          ELSE
            ! Successful step
            kflag_com = 0
            iredo = 0
            nst_com = nst_com + 1
            hu_com = h_com
            nqu_com = nq_com
            DO j = 1, l_com
              DO i = 1, n_com
                Yh(i,j) = Yh(i,j) + el_com(j) * Acor(i)
              END DO
            END DO
            ialth_com = ialth_com - 1

            IF( ialth_com /= 0 ) THEN
              IF( ialth_com <= 1 .AND. l_com /= lmax_com ) THEN
                Yh(1:n_com, lmax_com) = Acor(1:n_com)
              END IF
              r = 1._DP / tesco_com(2, nqu_com)
              Acor(1:n_com) = Acor(1:n_com) * r
              do_exit = .TRUE.
              EXIT corrector_outer
            ELSE
              ! Consider order change
              rhup = 0._DP
              IF( l_com /= lmax_com ) THEN
                Savf(1:n_com) = Acor(1:n_com) - Yh(1:n_com, lmax_com)
                dup = DVNRMS(n_com, Savf, Ewt) / tesco_com(3, nq_com)
                exup = 1._DP / (l_com + 1)
                rhup = 1._DP / (1.4_DP * dup**exup + 0.0000014_DP)
              END IF
              CALL compute_rh_factors()
              EXIT corrector_outer
            END IF
          END IF
        ELSE
          ! Not converged yet
          m = m + 1
          IF( m == 3 ) THEN
            ! 3 iterations failed
            IF( ipup_com /= 0 ) THEN
              ipup_com = miter_com
              CYCLE corrector_outer
            END IF
            ! Retract and reduce step
            tn_com = told
            ncf = ncf + 1
            rmax_com = 2._DP
            i1 = nqnyh_com + 1
            DO jb = 1, nq_com
              i1 = i1 - Nyh
              DO i = i1, nqnyh_com
                Yh1(i) = Yh1(i) - Yh1(i+Nyh)
              END DO
            END DO
            IF( ABS(h_com) <= hmin_com * 1.00001_DP ) THEN
              kflag_com = -2
              do_exit = .TRUE.
              EXIT corrector_outer
            ELSEIF( ncf /= 10 ) THEN
              rh = 0.25_DP
              ipup_com = miter_com
              iredo = 1
              rh = MAX(rh, hmin_com / ABS(h_com))
              do_rescale = .TRUE.
              EXIT corrector_outer
            ELSE
              kflag_com = -2
              do_exit = .TRUE.
              EXIT corrector_outer
            END IF
          ELSEIF( m >= 2 .AND. del > 2._DP * delp ) THEN
            ! Diverging
            IF( ipup_com /= 0 ) THEN
              ipup_com = miter_com
              CYCLE corrector_outer
            END IF
            tn_com = told
            ncf = ncf + 1
            rmax_com = 2._DP
            i1 = nqnyh_com + 1
            DO jb = 1, nq_com
              i1 = i1 - Nyh
              DO i = i1, nqnyh_com
                Yh1(i) = Yh1(i) - Yh1(i+Nyh)
              END DO
            END DO
            IF( ABS(h_com) <= hmin_com * 1.00001_DP ) THEN
              kflag_com = -2
              do_exit = .TRUE.
              EXIT corrector_outer
            ELSEIF( ncf /= 10 ) THEN
              rh = 0.25_DP
              ipup_com = miter_com
              iredo = 1
              rh = MAX(rh, hmin_com / ABS(h_com))
              do_rescale = .TRUE.
              EXIT corrector_outer
            ELSE
              kflag_com = -2
              do_exit = .TRUE.
              EXIT corrector_outer
            END IF
          ELSE
            ! Continue iteration
            delp = del
            CALL DF(tn_com, Y, Savf)
            nfe_com = nfe_com + 1
            CYCLE corrector_loop
          END IF
        END IF
      END DO corrector_loop

      do_corrector = .FALSE.
    END DO corrector_outer

  END DO main_loop

  ! Final exit
  hold_com = h_com
  jstart_com = 1

CONTAINS

  SUBROUTINE compute_rh_factors()
    ! Compute step size factors for order change decisions
    exsm = 1._DP / l_com
    rhsm = 1._DP / (1.2_DP * dsm**exsm + 0.0000012_DP)
    rhdn = 0._DP
    IF( nq_com /= 1 ) THEN
      ddn = DVNRMS(n_com, Yh(:,l_com), Ewt) / tesco_com(1, nq_com)
      exdn = 1._DP / nq_com
      rhdn = 1._DP / (1.3_DP * ddn**exdn + 0.0000013_DP)
    END IF

    change_order = .FALSE.

    IF( rhsm >= rhup ) THEN
      IF( rhsm >= rhdn ) THEN
        ! Keep same order
        newq = nq_com
        rh = rhsm
      ELSE
        ! Decrease order
        newq = nq_com - 1
        rh = rhdn
        IF( kflag_com < 0 .AND. rh > 1._DP ) rh = 1._DP
      END IF
    ELSEIF( rhup > rhdn ) THEN
      ! Increase order
      newq = l_com
      rh = rhup
      IF( rh < 1.1_DP ) THEN
        ialth_com = 3
        r = 1._DP / tesco_com(2, nqu_com)
        Acor(1:n_com) = Acor(1:n_com) * r
        do_exit = .TRUE.
        RETURN
      END IF
      r = el_com(l_com) / l_com
      Yh(1:n_com, newq+1) = Acor(1:n_com) * r
      change_order = .TRUE.
    ELSE
      ! Decrease order
      newq = nq_com - 1
      rh = rhdn
      IF( kflag_com < 0 .AND. rh > 1._DP ) rh = 1._DP
    END IF

    IF( kflag_com == 0 .AND. rh < 1.1_DP .AND. .NOT. change_order ) THEN
      ialth_com = 3
      r = 1._DP / tesco_com(2, nqu_com)
      Acor(1:n_com) = Acor(1:n_com) * r
      do_exit = .TRUE.
      RETURN
    END IF

    IF( kflag_com <= -2 ) rh = MIN(rh, 0.2_DP)

    IF( newq == nq_com ) THEN
      rh = MAX(rh, hmin_com / ABS(h_com))
      do_rescale = .TRUE.
    ELSE
      nq_com = newq
      l_com = nq_com + 1
      istate = 2
      update_coefficients = .TRUE.
    END IF
  END SUBROUTINE compute_rh_factors

END SUBROUTINE DSTOD
