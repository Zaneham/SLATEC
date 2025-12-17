!** ZSERI
PURE SUBROUTINE ZSERI(Z,Fnu,Kode,N,Y,Nz,Tol,Elim,Alim)
  !> Subsidiary to ZBESI and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CSERI-A, ZSERI-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z)>=0.0 BY
  !     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE
  !     REGION ABS(Z)<=2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
  !     NZ>0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
  !     DUE TO UNDERFLOW. NZ<0 MEANS UNDERFLOW OCCURRED, BUT THE
  !     CONDITION ABS(Z)<=2*SQRT(FNU+1) WAS VIOLATED AND THE
  !     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
  !
  !***
  ! **See also:**  ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, DGAMLN, ZABS, ZDIV, ZLOG, ZMLT, ZUCHK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZLOG to EXTERNAL statement.  (RWC)
!   251217  Eliminated GOTOs per MODERNISATION_GUIDE.md S1. (ZH)
!           Ref: ISO/IEC 1539-1:2018 S11.1.7.4.4 (DO construct with state)
!           Original: Amos, D.E. (SNL)
  USE service, ONLY : tiny_dp
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, ib, iflag, il, k, l, m, nn, nw, istate
  COMPLEX(DP) :: ak1, ck, coef, crsc, cz, hz, rz, s1, s2, w(2)
  REAL(DP) :: aa, acz, ak, arm, ascle, atol, az, dfnu, &
    fnup, rak1, rs, rtr1, s, ss, x
  !* FIRST EXECUTABLE STATEMENT  ZSERI
  !  State machine: 1=main_loop, 2=underflow, 3=recurrence, 4=scale_adjust, 5=early_exit
  Nz = 0
  az = ABS(Z)
  IF( az==0._DP ) THEN
    ! Early exit for az==0 (was GOTO 500)
    Y(1) = (0._DP,0._DP)
    IF( Fnu==0._DP ) Y(1) = (1._DP,0._DP)
    IF( N==1 ) RETURN
    Y = (0._DP,0._DP)
    RETURN
  END IF
  x = REAL(Z,DP)
  arm = 1.E+3_DP*tiny_dp
  rtr1 = SQRT(arm)
  crsc = CMPLX(1._DP,0._DP,DP)
  iflag = 0
  IF( az<arm ) THEN
    Nz = N
    IF( Fnu==0._DP ) Nz = Nz - 1
    ! Early exit for small az (was GOTO 500)
    Y(1) = (0._DP,0._DP)
    IF( Fnu==0._DP ) Y(1) = (1._DP,0._DP)
    IF( N==1 ) RETURN
    Y = (0._DP,0._DP)
    RETURN
  END IF
  hz = Z*CMPLX(0.5_DP,0._DP,DP)
  cz = (0._DP,0._DP)
  IF( az>rtr1 ) cz = hz*hz
  acz = ABS(cz)
  nn = N
  ck = LOG(hz)
  !
  istate = 1
  main_loop: DO
    SELECT CASE (istate)
    !
    CASE (1)  ! Main loop (was label 100)
      dfnu = Fnu + (nn-1)
      fnup = dfnu + 1._DP
      !-----------------------------------------------------------------------
      !     UNDERFLOW TEST
      !-----------------------------------------------------------------------
      ak1 = ck*CMPLX(dfnu,0._DP,DP)
      ak = LOG_GAMMA(fnup)
      ak1 = ak1 - CMPLX(ak,0._DP,DP)
      IF( Kode==2 ) ak1 = ak1 - CMPLX(x,0._DP,DP)
      rak1 = REAL(ak1,DP)
      IF( rak1>(-Elim) ) THEN
        IF( rak1<=(-Alim) ) THEN
          iflag = 1
          ss = 1._DP/Tol
          crsc = CMPLX(Tol,0._DP,DP)
          ascle = arm*ss
        END IF
        ak = AIMAG(ak1)
        aa = EXP(rak1)
        IF( iflag==1 ) aa = aa*ss
        coef = CMPLX(aa,0._DP,DP)*CMPLX(COS(ak),SIN(ak),DP)
        atol = Tol*acz/fnup
        il = MIN(2,nn)
        DO i = 1, il
          dfnu = Fnu + (nn-i)
          fnup = dfnu + 1._DP
          s1 = (1._DP,0._DP)
          IF( acz>=Tol*fnup ) THEN
            ak1 = (1._DP,0._DP)
            ak = fnup + 2._DP
            s = fnup
            aa = 2._DP
            DO
              rs = 1._DP/s
              ak1 = ak1*cz*CMPLX(rs,0._DP,DP)
              s1 = s1 + ak1
              s = s + ak
              ak = ak + 2._DP
              aa = aa*acz*rs
              IF( aa<=atol ) EXIT
            END DO
          END IF
          m = nn - i + 1
          s2 = s1*coef
          w(i) = s2
          IF( iflag/=0 ) THEN
            CALL ZUCHK(s2,nw,ascle,Tol)
            IF( nw/=0 ) THEN
              istate = 2  ! Go to underflow handling
              CYCLE main_loop
            END IF
          END IF
          Y(m) = s2*crsc
          IF( i/=il ) coef = coef*CMPLX(dfnu,0._DP,DP)/hz
        END DO
        IF( nn<=2 ) RETURN
        k = nn - 2
        ak = k
        rz = (2._DP,0._DP)/Z
        IF( iflag==1 ) THEN
          !-----------------------------------------------------------------------
          !     RECUR BACKWARD WITH SCALED VALUES
          !-----------------------------------------------------------------------
          s1 = w(1)
          s2 = w(2)
          DO l = 3, nn
            ck = s2
            s2 = s1 + CMPLX(ak+Fnu,0._DP,DP)*rz*s2
            s1 = ck
            ck = s2*crsc
            Y(k) = ck
            ak = ak - 1._DP
            k = k - 1
            IF( ABS(ck)>ascle ) THEN
              ! Scale adjustment (was GOTO 400)
              ib = l + 1
              IF( ib>nn ) RETURN
              istate = 3  ! Go to recurrence
              CYCLE main_loop
            END IF
          END DO
          RETURN
        ELSE
          ib = 3
          istate = 3  ! Go to recurrence
        END IF
      ELSE
        istate = 2  ! Go to underflow handling
      END IF
    !
    CASE (2)  ! Underflow handling (was label 200)
      Nz = Nz + 1
      Y(nn) = (0._DP,0._DP)
      IF( acz>dfnu ) THEN
        !-----------------------------------------------------------------------
        !     RETURN WITH NZ<0 IF ABS(Z*Z/4)>FNU+N-NZ-1 COMPLETE
        !     THE CALCULATION IN CBINU WITH N=N-ABS(NZ)
        !-----------------------------------------------------------------------
        Nz = -Nz
        RETURN
      ELSE
        nn = nn - 1
        IF( nn==0 ) RETURN
        istate = 1  ! Back to main loop
      END IF
    !
    CASE (3)  ! Recurrence (was label 300)
      DO i = ib, nn
        Y(k) = CMPLX(ak+Fnu,0._DP,DP)*rz*Y(k+1) + Y(k+2)
        ak = ak - 1._DP
        k = k - 1
      END DO
      RETURN
    !
    END SELECT
  END DO main_loop
  !
  RETURN
END SUBROUTINE ZSERI