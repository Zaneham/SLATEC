!** DXADJ
SUBROUTINE DXADJ(X,Ix,Ierror)
  !> To provide double-precision floating-point arithmetic
  !            with an extended exponent range.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  A3D
  !***
  ! **Type:**      DOUBLE PRECISION (XADJ-S, DXADJ-D)
  !***
  ! **Keywords:**  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
  !***
  ! **Author:**  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***
  ! **Description:**
  !     DOUBLE PRECISION X
  !     INTEGER IX
  !
  !                  TRANSFORMS (X,IX) SO THAT
  !                  RADIX**(-L) <= ABS(X) < RADIX**L.
  !                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
  !                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
  !                  THE NUMBER BASE OF DOUBLE-PRECISION ARITHMETIC.
  !
  !***
  ! **See also:**  DXSET
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XERMSG
  !***
  ! COMMON BLOCKS    DXBLK2

  !* REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.8 (IF construct)
  !           Original: Lozier, Smith (NBS)
  USE DXBLK ,ONLY: radixl_com, rad2l_com, l2_com, kmax_com

  INTEGER :: Ierror, Ix
  REAL(DP) :: X
  !
  !   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
  ! IS
  !     2*L <= KMAX
  !
  ! THIS CONDITION MUST BE MET BY APPROPRIATE CODING
  ! IN SUBROUTINE DXSET.
  !
  !* FIRST EXECUTABLE STATEMENT  DXADJ
  Ierror = 0
  IF( X==0._DP ) THEN
    Ix = 0
  ELSEIF( ABS(X)>=1._DP ) THEN
    IF( ABS(X)>=radixl_com ) THEN
      X = X/rad2l_com
      IF( Ix<=0 .OR. Ix<=kmax_com-l2_com ) THEN
        Ix = Ix + l2_com
        RETURN
      ELSE
        ! Overflow in adjustment (was GOTO 100)
        ERROR STOP 'DXADJ : overflow in auxiliary index'
        Ierror = 207
        RETURN
      END IF
    END IF
  ELSE
    IF( radixl_com*ABS(X)<1._DP ) THEN
      X = X*rad2l_com
      IF( Ix>=0 .OR. Ix>=-kmax_com+l2_com ) THEN
        Ix = Ix - l2_com
        RETURN
      ELSE
        ! Overflow in adjustment (was GOTO 100)
        ERROR STOP 'DXADJ : overflow in auxiliary index'
        Ierror = 207
        RETURN
      END IF
    END IF
  END IF
  ! Final overflow check (was label 200)
  IF( ABS(Ix)>kmax_com ) THEN
    ERROR STOP 'DXADJ : overflow in auxiliary index'
    Ierror = 207
  END IF
END SUBROUTINE DXADJ
