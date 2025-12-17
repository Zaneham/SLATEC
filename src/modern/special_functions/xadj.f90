!** XADJ
SUBROUTINE XADJ(X,Ix,Ierror)
  !> To provide single-precision floating-point arithmetic
  !            with an extended exponent range.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  A3D
  !***
  ! **Type:**      SINGLE PRECISION (XADJ-S, DXADJ-D)
  !***
  ! **Keywords:**  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC
  !***
  ! **Author:**  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***
  ! **Description:**
  !     REAL X
  !     INTEGER IX
  !
  !                  TRANSFORMS (X,IX) SO THAT
  !                  RADIX**(-L) <= ABS(X) < RADIX**L.
  !                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
  !                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
  !                  THE NUMBER BASE OF SINGLE-PRECISION ARITHMETIC.
  !
  !***
  ! **See also:**  XSET
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XERMSG
  !***
  ! COMMON BLOCKS    XBLK2

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
  USE XBLK ,ONLY: radixl_com, rad2l_com, l2_com, kmax_com

  INTEGER :: Ierror, Ix
  REAL(SP) :: X
  !
  !   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
  ! IS
  !     2*L <= KMAX
  !
  ! THIS CONDITION MUST BE MET BY APPROPRIATE CODING
  ! IN SUBROUTINE XSET.
  !
  !* FIRST EXECUTABLE STATEMENT  XADJ
  Ierror = 0
  IF( X==0._SP ) THEN
    Ix = 0
  ELSEIF( ABS(X)>=1._SP ) THEN
    IF( ABS(X)>=radixl_com ) THEN
      X = X/rad2l_com
      IF( Ix<=0 .OR. Ix<=kmax_com-l2_com ) THEN
        Ix = Ix + l2_com
        RETURN
      ELSE
        ! Overflow in adjustment (was GOTO 100)
        ERROR STOP 'XADJ : overflow in auxiliary index'
        Ierror = 107
        RETURN
      END IF
    END IF
  ELSE
    IF( radixl_com*ABS(X)<1._SP ) THEN
      X = X*rad2l_com
      IF( Ix>=0 .OR. Ix>=-kmax_com+l2_com ) THEN
        Ix = Ix - l2_com
        RETURN
      ELSE
        ! Overflow in adjustment (was GOTO 100)
        ERROR STOP 'XADJ : overflow in auxiliary index'
        Ierror = 107
        RETURN
      END IF
    END IF
  END IF
  ! Final overflow check (was label 200)
  IF( ABS(Ix)>kmax_com ) THEN
    ERROR STOP 'XADJ : overflow in auxiliary index'
    Ierror = 107
  END IF
END SUBROUTINE XADJ
