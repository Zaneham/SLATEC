!** DXADD
SUBROUTINE DXADD(X,Ix,Y,Iy,Z,Iz,Ierror)
  !> To provide single-precision floating-point arithmetic
  !            with an extended exponent range.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  A3D
  !***
  ! **Type:**      DOUBLE PRECISION (DXADD-S, DDXADD-D)
  !***
  ! **Keywords:**  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC
  !***
  ! **Author:**  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***
  ! **Description:**
  !     REAL X, Y, Z
  !     INTEGER IX, IY, IZ
  !
  !                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) =
  !                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED
  !                  BEFORE RETURNING. THE INPUT OPERANDS
  !                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR
  !                  PRINCIPAL PARTS MUST SATISFY
  !                  RADIX**(-2L)<=ABS(X)<=RADIX**(2L),
  !                  RADIX**(-2L)<=ABS(Y)<=RADIX**(2L).
  !
  !***
  ! **See also:**  XSET
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  DXADJ
  !***
  ! COMMON BLOCKS    XBLK2

  !* REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  !   2024xx  Modernized: eliminated 11 GOTOs using structured control flow
  USE XBLK ,ONLY: radixx_com, radixl_com, rad2l_com, l_com
  INTEGER :: Ierror, Ix, Iy, Iz
  REAL(DP) :: X, Y, Z
  INTEGER :: i, i1, i2, is, j
  REAL(DP) :: s, t
  !
  ! Control flow flags (replace GOTOs)
  LOGICAL :: do_label_300, do_label_400, return_s_is
  !
  !
  !   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE
  ! ARE
  !     (1) 1 < L <= 0.5*LOGR(0.5*DZERO)
  !
  !     (2) NRADPL < L <= KMAX/6
  !
  !     (3) KMAX <= (2**NBITS - 4*L - 1)/2
  !
  ! THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
  ! IN SUBROUTINE XSET.
  !
  !* FIRST EXECUTABLE STATEMENT  DXADD
  Ierror = 0
  IF( X==0._DP ) THEN
    Z = Y
    Iz = Iy
    CALL DXADJ(Z,Iz,Ierror)
    RETURN
  ELSEIF( Y/=0._DP ) THEN
    IF( Ix<0 .OR. Iy<0 ) THEN
      IF( Ix>=0 .OR. Iy>=0 ) THEN
        IF( ABS(Ix)>6*l_com .OR. ABS(Iy)>6*l_com ) THEN
          IF( Ix>=0 ) THEN
            Z = X
            Iz = Ix
          ELSE
            Z = Y
            Iz = Iy
          END IF
          CALL DXADJ(Z,Iz,Ierror)
          RETURN
        END IF
      END IF
    END IF
    i = Ix - Iy
    IF( i<0 ) THEN
      s = Y
      is = Iy
      t = X
    ELSEIF( i==0 ) THEN
      IF( ABS(X)>1._DP .AND. ABS(Y)>1._DP ) THEN
        s = X/radixl_com
        t = Y/radixl_com
        Z = s + t
        Iz = Ix + l_com
      ELSEIF( ABS(X)<1._DP .AND. ABS(Y)<1._DP ) THEN
        s = X*radixl_com
        t = Y*radixl_com
        Z = s + t
        Iz = Ix - l_com
      ELSE
        Z = X + Y
        Iz = Ix
      END IF
      CALL DXADJ(Z,Iz,Ierror)
      RETURN
    ELSE
      s = X
      is = Ix
      t = Y
    END IF
    !
    !  AT THIS POINT, THE ONE OF (X,IX) OR (Y,IY) THAT HAS THE
    ! LARGER AUXILIARY INDEX IS STORED IN (S,IS). THE PRINCIPAL
    ! PART OF THE OTHER INPUT IS STORED IN T.
    !
    i1 = ABS(i)/l_com
    i2 = MOD(ABS(i),l_com)

    ! Initialize control flow flags
    do_label_300 = .FALSE.
    do_label_400 = .FALSE.
    return_s_is = .FALSE.

    IF( ABS(t)>=radixl_com ) THEN
      j = i1 - 2
      IF( j>=0 ) THEN
        t = t*radixx_com**(-i2)/rad2l_com
        do_label_300 = .TRUE.
      END IF
    ELSEIF( ABS(t)>=1._DP ) THEN
      ! Was: GOTO 200
      j = i1 - 1
      IF( j>=0 ) THEN
        t = t*radixx_com**(-i2)/radixl_com
        do_label_300 = .TRUE.
      END IF
    ELSEIF( radixl_com*ABS(t)<1._DP ) THEN
      j = i1 + 1
      t = t*radixx_com**(l_com-i2)
      do_label_300 = .TRUE.
    END IF

    ! Was label 200 path (j<0 cases fall through to label 100)
    IF( .NOT. do_label_300 ) THEN
      j = i1 - 1
      IF( j>=0 ) THEN
        t = t*radixx_com**(-i2)/radixl_com
        do_label_300 = .TRUE.
      ELSE
        ! Was label 100
        j = i1
        t = t*radixx_com**(-i2)
        do_label_300 = .TRUE.
      END IF
    END IF

    !
    !  AT THIS POINT, SOME OR ALL OF THE DIFFERENCE IN THE
    ! AUXILIARY INDICES HAS BEEN USED TO EFFECT A LEFT SHIFT
    ! OF T.  THE SHIFTED VALUE OF T SATISFIES
    !
    !       RADIX**(-2*L) <= ABS(T) <= 1.0
    !
    ! AND, IF J=0, NO FURTHER SHIFTING REMAINS TO BE DONE.
    !
    ! Was label 300
    IF( do_label_300 ) THEN
      IF( j/=0 ) THEN
        IF( ABS(s)<radixl_com .AND. j<=3 ) THEN
          IF( ABS(s)>=1._DP ) THEN
            SELECT CASE (j)
              CASE (1)
                s = s*radixl_com
                do_label_400 = .TRUE.
              CASE (2,3)
                return_s_is = .TRUE.  ! Was: GOTO 350
              CASE DEFAULT
            END SELECT
          END IF
          IF( .NOT. do_label_400 .AND. .NOT. return_s_is ) THEN
            IF( radixl_com*ABS(s)>=1._DP ) THEN
              SELECT CASE (j)
                CASE (1)
                  s = s*radixl_com
                  do_label_400 = .TRUE.
                CASE (2)
                  s = s*radixl_com
                  s = s*radixl_com
                  do_label_400 = .TRUE.
                CASE (3)
                  return_s_is = .TRUE.  ! Was: GOTO 350
                CASE DEFAULT
                  ! Fall through to label 320
              END SELECT
            END IF
          END IF
          IF( .NOT. do_label_400 .AND. .NOT. return_s_is ) THEN
            ! Was label 320
            SELECT CASE (j)
              CASE (1)
                s = s*radixl_com
                do_label_400 = .TRUE.
              CASE (2)
                s = s*radixl_com
                s = s*radixl_com
                do_label_400 = .TRUE.
              CASE (3)
                s = s*radixl_com
                s = s*radixl_com
                s = s*radixl_com
                do_label_400 = .TRUE.
              CASE DEFAULT
                return_s_is = .TRUE.  ! Was: GOTO 350
            END SELECT
          END IF
        ELSE
          return_s_is = .TRUE.  ! Was: GOTO 350
        END IF
      ELSE
        do_label_400 = .TRUE.
      END IF

      ! Was label 350 - early return
      IF( return_s_is ) THEN
        Z = s
        Iz = is
        CALL DXADJ(Z,Iz,Ierror)
        RETURN
      END IF
    END IF

    !
    !   AT THIS POINT, THE REMAINING DIFFERENCE IN THE
    ! AUXILIARY INDICES HAS BEEN USED TO EFFECT A RIGHT SHIFT
    ! OF S.  IF THE SHIFTED VALUE OF S WOULD HAVE EXCEEDED
    ! RADIX**L, THEN (S,IS) IS RETURNED AS THE VALUE OF THE
    ! SUM.
    !
    ! Was label 400
    IF( do_label_400 ) THEN
      IF( ABS(s)>1._DP .AND. ABS(t)>1._DP ) THEN
        s = s/radixl_com
        t = t/radixl_com
        Z = s + t
        Iz = is - j*l_com + l_com
      ELSEIF( ABS(s)<1._DP .AND. ABS(t)<1._DP ) THEN
        s = s*radixl_com
        t = t*radixl_com
        Z = s + t
        Iz = is - j*l_com - l_com
      ELSE
        Z = s + t
        Iz = is - j*l_com
      END IF
      CALL DXADJ(Z,Iz,Ierror)
      RETURN
    END IF
  ELSE
    Z = X
    Iz = Ix
    CALL DXADJ(Z,Iz,Ierror)
    RETURN
  END IF
END SUBROUTINE DXADD
