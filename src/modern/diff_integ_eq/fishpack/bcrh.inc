!** BCRH
REAL(SP) PURE FUNCTION BCRH(Xll,Xrr,Iz,C,A,Bh,F,Sgn)
  !> Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (BCRH-S, BSRH-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  CBLKTR
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    CCBLK

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S8.1.6.2 (DO WHILE construct)
  !           Original: Swarztrauber, NCAR (FISHPACK)
  USE CCBLK, ONLY : cnv_com
  !
  INTERFACE
    REAL(SP) PURE FUNCTION F(X,Iz,C,A,Bh)
      IMPORT SP
      INTEGER, INTENT(IN) :: Iz
      REAL(SP), INTENT(IN) :: X, A(Iz), Bh(Iz), C(Iz)
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(IN) :: Iz
  REAL(SP), INTENT(IN) :: A(:), Bh(:), C(:)
  REAL(SP), INTENT(IN) :: Sgn, Xll, Xrr
  !
  REAL(SP) :: dx, x, xl, xr
  !* FIRST EXECUTABLE STATEMENT  BCRH
  xl = Xll
  xr = Xrr
  dx = 0.5_SP*ABS(xr-xl)
  ! Bisection iteration (was backward GOTO 100)
  DO WHILE( dx > cnv_com )
    x = 0.5_SP*(xl+xr)
    IF( Sgn*F(x,Iz,C,A,Bh) < 0 ) THEN
      xl = x
    ELSEIF( Sgn*F(x,Iz,C,A,Bh) == 0 ) THEN
      BCRH = 0.5_SP*(xl+xr)
      RETURN
    ELSE
      xr = x
    END IF
    dx = 0.5_SP*dx
  END DO
  BCRH = 0.5_SP*(xl+xr)
  !
  RETURN
END FUNCTION BCRH