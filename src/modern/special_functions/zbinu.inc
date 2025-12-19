!** ZBINU
PURE SUBROUTINE ZBINU(Z,Fnu,Kode,N,Cy,Nz,Rl,Fnul,Tol,Elim,Alim)
  !> Subsidiary to ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK and ZBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CBINU-A, ZBINU-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
  !
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBIRY
  !***
  ! **Routines called:**  ZABS, ZASYI, ZBUNI, ZMLRI, ZSERI, ZUOIK, ZWRSK
  !
  ! REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   2024xx  Modernized: eliminated 16 GOTOs using structured control flow

  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Fnul, Rl, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Cy(N)
  !
  INTEGER :: i, inw, nlast, nn, nui, nw
  REAL(DP) :: az, dfnu
  COMPLEX(DP) :: cw(2)
  ! Control flow
  logical :: skip_to_100, need_zbuni
  !* FIRST EXECUTABLE STATEMENT  ZBINU
  Nz = 0
  az = ABS(Z)
  nn = N
  dfnu = Fnu + (N-1)

  ! Determine if we should skip power series
  skip_to_100 = .false.
  IF( az>2._DP ) THEN
    IF( az*az*0.25_DP>dfnu+1._DP ) skip_to_100 = .true.
  END IF

  IF( .not. skip_to_100 ) THEN
    !-----------------------------------------------------------------------
    !     POWER SERIES
    !-----------------------------------------------------------------------
    CALL ZSERI(Z,Fnu,Kode,nn,Cy,nw,Tol,Elim,Alim)
    inw = ABS(nw)
    Nz = Nz + inw
    nn = nn - inw
    IF( nn==0 ) RETURN
    IF( nw>=0 ) RETURN  ! Success return
    dfnu = Fnu + (nn-1)
  END IF

  ! Label 100 equivalent
  IF( az>=Rl ) THEN
    IF( dfnu<=1._DP .OR. az+az>=dfnu*dfnu ) THEN
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR LARGE Z
      !-----------------------------------------------------------------------
      CALL ZASYI(Z,Fnu,Kode,nn,Cy,nw,Rl,Tol,Elim,Alim)
      IF( nw>=0 ) RETURN  ! Success
      Nz = -1
      IF( nw==(-2) ) Nz = -2
      RETURN
    END IF
  ELSEIF( dfnu<=1._DP ) THEN
    !-----------------------------------------------------------------------
    !     MILLER ALGORITHM NORMALIZED BY THE SERIES (label 400)
    !-----------------------------------------------------------------------
    CALL ZMLRI(Z,Fnu,Kode,nn,Cy,nw,Tol)
    IF( nw>=0 ) RETURN  ! Success
    Nz = -1
    IF( nw==(-2) ) Nz = -2
    RETURN
  END IF

  !-----------------------------------------------------------------------
  !     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
  !-----------------------------------------------------------------------
  ! Label 200 equivalent
  CALL ZUOIK(Z,Fnu,Kode,1,nn,Cy,nw,Tol,Elim,Alim)
  IF( nw<0 ) THEN
    Nz = -1
    IF( nw==(-2) ) Nz = -2
    RETURN
  END IF
  Nz = Nz + nw
  nn = nn - nw
  IF( nn==0 ) RETURN
  dfnu = Fnu + (nn-1)

  need_zbuni = (dfnu>Fnul .OR. az>Fnul)

  IF( need_zbuni ) THEN
    !-----------------------------------------------------------------------
    !     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
    !-----------------------------------------------------------------------
    ! Label 500 equivalent
    nui = INT( Fnul - dfnu ) + 1
    nui = MAX(nui,0)
    CALL ZBUNI(Z,Fnu,Kode,nn,Cy,nw,nui,nlast,Fnul,Tol,Elim,Alim)
    IF( nw<0 ) THEN
      Nz = -1
      IF( nw==(-2) ) Nz = -2
      RETURN
    END IF
    Nz = Nz + nw
    IF( nlast==0 ) RETURN  ! Success
    nn = nlast
    ! Fall through to label 300 equivalent
  END IF

  ! Label 300 equivalent
  IF( az>Rl ) THEN
    !-----------------------------------------------------------------------
    !     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
    !-----------------------------------------------------------------------
    !     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
    !-----------------------------------------------------------------------
    CALL ZUOIK(Z,Fnu,Kode,2,2,cw,nw,Tol,Elim,Alim)
    IF( nw>=0 ) THEN
      IF( nw>0 ) THEN
        Nz = -1
        IF( nw==(-2) ) Nz = -2
        RETURN
      END IF
      CALL ZWRSK(Z,Fnu,Kode,nn,Cy,nw,cw,Tol,Elim,Alim)
      IF( nw>=0 ) RETURN  ! Success
      Nz = -1
      IF( nw==(-2) ) Nz = -2
      RETURN
    ELSE
      Nz = nn
      DO i = 1, nn
        Cy(i) = (0._DP,0._DP)
      END DO
      RETURN
    END IF
  ELSE
    !-----------------------------------------------------------------------
    !     MILLER ALGORITHM NORMALIZED BY THE SERIES (label 400)
    !-----------------------------------------------------------------------
    CALL ZMLRI(Z,Fnu,Kode,nn,Cy,nw,Tol)
    IF( nw>=0 ) RETURN  ! Success
    Nz = -1
    IF( nw==(-2) ) Nz = -2
    RETURN
  END IF
  !
END SUBROUTINE ZBINU
