!** CBINU
PURE SUBROUTINE CBINU(Z,Fnu,Kode,N,Cy,Nz,Rl,Fnul,Tol,Elim,Alim)
  !> Subsidiary to ZAIRY, ZBESH, CBESI, ZBESJ, ZBESK and ZBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CBINU-A, CBINU-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
  !
  !***
  ! **See also:**  ZAIRY, ZBESH, CBESI, ZBESJ, ZBESK, ZBIRY
  !***
  ! **Routines called:**  ZABS, CASYI, CBUNI, CMLRI, CSERI, CUOIK, CWRSK
  !
  ! REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   2024xx  Modernized: eliminated 16 GOTOs using structured control flow

  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(SP), INTENT(IN) :: Alim, Elim, Fnu, Fnul, Rl, Tol
  COMPLEX(SP), INTENT(IN) :: Z
  COMPLEX(SP), INTENT(OUT) :: Cy(N)
  !
  INTEGER :: i, inw, nlast, nn, nui, nw
  REAL(SP) :: az, dfnu
  COMPLEX(SP) :: cw(2)
  ! Control flow
  logical :: skip_to_100, need_zbuni
  !* FIRST EXECUTABLE STATEMENT  CBINU
  Nz = 0
  az = ABS(Z)
  nn = N
  dfnu = Fnu + (N-1)

  ! Determine if we should skip power series
  skip_to_100 = .false.
  IF( az>2._SP ) THEN
    IF( az*az*0.25_SP>dfnu+1._SP ) skip_to_100 = .true.
  END IF

  IF( .not. skip_to_100 ) THEN
    !-----------------------------------------------------------------------
    !     POWER SERIES
    !-----------------------------------------------------------------------
    CALL CSERI(Z,Fnu,Kode,nn,Cy,nw,Tol,Elim,Alim)
    inw = ABS(nw)
    Nz = Nz + inw
    nn = nn - inw
    IF( nn==0 ) RETURN
    IF( nw>=0 ) RETURN  ! Success return
    dfnu = Fnu + (nn-1)
  END IF

  ! Label 100 equivalent
  IF( az>=Rl ) THEN
    IF( dfnu<=1._SP .OR. az+az>=dfnu*dfnu ) THEN
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR LARGE Z
      !-----------------------------------------------------------------------
      CALL CASYI(Z,Fnu,Kode,nn,Cy,nw,Rl,Tol,Elim,Alim)
      IF( nw>=0 ) RETURN  ! Success
      Nz = -1
      IF( nw==(-2) ) Nz = -2
      RETURN
    END IF
  ELSEIF( dfnu<=1._SP ) THEN
    !-----------------------------------------------------------------------
    !     MILLER ALGORITHM NORMALIZED BY THE SERIES (label 400)
    !-----------------------------------------------------------------------
    CALL CMLRI(Z,Fnu,Kode,nn,Cy,nw,Tol)
    IF( nw>=0 ) RETURN  ! Success
    Nz = -1
    IF( nw==(-2) ) Nz = -2
    RETURN
  END IF

  !-----------------------------------------------------------------------
  !     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
  !-----------------------------------------------------------------------
  ! Label 200 equivalent
  CALL CUOIK(Z,Fnu,Kode,1,nn,Cy,nw,Tol,Elim,Alim)
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
    CALL CBUNI(Z,Fnu,Kode,nn,Cy,nw,nui,nlast,Fnul,Tol,Elim,Alim)
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
    CALL CUOIK(Z,Fnu,Kode,2,2,cw,nw,Tol,Elim,Alim)
    IF( nw>=0 ) THEN
      IF( nw>0 ) THEN
        Nz = -1
        IF( nw==(-2) ) Nz = -2
        RETURN
      END IF
      CALL CWRSK(Z,Fnu,Kode,nn,Cy,nw,cw,Tol,Elim,Alim)
      IF( nw>=0 ) RETURN  ! Success
      Nz = -1
      IF( nw==(-2) ) Nz = -2
      RETURN
    ELSE
      Nz = nn
      DO i = 1, nn
        Cy(i) = (0._SP,0._SP)
      END DO
      RETURN
    END IF
  ELSE
    !-----------------------------------------------------------------------
    !     MILLER ALGORITHM NORMALIZED BY THE SERIES (label 400)
    !-----------------------------------------------------------------------
    CALL CMLRI(Z,Fnu,Kode,nn,Cy,nw,Tol)
    IF( nw>=0 ) RETURN  ! Success
    Nz = -1
    IF( nw==(-2) ) Nz = -2
    RETURN
  END IF
  !
END SUBROUTINE CBINU
