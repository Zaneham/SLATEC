!** PCHSP
PURE SUBROUTINE PCHSP(Ic,Vc,N,X,F,D,Incfd,Wk,Nwk,Ierr)
  !> Set derivatives needed to determine the Hermite representation of the cubic
  !  spline interpolant to given data, with specified boundary conditions.
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Category:**  E1A
  !***
  ! **Type:**      DOUBLE PRECISION (PCHSP-S, PCHSP-D)
  !***
  ! **Keywords:**  CUBIC HERMITE INTERPOLATION, PCHIP, PIECEWISE CUBIC INTERPOLATION,
  !                SPLINE INTERPOLATION
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             P.O. Box 808  (L-316)
  !             Livermore, CA  94550
  !             FTS 532-4275, (510) 422-4275
  !***
  ! **Description:**
  !
  !          PCHSP:   Piecewise Cubic Hermite Spline
  !
  !     Computes the Hermite representation of the cubic spline inter-
  !     polant to the data given in X and F satisfying the boundary
  !     conditions specified by IC and VC.
  !
  !***
  ! **References:**  Carl de Boor, A Practical Guide to Splines, Springer-
  !                 Verlag, New York, 1978, pp. 53-59.
  !***
  ! **Routines called:**  PCHDF, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   820503  DATE WRITTEN
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251220  Eliminated GOTOs using Jacob Williams pattern. (ZH)
  !           Ref: github.com/jacobwilliams/PCHIP

  INTEGER, INTENT(IN) :: Ic(2), N, Incfd, Nwk
  INTEGER, INTENT(OUT) :: Ierr
  REAL(SP), INTENT(IN) :: Vc(2), X(N), F(Incfd,N)
  REAL(SP), INTENT(OUT) :: D(Incfd,N), Wk(2,Nwk)
  !
  INTEGER :: ibeg, iend, indexx, j, nm1
  REAL(SP) :: g, stemp(3), xtemp(4)
  LOGICAL :: singular_error, x_not_increasing, skip_forward_pass
  !
  IF( N<2 ) THEN
    Ierr = -1
    ERROR STOP 'PCHSP : NUMBER OF DATA POINTS LESS THAN TWO'
    RETURN
  END IF
  !
  IF( Incfd<1 ) THEN
    Ierr = -2
    ERROR STOP 'PCHSP : INCREMENT LESS THAN ONE'
    RETURN
  END IF
  !
  x_not_increasing = .FALSE.
  DO j = 2, N
    IF( X(j)<=X(j-1) ) THEN
      x_not_increasing = .TRUE.
      EXIT
    END IF
  END DO
  !
  IF( x_not_increasing ) THEN
    Ierr = -3
    ERROR STOP 'PCHSP : X-ARRAY NOT STRICTLY INCREASING'
    RETURN
  END IF
  !
  ibeg = Ic(1)
  iend = Ic(2)
  Ierr = 0
  IF( (ibeg<0) .OR. (ibeg>4) ) Ierr = Ierr - 1
  IF( (iend<0) .OR. (iend>4) ) Ierr = Ierr - 2
  IF( Ierr<0 ) THEN
    Ierr = Ierr - 3
    ERROR STOP 'PCHSP : IC OUT OF RANGE'
    RETURN
  END IF
  !
  IF( Nwk<2*N ) THEN
    Ierr = -7
    ERROR STOP 'PCHSP : WORK ARRAY TOO SMALL'
    RETURN
  END IF
  !
  DO j = 2, N
    Wk(1,j) = X(j) - X(j-1)
    Wk(2,j) = (F(1,j)-F(1,j-1))/Wk(1,j)
  END DO
  !
  IF( ibeg>N ) ibeg = 0
  IF( iend>N ) iend = 0
  !
  IF( (ibeg==1) .OR. (ibeg==2) ) THEN
    D(1,1) = Vc(1)
  ELSEIF( ibeg>2 ) THEN
    DO j = 1, ibeg
      indexx = ibeg - j + 1
      xtemp(j) = X(indexx)
      IF( j<ibeg ) stemp(j) = Wk(2,indexx)
    END DO
    D(1,1) = PCHDF(ibeg,xtemp,stemp)
    ibeg = 1
  END IF
  !
  IF( (iend==1) .OR. (iend==2) ) THEN
    D(1,N) = Vc(2)
  ELSEIF( iend>2 ) THEN
    DO j = 1, iend
      indexx = N - iend + j
      xtemp(j) = X(indexx)
      IF( j<iend ) stemp(j) = Wk(2,indexx+1)
    END DO
    D(1,N) = PCHDF(iend,xtemp,stemp)
    iend = 1
  END IF
  !
  singular_error = .FALSE.
  skip_forward_pass = .FALSE.
  !
  IF( ibeg==0 ) THEN
    IF( N==2 ) THEN
      Wk(2,1) = 1._SP
      Wk(1,1) = 1._SP
      D(1,1) = 2._SP*Wk(2,2)
    ELSE
      Wk(2,1) = Wk(1,3)
      Wk(1,1) = Wk(1,2) + Wk(1,3)
      D(1,1) = ((Wk(1,2)+2._SP*Wk(1,1))*Wk(2,2)*Wk(1,3)+Wk(1,2)**2*Wk(2,3))/Wk(1,1)
    END IF
  ELSEIF( ibeg==1 ) THEN
    Wk(2,1) = 1._SP
    Wk(1,1) = 0._SP
  ELSE
    Wk(2,1) = 2._SP
    Wk(1,1) = 1._SP
    D(1,1) = 3._SP*Wk(2,2) - 0.5_SP*Wk(1,2)*D(1,1)
  END IF
  !
  nm1 = N - 1
  IF( nm1>1 ) THEN
    DO j = 2, nm1
      IF( Wk(2,j-1)==0._SP ) THEN
        singular_error = .TRUE.
        EXIT
      END IF
      g = -Wk(1,j+1)/Wk(2,j-1)
      D(1,j) = g*D(1,j-1) + 3._SP*(Wk(1,j)*Wk(2,j+1)+Wk(1,j+1)*Wk(2,j))
      Wk(2,j) = g*Wk(1,j-1) + 2._SP*(Wk(1,j)+Wk(1,j+1))
    END DO
  END IF
  !
  IF( .NOT. singular_error ) THEN
    IF( iend/=1 ) THEN
      IF( iend/=0 ) THEN
        D(1,N) = 3._SP*Wk(2,N) + 0.5_SP*Wk(1,N)*D(1,N)
        Wk(2,N) = 2._SP
        IF( Wk(2,N-1)==0._SP ) THEN
          singular_error = .TRUE.
        ELSE
          g = -1._SP/Wk(2,N-1)
        END IF
      ELSEIF( N==2 .AND. ibeg==0 ) THEN
        D(1,2) = Wk(2,2)
        skip_forward_pass = .TRUE.
      ELSEIF( (N==2) .OR. (N==3 .AND. ibeg==0) ) THEN
        D(1,N) = 2._SP*Wk(2,N)
        Wk(2,N) = 1._SP
        IF( Wk(2,N-1)==0._SP ) THEN
          singular_error = .TRUE.
        ELSE
          g = -1._SP/Wk(2,N-1)
        END IF
      ELSE
        g = Wk(1,N-1) + Wk(1,N)
        D(1,N) = ((Wk(1,N)+2._SP*g)*Wk(2,N)*Wk(1,N-1)+Wk(1,N)**2*(F(1,N-1)-F(1,N-2))/Wk(1,N-1))/g
        IF( Wk(2,N-1)==0._SP ) THEN
          singular_error = .TRUE.
        ELSE
          g = -g/Wk(2,N-1)
          Wk(2,N) = Wk(1,N-1)
        END IF
      END IF
      !
      IF( .NOT. singular_error .AND. .NOT. skip_forward_pass ) THEN
        Wk(2,N) = g*Wk(1,N-1) + Wk(2,N)
        IF( Wk(2,N)==0._SP ) THEN
          singular_error = .TRUE.
        ELSE
          D(1,N) = (g*D(1,N-1)+D(1,N))/Wk(2,N)
        END IF
      END IF
    END IF
  END IF
  !
  IF( .NOT. singular_error ) THEN
    DO j = nm1, 1, -1
      IF( Wk(2,j)==0._SP ) THEN
        singular_error = .TRUE.
        EXIT
      END IF
      D(1,j) = (D(1,j)-Wk(1,j)*D(1,j+1))/Wk(2,j)
    END DO
  END IF
  !
  IF( singular_error ) THEN
    Ierr = -8
    ERROR STOP 'PCHSP : SINGULAR LINEAR SYSTEM'
    RETURN
  END IF
  !
  RETURN
  !
END SUBROUTINE PCHSP
