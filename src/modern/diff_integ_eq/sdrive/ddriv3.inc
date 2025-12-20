!** DDRIV3
SUBROUTINE DDRIV3(N,T,Y,F,Nstate,Tout,Ntask,Nroot,Eps,Ewt,Ierror,Mint,&
    Miter,Impl,Ml,Mu,Mxord,Hmax,Work,Lenw,Iwork,Leniw,&
    JACOBN,FA,Nde,Mxstep,G,USERS,Ierflg)
  !> The function of DDRIV3 is to solve N ordinary differential equations of the
  !  form dY(I)/dT = F(Y(I),T), given the initial conditions Y(I) = YI.
  !  The program has options to allow the solution of both stiff and non-stiff
  !  differential equations.  Other important options are available.
  !  DDRIV3 uses double precision arithmetic.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Category:**  I1A2, I1A1B
  !***
  ! **Type:**      DOUBLE PRECISION (SDRIV3-S, DDRIV3-D, CDRIV3-C)
  !***
  ! **Keywords:**  DOUBLE PRECISION, GEAR'S METHOD, INITIAL VALUE PROBLEMS,
  !             ODE, ORDINARY DIFFERENTIAL EQUATIONS, SDRIVE, STIFF
  !***
  ! **Author:**  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***
  ! **Description:**
  !
  !  [Documentation preserved - see original file for full description]
  !
  !***
  ! **References:**  C. W. Gear, Numerical Initial Value Problems in
  !                 Ordinary Differential Equations, Prentice-Hall, 1971.
  !***
  ! **Routines called:**  D1MACH, DDNTP, DDSTP, DDZRO, DGBFA, DGBSL, DGEFA,
  !                    DGESL, DNRM2, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  !   251218  Refactored to eliminate all GOTO statements.
  !           Replaced labelled blocks with structured control flow.
  !           (C. Markwardt, modernisation project)
  USE service, ONLY : eps_dp, tiny_dp
  USE linpack, ONLY : DGBFA, DGEFA
  USE lapack, ONLY : DGBTRS, DGETRS
  !
  INTERFACE
    REAL(DP) PURE FUNCTION G(N,T,Y,Iroot)
      IMPORT DP
      INTEGER, INTENT(IN) :: N, Iroot
      REAL(DP), INTENT(IN) :: T, Y(N)
    END FUNCTION G
    PURE SUBROUTINE F(N,T,Y,Ydot)
      IMPORT DP
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: T, Y(:)
      REAL(DP), INTENT(OUT) :: Ydot(:)
    END SUBROUTINE F
    PURE SUBROUTINE JACOBN(N,T,Y,Dfdy,Matdim,Ml,Mu)
      IMPORT DP
      INTEGER, INTENT(IN) :: N, Matdim, Ml, Mu
      REAL(DP), INTENT(IN) :: T, Y(N)
      REAL(DP), INTENT(OUT) :: Dfdy(Matdim,N)
    END SUBROUTINE JACOBN
    PURE SUBROUTINE USERS(Y,Yh,Ywt,Save1,Save2,T,H,El,Impl,N,Nde,Iflag)
      IMPORT DP
      INTEGER, INTENT(IN) :: Impl, N, Nde, iflag
      REAL(DP), INTENT(IN) :: T, H, El
      REAL(DP), INTENT(IN) :: Y(N), Yh(N,13), Ywt(N)
      REAL(DP), INTENT(INOUT) :: Save1(N), Save2(N)
    END SUBROUTINE USERS
    PURE SUBROUTINE FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IMPORT DP
      INTEGER, INTENT(IN) :: N, Matdim, Ml, Mu, Nde
      REAL(DP), INTENT(IN) :: T, Y(N)
      REAL(DP), INTENT(INOUT) :: A(:,:)
    END SUBROUTINE FA
  END INTERFACE
  INTEGER, INTENT(IN) :: Ierror, Impl, Leniw, Lenw, Mint, Miter, Ml, Mu, Mxord, &
    Mxstep, N, Nde, Nroot, Ntask
  INTEGER, INTENT(INOUT) :: Nstate, Iwork(Leniw+N)
  INTEGER, INTENT(OUT) :: Ierflg
  REAL(DP), INTENT(IN) :: Hmax, Tout, Ewt(N)
  REAL(DP), INTENT(INOUT) :: Eps, T
  REAL(DP), INTENT(INOUT) :: Work(Lenw+Leniw), Y(N+1)
  !
  INTEGER :: i, ia, idfdy, ifac, iflag, ignow, imxerr, info, iroot, isave1, &
    isave2, itroot, iywt, j, jstate, jtroot, lenchk, liwchk, matdim, maxord, &
    ndecom, npar, nstepl
  REAL(DP) :: ae, big, glast, gnow, h, hsign, hused, re, sizee, summ, tlast, &
    troot, uround
  LOGICAL :: convrg
  REAL(DP), ALLOCATABLE :: a(:,:)
  CHARACTER(8) :: intgr1, intgr2
  CHARACTER(16) :: rl1, rl2
  REAL(DP), PARAMETER :: NROUND = 20._DP
  INTEGER, PARAMETER :: IAVGH = 1, IHUSED = 2, IAVGRD = 3, IEL = 4, IH = 160, &
    IHMAX = 161, IHOLD = 162, IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166, &
    ITOUT = 167, ITQ = 168, ITREND = 204, IMACH1 = 205, IMACH4 = 206, IYH = 251, &
    INDMXR = 1, INQUSE = 2, INSTEP = 3, INFE = 4, INJE = 5, INROOT = 6, ICNVRG = 7, &
    IJROOT = 8, IJTASK = 9, IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13, &
    INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17, IMTR = 18, IMXRDS = 19, &
    IMXORD = 20, INDPRT = 21, IJSTPL = 22, INDPVT = 51
  !
  ! Local variables for structured control flow
  LOGICAL :: do_return, do_return_500, fatal_error, singular_error
  LOGICAL :: ywt_done, skip_to_main_loop
  !
  !* FIRST EXECUTABLE STATEMENT  DDRIV3
  !
  ! Initialize control flags
  do_return = .FALSE.
  do_return_500 = .FALSE.
  fatal_error = .FALSE.
  singular_error = .FALSE.
  !
  IF( Nstate==12 ) THEN
    Ierflg = 999
    ERROR STOP 'DDRIV3 : Illegal input.  The value of NSTATE is 12 .'
    RETURN
  ELSEIF( Nstate<1 .OR. Nstate>12 ) THEN
    WRITE (intgr1,'(I8)') Nstate
    Ierflg = 26
    ERROR STOP 'DDRIV3 : Illegal input.  Improper value for NSTATE.'
    Nstate = 12
    RETURN
  END IF
  npar = N
  IF( Eps<0._DP ) THEN
    WRITE (rl1,'(D16.8)') Eps
    Ierflg = 27
    ERROR STOP 'DDRIV3 : Illegal input.  EPS is negative.'
    Nstate = 12
    RETURN
  END IF
  IF( N<=0 ) THEN
    WRITE (intgr1,'(I8)') N
    Ierflg = 22
    ERROR STOP 'DDRIV3 : Illegal input.  Number of equations is not positive.'
    Nstate = 12
    RETURN
  END IF
  IF( Mxord<=0 ) THEN
    WRITE (intgr1,'(I8)') Mxord
    Ierflg = 28
    ERROR STOP 'DDRIV3 : Illegal input.  Maximum order is not positive.'
    Nstate = 12
    RETURN
  END IF
  IF( Mint<1 .OR. Mint>3 ) THEN
    WRITE (intgr1,'(I8)') Mint
    Ierflg = 23
    ERROR STOP 'DDRIV3 : Illegal input.  Improper value for the integration method flag.'
    Nstate = 12
    RETURN
  ELSEIF( Miter<0 .OR. Miter>5 ) THEN
    WRITE (intgr1,'(I8)') Miter
    Ierflg = 24
    ERROR STOP 'DDRIV3 : Illegal input.  Improper value for MITER.'
    Nstate = 12
    RETURN
  ELSEIF( Impl<0 .OR. Impl>3 ) THEN
    WRITE (intgr1,'(I8)') Impl
    Ierflg = 25
    ERROR STOP 'DDRIV3 : Illegal input.  Improper value for IMPL.'
    Nstate = 12
    RETURN
  ELSEIF( Mint==3 .AND. (Miter==0 .OR. Miter==3 .OR. Impl/=0) ) THEN
    WRITE (intgr1,'(I8)') Miter
    WRITE (intgr2,'(I8)') Impl
    Ierflg = 29
    ERROR STOP 'DDRIV3 : Illegal input.  For MINT = 3, the value of MITER&
      & and/or IMPL is not allowed.'
    Nstate = 12
    RETURN
  ELSEIF( (Impl>=1 .AND. Impl<=3) .AND. Miter==0 ) THEN
    WRITE (intgr1,'(I8)') Impl
    Ierflg = 30
    ERROR STOP 'DDRIV3 : Illegal input.  For MITER = 0, the value of IMPL is not allowed.'
    Nstate = 12
    RETURN
  ELSEIF( (Impl==2 .OR. Impl==3) .AND. Mint==1 ) THEN
    WRITE (intgr1,'(I8)') Impl
    Ierflg = 31
    ERROR STOP 'DDRIV3 : Illegal input.  For MINT = 1, the value of IMPL is not allowed.'
    Nstate = 12
    RETURN
  END IF
  IF( Miter==0 .OR. Miter==3 ) THEN
    liwchk = INDPVT - 1
  ELSEIF( Miter==1 .OR. Miter==2 .OR. Miter==4 .OR. Miter==5 ) THEN
    liwchk = INDPVT + N - 1
  END IF
  IF( Leniw<liwchk ) THEN
    WRITE (intgr1,'(I8)') liwchk
    Ierflg = 33
    ERROR STOP 'DDRIV3 : Illegal input.  Insufficient storage&
      & allocated for the IWORK array.'
    Nstate = 12
    RETURN
  END IF
  !                                                Allocate the WORK array
  !                                         IYH is the index of YH in WORK
  IF( Mint==1 .OR. Mint==3 ) THEN
    maxord = MIN(Mxord,12)
  ELSEIF( Mint==2 ) THEN
    maxord = MIN(Mxord,5)
  END IF
  idfdy = IYH + (maxord+1)*N
  !                                             IDFDY is the index of DFDY
  !
  IF( Miter==0 .OR. Miter==3 ) THEN
    iywt = idfdy
  ELSEIF( Miter==1 .OR. Miter==2 ) THEN
    iywt = idfdy + N*N
  ELSEIF( Miter==4 .OR. Miter==5 ) THEN
    iywt = idfdy + (2*Ml+Mu+1)*N
  END IF
  !                                               IYWT is the index of YWT
  isave1 = iywt + N
  !                                           ISAVE1 is the index of SAVE1
  isave2 = isave1 + N
  !                                           ISAVE2 is the index of SAVE2
  ignow = isave2 + N
  !                                             IGNOW is the index of GNOW
  itroot = ignow + Nroot
  !                                           ITROOT is the index of TROOT
  ifac = itroot + Nroot
  !                                               IFAC is the index of FAC
  IF( Miter==2 .OR. Miter==5 .OR. Mint==3 ) THEN
    ia = ifac + N
  ELSE
    ia = ifac
  END IF
  !                                                   IA is the index of A
  IF( Impl==0 .OR. Miter==3 ) THEN
    lenchk = ia - 1
  ELSEIF( Impl==1 .AND. (Miter==1 .OR. Miter==2) ) THEN
    lenchk = ia - 1 + N*N
  ELSEIF( Impl==1 .AND. (Miter==4 .OR. Miter==5) ) THEN
    lenchk = ia - 1 + (2*Ml+Mu+1)*N
  ELSEIF( Impl==2 .AND. Miter/=3 ) THEN
    lenchk = ia - 1 + N
  ELSEIF( Impl==3 .AND. (Miter==1 .OR. Miter==2) ) THEN
    lenchk = ia - 1 + N*Nde
  ELSEIF( Impl==3 .AND. (Miter==4 .OR. Miter==5) ) THEN
    lenchk = ia - 1 + (2*Ml+Mu+1)*Nde
  END IF
  IF( Impl==0 .OR. Miter==3 ) THEN
    ALLOCATE( a(0,0) )
  ELSEIF( Impl==1 ) THEN
    ALLOCATE( a( (lenchk-ia+1)/N, N ) )
  ELSEIF( Impl==2 ) THEN
    ALLOCATE( a(N,1) )
  ELSEIF( Impl==3 ) THEN
    ALLOCATE( a( (lenchk-ia+1)/Nde, Nde ) )
  ENDIF
  IF( Lenw<lenchk ) THEN
    WRITE (intgr1,'(I8)') lenchk
    Ierflg = 32
    ERROR STOP 'DDRIV3 : Illegal input.  Insufficient storage&
      & allocated for the WORK array.'
    Nstate = 12
    RETURN
  END IF
  IF( Miter==0 .OR. Miter==3 ) THEN
    matdim = 1
  ELSEIF( Miter==1 .OR. Miter==2 ) THEN
    matdim = N
  ELSEIF( Miter==4 .OR. Miter==5 ) THEN
    matdim = 2*Ml + Mu + 1
  END IF
  IF( Impl==0 .OR. Impl==1 ) THEN
    ndecom = N
  ELSEIF( Impl==2 .OR. Impl==3 ) THEN
    ndecom = Nde
  END IF
  !
  skip_to_main_loop = .FALSE.
  !
  IF( Nstate==1 ) THEN
    !                                                  Initialize parameters
    IF( Mint==1 .OR. Mint==3 ) THEN
      Iwork(IMXORD) = MIN(Mxord,12)
    ELSEIF( Mint==2 ) THEN
      Iwork(IMXORD) = MIN(Mxord,5)
    END IF
    Iwork(IMXRDS) = Mxord
    IF( Mint==1 .OR. Mint==2 ) THEN
      Iwork(IMNT) = Mint
      Iwork(IMTR) = Miter
      Iwork(IMNTLD) = Mint
      Iwork(IMTRLD) = Miter
    ELSEIF( Mint==3 ) THEN
      Iwork(IMNT) = 1
      Iwork(IMTR) = 0
      Iwork(IMNTLD) = Iwork(IMNT)
      Iwork(IMTRLD) = Iwork(IMTR)
      Iwork(IMTRSV) = Miter
    END IF
    Work(IHMAX) = Hmax
    uround = eps_dp
    Work(IMACH4) = uround
    Work(IMACH1) = tiny_dp
    IF( Nroot/=0 ) THEN
      re = uround
      ae = Work(IMACH1)
    END IF
    h = (Tout-T)*(1._DP-4._DP*uround)
    h = SIGN(MIN(ABS(h),Hmax),h)
    Work(IH) = h
    hsign = SIGN(1._DP,h)
    Work(IHSIGN) = hsign
    Iwork(IJTASK) = 0
    Work(IAVGH) = 0._DP
    Work(IHUSED) = 0._DP
    Work(IAVGRD) = 0._DP
    Iwork(INDMXR) = 0
    Iwork(INQUSE) = 0
    Iwork(INSTEP) = 0
    Iwork(IJSTPL) = 0
    Iwork(INFE) = 0
    Iwork(INJE) = 0
    Iwork(INROOT) = 0
    Work(IT) = T
    Iwork(ICNVRG) = 0
    Iwork(INDPRT) = 0
    !                                                 Set initial conditions
    DO i = 1, N
      Work(i+IYH-1) = Y(i)
    END DO
    IF( T==Tout ) RETURN
    skip_to_main_loop = .TRUE.
  ELSE
    uround = Work(IMACH4)
    IF( Nroot/=0 ) THEN
      re = uround
      ae = Work(IMACH1)
    END IF
  END IF
  !
  IF( .NOT. skip_to_main_loop ) THEN
    !                                             On a continuation, check
    !                                             that output points have
    !                                             been or will be overtaken.
    IF( Iwork(ICNVRG)==1 ) THEN
      convrg = .TRUE.
    ELSE
      convrg = .FALSE.
    END IF
    T = Work(IT)
    h = Work(IH)
    hsign = Work(IHSIGN)
    IF( Iwork(IJTASK)/=0 ) THEN
      !
      !                                   IWORK(IJROOT) flags unreported
      !                                   roots, and is set to the value of
      !                                   NTASK when a root was last selected.
      !                                   It is set to zero when all roots
      !                                   have been reported.  IWORK(INROOT)
      !                                   contains the index and WORK(ITOUT)
      !                                   contains the value of the root last
      !                                   selected to be reported.
      !                                   IWORK(INRTLD) contains the value of
      !                                   NROOT and IWORK(INDTRT) contains
      !                                   the value of ITROOT when the array
      !                                   of roots was last calculated.
      IF( Nroot/=0 ) THEN
        IF( Iwork(IJROOT)>0 ) THEN
          !                                      TOUT has just been reported.
          !                                      If TROOT <= TOUT, report TROOT.
          IF( Nstate==5 ) THEN
            troot = T
            iroot = 0
            DO i = 1, Iwork(INRTLD)
              jtroot = i + Iwork(INDTRT) - 1
              IF( Work(jtroot)*hsign<=troot*hsign ) THEN
                !
                !                                              Check for multiple roots.
                !
                IF( Work(jtroot)==Work(ITOUT) .AND. i>Iwork(INROOT) ) THEN
                  iroot = i
                  troot = Work(jtroot)
                  EXIT
                END IF
                IF( Work(jtroot)*hsign>Work(ITOUT)*hsign ) THEN
                  iroot = i
                  troot = Work(jtroot)
                END IF
              END IF
            END DO
            Iwork(INROOT) = iroot
            Work(ITOUT) = troot
            Iwork(IJROOT) = Ntask
            IF( Ntask==1 ) THEN
              IF( iroot==0 ) THEN
                Iwork(IJROOT) = 0
              ELSEIF( Tout*hsign>=troot*hsign ) THEN
                CALL DDNTP(h,0,N,Iwork(INQ),T,troot,Work(IYH),Y)
                Nstate = 5
                T = troot
                Ierflg = 0
                do_return_500 = .TRUE.
              END IF
            ELSEIF( Ntask==2 .OR. Ntask==3 ) THEN
              !
              !                                     If there are no more roots, or the
              !                                     user has altered TOUT to be less
              !                                     than a root, set IJROOT to zero.
              !
              IF( iroot==0 .OR. (Tout*hsign<troot*hsign) ) THEN
                Iwork(IJROOT) = 0
              ELSE
                CALL DDNTP(h,0,N,Iwork(INQ),T,troot,Work(IYH),Y)
                Nstate = 5
                Ierflg = 0
                T = troot
                do_return_500 = .TRUE.
              END IF
            END IF
          ELSEIF( Tout*hsign>=Work(ITOUT)*hsign ) THEN
            troot = Work(ITOUT)
            CALL DDNTP(h,0,N,Iwork(INQ),T,troot,Work(IYH),Y)
            T = troot
            Nstate = 5
            Ierflg = 0
            do_return_500 = .TRUE.
            !                                         A root has just been reported.
            !                                         Select the next root.
          END IF
        END IF
      END IF
      !
      IF( .NOT. do_return_500 ) THEN
        IF( Ntask==1 ) THEN
          Nstate = 2
          IF( T*hsign>=Tout*hsign ) THEN
            CALL DDNTP(h,0,N,Iwork(INQ),T,Tout,Work(IYH),Y)
            T = Tout
            Ierflg = 0
            do_return_500 = .TRUE.
          END IF
        ELSEIF( Ntask==2 ) THEN
          !                                                      Check if TOUT has
          !                                                      been reset < T
          IF( T*hsign>Tout*hsign ) THEN
            WRITE (rl1,'(D16.8)') T
            WRITE (rl2,'(D16.8)') Tout
            Ierflg = 11
            Nstate = 11
            CALL DDNTP(h,0,N,Iwork(INQ),T,Tout,Work(IYH),Y)
            T = Tout
            do_return_500 = .TRUE.
          ELSEIF( ABS(Tout-T)<=NROUND*uround*MAX(ABS(T),ABS(Tout)) ) THEN
            !                                   Determine if TOUT has been overtaken
            T = Tout
            Nstate = 2
            Ierflg = 0
            do_return = .TRUE.
          ELSEIF( Nstate==5 ) THEN
            !                                             If there are no more roots
            !                                             to report, report T.
            Nstate = 2
            Ierflg = 0
            do_return = .TRUE.
          ELSE
            Nstate = 2
            !                                                       See if TOUT will
            !                                                       be overtaken.
            IF( (T+h)*hsign>Tout*hsign ) THEN
              h = Tout - T
              IF( (T+h)*hsign>Tout*hsign ) h = h*(1._DP-4._DP*uround)
              Work(IH) = h
              IF( h==0._DP ) THEN
                fatal_error = .TRUE.
              ELSE
                Iwork(IJTASK) = -1
              END IF
            END IF
          END IF
        ELSEIF( Ntask==3 ) THEN
          Nstate = 2
          IF( T*hsign>Tout*hsign ) THEN
            WRITE (rl1,'(D16.8)') T
            WRITE (rl2,'(D16.8)') Tout
            Ierflg = 11
            Nstate = 11
            CALL DDNTP(h,0,N,Iwork(INQ),T,Tout,Work(IYH),Y)
            T = Tout
            do_return_500 = .TRUE.
          ELSEIF( ABS(Tout-T)<=NROUND*uround*MAX(ABS(T),ABS(Tout)) ) THEN
            T = Tout
            Ierflg = 0
            do_return = .TRUE.
          ELSE
            IF( (T+h)*hsign>Tout*hsign ) THEN
              h = Tout - T
              IF( (T+h)*hsign>Tout*hsign ) h = h*(1._DP-4._DP*uround)
              Work(IH) = h
              IF( h==0._DP ) THEN
                fatal_error = .TRUE.
              ELSE
                Iwork(IJTASK) = -1
              END IF
            END IF
          END IF
        END IF
      END IF
      !
      IF( .NOT. (do_return .OR. do_return_500 .OR. fatal_error) ) THEN
        !                         Implement changes in MINT, MITER, and/or HMAX.
        IF( (Mint/=Iwork(IMNTLD) .OR. Miter/=Iwork(IMTRLD)) .AND. Mint/=3 .AND. &
          Iwork(IMNTLD)/=3 ) Iwork(IJTASK) = -1
        IF( Hmax/=Work(IHMAX) ) THEN
          h = SIGN(MIN(ABS(h),Hmax),h)
          IF( h/=Work(IH) ) THEN
            Iwork(IJTASK) = -1
            Work(IH) = h
          END IF
          Work(IHMAX) = Hmax
        END IF
      END IF
    END IF
  END IF
  !
  ! Handle fatal error from step size going to zero
  IF( fatal_error ) THEN
    WRITE (rl1,'(D16.8)') T
    Ierflg = 41
    ERROR STOP 'DDRIV3 : At T the attempted step size has gone to zero.&
      & Often this occurs if the problem setup is incorrect.'
    Nstate = 12
    RETURN
  END IF
  !
  ! Handle early returns
  IF( do_return ) THEN
    DO i = 1, N
      Y(i) = Work(i+IYH-1)
    END DO
    ! Fall through to return block 500
    do_return_500 = .TRUE.
  END IF
  !
  IF( do_return_500 ) THEN
    IF( Iwork(IJTASK)==0 ) RETURN
    big = 0._DP
    imxerr = 1
    DO i = 1, N
      sizee = ABS(Work(i+isave1-1)/Work(i+iywt-1))
      IF( big<sizee ) THEN
        big = sizee
        imxerr = i
      END IF
    END DO
    Iwork(INDMXR) = imxerr
    Work(IHUSED) = hused
    RETURN
  END IF
  !
  !-----------------------------------------------------------------------
  !     MAIN INTEGRATION LOOP
  !-----------------------------------------------------------------------
  !
  main_loop: DO
    nstepl = Iwork(INSTEP)
    DO i = 1, N
      Y(i) = Work(i+IYH-1)
    END DO
    IF( Nroot/=0 ) THEN
      DO i = 1, Nroot
        Work(i+ignow-1) = G(npar,T,Y,i)
        IF( npar==0 ) THEN
          Iwork(INROOT) = i
          Nstate = 7
          RETURN
        END IF
      END DO
    END IF
    !
    ! Set up YWT array based on IERROR
    ywt_done = .FALSE.
    !
    IF( Ierror==1 ) THEN
      DO i = 1, N
        Work(i+iywt-1) = 1._DP
      END DO
      ywt_done = .TRUE.
    ELSEIF( Ierror==5 ) THEN
      DO i = 1, N
        Work(i+iywt-1) = Ewt(i)
      END DO
      ywt_done = .TRUE.
    END IF
    !
    !                                       Reset YWT array.  Looping point.
    ywt_loop: DO WHILE( .NOT. ywt_done )
      IF( Ierror==2 ) THEN
        ! Check if any Y(i) is zero
        DO i = 1, N
          IF( Y(i)==0._DP ) THEN
            ! Handle zero Y case
            IF( Iwork(IJTASK)==0 ) THEN
              CALL F(npar,T,Y,Work(isave2:))
              IF( npar==0 ) THEN
                Nstate = 6
                RETURN
              END IF
              Iwork(INFE) = Iwork(INFE) + 1
              IF( Miter==3 .AND. Impl/=0 ) THEN
                iflag = 0
                CALL USERS(Y,Work(IYH),Work(iywt),Work(isave1),Work(isave2),T,h,&
                  Work(IEL),Impl,npar,ndecom,iflag)
                IF( iflag==-1 ) THEN
                  singular_error = .TRUE.
                  EXIT ywt_loop
                END IF
                IF( npar==0 ) THEN
                  Nstate = 10
                  RETURN
                END IF
              ELSEIF( Impl==1 ) THEN
                IF( Miter==1 .OR. Miter==2 ) THEN
                  CALL FA(npar,T,Y,a,matdim,Ml,Mu,ndecom)
                  IF( npar==0 ) THEN
                    Nstate = 9
                    RETURN
                  END IF
                  CALL DGEFA(a,matdim,N,Iwork(INDPVT),info)
                  IF( info/=0 ) THEN
                    singular_error = .TRUE.
                    EXIT ywt_loop
                  END IF
                  CALL DGETRS('N',N,1,a,matdim,Iwork(INDPVT),Work(isave2),N,info)
                ELSEIF( Miter==4 .OR. Miter==5 ) THEN
                  CALL FA(npar,T,Y,a(Ml+1:,:),matdim,Ml,Mu,ndecom)
                  IF( npar==0 ) THEN
                    Nstate = 9
                    RETURN
                  END IF
                  CALL DGBFA(a,matdim,N,Ml,Mu,Iwork(INDPVT),info)
                  IF( info/=0 ) THEN
                    singular_error = .TRUE.
                    EXIT ywt_loop
                  END IF
                  CALL DGBTRS('N',N,Ml,Mu,1,a,matdim,Iwork(INDPVT),Work(isave2),N,info)
                END IF
              ELSEIF( Impl==2 ) THEN
                CALL FA(npar,T,Y,a,matdim,Ml,Mu,ndecom)
                IF( npar==0 ) THEN
                  Nstate = 9
                  RETURN
                END IF
                DO j = 1, ndecom
                  IF( Work(j+ia-1)==0._DP ) THEN
                    singular_error = .TRUE.
                    EXIT ywt_loop
                  END IF
                  Work(j+isave2-1) = Work(j+isave2-1)/Work(j+ia-1)
                END DO
              ELSEIF( Impl==3 ) THEN
                IF( Miter==1 .OR. Miter==2 ) THEN
                  CALL FA(npar,T,Y,a,matdim,Ml,Mu,ndecom)
                  IF( npar==0 ) THEN
                    Nstate = 9
                    RETURN
                  END IF
                  CALL DGEFA(a,matdim,Nde,Iwork(INDPVT),info)
                  IF( info/=0 ) THEN
                    singular_error = .TRUE.
                    EXIT ywt_loop
                  END IF
                  CALL DGETRS('N',Nde,1,a,matdim,Iwork(INDPVT),Work(isave2),N,info)
                ELSEIF( Miter==4 .OR. Miter==5 ) THEN
                  CALL FA(npar,T,Y,a(Ml+1:,:),matdim,Ml,Mu,ndecom)
                  IF( npar==0 ) THEN
                    Nstate = 9
                    RETURN
                  END IF
                  CALL DGBFA(a,matdim,Nde,Ml,Mu,Iwork(INDPVT),info)
                  IF( info/=0 ) THEN
                    singular_error = .TRUE.
                    EXIT ywt_loop
                  END IF
                  CALL DGBTRS('N',Nde,Ml,Mu,1,a,matdim,Iwork(INDPVT),Work(isave2),N,info)
                END IF
              END IF
            END IF
            ! Set YWT for remaining components
            DO j = i, N
              IF( Y(j)/=0._DP ) THEN
                Work(j+iywt-1) = ABS(Y(j))
              ELSEIF( Iwork(IJTASK)==0 ) THEN
                Work(j+iywt-1) = ABS(h*Work(j+isave2-1))
              ELSE
                Work(j+iywt-1) = ABS(Work(j+IYH+N-1))
              END IF
              IF( Work(j+iywt-1)==0._DP ) Work(j+iywt-1) = uround
            END DO
            ywt_done = .TRUE.
            EXIT ywt_loop
          END IF
          Work(i+iywt-1) = ABS(Y(i))
        END DO
        ywt_done = .TRUE.
      ELSEIF( Ierror==3 ) THEN
        DO i = 1, N
          Work(i+iywt-1) = MAX(Ewt(1),ABS(Y(i)))
        END DO
        ywt_done = .TRUE.
      ELSEIF( Ierror==4 ) THEN
        DO i = 1, N
          Work(i+iywt-1) = MAX(Ewt(i),ABS(Y(i)))
        END DO
        ywt_done = .TRUE.
      ELSE
        ywt_done = .TRUE.
      END IF
    END DO ywt_loop
    !
    ! Handle singular matrix error
    IF( singular_error ) THEN
      WRITE (rl1,'(D16.8)') T
      Ierflg = 43
      ERROR STOP 'DDRIV3 : At T while solving A*YDOT = F, A is singular.'
      Nstate = 12
      RETURN
    END IF
    !
    ! Continue with weight scaling
    DO i = 1, N
      Work(i+isave2-1) = Y(i)/Work(i+iywt-1)
    END DO
    summ = NORM2(Work(isave2:isave2+N-1))/SQRT(REAL(N, DP))
    summ = MAX(1._DP,summ)
    IF( Eps<summ*uround ) THEN
      Eps = summ*uround*(1._DP+10._DP*uround)
      WRITE (rl1,'(D16.8)') T
      WRITE (rl2,'(D16.8)') Eps
      Ierflg = 4
      Nstate = 4
      ! Return via block 400
      DO i = 1, N
        Y(i) = Work(i+IYH-1)
      END DO
      EXIT main_loop
    END IF
    IF( ABS(h)>=uround*ABS(T) ) THEN
      Iwork(INDPRT) = 0
    ELSEIF( Iwork(INDPRT)==0 ) THEN
      WRITE (rl1,'(D16.8)') T
      WRITE (rl2,'(D16.8)') h
      Ierflg = 15
      Iwork(INDPRT) = 1
    END IF
    IF( Ntask/=2 ) THEN
      IF( (Iwork(INSTEP)-nstepl)==Mxstep ) THEN
        WRITE (rl1,'(D16.8)') T
        WRITE (intgr1,'(I8)') Mxstep
        WRITE (rl2,'(D16.8)') Tout
        Ierflg = 3
        Nstate = 3
        ! Return via block 400
        DO i = 1, N
          Y(i) = Work(i+IYH-1)
        END DO
        EXIT main_loop
      END IF
    END IF
    !
    CALL DDSTP(Eps,F,FA,Work(IHMAX),Impl,Ierror,JACOBN,matdim,Iwork(IMXORD),&
      Iwork(IMNT),Iwork(IMTR),Ml,Mu,npar,ndecom,Work(iywt),uround,&
      USERS,Work(IAVGH),Work(IAVGRD),Work(IH),hused,Iwork(IJTASK),&
      Iwork(IMNTLD),Iwork(IMTRLD),Iwork(INFE),Iwork(INJE),&
      Iwork(INQUSE),Iwork(INSTEP),Work(IT),Y,Work(IYH),a,&
      convrg,Work(idfdy),Work(IEL),Work(ifac),Work(IHOLD),&
      Iwork(INDPVT),jstate,Iwork(IJSTPL),Iwork(INQ),Iwork(INWAIT),&
      Work(IRC),Work(IRMAX),Work(isave1),Work(isave2),Work(ITQ),&
      Work(ITREND),Mint,Iwork(IMTRSV),Iwork(IMXRDS))
    T = Work(IT)
    h = Work(IH)
    IF( convrg ) THEN
      Iwork(ICNVRG) = 1
    ELSE
      Iwork(ICNVRG) = 0
    END IF
    !
    SELECT CASE (jstate)
      CASE (2)
        ! Step size went to zero
        WRITE (rl1,'(D16.8)') T
        Ierflg = 41
        ERROR STOP 'DDRIV3 : At T the attempted step size has gone to zero.&
          & Often this occurs if the problem setup is incorrect.'
        Nstate = 12
        RETURN
      CASE (3)
        WRITE (rl1,'(D16.8)') T
        Ierflg = 42
        ERROR STOP 'DDRIV3 : At T the step size has been reduced about 50 &
          & times without advancing the solution.&
          & Often this occurs if the problem setup is incorrect.'
        Nstate = 12
        RETURN
      CASE (4,5)
        ! Singular matrix
        WRITE (rl1,'(D16.8)') T
        Ierflg = 43
        ERROR STOP 'DDRIV3 : At T while solving A*YDOT = F, A is singular.'
        Nstate = 12
        RETURN
      CASE (6,7,8,9,10)
        Nstate = jstate
        RETURN
      CASE DEFAULT
        ! jstate = 1, successful step
        Iwork(IJTASK) = 1
        !                                 Determine if a root has been overtaken
        IF( Nroot/=0 ) THEN
          iroot = 0
          DO i = 1, Nroot
            glast = Work(i+ignow-1)
            gnow = G(npar,T,Y,i)
            IF( npar==0 ) THEN
              Iwork(INROOT) = i
              Nstate = 7
              RETURN
            END IF
            Work(i+ignow-1) = gnow
            IF( glast*gnow>0._DP ) THEN
              Work(i+itroot-1) = T + h
            ELSEIF( gnow==0._DP ) THEN
              Work(i+itroot-1) = T
              iroot = i
            ELSEIF( glast==0._DP ) THEN
              Work(i+itroot-1) = T + h
            ELSEIF( ABS(hused)>=uround*ABS(T) ) THEN
              tlast = T - hused
              iroot = i
              troot = T
              CALL DDZRO(ae,G,h,npar,Iwork(INQ),iroot,re,T,Work(IYH),uround,&
                troot,tlast,gnow,glast,Y)
              DO j = 1, N
                Y(j) = Work(IYH+j-1)
              END DO
              IF( npar==0 ) THEN
                Iwork(INROOT) = i
                Nstate = 7
                RETURN
              END IF
              Work(i+itroot-1) = troot
            ELSE
              Work(i+itroot-1) = T
              iroot = i
            END IF
          END DO
          IF( iroot==0 ) THEN
            Iwork(IJROOT) = 0
            !                                                  Select the first root
          ELSE
            Iwork(IJROOT) = Ntask
            Iwork(INRTLD) = Nroot
            Iwork(INDTRT) = itroot
            troot = T + h
            DO i = 1, Nroot
              IF( Work(i+itroot-1)*hsign<troot*hsign ) THEN
                troot = Work(i+itroot-1)
                iroot = i
              END IF
            END DO
            Iwork(INROOT) = iroot
            Work(ITOUT) = troot
            IF( troot*hsign<=Tout*hsign ) THEN
              CALL DDNTP(h,0,N,Iwork(INQ),T,troot,Work(IYH),Y)
              Nstate = 5
              T = troot
              Ierflg = 0
              EXIT main_loop
            END IF
          END IF
        END IF
        !                               Test for NTASK condition to be satisfied
        Nstate = 2
        IF( Ntask==1 ) THEN
          IF( T*hsign<Tout*hsign ) CYCLE main_loop
          CALL DDNTP(h,0,N,Iwork(INQ),T,Tout,Work(IYH),Y)
          T = Tout
          Ierflg = 0
          EXIT main_loop
          !                               TOUT is assumed to have been attained
          !                               exactly if T is within twenty roundoff
          !                               units of TOUT, relative to MAX(TOUT, T).
        ELSEIF( Ntask==2 ) THEN
          IF( ABS(Tout-T)<=NROUND*uround*MAX(ABS(T),ABS(Tout)) ) THEN
            T = Tout
          ELSEIF( (T+h)*hsign>Tout*hsign ) THEN
            h = Tout - T
            IF( (T+h)*hsign>Tout*hsign ) h = h*(1._DP-4._DP*uround)
            Work(IH) = h
            IF( h==0._DP ) THEN
              WRITE (rl1,'(D16.8)') T
              Ierflg = 41
              ERROR STOP 'DDRIV3 : At T the attempted step size has gone to zero.&
                & Often this occurs if the problem setup is incorrect.'
              Nstate = 12
              RETURN
            END IF
            Iwork(IJTASK) = -1
          END IF
        ELSEIF( Ntask==3 ) THEN
          IF( ABS(Tout-T)<=NROUND*uround*MAX(ABS(T),ABS(Tout)) ) THEN
            T = Tout
          ELSE
            IF( (T+h)*hsign>Tout*hsign ) THEN
              h = Tout - T
              IF( (T+h)*hsign>Tout*hsign ) h = h*(1._DP-4._DP*uround)
              Work(IH) = h
              IF( h==0._DP ) THEN
                WRITE (rl1,'(D16.8)') T
                Ierflg = 41
                ERROR STOP 'DDRIV3 : At T the attempted step size has gone to zero.&
                  & Often this occurs if the problem setup is incorrect.'
                Nstate = 12
                RETURN
              END IF
              Iwork(IJTASK) = -1
            END IF
            CYCLE main_loop
          END IF
        END IF
        Ierflg = 0
    END SELECT
    !
    ! Successful completion - exit main loop
    EXIT main_loop
    !
  END DO main_loop
  !
  !                                      All returns are made through this
  !                                      section.  IMXERR is determined.
  IF( Iwork(IJTASK)==0 ) RETURN
  big = 0._DP
  imxerr = 1
  DO i = 1, N
    !                                            SIZE = ABS(ERROR(I)/YWT(I))
    sizee = ABS(Work(i+isave1-1)/Work(i+iywt-1))
    IF( big<sizee ) THEN
      big = sizee
      imxerr = i
    END IF
  END DO
  Iwork(INDMXR) = imxerr
  Work(IHUSED) = hused
  !
END SUBROUTINE DDRIV3
