!** ZUNK1
PURE SUBROUTINE ZUNK1(Z,Fnu,Kode,Mr,N,Y,Nz,Tol,Elim,Alim)
  !> Subsidiary to ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CUNK1-A, ZUNK1-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
  !     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
  !     UNIFORM ASYMPTOTIC EXPANSION.
  !     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
  !     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
  !
  !***
  ! **See also:**  ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZS1S2, ZUCHK, ZUNIK
  !
  ! REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   2024xx  Modernized: eliminated 17 GOTOs using structured control flow
  USE service, ONLY : tiny_dp, huge_dp
  !
  INTEGER, INTENT(IN) :: Kode, Mr, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, ib, iflag, ifn, il, init(2), inu, iuf, k, kdflg, kflag, &
    kk, nw, j, ipard, initd, ic, m
  REAL(DP) :: ang, aphi, asc, ascle, bry(3), cpn, c2i, c2m, c2r, &
    fmr, fn, fnf, rs1, sgn, spn, x
  COMPLEX(DP) :: cfn, ck, crsc, cs, cscl, csgn, cspn, csr(3), css(3), &
    cwrk(16,3), cy(2), c1, c2, phi(2), rz, summ(2), s1, s2, &
    zeta1(2), zeta2(2), zr, phid, zeta1d, zeta2d, sumd
  REAL(DP), PARAMETER :: pi = 3.14159265358979324_DP
  ! Control flow flags
  logical :: underflow_occurred, exit_main_loop, do_forward_recur
  logical :: skip_param_reuse, exit_cont_loop, do_backward_recur
  !* FIRST EXECUTABLE STATEMENT  ZUNK1
  kdflg = 1
  Nz = 0
  !-----------------------------------------------------------------------
  !     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
  !     THE UNDERFLOW LIMIT
  !-----------------------------------------------------------------------
  cscl = CMPLX(1._DP/Tol,0._DP,DP)
  crsc = CMPLX(Tol,0._DP,DP)
  css(1) = cscl
  css(2) = (1._DP,0._DP)
  css(3) = crsc
  csr(1) = crsc
  csr(2) = (1._DP,0._DP)
  csr(3) = cscl
  bry(1) = 1.E+3_DP*tiny_dp/Tol
  bry(2) = 1._DP/bry(1)
  bry(3) = huge_dp
  x = REAL(Z,DP)
  zr = Z
  IF( x<0._DP ) zr = -Z
  j = 2
  exit_main_loop = .false.
  
  main_loop: DO i = 1, N
    !-----------------------------------------------------------------------
    !     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
    !-----------------------------------------------------------------------
    j = 3 - j
    fn = Fnu + (i-1)
    init(j) = 0
    CALL ZUNIK(zr,fn,2,0,Tol,init(j),phi(j),zeta1(j),zeta2(j),summ(j),cwrk(1,j))
    IF( Kode==1 ) THEN
      s1 = zeta1(j) - zeta2(j)
    ELSE
      cfn = CMPLX(fn,0._DP,DP)
      s1 = zeta1(j) - cfn*(cfn/(zr+zeta2(j)))
    END IF
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    rs1 = REAL(s1,DP)
    underflow_occurred = .false.
    
    IF( ABS(rs1)<=Elim ) THEN
      IF( kdflg==1 ) kflag = 2
      IF( ABS(rs1)>=Alim ) THEN
        !-----------------------------------------------------------------------
        !     REFINE  TEST AND SCALE
        !-----------------------------------------------------------------------
        aphi = ABS(phi(j))
        rs1 = rs1 + LOG(aphi)
        IF( ABS(rs1)>Elim ) underflow_occurred = .true.
        IF( .not. underflow_occurred ) THEN
          IF( kdflg==1 ) kflag = 1
          IF( rs1>=0._DP ) THEN
            IF( kdflg==1 ) kflag = 3
          END IF
        END IF
      END IF
      
      IF( .not. underflow_occurred ) THEN
        !-----------------------------------------------------------------------
        !     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
        !     EXPONENT EXTREMES
        !-----------------------------------------------------------------------
        s2 = phi(j)*summ(j)
        c2r = REAL(s1,DP)
        c2i = AIMAG(s1)
        c2m = EXP(c2r)*REAL(css(kflag),DP)
        s1 = CMPLX(c2m,0._DP,DP)*CMPLX(COS(c2i),SIN(c2i),DP)
        s2 = s2*s1
        IF( kflag==1 ) THEN
          CALL ZUCHK(s2,nw,bry(1),Tol)
          IF( nw/=0 ) underflow_occurred = .true.
        END IF
        
        IF( .not. underflow_occurred ) THEN
          cy(kdflg) = s2
          Y(i) = s2*csr(kflag)
          IF( kdflg==2 ) THEN
            exit_main_loop = .true.
            EXIT main_loop
          END IF
          kdflg = 2
          CYCLE main_loop
        END IF
      END IF
    ELSE
      underflow_occurred = .true.
    END IF
    
    ! Handle underflow (was label 50)
    IF( rs1>0._DP ) THEN
      Nz = -1
      RETURN
    END IF
    !-----------------------------------------------------------------------
    !     FOR X<0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
    !-----------------------------------------------------------------------
    IF( x<0._DP ) THEN
      Nz = -1
      RETURN
    END IF
    kdflg = 1
    Y(i) = (0._DP,0._DP)
    Nz = Nz + 1
    IF( i/=1 ) THEN
      IF( Y(i-1)/=(0._DP,0._DP) ) THEN
        Y(i-1) = (0._DP,0._DP)
        Nz = Nz + 1
      END IF
    END IF
  END DO main_loop
  
  IF( .not. exit_main_loop ) i = N
  
  ! Label 100 equivalent - setup for forward recurrence
  rz = CMPLX(2._DP,0._DP,DP)/zr
  ck = CMPLX(fn,0._DP,DP)*rz
  ib = i + 1
  do_forward_recur = (N >= ib)

  IF( do_forward_recur ) THEN
    !-----------------------------------------------------------------------
    !  TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO ON UNDERFLOW
    !-----------------------------------------------------------------------
    fn = Fnu + (N-1)
    ipard = 1
    IF( Mr/=0 ) ipard = 0
    initd = 0
    CALL ZUNIK(zr,fn,2,ipard,Tol,initd,phid,zeta1d,zeta2d,sumd,cwrk(1,3))
    IF( Kode==1 ) THEN
      s1 = zeta1d - zeta2d
    ELSE
      cfn = CMPLX(fn,0._DP,DP)
      s1 = zeta1d - cfn*(cfn/(zr+zeta2d))
    END IF
    rs1 = REAL(s1,DP)
    do_forward_recur = .false.  ! Assume we will not do forward recurrence
    
    IF( ABS(rs1)<=Elim ) THEN
      IF( ABS(rs1)<Alim ) THEN
        do_forward_recur = .true.
      ELSE
        !-----------------------------------------------------------------------
        !     REFINE ESTIMATE AND TEST
        !-----------------------------------------------------------------------
        aphi = ABS(phid)
        rs1 = rs1 + LOG(aphi)
        IF( ABS(rs1)<Elim ) do_forward_recur = .true.
      END IF
    END IF
    
    IF( .not. do_forward_recur ) THEN
      IF( rs1>0._DP ) THEN
        Nz = -1
        RETURN
      END IF
      !-----------------------------------------------------------------------
      !     FOR X<0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
      !-----------------------------------------------------------------------
      IF( x<0._DP ) THEN
        Nz = -1
        RETURN
      END IF
      Nz = N
      DO i = 1, N
        Y(i) = (0._DP,0._DP)
      END DO
      RETURN
    END IF
    
    !-----------------------------------------------------------------------
    !     RECUR FORWARD FOR REMAINDER OF THE SEQUENCE (was label 200)
    !-----------------------------------------------------------------------
    s1 = cy(1)
    s2 = cy(2)
    c1 = csr(kflag)
    ascle = bry(kflag)
    DO i = ib, N
      c2 = s2
      s2 = ck*s2 + s1
      s1 = c2
      ck = ck + rz
      c2 = s2*c1
      Y(i) = c2
      IF( kflag<3 ) THEN
        c2r = REAL(c2,DP)
        c2i = AIMAG(c2)
        c2r = ABS(c2r)
        c2i = ABS(c2i)
        c2m = MAX(c2r,c2i)
        IF( c2m>ascle ) THEN
          kflag = kflag + 1
          ascle = bry(kflag)
          s1 = s1*c1
          s2 = c2
          s1 = s1*css(kflag)
          s2 = s2*css(kflag)
          c1 = csr(kflag)
        END IF
      END IF
    END DO
  END IF
  
  ! Label 300 equivalent - Analytic continuation section
  IF( Mr==0 ) RETURN
  !-----------------------------------------------------------------------
  !     ANALYTIC CONTINUATION FOR RE(Z)<0.0E0
  !-----------------------------------------------------------------------
  Nz = 0
  fmr = Mr
  sgn = -SIGN(pi,fmr)
  !-----------------------------------------------------------------------
  !     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
  !-----------------------------------------------------------------------
  csgn = CMPLX(0._DP,sgn,DP)
  inu = INT( Fnu )
  fnf = Fnu - inu
  ifn = inu + N - 1
  ang = fnf*sgn
  cpn = COS(ang)
  spn = SIN(ang)
  cspn = CMPLX(cpn,spn,DP)
  IF( MOD(ifn,2)==1 ) cspn = -cspn
  asc = bry(1)
  kk = N
  iuf = 0
  kdflg = 1
  ib = ib - 1
  ic = ib - 1
  do_backward_recur = .false.
  
  cont_loop: DO k = 1, N
    fn = Fnu + (kk-1)
    !-----------------------------------------------------------------------
    !     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
    !     FUNCTION ABOVE
    !-----------------------------------------------------------------------
    skip_param_reuse = .false.
    m = 3
    IF( N>2 ) THEN
      IF( (kk==N) .AND. (ib<N) ) THEN
        skip_param_reuse = .true.
      ELSE IF( (kk/=ib) .AND. (kk/=ic) ) THEN
        initd = 0
        skip_param_reuse = .true.
      END IF
    END IF
    
    IF( .not. skip_param_reuse ) THEN
      initd = init(j)
      phid = phi(j)
      zeta1d = zeta1(j)
      zeta2d = zeta2(j)
      sumd = summ(j)
      m = j
      j = 3 - j
    END IF
    
    ! Label 350 equivalent - call ZUNIK
    CALL ZUNIK(zr,fn,1,0,Tol,initd,phid,zeta1d,zeta2d,sumd,cwrk(1,m))
    IF( Kode==1 ) THEN
      s1 = -zeta1d + zeta2d
    ELSE
      cfn = CMPLX(fn,0._DP,DP)
      s1 = -zeta1d + cfn*(cfn/(zr+zeta2d))
    END IF
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    rs1 = REAL(s1,DP)
    underflow_occurred = .false.
    
    IF( ABS(rs1)<=Elim ) THEN
      IF( kdflg==1 ) iflag = 2
      IF( ABS(rs1)>=Alim ) THEN
        !-----------------------------------------------------------------------
        !     REFINE  TEST AND SCALE
        !-----------------------------------------------------------------------
        aphi = ABS(phid)
        rs1 = rs1 + LOG(aphi)
        IF( ABS(rs1)>Elim ) underflow_occurred = .true.
        IF( .not. underflow_occurred ) THEN
          IF( kdflg==1 ) iflag = 1
          IF( rs1>=0._DP ) THEN
            IF( kdflg==1 ) iflag = 3
          END IF
        END IF
      END IF
      
      IF( .not. underflow_occurred ) THEN
        s2 = csgn*phid*sumd
        c2r = REAL(s1,DP)
        c2i = AIMAG(s1)
        c2m = EXP(c2r)*REAL(css(iflag),DP)
        s1 = CMPLX(c2m,0._DP,DP)*CMPLX(COS(c2i),SIN(c2i),DP)
        s2 = s2*s1
        IF( iflag==1 ) THEN
          CALL ZUCHK(s2,nw,bry(1),Tol)
          IF( nw/=0 ) s2 = CMPLX(0._DP,0._DP,DP)
        END IF
      END IF
    ELSE
      underflow_occurred = .true.
    END IF
    
    ! Handle underflow case (was label 450)
    IF( underflow_occurred ) THEN
      IF( rs1>0._DP ) THEN
        Nz = -1
        RETURN
      END IF
      s2 = (0._DP,0._DP)
    END IF
    
    ! Label 400 equivalent - store result and combine I and K functions
    cy(kdflg) = s2
    c2 = s2
    s2 = s2*csr(iflag)
    !-----------------------------------------------------------------------
    !     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
    !-----------------------------------------------------------------------
    s1 = Y(kk)
    IF( Kode/=1 ) THEN
      CALL ZS1S2(zr,s1,s2,nw,asc,Alim,iuf)
      Nz = Nz + nw
    END IF
    Y(kk) = s1*cspn + s2
    kk = kk - 1
    cspn = -cspn
    IF( c2/=(0._DP,0._DP) ) THEN
      IF( kdflg==2 ) THEN
        do_backward_recur = .true.
        EXIT cont_loop
      END IF
      kdflg = 2
    ELSE
      kdflg = 1
    END IF
  END DO cont_loop
  
  IF( .not. do_backward_recur ) k = N
  
  ! Label 500 equivalent - backward recurrence
  il = N - k
  IF( il==0 ) RETURN
  !-----------------------------------------------------------------------
  !     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
  !     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
  !     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
  !-----------------------------------------------------------------------
  s1 = cy(1)
  s2 = cy(2)
  cs = csr(iflag)
  ascle = bry(iflag)
  fn = (inu+il)
  DO i = 1, il
    c2 = s2
    s2 = s1 + CMPLX(fn+fnf,0._DP,DP)*rz*s2
    s1 = c2
    fn = fn - 1._DP
    c2 = s2*cs
    ck = c2
    c1 = Y(kk)
    IF( Kode/=1 ) THEN
      CALL ZS1S2(zr,c1,c2,nw,asc,Alim,iuf)
      Nz = Nz + nw
    END IF
    Y(kk) = c1*cspn + c2
    kk = kk - 1
    cspn = -cspn
    IF( iflag<3 ) THEN
      c2r = REAL(ck,DP)
      c2i = AIMAG(ck)
      c2r = ABS(c2r)
      c2i = ABS(c2i)
      c2m = MAX(c2r,c2i)
      IF( c2m>ascle ) THEN
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1*cs
        s2 = ck
        s1 = s1*css(iflag)
        s2 = s2*css(iflag)
        cs = csr(iflag)
      END IF
    END IF
  END DO
  !
END SUBROUTINE ZUNK1
