!** POLFIT
PURE SUBROUTINE POLFIT(N,X,Y,W,Maxdeg,Ndeg,Eps,R,Ierr,A)
  !> Fit discrete data in a least squares sense by polynomials in one variable.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  K1A1A2
  !***
  ! **Type:**      SINGLE PRECISION (POLFIT-S, POLFIT-D)
  !***
  ! **Keywords:**  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
  !***
  ! **Author:**  Shampine, L. F., (SNLA)
  !           Davenport, S. M., (SNLA)
  !           Huddleston, R. E., (SNLL)
  !***
  ! **Description:**
  !
  !     Abstract
  !
  !     Given a collection of points X(I) and a set of values Y(I) which
  !     correspond to some function or measurement at each of the X(I),
  !     subroutine  POLFIT  computes the weighted least-squares polynomial
  !     fits of all degrees up to some degree either specified by the user
  !     or determined by the routine.  The fits thus obtained are in
  !     orthogonal polynomial form.  Subroutine  P1VLU  may then be
  !     called to evaluate the fitted polynomials and any of their
  !     derivatives at any point.  The subroutine  PCOEF  may be used to
  !     express the polynomial fits as powers of (X-C) for any specified
  !     point C.
  !
  !     The parameters for  POLFIT  are
  !
  !     Input -- All TYPE REAL variables are SINGLE PRECISION
  !         N -      the number of data points.  The arrays X, Y and W
  !                  must be dimensioned at least  N  (N >= 1).
  !         X -      array of values of the independent variable.  These
  !                  values may appear in any order and need not all be
  !                  distinct.
  !         Y -      array of corresponding function values.
  !         W -      array of positive values to be used as weights.  If
  !                  W(1) is negative,  POLFIT  will set all the weights
  !                  to 1.0, which means unweighted least squares error
  !                  will be minimized.  To minimize relative error, the
  !                  user should set the weights to:  W(I) = 1.0/Y(I)**2,
  !                  I = 1,...,N .
  !         MAXDEG - maximum degree to be allowed for polynomial fit.
  !                  MAXDEG  may be any non-negative integer less than  N.
  !                  Note -- MAXDEG  cannot be equal to  N-1  when a
  !                  statistical test is to be used for degree selection,
  !                  i.e., when input value of  EPS  is negative.
  !         EPS -    specifies the criterion to be used in determining
  !                  the degree of fit to be computed.
  !                  (1)  If  EPS  is input negative,  POLFIT  chooses the
  !                       degree based on a statistical F test of
  !                       significance.  One of three possible
  !                       significance levels will be used:  .01, .05 or
  !                       .10.  If  EPS=-1.0, the routine will
  !                       automatically select one of these levels based
  !                       on the number of data points and the maximum
  !                       degree to be considered.  If  EPS  is input as
  !                       -.01, -.05, or -.10, a significance level of
  !                       .01, .05, or .10, respectively, will be used.
  !                  (2)  If  EPS  is set to 0.,  POLFIT  computes the
  !                       polynomials of degrees 0 through  MAXDEG .
  !                  (3)  If  EPS  is input positive,  EPS  is the RMS
  !                       error tolerance which must be satisfied by the
  !                       fitted polynomial.  POLFIT  will increase the
  !                       degree of fit until this criterion is met or
  !                       until the maximum degree is reached.
  !
  !     Output -- All TYPE REAL variables are SINGLE PRECISION
  !         NDEG -   degree of the highest degree fit computed.
  !         EPS -    RMS error of the polynomial of degree  NDEG .
  !         R -      vector of dimension at least NDEG containing values
  !                  of the fit of degree  NDEG  at each of the  X(I) .
  !                  Except when the statistical test is used, these
  !                  values are more accurate than results from subroutine
  !                  P1VLU  normally are.
  !         IERR -   error flag with the following possible values.
  !             1 -- indicates normal execution, i.e., either
  !                  (1)  the input value of  EPS  was negative, and the
  !                       computed polynomial fit of degree  NDEG
  !                       satisfies the specified F test, or
  !                  (2)  the input value of  EPS  was 0., and the fits of
  !                       all degrees up to  MAXDEG  are complete, or
  !                  (3)  the input value of  EPS  was positive, and the
  !                       polynomial of degree  NDEG  satisfies the RMS
  !                       error requirement.
  !             2 -- invalid input parameter.  At least one of the input
  !                  parameters has an illegal value and must be corrected
  !                  before  POLFIT  can proceed.  Valid input results
  !                  when the following restrictions are observed
  !                       N >= 1
  !                       0 <= MAXDEG <= N-1  for  EPS >= 0.
  !                       0 <= MAXDEG <= N-2  for  EPS < 0.
  !                       W(1)=-1.0  or  W(I) > 0., I=1,...,N .
  !             3 -- cannot satisfy the RMS error requirement with a
  !                  polynomial of degree no greater than  MAXDEG .  Best
  !                  fit found is of degree  MAXDEG .
  !             4 -- cannot satisfy the test for significance using
  !                  current value of  MAXDEG .  Statistically, the
  !                  best fit found is of order  NORD .  (In this case,
  !                  NDEG will have one of the values:  MAXDEG-2,
  !                  MAXDEG-1, or MAXDEG).  Using a higher value of
  !                  MAXDEG  may result in passing the test.
  !         A -      work and output array having at least 3N+3MAXDEG+3
  !                  locations
  !
  !     Note - POLFIT  calculates all fits of degrees up to and including
  !            NDEG .  Any or all of these fits can be evaluated or
  !            expressed as powers of (X-C) using  P1VLU  and  PCOEF
  !            after just one call to  POLFIT .
  !
  !***
  ! **References:**  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
  !                 Curve fitting by polynomials in one variable, Report
  !                 SLA-74-0270, Sandia Laboratories, June 1974.
  !***
  ! **Routines called:**  P1VLU, XERMSG
  !* REVISION HISTORY  (YYMMDD)
  !   740601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900911  Added variable YP to SINGLE PRECISION declaration.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   920527  Corrected erroneous statements in DESCRIPTION.  (WRB)
  !   2024xx  Eliminated GOTOs using structured control flow.  (modernization)
  !
  ! MODERNIZATION NOTES:
  ! --------------------
  ! The original code used GOTOs for three distinct purposes:
  !   1. Input validation jumps to error return (label 700)
  !   2. Mode selection based on EPS sign (labels 200, 300, 400)
  !   3. Common exit handling (labels 500-900)
  !
  ! This version uses:
  !   - Early RETURN for input validation errors
  !   - INTEGER mode variable with SELECT CASE for EPS modes
  !   - Logical flags (done, f_test_passed) for loop control
  !   - common_exit: BLOCK for shared cleanup code
  !
  ! Algorithm: Forsythe's method for orthogonal polynomial fitting
  ! using three-term recurrence. See SLA-74-0270 for details.
  !
  INTEGER, INTENT(IN) :: Maxdeg, N
  INTEGER, INTENT(OUT) :: Ierr, Ndeg
  REAL(SP), INTENT(INOUT) :: Eps
  REAL(SP), INTENT(IN) :: X(N), Y(N)
  REAL(SP), INTENT(INOUT) :: A(3*(N+Maxdeg+1)), W(N)
  REAL(SP), INTENT(OUT) :: R(:)
  !
  INTEGER :: i, idegf, j, jp1, jpas, k1, k1pj, k2, k2pj, k3, k3pi, k4, k4pi, k5, &
    k5pi, ksig, m, mop1, nder, nfail
  REAL(SP) :: temd1, temd2, degf, den, etst, f, fcrit, sig, sigj, sigjm1, sigpas, &
    temp, xm, yp(1), w1, w11
  REAL(SP), PARAMETER :: co(4,3)= RESHAPE( [ &
    -13.086850_SP, -2.4648165_SP, -3.3846535_SP, -1.2973162_SP, &
    -3.3381146_SP, -1.7812271_SP, -3.2578406_SP, -1.6589279_SP, &
    -1.6282703_SP, -1.3152745_SP, -3.2640179_SP, -1.9829776_SP ], [4,3] )
  ! Mode constants for degree selection criterion
  INTEGER, PARAMETER :: MODE_FTEST = 1    ! EPS < 0: Statistical F-test
  INTEGER, PARAMETER :: MODE_ALLDEG = 2   ! EPS = 0: Compute all degrees
  INTEGER, PARAMETER :: MODE_RMS = 3      ! EPS > 0: RMS error criterion
  INTEGER :: mode
  LOGICAL :: done, f_test_passed
  !* FIRST EXECUTABLE STATEMENT  POLFIT
  m = ABS(N)
  yp = 0._SP

  !---------------------------------------------------------------------------
  ! INPUT VALIDATION
  ! Check all input parameters and return with Ierr=2 if any are invalid.
  !---------------------------------------------------------------------------
  IF( m==0 ) THEN
    Ierr = 2
    ERROR STOP 'POLFIT : INVALID INPUT PARAMETER.'
    RETURN
  END IF
  IF( Maxdeg<0 ) THEN
    Ierr = 2
    ERROR STOP 'POLFIT : INVALID INPUT PARAMETER.'
    RETURN
  END IF
  A(1) = Maxdeg
  mop1 = Maxdeg + 1
  IF( m<mop1 ) THEN
    Ierr = 2
    ERROR STOP 'POLFIT : INVALID INPUT PARAMETER.'
    RETURN
  END IF
  IF( Eps<0._SP .AND. m==mop1 ) THEN
    Ierr = 2
    ERROR STOP 'POLFIT : INVALID INPUT PARAMETER.'
    RETURN
  END IF
  xm = m
  etst = Eps*Eps*xm

  ! Validate or initialize weights
  IF( W(1)<0._SP ) THEN
    ! Negative W(1) signals uniform weights
    DO i = 1, m
      W(i) = 1._SP
    END DO
  ELSE
    ! Check all weights are positive
    DO i = 1, m
      IF( W(i)<=0._SP ) THEN
        Ierr = 2
        ERROR STOP 'POLFIT : INVALID INPUT PARAMETER.'
        RETURN
      END IF
    END DO
  END IF

  !---------------------------------------------------------------------------
  ! DETERMINE FITTING MODE AND SIGNIFICANCE LEVEL
  ! Mode is based on the sign of EPS:
  !   EPS < 0: Use F-test for degree selection (MODE_FTEST)
  !   EPS = 0: Compute all degrees up to MAXDEG (MODE_ALLDEG)
  !   EPS > 0: Stop when RMS error < EPS (MODE_RMS)
  !---------------------------------------------------------------------------
  IF( Eps<0._SP ) THEN
    mode = MODE_FTEST
    ! Determine significance level index for F-test
    IF( Eps>(-.55_SP) ) THEN
      ksig = 1
      IF( Eps<(-.03_SP) ) ksig = 2
      IF( Eps<(-.07_SP) ) ksig = 3
    ELSE
      idegf = m - Maxdeg - 1
      ksig = 1
      IF( idegf<10 ) ksig = 2
      IF( idegf<5 ) ksig = 3
    END IF
  ELSEIF( Eps==0._SP ) THEN
    mode = MODE_ALLDEG
  ELSE
    mode = MODE_RMS
  END IF

  !---------------------------------------------------------------------------
  ! INITIALIZE INDEXES AND COEFFICIENTS
  ! Set up array partitioning for orthogonal polynomial computation.
  ! A is partitioned as: [alpha coeffs | beta coeffs | c coeffs | P_j | P_{j-1}]
  !---------------------------------------------------------------------------
  k1 = Maxdeg + 1
  k2 = k1 + Maxdeg
  k3 = k2 + Maxdeg + 2
  k4 = k3 + m
  k5 = k4 + m
  DO i = 2, k4
    A(i) = 0._SP
  END DO
  w11 = 0._SP

  IF( N<0 ) THEN
    ! Constrained case (N passed as negative)
    DO i = 1, m
      k4pi = k4 + i
      w11 = w11 + W(i)*A(k4pi)**2
    END DO
  ELSE
    ! Unconstrained case (normal)
    DO i = 1, m
      k4pi = k4 + i
      A(k4pi) = 1._SP
      w11 = w11 + W(i)
    END DO
  END IF

  !---------------------------------------------------------------------------
  ! COMPUTE FIT OF DEGREE ZERO
  ! This is the weighted mean: sum(w*y*P0) / sum(w*P0^2)
  !---------------------------------------------------------------------------
  temd1 = 0._SP
  DO i = 1, m
    k4pi = k4 + i
    temd1 = temd1 + W(i)*Y(i)*A(k4pi)
  END DO
  temd1 = temd1/w11
  A(k2+1) = temd1
  sigj = 0._SP
  DO i = 1, m
    k4pi = k4 + i
    k5pi = k5 + i
    temd2 = temd1*A(k4pi)
    R(i) = temd2
    A(k5pi) = temd2 - R(i)
    sigj = sigj + W(i)*((Y(i)-R(i))-A(k5pi))**2
  END DO
  j = 0

  !---------------------------------------------------------------------------
  ! CHECK IF DEGREE 0 SATISFIES CRITERION
  ! Before entering the main loop, check if we're already done.
  !---------------------------------------------------------------------------
  done = .FALSE.
  nfail = 0
  jpas = 0
  sigpas = sigj

  SELECT CASE( mode )
  CASE( MODE_FTEST )
    ! F-test mode: degree 0 always "passes" as initial reference
    f_test_passed = .TRUE.
    sigpas = sigj
    jpas = 0
    IF( Maxdeg==0 ) THEN
      ! Can't do F-test with only degree 0
      Ierr = 4
      Ndeg = jpas
      sig = sigpas
      done = .TRUE.
    END IF
  CASE( MODE_ALLDEG )
    ! Compute all degrees: continue unless MAXDEG=0
    IF( Maxdeg==0 ) THEN
      Ierr = 1
      Ndeg = 0
      sig = sigj
      done = .TRUE.
    END IF
  CASE( MODE_RMS )
    ! RMS mode: check if degree 0 already satisfies criterion
    IF( sigj<=etst ) THEN
      Ierr = 1
      Ndeg = 0
      sig = sigj
      done = .TRUE.
    ELSEIF( Maxdeg==0 ) THEN
      Ierr = 3
      Ndeg = 0
      sig = sigj
      done = .TRUE.
    END IF
  END SELECT

  !---------------------------------------------------------------------------
  ! MAIN DEGREE INCREMENT LOOP
  ! Increment polynomial degree using three-term recurrence for orthogonal
  ! polynomials. Continue until stopping criterion is met or MAXDEG reached.
  !---------------------------------------------------------------------------
  degree_loop: DO WHILE( .NOT. done )
    ! Increment degree
    j = j + 1
    jp1 = j + 1
    k1pj = k1 + j
    k2pj = k2 + j
    sigjm1 = sigj

    ! Compute new B coefficient (beta in three-term recurrence)
    ! Skip for j=1 since P_{-1} = 0
    IF( j>1 ) A(k1pj) = w11/w1

    ! Compute new A coefficient (alpha in three-term recurrence)
    ! alpha_j = sum(x * w * P_j^2) / sum(w * P_j^2)
    temd1 = 0._SP
    DO i = 1, m
      k4pi = k4 + i
      temd2 = A(k4pi)
      temd1 = temd1 + X(i)*W(i)*temd2*temd2
    END DO
    A(jp1) = temd1/w11

    ! Evaluate orthogonal polynomial P_{j+1} at all data points
    ! P_{j+1}(x) = (x - alpha_j)*P_j(x) - beta_j*P_{j-1}(x)
    w1 = w11
    w11 = 0._SP
    DO i = 1, m
      k3pi = k3 + i
      k4pi = k4 + i
      temp = A(k3pi)
      A(k3pi) = A(k4pi)
      A(k4pi) = (X(i)-A(jp1))*A(k3pi) - A(k1pj)*temp
      w11 = w11 + W(i)*A(k4pi)**2
    END DO

    ! Get new orthogonal polynomial coefficient using partial single precision
    ! c_j = sum(w * (y - R) * P_j) / sum(w * P_j^2)
    temd1 = 0._SP
    DO i = 1, m
      k4pi = k4 + i
      k5pi = k5 + i
      temd2 = W(i)*((Y(i)-R(i))-A(k5pi))*A(k4pi)
      temd1 = temd1 + temd2
    END DO
    temd1 = temd1/w11
    A(k2pj+1) = temd1

    ! Update polynomial evaluations and accumulate sum of squared errors
    ! R(i) stores the fit value, A(k5pi) stores the low-order bits
    sigj = 0._SP
    DO i = 1, m
      k4pi = k4 + i
      k5pi = k5 + i
      temd2 = R(i) + A(k5pi) + temd1*A(k4pi)
      R(i) = temd2
      A(k5pi) = temd2 - R(i)
      sigj = sigj + W(i)*((Y(i)-R(i))-A(k5pi))**2
    END DO

    !-------------------------------------------------------------------------
    ! CHECK DEGREE SELECTION CRITERION
    !-------------------------------------------------------------------------
    SELECT CASE( mode )

    CASE( MODE_FTEST )
      ! Statistical F-test for significance
      IF( sigj==0._SP ) THEN
        ! Perfect fit achieved
        Ierr = 1
        Ndeg = jpas
        sig = sigpas
        done = .TRUE.
        CYCLE degree_loop
      END IF

      ! Compute F statistic and critical value
      degf = m - j - 1
      den = (co(4,ksig)*degf+1._SP)*degf
      fcrit = (((co(3,ksig)*degf)+co(2,ksig))*degf+co(1,ksig))/den
      fcrit = fcrit*fcrit
      f = (sigjm1-sigj)*degf/sigj

      IF( f>=fcrit ) THEN
        ! Polynomial of degree j passes F-test
        sigpas = sigj
        jpas = j
        nfail = 0
        IF( Maxdeg==j ) THEN
          ! Reached MAXDEG, cannot continue test
          Ierr = 4
          Ndeg = jpas
          sig = sigpas
          done = .TRUE.
        END IF
      ELSE
        ! Polynomial of degree j fails F-test
        nfail = nfail + 1
        IF( nfail>=3 ) THEN
          ! Three successive failures: statistically best degree found
          Ierr = 1
          Ndeg = jpas
          sig = sigpas
          done = .TRUE.
        ELSEIF( Maxdeg==j ) THEN
          ! Reached MAXDEG with failure
          Ierr = 4
          Ndeg = jpas
          sig = sigpas
          done = .TRUE.
        END IF
      END IF

    CASE( MODE_ALLDEG )
      ! Compute all degrees up to MAXDEG
      IF( Maxdeg==j ) THEN
        Ierr = 1
        Ndeg = j
        sig = sigj
        done = .TRUE.
      END IF

    CASE( MODE_RMS )
      ! RMS error criterion
      IF( sigj<=etst ) THEN
        ! RMS error satisfied
        Ierr = 1
        Ndeg = j
        sig = sigj
        done = .TRUE.
      ELSEIF( Maxdeg==j ) THEN
        ! Reached MAXDEG without satisfying criterion
        Ierr = 3
        Ndeg = Maxdeg
        sig = sigj
        done = .TRUE.
      END IF

    END SELECT
  END DO degree_loop

  !---------------------------------------------------------------------------
  ! COMMON EXIT PROCESSING
  ! Store final degree and optionally re-evaluate polynomial at data points.
  !---------------------------------------------------------------------------
  A(k3) = Ndeg

  ! When statistical test was used and we didn't reach MAXDEG,
  ! re-evaluate the best polynomial at all data points
  IF( Eps<0.0_SP .AND. Ndeg/=Maxdeg ) THEN
    nder = 0
    DO i = 1, m
      CALL P1VLU(Ndeg,nder,X(i),R(i),yp,A)
    END DO
  END IF

  ! Convert sum of squared errors to RMS error
  Eps = SQRT(sig/xm)

END SUBROUTINE POLFIT
