!** DQAWOE
PURE SUBROUTINE DQAWOE(F,A,B,Omega,Integr,Epsabs,Epsrel,Limit,Icall,Maxp1,Result,&
    Abserr,Neval,Ier,Last,Alist,Blist,Rlist,Elist,Iord,Nnlog,Momcom,Chebmo)
  USE service, ONLY : tiny_dp, huge_dp, eps_dp
  INTERFACE
    REAL(DP) PURE FUNCTION F(X)
      IMPORT DP
      REAL(DP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(IN) :: Icall, Integr, Limit, Maxp1
  INTEGER, INTENT(INOUT) :: Momcom
  INTEGER, INTENT(OUT) :: Ier, Last, Neval
  INTEGER, INTENT(OUT) :: Iord(Limit), Nnlog(Limit)
  REAL(DP), INTENT(IN) :: A, B, Epsabs, Epsrel, Omega
  REAL(DP), INTENT(OUT) :: Abserr, Result
  REAL(DP), INTENT(INOUT) :: Chebmo(Maxp1,25)
  REAL(DP), INTENT(OUT) :: Alist(Limit), Blist(Limit), Elist(Limit), Rlist(Limit)
  INTEGER :: id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, maxerr, &
    nev, nres, nrmax, nrmom, numrl2
  REAL(DP) :: abseps, area, area1, area12, area2, a1, a2, b1, b2, correc, defab1, &
    defab2, defabs, domega, dres, epmach, erlarg, erlast, errbnd, errmax, error1, &
    erro12, error2, errsum, ertest, oflow, resabs, reseps, res3la(3), rlist2(52), &
    small, uflow, width
  LOGICAL :: extrap, noext, extall, done, test
  epmach = eps_dp
  Ier = 0; Neval = 0; Last = 0; Result = 0._DP; Abserr = 0._DP
  Alist(1) = A; Blist(1) = B; Rlist(1) = 0._DP; Elist(1) = 0._DP
  Iord(1) = 0; Nnlog(1) = 0
  IF( (Integr/=1 .AND. Integr/=2) .OR. &
    (Epsabs<=0._DP .AND. Epsrel<MAX(0.5E+02_DP*epmach,0.5E-28_DP)) .OR. &
    Icall<1 .OR. Maxp1<1 ) Ier = 6
  IF( Ier==6 ) RETURN
  done = .FALSE.
  domega = ABS(Omega); nrmom = 0
  IF( Icall<=1 ) Momcom = 0
  CALL DQC25F(F,A,B,domega,Integr,nrmom,Maxp1,0,Result,Abserr,Neval,defabs,resabs,Momcom,Chebmo)
  dres = ABS(Result); errbnd = MAX(Epsabs,Epsrel*dres)
  Rlist(1) = Result; Elist(1) = Abserr; Iord(1) = 1
  IF( Abserr<=100._DP*epmach*defabs .AND. Abserr>errbnd ) Ier = 2
  IF( Limit==1 ) Ier = 1
  IF( Ier/=0 .OR. Abserr<=errbnd ) THEN
    IF( Integr==2 .AND. Omega<0._DP ) Result = -Result
    RETURN
  END IF
  uflow = tiny_dp; oflow = huge_dp; errmax = Abserr; maxerr = 1
  area = Result; errsum = Abserr; Abserr = oflow; nrmax = 1
  extrap = .FALSE.; noext = .FALSE.; ierro = 0
  iroff1 = 0; iroff2 = 0; iroff3 = 0; ktmin = 0
  small = ABS(B-A)*0.75_DP; nres = 0; numrl2 = 0; extall = .FALSE.
  IF( 0.5_DP*ABS(B-A)*domega<=2._DP ) THEN
    numrl2 = 1; extall = .TRUE.; rlist2(1) = Result
  END IF
  IF( 0.25_DP*ABS(B-A)*domega<=2._DP ) extall = .TRUE.
  ksgn = -1
  IF( dres>=(1._DP-0.5E+02_DP*epmach)*defabs ) ksgn = 1
  main_loop: DO Last = 2, Limit
    nrmom = Nnlog(maxerr) + 1
    a1 = Alist(maxerr); b1 = 0.5_DP*(Alist(maxerr)+Blist(maxerr))
    a2 = b1; b2 = Blist(maxerr); erlast = errmax
    CALL DQC25F(F,a1,b1,domega,Integr,nrmom,Maxp1,0,area1,error1,nev,resabs,defab1,Momcom,Chebmo)
    Neval = Neval + nev
    CALL DQC25F(F,a2,b2,domega,Integr,nrmom,Maxp1,1,area2,error2,nev,resabs,defab2,Momcom,Chebmo)
    Neval = Neval + nev
    area12 = area1 + area2; erro12 = error1 + error2
    errsum = errsum + erro12 - errmax; area = area + area12 - Rlist(maxerr)
    IF( defab1/=error1 .AND. defab2/=error2 ) THEN
      IF( ABS(Rlist(maxerr)-area12)<=0.1D-04*ABS(area12) .AND. erro12>=0.99_DP*errmax ) THEN
        IF( extrap ) iroff2 = iroff2 + 1
        IF( .NOT. extrap ) iroff1 = iroff1 + 1
      END IF
      IF( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
    END IF
    Rlist(maxerr) = area1; Rlist(Last) = area2
    Nnlog(maxerr) = nrmom; Nnlog(Last) = nrmom
    errbnd = MAX(Epsabs,Epsrel*ABS(area))
    IF( iroff1+iroff2>=10 .OR. iroff3>=20 ) Ier = 2
    IF( iroff2>=5 ) ierro = 3
    IF( Last==Limit ) Ier = 1
    IF( MAX(ABS(a1),ABS(b2))<=(1._DP+100._DP*epmach)*(ABS(a2)+0.1D+04*uflow) ) Ier = 4
    IF( error2>error1 ) THEN
      Alist(maxerr) = a2; Alist(Last) = a1; Blist(Last) = b1
      Rlist(maxerr) = area2; Rlist(Last) = area1
      Elist(maxerr) = error2; Elist(Last) = error1
    ELSE
      Alist(Last) = a2; Blist(maxerr) = b1; Blist(Last) = b2
      Elist(maxerr) = error1; Elist(Last) = error2
    END IF
    CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
    IF( errsum<=errbnd ) THEN; done = .TRUE.; EXIT main_loop; END IF
    IF( Ier/=0 ) EXIT main_loop
    IF( Last==2 .AND. extall ) THEN
      small = small*0.5_DP; numrl2 = numrl2 + 1; rlist2(numrl2) = area
    ELSE
      IF( noext ) CYCLE main_loop
      test = .TRUE.
      IF( extall ) THEN
        erlarg = erlarg - erlast
        IF( ABS(b1-a1)>small ) erlarg = erlarg + erro12
        IF( extrap ) test = .FALSE.
      END IF
      IF( test ) THEN
        width = ABS(Blist(maxerr)-Alist(maxerr))
        IF( width>small ) CYCLE main_loop
        IF( extall ) THEN
          extrap = .TRUE.; nrmax = 2
        ELSE
          small = small*0.5_DP
          IF( 0.25_DP*width*domega>2._DP ) CYCLE main_loop
          extall = .TRUE.; ertest = errbnd; erlarg = errsum
          CYCLE main_loop
        END IF
      END IF
      IF( ierro/=3 .AND. erlarg>ertest ) THEN
        jupbnd = Last
        IF( Last>(Limit/2+2) ) jupbnd = Limit + 3 - Last
        id = nrmax
        DO k = id, jupbnd
          maxerr = Iord(nrmax); errmax = Elist(maxerr)
          IF( ABS(Blist(maxerr)-Alist(maxerr))>small ) CYCLE main_loop
          nrmax = nrmax + 1
        END DO
      END IF
      numrl2 = numrl2 + 1; rlist2(numrl2) = area
      IF( numrl2>=3 ) THEN
        CALL DQELG(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin + 1
        IF( ktmin>5 .AND. Abserr<0.1D-02*errsum ) Ier = 5
        IF( abseps<Abserr ) THEN
          ktmin = 0; Abserr = abseps; Result = reseps
          correc = erlarg; ertest = MAX(Epsabs,Epsrel*ABS(reseps))
          IF( Abserr<=ertest ) EXIT main_loop
        END IF
        IF( numrl2==1 ) noext = .TRUE.
        IF( Ier==5 ) EXIT main_loop
      END IF
      maxerr = Iord(1); errmax = Elist(maxerr); nrmax = 1
      extrap = .FALSE.; small = small*0.5_DP; erlarg = errsum
      CYCLE main_loop
    END IF
    ertest = errbnd; erlarg = errsum
  END DO main_loop
  final_result: BLOCK
    IF( done ) EXIT final_result
    IF( Abserr/=oflow .AND. nres/=0 ) THEN
      IF( Ier+ierro/=0 ) THEN
        IF( ierro==3 ) Abserr = Abserr + correc
        IF( Ier==0 ) Ier = 3
        IF( Result==0._DP .OR. area==0._DP ) THEN
          IF( Abserr>errsum ) EXIT final_result
          IF( area==0._DP ) THEN
            IF( Ier>2 ) Ier = Ier - 1
            IF( Integr==2 .AND. Omega<0._DP ) Result = -Result
            RETURN
          END IF
        ELSEIF( Abserr/ABS(Result)>errsum/ABS(area) ) THEN
          EXIT final_result
        END IF
      END IF
      IF( ksgn/=(-1) .OR. MAX(ABS(Result),ABS(area))>defabs*0.1D-01 ) THEN
        IF( 0.1D-01>(Result/area) .OR. (Result/area)>0.1D+03 .OR. errsum>=ABS(area) ) Ier = 6
      END IF
      IF( Ier>2 ) Ier = Ier - 1
      IF( Integr==2 .AND. Omega<0._DP ) Result = -Result
      RETURN
    END IF
  END BLOCK final_result
  Result = 0._DP
  DO k = 1, Last
    Result = Result + Rlist(k)
  END DO
  Abserr = errsum
  IF( Ier>2 ) Ier = Ier - 1
  IF( Integr==2 .AND. Omega<0._DP ) Result = -Result
END SUBROUTINE DQAWOE
