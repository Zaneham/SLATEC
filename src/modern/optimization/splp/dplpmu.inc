!** DPLPMU
SUBROUTINE DPLPMU(Mrelas,Nvars,Lmx,Lbm,Nredc,Info,Ienter,Ileave,Npp,&
    Jstrt,Ibasis,Imat,Ibrc,Ipr,Iwr,Ind,Ibb,Anorm,Eps,Uu,Gg,Rprnrm,Erdnrm,Dulnrm,&
    Theta,Costsc,Xlamda,Rhsnrm,Amat,Basmat,Csc,Wr,Rprim,Ww,Bu,Bl,Rhs,Erd,Erp,Rz,&
    Rg,Colnrm,Costs,Primal,Duals,Singlr,Redbas,Zerolv,Stpedg)
  !> Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SPLPMU-S, DPLPMU-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
  !     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
  !
  !     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
  !     /REAL (12 BLANKS)/DOUBLE PRECISION/,
  !     /SASUM/DASUM/,/SCOPY/DCOPY/,/SDOT/DDOT/,
  !     /.E0/.D0/
  !
  !     THIS SUBPROGRAM IS FROM THE DSPLP( ) PACKAGE.  IT PERFORMS THE
  !     TASKS OF UPDATING THE PRIMAL SOLUTION, EDGE WEIGHTS, REDUCED
  !     COSTS, AND MATRIX DECOMPOSITION.
  !     IT IS THE MAIN PART OF THE PROCEDURE (MAKE MOVE AND UPDATE).
  !
  !     REVISED 821122-1100
  !     REVISED YYMMDD
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  DASUM, DCOPY, DDOT, DPLPDM, DPNNZR, DPRWPG, IDLOC,
  !                    LA05BD, LA05CD, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   890606  Changed references from IPLOC to IDLOC.  (WRB)
  !   890606  Removed unused COMMON block LA05DD.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !
  ! Modernisation notes (202412):
  !   Eliminated GOTOs using internal subroutines and logical flags.
  !   Original F77 used computed GOTO pattern with npr001/npr003.
  !   Extracted "procedures" 1800/1900 as internal subroutines.

  INTEGER, INTENT(IN) :: Ienter, Ileave, Lbm, Lmx, Mrelas, Npp, Nvars
  INTEGER, INTENT(INOUT) :: Jstrt
  INTEGER, INTENT(OUT) :: Info, Nredc
  REAL(DP), INTENT(IN) :: Costsc, Erdnrm, Eps, Xlamda
  REAL(DP), INTENT(INOUT) :: Dulnrm, Gg, Rhsnrm, Rprnrm, Theta, Uu
  REAL(DP), INTENT(OUT) :: Anorm
  LOGICAL, INTENT(IN) :: Stpedg, Zerolv
  LOGICAL, INTENT(OUT) :: Redbas, Singlr
  INTEGER, INTENT(IN) :: Imat(Lmx), Ind(Nvars+Mrelas)
  INTEGER, INTENT(INOUT) :: Ibasis(Nvars+Mrelas), Ibb(Nvars+Mrelas), Ibrc(Lbm,2), &
    Ipr(2*Mrelas), Iwr(8*Mrelas)
  REAL(DP), INTENT(IN) :: Amat(Lmx), Csc(Nvars), Bu(Nvars+Mrelas), Bl(Nvars+Mrelas), &
    Erp(Mrelas), Costs(Nvars), Colnrm(Nvars)
  REAL(DP), INTENT(INOUT) :: Basmat(Lbm), Primal(Nvars+Mrelas), Rg(Nvars+Mrelas), &
    Rhs(Mrelas), Rprim(Mrelas), Rz(Nvars+Mrelas), Ww(Mrelas)
  REAL(DP), INTENT(OUT) :: Duals(Nvars+Mrelas), Erd(Mrelas), Wr(Mrelas)
  INTEGER :: i, ibas, ihi, il1, ilow, ipage, iplace, iu1, j, k, lpg, &
    nerr, nnegrc
  REAL(DP) :: aij, alpha, gama, gq, rzj, scalr, wp, cnorm, rcost
  LOGICAL :: pagepl, trans
  !
  !* FIRST EXECUTABLE STATEMENT  DPLPMU
  lpg = Lmx - (Nvars+4)
  !
  !     UPDATE THE PRIMAL SOLUTION WITH A MULTIPLE OF THE SEARCH
  !     DIRECTION.
  DO i = 1, Mrelas
    Rprim(i) = Rprim(i) - Theta*Ww(i)
  END DO
  !
  !     IF EJECTED VARIABLE IS LEAVING AT AN UPPER BOUND,  THEN
  !     TRANSLATE RIGHT HAND SIDE.
  IF( Ileave<0 ) THEN
    ibas = Ibasis(ABS(Ileave))
    scalr = Rprim(ABS(Ileave))
    CALL translate_rhs()
    Ibb(ibas) = ABS(Ibb(ibas)) + 1
  END IF
  !
  !     IF ENTERING VARIABLE IS RESTRICTED TO ITS UPPER BOUND, TRANSLATE
  !     RIGHT HAND SIDE.  IF THE VARIABLE DECREASED FROM ITS UPPER
  !     BOUND, A SIGN CHANGE IS REQUIRED IN THE TRANSLATION.
  IF( Ienter/=Ileave ) THEN
    ibas = Ibasis(Ienter)
    !
    !     IF ENTERING VARIABLE IS DECREASING FROM ITS UPPER BOUND,
    !     COMPLEMENT ITS PRIMAL VALUE.
    IF( Ind(ibas)==3 .AND. MOD(Ibb(ibas),2)==0 ) THEN
      scalr = -(Bu(ibas)-Bl(ibas))
      IF( ibas<=Nvars ) scalr = scalr/Csc(ibas)
      CALL translate_rhs()
      Theta = -scalr - Theta
      Ibb(ibas) = Ibb(ibas) + 1
    END IF
    Rprim(ABS(Ileave)) = Theta
    Ibb(ibas) = -ABS(Ibb(ibas))
    i = Ibasis(ABS(Ileave))
    Ibb(i) = ABS(Ibb(i))
    IF( Primal(ABS(Ileave)+Nvars)>0._DP ) Ibb(i) = Ibb(i) + 1
  ELSE
    ibas = Ibasis(Ienter)
    scalr = Theta
    IF( MOD(Ibb(ibas),2)==0 ) scalr = -scalr
    CALL translate_rhs()
    Ibb(ibas) = Ibb(ibas) + 1
  END IF
  !
  !     INTERCHANGE COLUMN POINTERS TO NOTE EXCHANGE OF COLUMNS.
  ibas = Ibasis(Ienter)
  Ibasis(Ienter) = Ibasis(ABS(Ileave))
  Ibasis(ABS(Ileave)) = ibas
  !
  !     IF VARIABLE WAS EXCHANGED AT A ZERO LEVEL, MARK IT SO THAT
  !     IT CAN'T BE BROUGHT BACK IN.  THIS IS TO HELP PREVENT CYCLING.
  IF( Zerolv ) Ibasis(Ienter) = -ABS(Ibasis(Ienter))
  Rprnrm = MAX(Rprnrm,SUM(ABS(Rprim(1:Mrelas))) )
  !
  !     SEE IF VARIABLES THAT WERE CLASSIFIED AS INFEASIBLE HAVE NOW
  !     BECOME FEASIBLE.  THIS MAY REQUIRED TRANSLATING UPPER BOUNDED
  !     VARIABLES.
  k = 1
  feasibility_check: DO WHILE( k<=Mrelas )
    IF( Primal(k+Nvars)==0._DP .OR. ABS(Rprim(k))>Rprnrm*Erp(k) ) THEN
      k = k + 1
    ELSE
      IF( Primal(k+Nvars)>0._DP ) THEN
        ibas = Ibasis(k)
        scalr = -(Bu(ibas)-Bl(ibas))
        IF( ibas<=Nvars ) scalr = scalr/Csc(ibas)
        CALL translate_rhs()
        Rprim(k) = -scalr
        Rprnrm = Rprnrm - scalr
      END IF
      Primal(k+Nvars) = 0._DP
      k = k + 1
    END IF
  END DO feasibility_check
  !
  !     UPDATE REDUCED COSTS, EDGE WEIGHTS, AND MATRIX DECOMPOSITION.
  IF( Ienter==Ileave ) THEN
    !
    !     THIS IS NECESSARY ONLY FOR PRINTING OF INTERMEDIATE RESULTS.
    CALL compute_new_duals()
    RETURN
  END IF
  !
  !     THE INCOMING VARIABLE IS ALWAYS CLASSIFIED AS FEASIBLE.
  Primal(ABS(Ileave)+Nvars) = 0._DP
  !
  wp = Ww(ABS(Ileave))
  gq = NORM2(Ww(1:Mrelas))**2 + 1._DP
  !
  !     COMPUTE INVERSE (TRANSPOSE) TIMES SEARCH DIRECTION.
  trans = .TRUE.
  CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
  !
  !     UPDATE THE MATRIX DECOMPOSITION.  COL. ABS(ILEAVE) IS LEAVING.
  !     THE ARRAY DUALS(*) CONTAINS INTERMEDIATE RESULTS FOR THE
  !     INCOMING COLUMN.
  CALL LA05CD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Duals,Gg,Uu,ABS(Ileave))
  Redbas = .FALSE.
  IF( Gg<0._DP ) THEN
    !
    !     REDECOMPOSE BASIS MATRIX WHEN AN ERROR RETURN FROM
    !     LA05CD( ) IS NOTED.  THIS WILL PROBABLY BE DUE TO
    !     SPACE BEING EXHAUSTED, GG=-7.
    CALL DPLPDM(Mrelas,Nvars,Lbm,Nredc,Info,Ibasis,Imat,Ibrc,Ipr,&
      Iwr,Ind,Anorm,Eps,Uu,Gg,Amat,Basmat,Csc,Wr,Singlr,Redbas)
    IF( Singlr ) THEN
      nerr = 26
      Info = -nerr
      ERROR STOP 'DPLPMU : IN DSPLP, MOVED TO A SINGULAR POINT.  THIS SHOULD NOT HAPPEN.'
      RETURN
    END IF
    !     PROCEDURE (COMPUTE NEW PRIMAL)
    !
    !     COPY RHS INTO WW(*), SOLVE SYSTEM.
    Ww(1:Mrelas) = Rhs(1:Mrelas)
    trans = .FALSE.
    CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
    Rprim(1:Mrelas) = Ww(1:Mrelas)
    Rprnrm = SUM(ABS(Rprim(1:Mrelas)))
  END IF
  !
  !     IF STEEPEST EDGE PRICING IS USED, UPDATE REDUCED COSTS
  !     AND EDGE WEIGHTS.
  IF( Stpedg ) THEN
    !
    !     COMPUTE COL. ABS(ILEAVE) OF THE NEW INVERSE (TRANSPOSE) MATRIX
    !     HERE ABS(ILEAVE) POINTS TO THE EJECTED COLUMN.
    !     USE ERD(*) FOR TEMP. STORAGE.
    Erd(1:Mrelas) = 0._DP
    Erd(ABS(Ileave)) = 1._DP
    trans = .TRUE.
    CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Erd,trans)
    !
    !     COMPUTE UPDATED DUAL VARIABLES IN DUALS(*).
    CALL compute_new_duals()
    !
    !     COMPUTE THE DOT PRODUCT OF COL. J OF THE NEW INVERSE (TRANSPOSE)
    !     WITH EACH NON-BASIC COLUMN.  ALSO COMPUTE THE DOT PRODUCT OF THE
    !     INVERSE (TRANSPOSE) OF NON-UPDATED MATRIX (TIMES) THE
    !     SEARCH DIRECTION WITH EACH NON-BASIC COLUMN.
    !     RECOMPUTE REDUCED COSTS.
    pagepl = .TRUE.
    Rz(1:Nvars+Mrelas) = 0._DP
    nnegrc = 0
    j = Jstrt
    steepest_edge: DO
      IF( Ibb(j)<=0 ) THEN
        pagepl = .TRUE.
        Rg(j) = 1._DP
        !
        !     NONBASIC INDEPENDENT VARIABLES (COLUMN IN SPARSE MATRIX STORAGE)
      ELSEIF( j>Nvars ) THEN
        pagepl = .TRUE.
        scalr = -1._DP
        IF( Ind(j)==2 ) scalr = 1._DP
        i = j - Nvars
        alpha = scalr*Erd(i)
        Rz(j) = -scalr*Duals(i)
        gama = scalr*Ww(i)
        Rg(j) = MAX(Rg(j)-2._DP*alpha*gama+alpha**2*gq,1._DP+alpha**2)
      ELSE
        rzj = Costs(j)*Costsc
        alpha = 0._DP
        gama = 0._DP
        !
        !     COMPUTE THE DOT PRODUCT OF THE SPARSE MATRIX NONBASIC COLUMNS
        !     WITH THREE VECTORS INVOLVED IN THE UPDATING STEP.
        IF( j/=1 ) THEN
          ilow = Imat(j+3) + 1
        ELSE
          ilow = Nvars + 5
        END IF
        IF( .NOT. pagepl ) THEN
          il1 = ihi + 1
        ELSE
          il1 = IDLOC(ilow,Imat)
          IF( il1>=Lmx-1 ) THEN
            ilow = ilow + 2
            il1 = IDLOC(ilow,Imat)
          END IF
          ipage = ABS(Imat(Lmx-1))
        END IF
        ihi = Imat(j+4) - (ilow-il1)
        DO
          iu1 = MIN(Lmx-2,ihi)
          IF( il1>iu1 ) EXIT
          DO i = il1, iu1
            rzj = rzj - Amat(i)*Duals(Imat(i))
            alpha = alpha + Amat(i)*Erd(Imat(i))
            gama = gama + Amat(i)*Ww(Imat(i))
          END DO
          IF( ihi<=Lmx-2 ) EXIT
          il1 = Nvars + 5
          ihi = ihi - lpg
        END DO
        pagepl = ihi==(Lmx-2)
        Rz(j) = rzj*Csc(j)
        alpha = alpha*Csc(j)
        gama = gama*Csc(j)
        !
        !     NONBASIC DEPENDENT VARIABLES (COLUMNS DEFINED IMPLICITLY)
        Rg(j) = MAX(Rg(j)-2._DP*alpha*gama+alpha**2*gq,1._DP+alpha**2)
      END IF
      !
      rcost = Rz(j)
      IF( MOD(Ibb(j),2)==0 ) rcost = -rcost
      IF( Ind(j)==3 ) THEN
        IF( Bu(j)==Bl(j) ) rcost = 0._DP
      END IF
      IF( Ind(j)==4 ) rcost = -ABS(rcost)
      cnorm = 1._DP
      IF( j<=Nvars ) cnorm = Colnrm(j)
      IF( rcost+Erdnrm*Dulnrm*cnorm<0._DP ) nnegrc = nnegrc + 1
      j = MOD(j,Mrelas+Nvars) + 1
      IF( nnegrc>=Npp .OR. j==Jstrt ) EXIT steepest_edge
    END DO steepest_edge
    Jstrt = j
    !
    !     UPDATE THE EDGE WEIGHT FOR THE EJECTED VARIABLE.
    Rg(ABS(Ibasis(Ienter))) = gq/wp**2
  ELSE
    !
    !     IF MINIMUM REDUCED COST (DANTZIG) PRICING IS USED,
    !     CALCULATE THE NEW REDUCED COSTS.
    !
    !     COMPUTE THE UPDATED DUALS IN DUALS(*).
    CALL compute_new_duals()
    !
    Rz(1:Nvars+Mrelas) = 0._DP
    nnegrc = 0
    j = Jstrt
    pagepl = .TRUE.
    !
    dantzig_pricing: DO
      IF( Ibb(j)<=0 ) THEN
        pagepl = .TRUE.
        !
        !     NONBASIC INDEPENDENT VARIABLE (COLUMN IN SPARSE MATRIX STORAGE)
      ELSEIF( j>Nvars ) THEN
        pagepl = .TRUE.
        scalr = -1._DP
        IF( Ind(j)==2 ) scalr = 1._DP
        i = j - Nvars
        Rz(j) = -scalr*Duals(i)
      ELSE
        Rz(j) = Costs(j)*Costsc
        IF( j/=1 ) THEN
          ilow = Imat(j+3) + 1
        ELSE
          ilow = Nvars + 5
        END IF
        IF( .NOT. pagepl ) THEN
          il1 = ihi + 1
        ELSE
          il1 = IDLOC(ilow,Imat)
          IF( il1>=Lmx-1 ) THEN
            ilow = ilow + 2
            il1 = IDLOC(ilow,Imat)
          END IF
          ipage = ABS(Imat(Lmx-1))
        END IF
        ihi = Imat(j+4) - (ilow-il1)
        dot_product_loop: DO
          iu1 = MIN(Lmx-2,ihi)
          IF( iu1>=il1 .AND. MOD(iu1-il1,2)==0 ) THEN
            Rz(j) = Rz(j) - Amat(il1)*Duals(Imat(il1))
            il1 = il1 + 1
          END IF
          IF( il1<=iu1 ) THEN
            !
            !     UNROLL THE DOT PRODUCT LOOP TO A DEPTH OF TWO.  (THIS IS DONE
            !     FOR INCREASED EFFICIENCY).
            DO i = il1, iu1, 2
              Rz(j) = Rz(j) - Amat(i)*Duals(Imat(i)) - Amat(i+1)*Duals(Imat(i+1))
            END DO
            IF( ihi>Lmx-2 ) THEN
              il1 = Nvars + 5
              ihi = ihi - lpg
              CYCLE dot_product_loop
            END IF
          END IF
          EXIT dot_product_loop
        END DO dot_product_loop
        pagepl = ihi==(Lmx-2)
        !
        !     NONBASIC DEPENDENT VARIABLES (COLUMNS DEFINED IMPLICITLY)
        Rz(j) = Rz(j)*Csc(j)
      END IF
      !
      rcost = Rz(j)
      IF( MOD(Ibb(j),2)==0 ) rcost = -rcost
      IF( Ind(j)==3 ) THEN
        IF( Bu(j)==Bl(j) ) rcost = 0._DP
      END IF
      IF( Ind(j)==4 ) rcost = -ABS(rcost)
      cnorm = 1._DP
      IF( j<=Nvars ) cnorm = Colnrm(j)
      IF( rcost+Erdnrm*Dulnrm*cnorm<0._DP ) nnegrc = nnegrc + 1
      j = MOD(j,Mrelas+Nvars) + 1
      IF( nnegrc>=Npp .OR. j==Jstrt ) EXIT dantzig_pricing
    END DO dantzig_pricing
    Jstrt = j
  END IF
  !
CONTAINS
  !
  !     PROCEDURE (TRANSLATE RIGHT HAND SIDE)
  !
  !     PERFORM THE TRANSLATION ON THE RIGHT-HAND SIDE.
  SUBROUTINE translate_rhs()
    IF( ibas>Nvars ) THEN
      i = ibas - Nvars
      IF( Ind(ibas)/=2 ) THEN
        Rhs(i) = Rhs(i) + scalr
      ELSE
        Rhs(i) = Rhs(i) - scalr
      END IF
    ELSE
      i = 0
      DO
        CALL DPNNZR(i,aij,iplace,Amat,Imat,ibas)
        IF( i<=0 ) EXIT
        Rhs(i) = Rhs(i) - scalr*aij*Csc(ibas)
      END DO
    END IF
    Rhsnrm = MAX(Rhsnrm,SUM(ABS(Rhs(1:Mrelas))) )
  END SUBROUTINE translate_rhs
  !
  !     PROCEDURE (COMPUTE NEW DUALS)
  !
  !     SOLVE FOR DUAL VARIABLES. FIRST COPY COSTS INTO DUALS(*).
  SUBROUTINE compute_new_duals()
    INTEGER :: ii, jj
    DO ii = 1, Mrelas
      jj = Ibasis(ii)
      IF( jj>Nvars ) THEN
        Duals(ii) = Xlamda*Primal(ii+Nvars)
      ELSE
        Duals(ii) = Costsc*Costs(jj)*Csc(jj) + Xlamda*Primal(ii+Nvars)
      END IF
    END DO
    !
    trans = .TRUE.
    CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Duals,trans)
    Dulnrm = SUM(ABS(Duals(1:Mrelas)))
  END SUBROUTINE compute_new_duals
  !
END SUBROUTINE DPLPMU
