!** SPLPMN
SUBROUTINE SPLPMN(USRMAT,Mrelas,Nvars,Costs,Prgopt,Dattrv,Bl,Bu,Ind,Info,&
    Primal,Duals,Amat,Csc,Colnrm,Erd,Erp,Basmat,Wr,Rz,Rg,&
    Rprim,Rhs,Ww,Lmx,Lbm,Ibasis,Ibb,Imat,Ibrc,Ipr,Iwr)
  !> Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SPLPMN-S, DPLPMN-D)
  !***
  ! **Author:**  Hanson, R.J. & Hiebert, K.L. (SNLA)
  !***
  ! **Description:**
  !
  !     MARVEL OPTION(S).. OUTPUT=YES/NO TO ELIMINATE PRINTED OUTPUT.
  !     THIS DOES NOT APPLY TO THE CALLS TO THE ERROR PROCESSOR.
  !
  !     MAIN SUBROUTINE FOR SPLP PACKAGE.
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  IVOUT, LA05BS, PINITM, PNNZRS, PRWPGE, SASUM,
  !                    SCLOSM, SCOPY, SDOT, SPINCW, SPINIT, SPLPCE,
  !                    SPLPDM, SPLPFE, SPLPFL, SPLPMU, SPLPUP, SPOPT,
  !                    SVOUT, XERMSG
  !***
  ! COMMON BLOCKS    LA05DS

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Corrected references to XERRWV.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   251218  Eliminated GOTOs per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.7.4.3, S15.6.2.5
  !           Converted pseudo-coroutines to internal subroutines.
  !           Original algorithm: SAND81-0297 (Hanson & Hiebert, 1981)
  USE service, ONLY : IVOUT, SVOUT
  USE LA05DS, ONLY : lp_com

  INTERFACE
    PURE SUBROUTINE USRMAT(I,J,Aij,Indcat,Dattrv,Iflag)
      IMPORT SP
      INTEGER, INTENT(OUT) :: I, J, Indcat
      INTEGER, INTENT(INOUT) :: Iflag(4)
      REAL(SP), INTENT(IN) :: Dattrv(:)
      REAL(SP), INTENT(OUT) :: Aij
    END SUBROUTINE USRMAT
  END INTERFACE
  INTEGER, INTENT(IN) :: Lbm, Lmx, Mrelas, Nvars
  INTEGER, INTENT(OUT) :: Info
  INTEGER, INTENT(INOUT) :: Ibasis(Nvars+Mrelas), Ibrc(Lbm,2), Ind(Nvars+Mrelas), &
    Ipr(2*Mrelas), Iwr(8*Mrelas)
  INTEGER, INTENT(OUT) :: Ibb(Nvars+Mrelas), Imat(Lmx)
  REAL(SP), INTENT(IN) :: Costs(Nvars), Dattrv(:), Prgopt(:)
  REAL(SP), INTENT(INOUT) :: Basmat(Lbm), Bl(Nvars+Mrelas), Bu(Nvars+Mrelas), &
    Erp(Mrelas), Rprim(Mrelas), Ww(Mrelas)
  REAL(SP), INTENT(OUT) :: Amat(Lmx), Colnrm(Nvars), Csc(Nvars), Duals(Nvars+Mrelas), &
    Erd(Mrelas), Primal(Nvars+Mrelas), Rg(Nvars+Mrelas), Rhs(Mrelas), Rz(Nvars+Mrelas), &
    Wr(Mrelas)

  ! Local variables
  INTEGER :: i, ibas, ienter, ileave, iopt, ipage, iplace, itlp, j, jstrt, k, &
    key, lpg, lpr, lpr1, nerr, np, nparm, nredc, ntries, idum(01)
  INTEGER, TARGET :: intopt(08)
  INTEGER, POINTER :: idg, ipagef, isave, mxitlp, kprint, itbrc, npp, lprg
  REAL(SP) :: aij, anorm, dirnrm, dulnrm, erdnrm, factor, gg, resnrm, rhsnrm, &
    rprnrm, rzj, scalr, scosts, sizee, theta, upbnd, uu, xlamda, xval, rdum(01)
  REAL(SP), TARGET :: ropt(07)
  REAL(SP), POINTER :: eps, asmall, abig, costsc, tolls, tune, tolabs

  LOGICAL, TARGET :: lopt(8)
  LOGICAL, POINTER :: colscp, savedt, contin, cstscp, minprb, sizeup, stpedg, usrbas
  LOGICAL :: unbnd, feas, finite, found, redbas, singlr, trans, zerolv
  CHARACTER(8) :: xern1, xern2

  ! Set up pointer associations for options
  contin => lopt(1)
  usrbas => lopt(2)
  sizeup => lopt(3)
  savedt => lopt(4)
  colscp => lopt(5)
  cstscp => lopt(6)
  minprb => lopt(7)
  stpedg => lopt(8)
  idg => intopt(1)
  ipagef => intopt(2)
  isave => intopt(3)
  mxitlp => intopt(4)
  kprint => intopt(5)
  itbrc => intopt(6)
  npp => intopt(7)
  lprg => intopt(8)
  eps => ropt(1)
  asmall => ropt(2)
  abig => ropt(3)
  costsc => ropt(4)
  tolls => ropt(5)
  tune => ropt(6)
  tolabs => ropt(7)

  !* FIRST EXECUTABLE STATEMENT  SPLPMN
  lp_com = 0
  factor = 0.01_SP
  lpg = Lmx - (Nvars+4)
  iopt = 1
  Info = 0
  unbnd = .FALSE.
  jstrt = 1

  ! Process user options
  CALL SPOPT(Prgopt,Mrelas,Nvars,Info,Csc,Ibasis,ropt,intopt,lopt)
  IF( Info<0 ) THEN
    CALL do_save_and_return()
    RETURN
  END IF

  IF( .NOT. contin ) THEN
    ! Initialize sparse data matrix
    CALL PINITM(Mrelas,Nvars,Amat,Imat,Lmx,ipagef)
  ELSE
    ! Retrieve saved data from file
    CALL retrieve_saved_data()
  END IF

  ! Main algorithm loop - replaces the GOTO state machine
  main_algorithm: DO
    ! Update matrix data and check bounds
    CALL SPLPUP(USRMAT,Mrelas,Nvars,Dattrv,Bl,Bu,Ind,Info,Amat,Imat,&
      sizeup,asmall,abig)
    IF( Info<0 ) EXIT main_algorithm

    ! Print prologue if requested
    IF( kprint>=1 ) CALL print_prologue()

    ! Initialize: scale data, normalize bounds, form column check sums
    CALL SPINIT(Mrelas,Nvars,Costs,Bl,Bu,Ind,Primal,Amat,Csc,costsc,&
      Colnrm,xlamda,anorm,Rhs,rhsnrm,Ibasis,Ibb,Imat,lopt)
    IF( Info<0 ) EXIT main_algorithm

    nredc = 0

    ! Decompose initial basis matrix
    CALL decompose_basis(singlr)
    IF( Info<0 ) EXIT main_algorithm

    IF( singlr ) THEN
      nerr = 23
      Info = -nerr
      ERROR STOP 'SPLPMN : IN SPLP,  A SINGULAR INITIAL BASIS WAS ENCOUNTERED.'
      EXIT main_algorithm
    END IF

    ! Compute error in dual and primal systems
    ntries = 1
    CALL compute_error(ntries, singlr)
    IF( singlr ) THEN
      ! Try redecomposition
      CALL decompose_basis(singlr)
      IF( Info<0 ) EXIT main_algorithm
      ntries = 2
      CALL compute_error(ntries, singlr)
      IF( singlr ) THEN
        nerr = 26
        Info = -nerr
        PRINT*,'SPLPMN : IN SPLP, MOVED TO A SINGULAR POINT.'
        EXIT main_algorithm
      END IF
    END IF

    ! Compute new primal
    CALL compute_primal()

    ! Classify variables
    CALL classify_variables()

    ! Check user-provided basis feasibility
    IF( usrbas ) THEN
      CALL check_feasibility(feas)
      IF( .NOT. feas ) THEN
        nerr = 24
        Info = -nerr
        ERROR STOP 'SPLPMN : IN SPLP, AN INFEASIBLE INITIAL BASIS WAS ENCOUNTERED.'
        EXIT main_algorithm
      END IF
    END IF

    itlp = 0

    ! Penalty method phase
    CALL perform_simplex_phase(xlamda, scosts, feas, unbnd, itlp)
    IF( Info<0 ) EXIT main_algorithm

    ! Rescale and rearrange variables
    CALL rescale_and_rearrange(feas, unbnd)

    ! Set final info value
    IF( feas .AND. (.NOT. unbnd) ) THEN
      Info = 1
    ELSEIF( (.NOT. feas) .AND. (.NOT. unbnd) ) THEN
      nerr = 1
      Info = -nerr
      ERROR STOP 'SPLPMN : IN SPLP, THE PROBLEM APPEARS TO BE INFEASIBLE'
    ELSEIF( feas .AND. unbnd ) THEN
      nerr = 2
      Info = -nerr
      ERROR STOP 'SPLPMN : IN SPLP, THE PROBLEM APPEARS TO HAVE NO FINITE SOLUTION.'
    ELSEIF( (.NOT. feas) .AND. unbnd ) THEN
      nerr = 3
      Info = -nerr
      ERROR STOP 'SPLPMN : IN SPLP, THE PROBLEM APPEARS TO BE INFEASIBLE AND NO FINITE SOLN.'
    END IF

    ! Handle infeasibility markers
    IF( Info==(-1) .OR. Info==(-3) ) CALL mark_infeasible_vars()

    ! Handle unboundedness markers
    IF( Info==(-2) .OR. Info==(-3) ) CALL mark_unbounded_vars()

    ! Print summary if requested
    IF( kprint>=1 ) CALL print_summary()

    EXIT main_algorithm
  END DO main_algorithm

  ! Save data and return
  CALL do_save_and_return()
  RETURN

CONTAINS

  !---------------------------------------------------------------------------
  ! Internal subroutine: Retrieve saved data from file
  !---------------------------------------------------------------------------
  SUBROUTINE retrieve_saved_data()
    lpr = Nvars + 4
    REWIND isave
    READ (isave) (Amat(i),i=1,lpr), (Imat(i),i=1,lpr)
    key = 2
    ipage = 1

    retrieve_loop: DO
      lpr1 = lpr + 1
      READ (isave) (Amat(i),i=lpr1,Lmx), (Imat(i),i=lpr1,Lmx)
      np = Imat(Lmx-1)
      ipage = ipage + 1
      IF( np < 0 ) EXIT retrieve_loop
    END DO retrieve_loop

    nparm = Nvars + Mrelas
    READ (isave) (Ibasis(i),i=1,nparm)
    REWIND isave
  END SUBROUTINE retrieve_saved_data

  !---------------------------------------------------------------------------
  ! Internal subroutine: Save data to file
  !---------------------------------------------------------------------------
  SUBROUTINE save_data_to_file()
    lpr1 = lpr + 1

    save_loop: DO
      WRITE (isave) (Amat(i),i=lpr1,Lmx), (Imat(i),i=lpr1,Lmx)
      np = Imat(Lmx-1)
      ipage = ipage + 1
      IF( np < 0 ) EXIT save_loop
    END DO save_loop

    nparm = Nvars + Mrelas
    WRITE (isave) (Ibasis(i),i=1,nparm)
    ENDFILE isave
  END SUBROUTINE save_data_to_file

  !---------------------------------------------------------------------------
  ! Internal subroutine: Decompose basis matrix
  !---------------------------------------------------------------------------
  SUBROUTINE decompose_basis(is_singular)
    LOGICAL, INTENT(OUT) :: is_singular

    IF( kprint>=2 ) CALL IVOUT(Mrelas,Ibasis,&
      '('' SUBSCRIPTS OF BASIC VARIABLES DURING REDECOMPOSITION'')',idg)

    uu = 0.1_SP
    CALL SPLPDM(Mrelas,Nvars,Lbm,nredc,Info,Ibasis,Imat,Ibrc,Ipr,Iwr,&
      Ind,anorm,eps,uu,gg,Amat,Basmat,Csc,Wr,is_singular,redbas)
  END SUBROUTINE decompose_basis

  !---------------------------------------------------------------------------
  ! Internal subroutine: Compute error in dual and primal systems
  !---------------------------------------------------------------------------
  SUBROUTINE compute_error(num_tries, is_singular)
    INTEGER, INTENT(INOUT) :: num_tries
    LOGICAL, INTENT(OUT) :: is_singular

    IF( num_tries <= 2 ) THEN
      CALL SPLPCE(Mrelas,Nvars,Lmx,Lbm,itlp,itbrc,Ibasis,Imat,Ibrc,Ipr,Iwr,&
        Ind,Ibb,erdnrm,eps,tune,gg,Amat,Basmat,Csc,Wr,Ww,Primal,Erd,Erp,is_singular,redbas)

      IF( .NOT. is_singular ) THEN
        IF( kprint>=3 ) THEN
          CALL SVOUT(Mrelas,Erp,'('' EST. ERROR IN PRIMAL COMPS.'')',idg)
          CALL SVOUT(Mrelas,Erd,'('' EST. ERROR IN DUAL COMPS.'')',idg)
        END IF
      END IF
    ELSE
      is_singular = .TRUE.
    END IF
  END SUBROUTINE compute_error

  !---------------------------------------------------------------------------
  ! Internal subroutine: Classify variables
  !---------------------------------------------------------------------------
  SUBROUTINE classify_variables()
    Primal(Nvars+1:Nvars+Mrelas) = 0._SP

    DO i = 1, Mrelas
      j = Ibasis(i)
      IF( Ind(j)/=4 ) THEN
        IF( Rprim(i)<0._SP ) THEN
          Primal(i+Nvars) = -1._SP
        ELSEIF( Ind(j)==3 ) THEN
          upbnd = Bu(j) - Bl(j)
          IF( j<=Nvars ) upbnd = upbnd/Csc(j)
          IF( Rprim(i)>upbnd ) THEN
            Rprim(i) = Rprim(i) - upbnd
            IF( j>Nvars ) THEN
              Rhs(j-Nvars) = Rhs(j-Nvars) + upbnd
            ELSE
              k = 0
              DO
                CALL PNNZRS(k,aij,iplace,Amat,Imat,j)
                IF( k<=0 ) EXIT
                Rhs(k) = Rhs(k) - upbnd*aij*Csc(j)
              END DO
            END IF
            Primal(i+Nvars) = 1._SP
          END IF
        END IF
      END IF
    END DO
  END SUBROUTINE classify_variables

  !---------------------------------------------------------------------------
  ! Internal subroutine: Compute right hand side
  !---------------------------------------------------------------------------
  SUBROUTINE compute_rhs()
    Rhs(1:Mrelas) = 0._SP

    j = 1
    DO WHILE( j <= Nvars + Mrelas )
      SELECT CASE( Ind(j) )
        CASE (2)
          scalr = -Bu(j)
        CASE (3)
          scalr = -Bl(j)
        CASE (4)
          scalr = 0._SP
        CASE DEFAULT
          scalr = -Bl(j)
      END SELECT

      IF( scalr==0._SP ) THEN
        j = j + 1
      ELSEIF( j>Nvars ) THEN
        Rhs(j-Nvars) = Rhs(j-Nvars) - scalr
        j = j + 1
      ELSE
        i = 0
        DO
          CALL PNNZRS(i,aij,iplace,Amat,Imat,j)
          IF( i>0 ) THEN
            Rhs(i) = Rhs(i) + aij*scalr
          ELSE
            j = j + 1
            EXIT
          END IF
        END DO
      END IF
    END DO

    j = 1
    DO WHILE( j <= Nvars + Mrelas )
      scalr = 0._SP
      IF( Ind(j)==3 .AND. MOD(Ibb(j),2)==0 ) scalr = Bu(j) - Bl(j)
      IF( scalr==0._SP ) THEN
        j = j + 1
      ELSEIF( j>Nvars ) THEN
        Rhs(j-Nvars) = Rhs(j-Nvars) + scalr
        j = j + 1
      ELSE
        i = 0
        DO
          CALL PNNZRS(i,aij,iplace,Amat,Imat,j)
          IF( i>0 ) THEN
            Rhs(i) = Rhs(i) - aij*scalr
          ELSE
            j = j + 1
            EXIT
          END IF
        END DO
      END IF
    END DO
  END SUBROUTINE compute_rhs

  !---------------------------------------------------------------------------
  ! Internal subroutine: Check feasibility
  !---------------------------------------------------------------------------
  SUBROUTINE check_feasibility(is_feasible)
    LOGICAL, INTENT(OUT) :: is_feasible

    Ww(1:Mrelas) = Rhs(1:Mrelas)

    DO j = 1, Mrelas
      ibas = Ibasis(j)
      xval = Rprim(j)

      IF( Ind(ibas)<=3 ) xval = MAX(0._SP,xval)

      IF( Ind(ibas)==3 ) THEN
        upbnd = Bu(ibas) - Bl(ibas)
        IF( ibas<=Nvars ) upbnd = upbnd/Csc(ibas)
        xval = MIN(upbnd,xval)
      END IF

      IF( xval/=0._SP ) THEN
        IF( ibas>Nvars ) THEN
          IF( Ind(ibas)/=2 ) THEN
            Ww(ibas-Nvars) = Ww(ibas-Nvars) + xval
          ELSE
            Ww(ibas-Nvars) = Ww(ibas-Nvars) - xval
          END IF
        ELSE
          i = 0
          DO
            CALL PNNZRS(i,aij,iplace,Amat,Imat,ibas)
            IF( i>0 ) THEN
              Ww(i) = Ww(i) - xval*aij*Csc(ibas)
            ELSE
              EXIT
            END IF
          END DO
        END IF
      END IF
    END DO

    resnrm = SUM(ABS(Ww(1:Mrelas)))
    is_feasible = resnrm<=tolls*(rprnrm*anorm+rhsnrm)

    IF( .NOT. is_feasible ) is_feasible = resnrm<=tolabs
    IF( is_feasible ) Primal(Nvars+1:Nvars+Mrelas) = 0._SP
  END SUBROUTINE check_feasibility

  !---------------------------------------------------------------------------
  ! Internal subroutine: Compute new primal
  !---------------------------------------------------------------------------
  SUBROUTINE compute_primal()
    Ww(1:Mrelas) = Rhs(1:Mrelas)
    trans = .FALSE.
    CALL LA05BS(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,gg,Ww,trans)
    Rprim(1:Mrelas) = Ww(1:Mrelas)
    rprnrm = SUM(ABS(Rprim(1:Mrelas)))
  END SUBROUTINE compute_primal

  !---------------------------------------------------------------------------
  ! Internal subroutine: Compute new duals
  !---------------------------------------------------------------------------
  SUBROUTINE compute_duals()
    DO i = 1, Mrelas
      j = Ibasis(i)
      IF( j>Nvars ) THEN
        Duals(i) = xlamda*Primal(i+Nvars)
      ELSE
        Duals(i) = costsc*Costs(j)*Csc(j) + xlamda*Primal(i+Nvars)
      END IF
    END DO

    trans = .TRUE.
    CALL LA05BS(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,gg,Duals,trans)
    dulnrm = SUM(ABS(Duals(1:Mrelas)))
  END SUBROUTINE compute_duals

  !---------------------------------------------------------------------------
  ! Internal subroutine: Initialize reduced costs
  !---------------------------------------------------------------------------
  SUBROUTINE init_reduced_costs()
    CALL SPINCW(Mrelas,Nvars,Lmx,Lbm,npp,jstrt,Imat,Ibrc,Ipr,Iwr,Ind,&
      Ibb,costsc,gg,erdnrm,dulnrm,Amat,Basmat,Csc,Wr,Ww,Rz,Rg,Costs,Colnrm,Duals,stpedg)
  END SUBROUTINE init_reduced_costs

  !---------------------------------------------------------------------------
  ! Internal subroutine: Find entering variable
  !---------------------------------------------------------------------------
  SUBROUTINE find_entering_variable(var_found)
    LOGICAL, INTENT(OUT) :: var_found

    CALL SPLPFE(Mrelas,Nvars,Lmx,Lbm,ienter,Ibasis,Imat,Ibrc,Ipr,Iwr,Ind,Ibb,&
      erdnrm,eps,gg,dulnrm,dirnrm,Amat,Basmat,Csc,Wr,Ww,Bl,Bu,Rz,Rg,Colnrm,Duals,var_found)
  END SUBROUTINE find_entering_variable

  !---------------------------------------------------------------------------
  ! Internal subroutine: Perform simplex phase (penalty, phase-1, or phase-2)
  !---------------------------------------------------------------------------
  SUBROUTINE perform_simplex_phase(lambda, saved_costsc, is_feasible, is_unbounded, iter_count)
    REAL(SP), INTENT(INOUT) :: lambda, saved_costsc
    LOGICAL, INTENT(INOUT) :: is_feasible, is_unbounded
    INTEGER, INTENT(INOUT) :: iter_count

    LOGICAL :: phase1_needed, phase2_started, any_infeasible
    INTEGER :: phase

    phase1_needed = .FALSE.
    phase2_started = .FALSE.
    phase = 0  ! 0=penalty, 1=phase-1, 2=phase-2

    simplex_phases: DO
      ! Compute duals for current phase
      CALL compute_duals()

      ! Main simplex iteration loop
      simplex_loop: DO
        ! Initialize reduced costs
        CALL init_reduced_costs()

        IF( kprint>2 ) THEN
          CALL SVOUT(Mrelas,Duals,'('' BASIC (INTERNAL) DUAL SOLN.'')',idg)
          CALL SVOUT(Nvars+Mrelas,Rz,'('' REDUCED COSTS'')',idg)
        END IF

        ! Find variable to enter basis
        CALL find_entering_variable(found)

        IF( .NOT. found ) THEN
          ! Check if redecomposition needed
          IF( .NOT. redbas ) THEN
            CALL decompose_basis(singlr)
            IF( Info<0 ) RETURN
            ntries = 1
            CALL compute_error(ntries, singlr)
            IF( singlr ) THEN
              nerr = 26
              Info = -nerr
              PRINT*,'SPLPMN : IN SPLP, MOVED TO A SINGULAR POINT.'
              RETURN
            END IF
            CALL compute_primal()
            CALL compute_duals()
            CALL init_reduced_costs()

            ! Erase non-cycling markers
            DO i = Mrelas + 1, Mrelas + Nvars
              Ibasis(i) = ABS(Ibasis(i))
            END DO

            CALL find_entering_variable(found)
          END IF

          IF( .NOT. found ) EXIT simplex_loop
        END IF

        IF( kprint>=3 ) CALL SVOUT(Mrelas,Ww,'('' SEARCH DIRECTION'')',idg)

        ! Choose variable to leave basis
        CALL SPLPFL(Mrelas,Nvars,ienter,ileave,Ibasis,Ind,theta,dirnrm,&
          rprnrm,Csc,Ww,Bl,Bu,Erp,Rprim,Primal,finite,zerolv)

        IF( .NOT. finite ) THEN
          is_unbounded = .TRUE.
          Ibb(Ibasis(ienter)) = 0
        ELSE
          ! Make move and update
          CALL SPLPMU(Mrelas,Nvars,Lmx,Lbm,nredc,Info,ienter,ileave,npp,&
            jstrt,Ibasis,Imat,Ibrc,Ipr,Iwr,Ind,Ibb,anorm,eps,uu,gg,&
            rprnrm,erdnrm,dulnrm,theta,costsc,xlamda,rhsnrm,Amat,&
            Basmat,Csc,Wr,Rprim,Ww,Bu,Bl,Rhs,Erd,Erp,Rz,Rg,Colnrm,&
            Costs,Primal,Duals,singlr,redbas,zerolv,stpedg)
          IF( Info==(-26) ) RETURN

          IF( kprint>=2 ) CALL print_iteration_summary()

          ! Compute error periodically
          ntries = 1
          CALL compute_error(ntries, singlr)
          IF( singlr ) THEN
            CALL decompose_basis(singlr)
            IF( Info<0 ) RETURN
            ntries = 2
            CALL compute_error(ntries, singlr)
            IF( singlr ) THEN
              nerr = 26
              Info = -nerr
              PRINT*,'SPLPMN : IN SPLP, MOVED TO A SINGULAR POINT.'
              RETURN
            END IF
          END IF
        END IF

        iter_count = iter_count + 1

        ! Check max iterations
        IF( iter_count > mxitlp ) THEN
          nerr = 25
          IF( kprint>=1 ) CALL print_summary()
          idum(1) = 0
          IF( savedt ) idum(1) = isave
          WRITE (xern1,'(I8)') mxitlp
          WRITE (xern2,'(I8)') idum(1)
          PRINT*,'SPLPMN','IN SPLP, MAX ITERATIONS = '//xern1//&
            ' TAKEN.  UP-TO-DATE RESULTS SAVED ON FILE NO. '//xern2
          Info = -nerr
          RETURN
        END IF
      END DO simplex_loop

      ! After simplex loop exits (optimality reached for current phase)
      CALL compute_rhs()
      CALL compute_primal()
      CALL check_feasibility(is_feasible)

      SELECT CASE (phase)
        CASE (0)  ! Penalty phase
          IF( is_feasible ) THEN
            ! Check if any basic vars still infeasible
            any_infeasible = .FALSE.
            DO i = 1, Mrelas
              IF( Primal(i+Nvars)/=0._SP ) THEN
                any_infeasible = .TRUE.
                EXIT
              END IF
            END DO

            IF( any_infeasible ) THEN
              lambda = 0._SP
              phase = 2
            ELSE
              EXIT simplex_phases  ! Done!
            END IF
          ELSE
            ! Need phase-1
            IF( kprint>=2 ) CALL IVOUT(0,idum,'('' ENTER STANDARD PHASE-1'')',idg)
            saved_costsc = costsc
            costsc = 0._SP
            CALL classify_variables()
            phase = 1
          END IF

        CASE (1)  ! Phase-1
          CALL compute_rhs()
          CALL compute_primal()
          CALL check_feasibility(is_feasible)

          IF( is_feasible ) THEN
            IF( kprint>1 ) CALL IVOUT(0,idum,'('' ENTER STANDARD PHASE-2'')',idg)
            lambda = 0._SP
            costsc = saved_costsc
            phase = 2
          ELSE
            EXIT simplex_phases  ! Infeasible
          END IF

        CASE (2)  ! Phase-2
          EXIT simplex_phases  ! Done
      END SELECT

      ! Update duals for next phase
      CALL compute_duals()
    END DO simplex_phases
  END SUBROUTINE perform_simplex_phase

  !---------------------------------------------------------------------------
  ! Internal subroutine: Rescale and rearrange variables
  !---------------------------------------------------------------------------
  SUBROUTINE rescale_and_rearrange(is_feasible, is_unbounded)
    LOGICAL, INTENT(IN) :: is_feasible, is_unbounded

    CALL compute_duals()

    IF( costsc/=0._SP ) THEN
      DO i = 1, Mrelas
        Duals(i) = Duals(i)/costsc
      END DO
    END IF

    CALL compute_primal()

    ! Reapply column scaling to primal
    DO i = 1, Mrelas
      j = Ibasis(i)
      IF( j<=Nvars ) THEN
        scalr = Csc(j)
        IF( Ind(j)==2 ) scalr = -scalr
        Rprim(i) = Rprim(i)*scalr
      END IF
    END DO

    ! Replace translated basic variables into array Primal
    Primal(1:Nvars+Mrelas) = 0._SP
    DO j = 1, Nvars + Mrelas
      ibas = ABS(Ibasis(j))
      xval = 0._SP
      IF( j<=Mrelas ) xval = Rprim(j)
      IF( Ind(ibas)==1 ) xval = xval + Bl(ibas)
      IF( Ind(ibas)==2 ) xval = Bu(ibas) - xval
      IF( Ind(ibas)==3 ) THEN
        IF( MOD(Ibb(ibas),2)==0 ) xval = Bu(ibas) - Bl(ibas) - xval
        xval = xval + Bl(ibas)
      END IF
      Primal(ibas) = xval
    END DO

    ! Compute duals for independent variables with bounds
    DO j = 1, Nvars
      rzj = 0._SP
      IF( Ibb(j)>0._SP .AND. Ind(j)/=4 ) THEN
        rzj = Costs(j)
        i = 0
        DO
          CALL PNNZRS(i,aij,iplace,Amat,Imat,j)
          IF( i<=0 ) EXIT
          rzj = rzj - aij*Duals(i)
        END DO
      END IF
      Duals(Mrelas+j) = rzj
    END DO
  END SUBROUTINE rescale_and_rearrange

  !---------------------------------------------------------------------------
  ! Internal subroutine: Mark infeasible variables
  !---------------------------------------------------------------------------
  SUBROUTINE mark_infeasible_vars()
    sizee = SUM(ABS(Primal(1:Nvars)))*anorm
    sizee = sizee/SUM(ABS(Csc(1:Nvars)))
    sizee = sizee + SUM(ABS(Primal(Nvars+1:Nvars+Mrelas)))

    DO i = 1, Nvars + Mrelas
      SELECT CASE( Ind(i) )
        CASE (1)
          IF( sizee+ABS(Primal(i)-Bl(i))*factor/=sizee ) THEN
            IF( Primal(i)<=Bl(i) ) Ind(i) = -4
          END IF
        CASE (2)
          IF( sizee+ABS(Primal(i)-Bu(i))*factor/=sizee ) THEN
            IF( Primal(i)>=Bu(i) ) Ind(i) = -4
          END IF
        CASE (3)
          IF( sizee+ABS(Primal(i)-Bl(i))*factor/=sizee ) THEN
            IF( Primal(i)<Bl(i) ) THEN
              Ind(i) = -4
            ELSEIF( sizee+ABS(Primal(i)-Bu(i))*factor/=sizee ) THEN
              IF( Primal(i)>Bu(i) ) Ind(i) = -4
            END IF
          END IF
        CASE (4)
          ! Free variable, no bounds to check
      END SELECT
    END DO
  END SUBROUTINE mark_infeasible_vars

  !---------------------------------------------------------------------------
  ! Internal subroutine: Mark unbounded variables
  !---------------------------------------------------------------------------
  SUBROUTINE mark_unbounded_vars()
    DO j = 1, Nvars
      IF( Ibb(j)==0 ) THEN
        SELECT CASE( Ind(j) )
          CASE (1)
            Bu(j) = Bl(j)
            Ind(j) = -3
          CASE (2)
            Bl(j) = Bu(j)
            Ind(j) = -3
          CASE (3)
            ! Has both bounds, keep as is
          CASE (4)
            Bl(j) = 0._SP
            Bu(j) = 0._SP
            Ind(j) = -3
        END SELECT
      END IF
    END DO
  END SUBROUTINE mark_unbounded_vars

  !---------------------------------------------------------------------------
  ! Internal subroutine: Print prologue
  !---------------------------------------------------------------------------
  SUBROUTINE print_prologue()
    idum(1) = Mrelas
    CALL IVOUT(1,idum,'(''1NUM. OF DEPENDENT VARS., MRELAS'')',idg)
    idum(1) = Nvars
    CALL IVOUT(1,idum,'('' NUM. OF INDEPENDENT VARS., NVARS'')',idg)
    CALL IVOUT(1,idum,'('' DIMENSION OF COSTS(*)='')',idg)
    idum(1) = Nvars + Mrelas
    CALL IVOUT(1,idum,&
      '('' DIMENSIONS OF BL(*),BU(*),IND(*)''/'' PRIMAL(*),DUALS(*) ='')',idg)
    CALL IVOUT(1,idum,'('' DIMENSION OF IBASIS(*)='')',idg)
    idum(1) = lprg + 1
    CALL IVOUT(1,idum,'('' DIMENSION OF PRGOPT(*)='')',idg)
    CALL IVOUT(0,idum,'('' 1-NVARS=INDEPENDENT VARIABLE INDICES.''/&
      &'' (NVARS+1)-(NVARS+MRELAS)=DEPENDENT VARIABLE INDICES.''/&
      &'' CONSTRAINT INDICATORS ARE 1-4 AND MEAN'')',idg)
    CALL IVOUT(0,idum,'('' 1=VARIABLE HAS ONLY LOWER BOUND.''/&
      &'' 2=VARIABLE HAS ONLY UPPER BOUND.''/&
      &'' 3=VARIABLE HAS BOTH BOUNDS.''/&
      &'' 4=VARIABLE HAS NO BOUNDS, IT IS FREE.'')',idg)
    CALL SVOUT(Nvars,Costs,'('' ARRAY OF COSTS'')',idg)
    CALL IVOUT(Nvars+Mrelas,Ind,'('' CONSTRAINT INDICATORS'')',idg)
    CALL SVOUT(Nvars+Mrelas,Bl,&
      '('' LOWER BOUNDS FOR VARIABLES  (IGNORE UNUSED ENTRIES.)'')',idg)
    CALL SVOUT(Nvars+Mrelas,Bu,&
      '('' UPPER BOUNDS FOR VARIABLES  (IGNORE UNUSED ENTRIES.)'')',idg)
    IF( kprint>=2 ) THEN
      CALL IVOUT(0,idum,&
        '(''0NON-BASIC INDICES THAT ARE NEGATIVE SHOW VARIABLES EXCHANGED AT A ZERO''/&
        &'' STEP LENGTH'')',idg)
      CALL IVOUT(0,idum,'('' WHEN COL. NO. LEAVING=COL. NO. ENTERING, THE ENTERING &
        &VARIABLE MOVED''/'' TO ITS BOUND.  IT REMAINS NON-BASIC.''/&
        &'' WHEN COL. NO. OF BASIS EXCHANGED IS NEGATIVE, THE LEAVING''/&
        &'' VARIABLE IS AT ITS UPPER BOUND.'')',idg)
    END IF
  END SUBROUTINE print_prologue

  !---------------------------------------------------------------------------
  ! Internal subroutine: Print iteration summary
  !---------------------------------------------------------------------------
  SUBROUTINE print_iteration_summary()
    idum(1) = itlp + 1
    CALL IVOUT(1,idum,'(''0ITERATION NUMBER'')',idg)
    idum(1) = Ibasis(ABS(ileave))
    CALL IVOUT(1,idum,'('' INDEX OF VARIABLE ENTERING THE BASIS'')',idg)
    idum(1) = ileave
    CALL IVOUT(1,idum,'('' COLUMN OF THE BASIS EXCHANGED'')',idg)
    idum(1) = Ibasis(ienter)
    CALL IVOUT(1,idum,'('' INDEX OF VARIABLE LEAVING THE BASIS'')',idg)
    rdum(1) = theta
    CALL SVOUT(1,rdum,'('' LENGTH OF THE EXCHANGE STEP'')',idg)
    IF( kprint>=3 ) THEN
      CALL SVOUT(Mrelas,Rprim,'('' BASIC (INTERNAL) PRIMAL SOLN.'')',idg)
      CALL IVOUT(Nvars+Mrelas,Ibasis,&
        '('' VARIABLE INDICES IN POSITIONS 1-MRELAS ARE BASIC.'')',idg)
      CALL IVOUT(Nvars+Mrelas,Ibb,'('' IBB ARRAY'')',idg)
      CALL SVOUT(Mrelas,Rhs,'('' TRANSLATED RHS'')',idg)
      CALL SVOUT(Mrelas,Duals,'('' BASIC (INTERNAL) DUAL SOLN.'')',idg)
    END IF
  END SUBROUTINE print_iteration_summary

  !---------------------------------------------------------------------------
  ! Internal subroutine: Print summary
  !---------------------------------------------------------------------------
  SUBROUTINE print_summary()
    idum(1) = Info
    CALL IVOUT(1,idum,'('' THE OUTPUT VALUE OF INFO IS'')',idg)
    IF( .NOT. minprb ) THEN
      CALL IVOUT(0,idum,'('' THIS IS A MAXIMIZATION PROBLEM.'')',idg)
    ELSE
      CALL IVOUT(0,idum,'('' THIS IS A MINIMIZATION PROBLEM.'')',idg)
    END IF
    IF( .NOT. stpedg ) THEN
      CALL IVOUT(0,idum,'('' MINIMUM REDUCED COST PRICING WAS USED.'')',idg)
    ELSE
      CALL IVOUT(0,idum,'('' STEEPEST EDGE PRICING WAS USED.'')',idg)
    END IF
    rdum(1) = DOT_PRODUCT(Costs(1:Nvars),Primal(1:Nvars))
    CALL SVOUT(1,rdum,'('' OUTPUT VALUE OF THE OBJECTIVE FUNCTION'')',idg)
    CALL SVOUT(Nvars+Mrelas,Primal,&
      '('' THE OUTPUT INDEPENDENT AND DEPENDENT VARIABLES'')',idg)
    CALL SVOUT(Mrelas+Nvars,Duals,'('' THE OUTPUT DUAL VARIABLES'')',idg)
    CALL IVOUT(Nvars+Mrelas,Ibasis,&
      '('' VARIABLE INDICES IN POSITIONS 1-MRELAS ARE BASIC.'')',idg)
    idum(1) = itlp
    CALL IVOUT(1,idum,'('' NO. OF ITERATIONS'')',idg)
    idum(1) = nredc
    CALL IVOUT(1,idum,'('' NO. OF FULL REDECOMPS'')',idg)
  END SUBROUTINE print_summary

  !---------------------------------------------------------------------------
  ! Internal subroutine: Save data and return
  !---------------------------------------------------------------------------
  SUBROUTINE do_save_and_return()
    IF( savedt ) THEN
      key = INT( Amat(4) )
      Amat(4) = 0._SP
      lpr = Nvars + 4
      WRITE (isave) (Amat(i),i=1,lpr), (Imat(i),i=1,lpr)
      Amat(4) = key
      ipage = 1
      key = 1
      CALL save_data_to_file()
    END IF
  END SUBROUTINE do_save_and_return

END SUBROUTINE SPLPMN
