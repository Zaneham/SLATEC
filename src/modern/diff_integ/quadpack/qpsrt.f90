!** QPSRT
PURE SUBROUTINE QPSRT(Limit,Last,Maxerr,Ermax,Elist,Iord,Nrmax)
  !> Subsidiary to QAGE, QAGIE, QAGPE, QAGSE, QAWCE, QAWOE and QAWSE
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (QPSRT-S, DQPSRT-D)
  !***
  ! **Keywords:**  SEQUENTIAL SORTING
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  ! 1.        QPSRT
  !           Ordering Routine
  !              Standard FORTRAN Subroutine
  !              REAL Version
  !
  ! 2.        PURPOSE
  !              This routine maintains the descending ordering
  !              in the list of the local error estimates resulting from
  !              the interval subdivision process. At each call two error
  !              estimates are inserted using the sequential search
  !              method, top-down for the largest error estimate
  !              and bottom-up for the smallest error estimate.
  !
  ! 3.        CALLING SEQUENCE
  !              CALL QPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
  !
  !           PARAMETERS (MEANING AT OUTPUT)
  !              LIMIT  - INTEGER
  !                       Maximum number of error estimates the list can contain
  !
  !              LAST   - INTEGER
  !                       Number of error estimates currently in the list
  !
  !              MAXERR - INTEGER
  !                       MAXERR points to the NRMAX-th largest error
  !                       estimate currently in the list
  !
  !              ERMAX  - REAL
  !                       NRMAX-th largest error estimate
  !                       ERMAX = ELIST(MAXERR)
  !
  !              ELIST  - REAL
  !                       Vector of dimension LAST containing the error estimates
  !
  !              IORD   - INTEGER
  !                       Vector of dimension LAST, the first K
  !                       elements of which contain pointers
  !                       to the error estimates, such that
  !                       ELIST(IORD(1)),..., ELIST(IORD(K))
  !                       form a decreasing sequence, with
  !                       K = LAST if LAST<=(LIMIT/2+2), and
  !                       K = LIMIT+1-LAST otherwise
  !
  !              NRMAX  - INTEGER
  !                       MAXERR = IORD(NRMAX)
  !
  !***
  ! **See also:**  QAGE, QAGIE, QAGPE, QAGSE, QAWCE, QAWOE, QAWSE
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200/300 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.7.4.4, S11.1.12 (DO, EXIT)
  !           Original: Piessens, de Doncker (K. U. Leuven)

  !
  INTEGER, INTENT(IN) :: Last, Limit
  INTEGER, INTENT(INOUT) :: Nrmax, Iord(Last), Maxerr
  REAL(SP), INTENT(IN) :: Elist(Last)
  REAL(SP), INTENT(OUT) :: Ermax
  !
  INTEGER :: i, ibeg, ido, isucc, j, jbnd, jupbn, k
  REAL(SP) :: errmax, errmin
  !
  !           CHECK WHETHER THE LIST CONTAINS MORE THAN
  !           TWO ERROR ESTIMATES.
  !
  !* FIRST EXECUTABLE STATEMENT  QPSRT
  IF( Last>2 ) THEN
    !
    !           THIS PART OF THE ROUTINE IS ONLY EXECUTED
    !           IF, DUE TO A DIFFICULT INTEGRAND, SUBDIVISION
    !           INCREASED THE ERROR ESTIMATE. IN THE NORMAL CASE
    !           THE INSERT PROCEDURE SHOULD START AFTER THE
    !           NRMAX-TH LARGEST ERROR ESTIMATE.
    !
    errmax = Elist(Maxerr)
    IF( Nrmax/=1 ) THEN
      ido = Nrmax - 1
      DO i = 1, ido
        isucc = Iord(Nrmax-1)
        !- **JUMP OUT OF DO-LOOP
        IF( errmax<=Elist(isucc) ) EXIT
        Iord(Nrmax) = isucc
        Nrmax = Nrmax - 1
      END DO
    END IF
    !
    !           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO
    !           BE MAINTAINED IN DESCENDING ORDER. THIS NUMBER
    !           DEPENDS ON THE NUMBER OF SUBDIVISIONS STILL
    !           ALLOWED.
    !
    jupbn = Last
    IF( Last>(Limit/2+2) ) jupbn = Limit + 3 - Last
    errmin = Elist(Last)
    !
    !           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
    !           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
    !
    jbnd = jupbn - 1
    ibeg = Nrmax + 1
    IF( ibeg<=jbnd ) THEN
      search_errmax: DO i = ibeg, jbnd
        isucc = Iord(i)
        IF( errmax>=Elist(isucc) ) EXIT search_errmax  ! Found insertion point (was GOTO 100)
        Iord(i-1) = isucc
      END DO search_errmax
      IF( i<=jbnd ) THEN  ! Found early (was label 100)
        Iord(i-1) = Maxerr
        k = jbnd
        search_errmin: DO j = i, jbnd
          isucc = Iord(k)
          IF( errmin<Elist(isucc) ) EXIT search_errmin  ! Found (was GOTO 200)
          Iord(k+1) = isucc
          k = k - 1
        END DO search_errmin
        IF( j<=jbnd ) THEN
          Iord(k+1) = Last  ! (Was label 200)
        ELSE
          Iord(i) = Last  ! (Was after second loop)
        END IF
      ELSE  ! First loop completed
        Iord(jbnd) = Maxerr
        Iord(jupbn) = Last
      END IF
    END IF
  ELSE
    Iord(1) = 1
    Iord(2) = 2
  END IF
  !
  !           SET MAXERR AND ERMAX.
  !
  ! Set MAXERR and ERMAX (label 300 removed)
  Maxerr = Iord(Nrmax)
  Ermax = Elist(Maxerr)
  !
END SUBROUTINE QPSRT