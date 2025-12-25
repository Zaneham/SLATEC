!** LA05CD
SUBROUTINE LA05CD(A,Ind,Ia,N,Ip,Iw,W,G,U,Mm)
  !> Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (LA05CS-D, LA05CD-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
  !     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
  !     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
  !     THE FINAL LETTER =D= IN THE NAMES USED HERE.
  !     REVISED SEP. 13, 1979.
  !
  !     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
  !     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
  !     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
  !     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
  !     DSPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  LA05ED, XERMSG, XSETUN
  !***
  ! COMMON BLOCKS    LA05DD

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900402  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920410  Corrected second dimension on IW declaration.  (WRB)
  !   920422  Changed upper limit on DO from LAST to LAST-1.  (WRB)
  USE LA05DD, ONLY : lp_com, lcol_com, lenl_com, lenu_com, ncp_com, lrow_com, small_com

  INTEGER, INTENT(IN) :: Ia, Mm, N
  INTEGER, INTENT(INOUT) :: Ind(Ia,2), Iw(N,8), Ip(N,2)
  REAL(DP), INTENT(IN) :: U
  REAL(DP), INTENT(INOUT) :: G, A(:), W(:)
  INTEGER :: i, ii, ij, im, in, ins, ipp, ir, is, j, jm, jns, jp, k, kj, kk, kl, &
    km, knp, kp, kpl, kq, kr, krl, ks, l, last, last1, last2, m, m1, mcp, nz
  REAL(DP) :: am, au
  CHARACTER(8) :: xern1
  INTEGER :: error_code
  LOGICAL :: skip_row_create, found_element, skip_interchange, break_singleton
  LOGICAL :: fill_placed
  !* FIRST EXECUTABLE STATEMENT  LA05CD
  !
  ! Modernisation notes (202412):
  !   Eliminated GOTOs using BLOCK/EXIT and logical flags.
  !   Original F77 used labels 20/40/50/100/110/120/122/124/200/300.
  !   Algorithm unchanged - sparse LU update (Harwell LA05C).
  error_code = 0
  !
  IF( G<0._DP ) THEN
    !
    IF( lp_com>0 ) ERROR STOP 'LA05CD : EARLIER ENTRY GAVE ERROR RETURN.'
    G = -8._DP
    RETURN
  END IF
  !
  main_body: BLOCK
    jm = Mm
    ! MCP LIMITS THE VALUE OF NCP PERMITTED BEFORE AN ERROR RETURN RESULTS.
    mcp = ncp_com + 20
    ! REMOVE OLD COLUMN
    lenu_com = lenu_com - Iw(jm,2)
    kp = Ip(jm,2)
    im = Ind(kp,1)
    kl = kp + Iw(jm,2) - 1
    Iw(jm,2) = 0
    DO k = kp, kl
      i = Ind(k,1)
      Ind(k,1) = 0
      kr = Ip(i,1)
      nz = Iw(i,1) - 1
      Iw(i,1) = nz
      krl = kr + nz
      DO km = kr, krl
        IF( Ind(km,2)==jm ) EXIT
      END DO
      A(km) = A(krl)
      Ind(km,2) = Ind(krl,2)
      Ind(krl,2) = 0
    END DO
    !
    ! INSERT NEW COLUMN
    insert_col: DO ii = 1, N
      i = Iw(ii,3)
      IF( i==im ) m = ii
      IF( ABS(W(i))<=small_com ) THEN
        W(i) = 0._DP
        CYCLE insert_col
      END IF
      lenu_com = lenu_com + 1
      last = ii
      IF( lcol_com+lenl_com>=Ia ) THEN
        ! COMPRESS COLUMN FILE IF NECESSARY.
        IF( ncp_com>=mcp .OR. lenl_com+lenu_com>=Ia ) THEN
          error_code = 300
          EXIT main_body
        END IF
        CALL LA05ED(A,Ind,Ip(:,2),N,Iw(1,2),.FALSE.)
      END IF
      lcol_com = lcol_com + 1
      nz = Iw(jm,2)
      IF( nz==0 ) Ip(jm,2) = lcol_com
      Iw(jm,2) = nz + 1
      Ind(lcol_com,1) = i
      nz = Iw(i,1)
      kpl = Ip(i,1) + nz
      skip_row_create = .FALSE.
      IF( kpl<=lrow_com ) THEN
        IF( Ind(kpl,2)==0 ) skip_row_create = .TRUE.
      END IF
      IF( .NOT. skip_row_create ) THEN
        ! NEW ENTRY HAS TO BE CREATED.
        IF( lenl_com+lrow_com+nz>=Ia ) THEN
          IF( ncp_com>=mcp .OR. lenl_com+lenu_com+nz>=Ia ) THEN
            error_code = 300
            EXIT main_body
          END IF
          ! COMPRESS ROW FILE IF NECESSARY.
          CALL LA05ED(A,Ind(:,2),Ip,N,Iw,.TRUE.)
        END IF
        kp = Ip(i,1)
        Ip(i,1) = lrow_com + 1
        IF( nz/=0 ) THEN
          kpl = kp + nz - 1
          DO k = kp, kpl
            lrow_com = lrow_com + 1
            A(lrow_com) = A(k)
            Ind(lrow_com,2) = Ind(k,2)
            Ind(k,2) = 0
          END DO
        END IF
        lrow_com = lrow_com + 1
        kpl = lrow_com
      END IF
      ! PLACE NEW ELEMENT AT END OF ROW.
      Iw(i,1) = nz + 1
      A(kpl) = W(i)
      Ind(kpl,2) = jm
      W(i) = 0._DP
    END DO insert_col
    IF( Iw(im,1)==0 .OR. Iw(jm,2)==0 .OR. m>last ) THEN
      error_code = 200
      EXIT main_body
    END IF
    !
    ! FIND COLUMN SINGLETONS, OTHER THAN THE SPIKE. NON-SINGLETONS ARE
    !     MARKED WITH W(J)=1. ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED
    !     FOR WORKSPACE.
    ins = m
    m1 = m
    W(jm) = 1._DP
    DO ii = m, last
      i = Iw(ii,3)
      j = Iw(ii,4)
      IF( W(j)==0._DP ) THEN
        ! PLACE SINGLETONS IN NEW POSITION.
        Iw(m1,3) = i
        m1 = m1 + 1
      ELSE
        kp = Ip(i,1)
        kl = kp + Iw(i,1) - 1
        DO k = kp, kl
          j = Ind(k,2)
          W(j) = 1._DP
        END DO
        Iw(ins,4) = i
        ins = ins + 1
      END IF
    END DO
    ! PLACE NON-SINGLETONS IN NEW POSITION.
    ij = m + 1
    DO ii = m1, last - 1
      Iw(ii,3) = Iw(ij,4)
      ij = ij + 1
    END DO
    ! PLACE SPIKE AT END.
    Iw(last,3) = im
    !
    ! FIND ROW SINGLETONS, APART FROM SPIKE ROW. NON-SINGLETONS ARE MARKED
    !     WITH W(I)=2. AGAIN ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED
    !     FOR WORKSPACE.
    last1 = last
    jns = last
    W(im) = 2._DP
    j = jm
    DO ij = m1, last
      ii = last + m1 - ij
      i = Iw(ii,3)
      IF( W(i)/=2._DP ) THEN
        Iw(last1,3) = i
        last1 = last1 - 1
      ELSE
        k = Ip(i,1)
        IF( ii/=last ) j = Ind(k,2)
        kp = Ip(j,2)
        kl = kp + Iw(j,2) - 1
        Iw(jns,4) = i
        jns = jns - 1
        DO k = kp, kl
          i = Ind(k,1)
          W(i) = 2._DP
        END DO
      END IF
    END DO
    DO ii = m1, last1
      jns = jns + 1
      i = Iw(jns,4)
      W(i) = 3._DP
      Iw(ii,3) = i
    END DO
    !
    ! DEAL WITH SINGLETON SPIKE COLUMN. NOTE THAT BUMP ROWS ARE MARKED BY
    !    W(I)=3.
    break_singleton = .FALSE.
    singleton_loop: DO ii = m1, last1
      kp = Ip(jm,2)
      kl = kp + Iw(jm,2) - 1
      is = 0
      DO k = kp, kl
        l = Ind(k,1)
        IF( W(l)==3._DP ) THEN
          IF( is/=0 ) THEN
            break_singleton = .TRUE.
            EXIT singleton_loop
          END IF
          i = l
          knp = k
          is = 1
        END IF
      END DO
      IF( is==0 ) THEN
        error_code = 200
        EXIT main_body
      END IF
      ! MAKE A(I,JM) A PIVOT.
      Ind(knp,1) = Ind(kp,1)
      Ind(kp,1) = i
      kp = Ip(i,1)
      DO k = kp, Ia
        IF( Ind(k,2)==jm ) EXIT
      END DO
      am = A(kp)
      A(kp) = A(k)
      A(k) = am
      Ind(k,2) = Ind(kp,2)
      Ind(kp,2) = jm
      jm = Ind(k,2)
      Iw(ii,4) = i
      W(i) = 2._DP
    END DO singleton_loop
    IF( break_singleton ) THEN
      ! When second singleton found, copy remaining data
      in = m1
      DO ij = ii, last1
        Iw(ij,4) = Iw(in,3)
        in = in + 1
      END DO
    ELSE
      ii = last1
    END IF
    !
    last2 = last1 - 1
  IF( m1/=last1 ) THEN
    DO i = m1, last2
      Iw(i,3) = Iw(i,4)
    END DO
    m1 = ii
    IF( m1/=last1 ) THEN
      !
      ! CLEAR W
      DO i = 1, N
        W(i) = 0._DP
      END DO
      !
      ! PERFORM ELIMINATION
      ir = Iw(last1,3)
      elim_loop: DO ii = m1, last1
        ipp = Iw(ii,3)
        kp = Ip(ipp,1)
        kr = Ip(ir,1)
        jp = Ind(kp,2)
        IF( ii==last1 ) jp = jm
        ! SEARCH NON-PIVOT ROW FOR ELEMENT TO BE ELIMINATED.
        !  AND BRING IT TO FRONT OF ITS ROW
        krl = kr + Iw(ir,1) - 1
        found_element = .FALSE.
        DO knp = kr, krl
          IF( jp==Ind(knp,2) ) THEN
            found_element = .TRUE.
            EXIT
          END IF
        END DO
        IF( .NOT. found_element ) THEN
          IF( ii==last1 ) THEN
            error_code = 200
            EXIT main_body
          END IF
          CYCLE elim_loop
        END IF
        ! BRING ELEMENT TO BE ELIMINATED TO FRONT OF ITS ROW.
        am = A(knp)
        A(knp) = A(kr)
        A(kr) = am
        Ind(knp,2) = Ind(kr,2)
        Ind(kr,2) = jp
        skip_interchange = .FALSE.
        IF( ii/=last1 ) THEN
          IF( ABS(A(kp))>=U*ABS(am) ) THEN
            IF( ABS(am)<U*ABS(A(kp)) ) THEN
              skip_interchange = .TRUE.
            ELSEIF( Iw(ipp,1)<=Iw(ir,1) ) THEN
              skip_interchange = .TRUE.
            END IF
          END IF
        END IF
        IF( .NOT. skip_interchange ) THEN
          ! PERFORM INTERCHANGE
          Iw(last1,3) = ipp
          Iw(ii,3) = ir
          ir = ipp
          ipp = Iw(ii,3)
          k = kr
          kr = kp
          kp = k
          kj = Ip(jp,2)
          DO k = kj, Ia
            IF( Ind(k,1)==ipp ) EXIT
          END DO
          Ind(k,1) = Ind(kj,1)
          Ind(kj,1) = ipp
        END IF
        IF( A(kp)==0._DP ) THEN
          error_code = 200
          EXIT main_body
        END IF
        IF( ii/=last1 ) THEN
          am = -A(kr)/A(kp)
          ! COMPRESS ROW FILE UNLESS IT IS CERTAIN THAT THERE IS ROOM FOR NEW ROW.
          IF( lrow_com+Iw(ir,1)+Iw(ipp,1)+lenl_com>Ia ) THEN
            IF( ncp_com>=mcp .OR. lenu_com+Iw(ir,1)+Iw(ipp,1)+lenl_com>Ia ) THEN
              error_code = 300
              EXIT main_body
            END IF
            CALL LA05ED(A,Ind(1,2),Ip(:,1),N,Iw,.TRUE.)
            kp = Ip(ipp,1)
            kr = Ip(ir,1)
          END IF
          krl = kr + Iw(ir,1) - 1
          kq = kp + 1
          kpl = kp + Iw(ipp,1) - 1
          ! PLACE PIVOT ROW (EXCLUDING PIVOT ITSELF) IN W.
          IF( kq<=kpl ) THEN
            DO k = kq, kpl
              j = Ind(k,2)
              W(j) = A(k)
            END DO
          END IF
          Ip(ir,1) = lrow_com + 1
          !
          ! TRANSFER MODIFIED ELEMENTS.
          Ind(kr,2) = 0
          kr = kr + 1
          IF( kr<=krl ) THEN
            DO ks = kr, krl
              j = Ind(ks,2)
              au = A(ks) + am*W(j)
              Ind(ks,2) = 0
              ! IF ELEMENT IS VERY SMALL REMOVE IT FROM U.
              IF( ABS(au)<=small_com ) THEN
                lenu_com = lenu_com - 1
                ! REMOVE ELEMENT FROM COL FILE.
                k = Ip(j,2)
                kl = k + Iw(j,2) - 1
                Iw(j,2) = kl - k
                DO kk = k, kl
                  IF( Ind(kk,1)==ir ) EXIT
                END DO
                Ind(kk,1) = Ind(kl,1)
                Ind(kl,1) = 0
              ELSE
                G = MAX(G,ABS(au))
                lrow_com = lrow_com + 1
                A(lrow_com) = au
                Ind(lrow_com,2) = j
              END IF
              W(j) = 0._DP
            END DO
          END IF
          !
          ! SCAN PIVOT ROW FOR FILLS.
          IF( kq<=kpl ) THEN
            fill_scan: DO ks = kq, kpl
              j = Ind(ks,2)
              au = am*W(j)
              IF( ABS(au)<=small_com ) THEN
                W(j) = 0._DP
                CYCLE fill_scan
              END IF
              lrow_com = lrow_com + 1
              A(lrow_com) = au
              Ind(lrow_com,2) = j
              lenu_com = lenu_com + 1
              !
              ! CREATE FILL IN COLUMN FILE.
              nz = Iw(j,2)
              k = Ip(j,2)
              kl = k + nz - 1
              fill_placed = .FALSE.
              ! IF POSSIBLE PLACE NEW ELEMENT AT END OF PRESENT ENTRY.
              IF( kl/=lcol_com ) THEN
                IF( Ind(kl+1,1)==0 ) THEN
                  Ind(kl+1,1) = ir
                  fill_placed = .TRUE.
                END IF
              ELSEIF( lcol_com+lenl_com<Ia ) THEN
                lcol_com = lcol_com + 1
                Ind(kl+1,1) = ir
                fill_placed = .TRUE.
              END IF
              IF( .NOT. fill_placed ) THEN
                ! NEW ENTRY HAS TO BE CREATED.
                IF( lcol_com+lenl_com+nz+1>=Ia ) THEN
                  ! COMPRESS COLUMN FILE IF THERE IS NOT ROOM FOR NEW ENTRY.
                  IF( ncp_com>=mcp .OR. lenu_com+lenl_com+nz+1>=Ia ) THEN
                    error_code = 300
                    EXIT main_body
                  END IF
                  CALL LA05ED(A,Ind,Ip(:,2),N,Iw(1,2),.FALSE.)
                  k = Ip(j,2)
                  kl = k + nz - 1
                END IF
                ! TRANSFER OLD ENTRY INTO NEW.
                Ip(j,2) = lcol_com + 1
                DO kk = k, kl
                  lcol_com = lcol_com + 1
                  Ind(lcol_com,1) = Ind(kk,1)
                  Ind(kk,1) = 0
                END DO
                ! ADD NEW ELEMENT.
                lcol_com = lcol_com + 1
                Ind(lcol_com,1) = ir
              END IF
              G = MAX(G,ABS(au))
              Iw(j,2) = nz + 1
              W(j) = 0._DP
            END DO fill_scan
          END IF
          Iw(ir,1) = lrow_com + 1 - Ip(ir,1)
          !
          ! STORE MULTIPLIER
          IF( lenl_com+lcol_com+1>Ia ) THEN
            ! COMPRESS COL FILE IF NECESSARY.
            IF( ncp_com>=mcp ) THEN
              error_code = 300
              EXIT main_body
            END IF
            CALL LA05ED(A,Ind,Ip(:,2),N,Iw(1,2),.FALSE.)
          END IF
          k = Ia - lenl_com
          lenl_com = lenl_com + 1
          A(k) = am
          Ind(k,1) = ipp
          Ind(k,2) = ir
          ! CREATE BLANK IN PIVOTAL COLUMN.
          kp = Ip(jp,2)
          nz = Iw(jp,2) - 1
          kl = kp + nz
          DO k = kp, kl
            IF( Ind(k,1)==ir ) EXIT
          END DO
          Ind(k,1) = Ind(kl,1)
          Iw(jp,2) = nz
          Ind(kl,1) = 0
          lenu_com = lenu_com - 1
        END IF
      END DO elim_loop
    END IF
  END IF
  !
  ! CONSTRUCT COLUMN PERMUTATION AND STORE IT IN IW(.,4)
  DO ii = m, last
    i = Iw(ii,3)
    k = Ip(i,1)
    j = Ind(k,2)
    Iw(ii,4) = j
  END DO
  !
  END BLOCK main_body
  !
  ! ERROR HANDLING
  !
  IF( error_code/=0 ) THEN
    SELECT CASE( error_code )
    CASE( 200 )
      IF( lp_com>0 ) THEN
        WRITE (xern1,'(I8)') Mm
        ERROR STOP 'LA05CD : SINGULAR MATRIX AFTER REPLACEMENT OF COLUMN'
      END IF
      G = -6._DP
    CASE( 300 )
      IF( lp_com>0 ) ERROR STOP 'LA05CD : LENGTHS OF ARRAYS A(*) AND IND(*,2) ARE TOO SMALL.'
      G = -7._DP
    END SELECT
  END IF
  !
END SUBROUTINE LA05CD
