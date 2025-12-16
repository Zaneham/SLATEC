#!/usr/bin/env python3
"""Fix GOTO in dreort.f90 and reort.f90 - BVSUP reorthonormalization"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/bvsup/dreort.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/bvsup/reort.f90', 'SP'),
]

for filepath, prec in files:
    with open(filepath, 'r') as f:
        content = f.read()

    # Add revision history
    old_rev = '!   910722  Updated AUTHOR section.  (ALS)'
    new_rev = '''!   910722  Updated AUTHOR section.  (ALS)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 50 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT with flag)
  !           Original: Watts, SNLA'''

    content = content.replace(old_rev, new_rev)

    # The GOTO 50 pattern: skip to orthonormalization section when S(ijk) > 1.0E20
    # We'll use a logical flag to handle the early orthonormalization trigger
    if prec == 'DP':
        # Add flag declaration
        old_decl = '''INTEGER :: ijk, j, k, kk, l, mflag, nfcp
  REAL(DP) :: dnd, dndt, dx, srp, vnorm, wcnd, ypnm'''
        new_decl = '''INTEGER :: ijk, j, k, kk, l, mflag, nfcp
  REAL(DP) :: dnd, dndt, dx, srp, vnorm, wcnd, ypnm
  LOGICAL :: need_orthonorm  ! Flag to trigger orthonormalization'''
        content = content.replace(old_decl, new_decl)

        # Replace the DO loop with GOTO 50 to use EXIT and flag
        old_loop = '''        DO ijk = 1, nfcp
          !              ......EXIT
          IF( S(ijk)>1.0D20 ) GOTO 50
        END DO'''
        new_loop = '''        need_orthonorm = .FALSE.
        DO ijk = 1, nfcp
          IF( S(ijk)>1.0D20 ) THEN
            need_orthonorm = .TRUE.
            EXIT
          END IF
        END DO
        IF( need_orthonorm ) GOTO 50  ! Temporary - will refactor below'''
        # Actually let's restructure this better - wrap the subsequent code in IF

        # The structure is complex: if any S(ijk) > 1.0D20, skip to orthonormalization
        # Otherwise continue with the extrapolation code
        # Let's use a cleaner approach: flag set if NO orthonorm needed
        old_block = '''        DO ijk = 1, nfcp
          !              ......EXIT
          IF( S(ijk)>1.0D20 ) GOTO 50
        END DO
        !
        !                 USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE
        !                 NORM DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION
        !                 CHECKPOINT.  OTHER CONTROLS ON THE NUMBER OF STEPS TO
        !                 THE NEXT CHECKPOINT ARE ADDED FOR SAFETY PURPOSES.
        !
        nswot_com = knswot_com
        knswot_com = 0
        lotjp_com = 0
        wcnd = LOG10(wcnd)
        IF( wcnd>tnd_com+3._DP ) nswot_com = 2*nswot_com
        IF( wcnd<pwcnd_com ) THEN
          dx = x_com - px_com
          dnd = pwcnd_com - wcnd
          IF( dnd>=4 ) nswot_com = nswot_com/2
          dndt = wcnd - tnd_com
          IF( ABS(dx*dndt)<=dnd*ABS(xend_com-x_com) ) THEN
            xot_com = x_com + dx*dndt/dnd
            nswot_com = MIN(mnswot_com,nswot_com)
            pwcnd_com = wcnd
            !           ......EXIT
            px_com = x_com
          ELSE
            xot_com = xend_com
            nswot_com = MIN(mnswot_com,nswot_com)
            pwcnd_com = wcnd
            px_com = x_com
          END IF
        ELSE
          xot_com = xend_com
          nswot_com = MIN(mnswot_com,nswot_com)
          pwcnd_com = wcnd
          px_com = x_com
        END IF
        RETURN
      END IF
    END IF
    !
    !              *********************************************************
    !
    !              ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE
    !              HOMOGENEOUS SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
    !
    50  nswot_com = 1'''

        new_block = '''        need_orthonorm = .FALSE.
        DO ijk = 1, nfcp
          IF( S(ijk)>1.0D20 ) THEN
            need_orthonorm = .TRUE.
            EXIT
          END IF
        END DO
        !
        IF( .NOT. need_orthonorm ) THEN
          !                 USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE
          !                 NORM DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION
          !                 CHECKPOINT.  OTHER CONTROLS ON THE NUMBER OF STEPS TO
          !                 THE NEXT CHECKPOINT ARE ADDED FOR SAFETY PURPOSES.
          !
          nswot_com = knswot_com
          knswot_com = 0
          lotjp_com = 0
          wcnd = LOG10(wcnd)
          IF( wcnd>tnd_com+3._DP ) nswot_com = 2*nswot_com
          IF( wcnd<pwcnd_com ) THEN
            dx = x_com - px_com
            dnd = pwcnd_com - wcnd
            IF( dnd>=4 ) nswot_com = nswot_com/2
            dndt = wcnd - tnd_com
            IF( ABS(dx*dndt)<=dnd*ABS(xend_com-x_com) ) THEN
              xot_com = x_com + dx*dndt/dnd
              nswot_com = MIN(mnswot_com,nswot_com)
              pwcnd_com = wcnd
              px_com = x_com
            ELSE
              xot_com = xend_com
              nswot_com = MIN(mnswot_com,nswot_com)
              pwcnd_com = wcnd
              px_com = x_com
            END IF
          ELSE
            xot_com = xend_com
            nswot_com = MIN(mnswot_com,nswot_com)
            pwcnd_com = wcnd
            px_com = x_com
          END IF
          RETURN
        END IF
      END IF
    END IF
    !
    !              *********************************************************
    !
    !              ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE
    !              HOMOGENEOUS SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
    !              (Label 50 removed - reached via need_orthonorm flag or fallthrough)
    !
    nswot_com = 1'''
        content = content.replace(old_block, new_block)

    else:  # SP version
        # Add flag declaration
        old_decl = '''INTEGER :: nfcp,ijk, j, k, kk, l, mflag
  REAL(SP) :: dnd, dndt, dx, srp, vnorm, wcnd, ypnm'''
        new_decl = '''INTEGER :: nfcp,ijk, j, k, kk, l, mflag
  REAL(SP) :: dnd, dndt, dx, srp, vnorm, wcnd, ypnm
  LOGICAL :: need_orthonorm  ! Flag to trigger orthonormalization'''
        content = content.replace(old_decl, new_decl)

        # Replace the block in SP version
        old_block = '''        DO ijk = 1, nfcp
          IF( S(ijk)>1.0E+20 ) GOTO 50
        END DO
        !
        !     USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE NORM
        !     DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION CHECKPOINT.
        !     OTHER CONTROLS ON THE NUMBER OF STEPS TO THE NEXT CHECKPOINT
        !     ARE ADDED FOR SAFETY PURPOSES.
        !
        nswot_com = knswot_com
        knswot_com = 0
        lotjp_com = 0
        wcnd = LOG10(wcnd)
        IF( wcnd>tnd_com+3. ) nswot_com = 2*nswot_com
        IF( wcnd>=pwcnd_com ) THEN
          xot_com = xend_com
        ELSE
          dx = x_com - px_com
          dnd = pwcnd_com - wcnd
          IF( dnd>=4 ) nswot_com = nswot_com/2
          dndt = wcnd - tnd_com
          IF( ABS(dx*dndt)>dnd*ABS(xend_com-x_com) ) THEN
            xot_com = xend_com
          ELSE
            xot_com = x_com + dx*dndt/dnd
          END IF
        END IF
        nswot_com = MIN(mnswot_com,nswot_com)
        pwcnd_com = wcnd
        px_com = x_com
        RETURN
      END IF
    END IF
    !
    !     ****************************************
    !
    !     ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE HOMOGENEOUS
    !     SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
    !
    50  nswot_com = 1'''

        new_block = '''        need_orthonorm = .FALSE.
        DO ijk = 1, nfcp
          IF( S(ijk)>1.0E+20 ) THEN
            need_orthonorm = .TRUE.
            EXIT
          END IF
        END DO
        !
        IF( .NOT. need_orthonorm ) THEN
          !     USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE NORM
          !     DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION CHECKPOINT.
          !     OTHER CONTROLS ON THE NUMBER OF STEPS TO THE NEXT CHECKPOINT
          !     ARE ADDED FOR SAFETY PURPOSES.
          !
          nswot_com = knswot_com
          knswot_com = 0
          lotjp_com = 0
          wcnd = LOG10(wcnd)
          IF( wcnd>tnd_com+3. ) nswot_com = 2*nswot_com
          IF( wcnd>=pwcnd_com ) THEN
            xot_com = xend_com
          ELSE
            dx = x_com - px_com
            dnd = pwcnd_com - wcnd
            IF( dnd>=4 ) nswot_com = nswot_com/2
            dndt = wcnd - tnd_com
            IF( ABS(dx*dndt)>dnd*ABS(xend_com-x_com) ) THEN
              xot_com = xend_com
            ELSE
              xot_com = x_com + dx*dndt/dnd
            END IF
          END IF
          nswot_com = MIN(mnswot_com,nswot_com)
          pwcnd_com = wcnd
          px_com = x_com
          RETURN
        END IF
      END IF
    END IF
    !
    !     ****************************************
    !
    !     ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE HOMOGENEOUS
    !     SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
    !     (Label 50 removed - reached via need_orthonorm flag or fallthrough)
    !
    nswot_com = 1'''
        content = content.replace(old_block, new_block)

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 50 eliminated (skip-to-section -> flag)!')
