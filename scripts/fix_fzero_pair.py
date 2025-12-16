#!/usr/bin/env python
"""Fix GOTO 20 in dfzero.f90 and fzero.f90 - Dekker root finding.

Pattern: Skip-ahead GOTO after bisection step
Replace with: Flag to skip increment selection

Reference: ISO/IEC 1539-1:2018 S11.1.8 (IF-THEN-ELSE construct)
Original: Shampine & Watts, SNLA (1970); T.J. Dekker, 1969
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/nonlin_eq/dfzero.f90', 'DP', 'DFZERO'),
    ('C:/dev/slatec-modern/src/modern/nonlin_eq/fzero.f90', 'SP', 'FZERO'),
]

for filepath, prec, name in files:
    if not os.path.exists(filepath):
        print(f'{os.path.basename(filepath)}: NOT FOUND, skipping')
        continue

    with open(filepath, 'r') as f:
        content = f.read()

    # Skip if already fixed
    if 'Eliminated GOTO' in content:
        print(f'{os.path.basename(filepath)}: already fixed, skipping')
        continue

    # Skip if no GOTO 20
    if 'GOTO 20' not in content:
        print(f'{os.path.basename(filepath)}: no GOTO 20, skipping')
        continue

    # Add revision history
    old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
    new_rev = """!   920501  Reformatted the REFERENCES section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 20 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.8 (IF construct)
  !           Original: Shampine & Watts, SNLA; Dekker, 1969"""

    content = content.replace(old_rev, new_rev)

    # Add flag variable declaration
    if prec == 'DP':
        old_decl = '  INTEGER :: ic, kount\n  REAL(DP) :: a, acbs'
        new_decl = '  INTEGER :: ic, kount\n  LOGICAL :: used_bisection\n  REAL(DP) :: a, acbs'
    else:
        old_decl = '  INTEGER :: ic, kount\n  REAL(SP) :: a, acbs'
        new_decl = '  INTEGER :: ic, kount\n  LOGICAL :: used_bisection\n  REAL(SP) :: a, acbs'

    content = content.replace(old_decl, new_decl)

    # Set up the pattern replacement
    # Original: bisection with GOTO 20, then increment selection, then label 20
    if prec == 'DP':
        old_pattern = """        IF( ic>=4 ) THEN
          IF( 8._DP*acmb>=acbs ) THEN
            !
            !   Use bisection (C+B)/2.
            !
            B = B + cmb
            GOTO 20
          ELSE
            ic = 0
            acbs = acmb
          END IF
        END IF
        !
        !   Test for too small a change.
        !
        IF( p<=ABS(q)*tol ) THEN
          !
          !   Increment by TOLerance.
          !
          B = B + SIGN(tol,cmb)
          !
          !   Root ought to be between B and (C+B)/2.
          !
        ELSEIF( p>=cmb*q ) THEN
          B = B + cmb
        ELSE
          !
          !   Use secant rule.
          !
          B = B + p/q
        END IF
      END IF
      !
      !   Have completed computation for new iterate B.
      !
      20  t = B"""

        new_pattern = """        used_bisection = .FALSE.
        IF( ic>=4 ) THEN
          IF( 8._DP*acmb>=acbs ) THEN
            !
            !   Use bisection (C+B)/2.
            !
            B = B + cmb
            used_bisection = .TRUE.
          ELSE
            ic = 0
            acbs = acmb
          END IF
        END IF
        !
        !   Test for too small a change (skip if bisection used).
        !
        IF( .NOT. used_bisection ) THEN
          IF( p<=ABS(q)*tol ) THEN
            !
            !   Increment by TOLerance.
            !
            B = B + SIGN(tol,cmb)
            !
            !   Root ought to be between B and (C+B)/2.
            !
          ELSEIF( p>=cmb*q ) THEN
            B = B + cmb
          ELSE
            !
            !   Use secant rule.
            !
            B = B + p/q
          END IF
        END IF
      END IF
      !
      !   Have completed computation for new iterate B.
      !   (Label 20 removed - bisection handled by flag)
      !
      t = B"""

    else:  # SP version
        old_pattern = """        IF( ic>=4 ) THEN
          IF( 8.0E0*acmb>=acbs ) THEN
            !
            !   Use bisection (C+B)/2.
            !
            B = B + cmb
            GOTO 20
          ELSE
            ic = 0
            acbs = acmb
          END IF
        END IF
        !
        !   Test for too small a change.
        !
        IF( p<=ABS(q)*tol ) THEN
          !
          !   Increment by TOLerance.
          !
          B = B + SIGN(tol,cmb)
          !
          !   Root ought to be between B and (C+B)/2.
          !
        ELSEIF( p>=cmb*q ) THEN
          B = B + cmb
        ELSE
          !
          !   Use secant rule.
          !
          B = B + p/q
        END IF
      END IF
      !
      !   Have completed computation for new iterate B.
      !
      20  t = B"""

        new_pattern = """        used_bisection = .FALSE.
        IF( ic>=4 ) THEN
          IF( 8.0E0*acmb>=acbs ) THEN
            !
            !   Use bisection (C+B)/2.
            !
            B = B + cmb
            used_bisection = .TRUE.
          ELSE
            ic = 0
            acbs = acmb
          END IF
        END IF
        !
        !   Test for too small a change (skip if bisection used).
        !
        IF( .NOT. used_bisection ) THEN
          IF( p<=ABS(q)*tol ) THEN
            !
            !   Increment by TOLerance.
            !
            B = B + SIGN(tol,cmb)
            !
            !   Root ought to be between B and (C+B)/2.
            !
          ELSEIF( p>=cmb*q ) THEN
            B = B + cmb
          ELSE
            !
            !   Use secant rule.
            !
            B = B + p/q
          END IF
        END IF
      END IF
      !
      !   Have completed computation for new iterate B.
      !   (Label 20 removed - bisection handled by flag)
      !
      t = B"""

    content = content.replace(old_pattern, new_pattern)

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 20 eliminated (bisection flag)!')
