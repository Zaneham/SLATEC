#!/usr/bin/env python
"""Fix GOTO 100/200 in ddzro.f90, sdzro.f90, cdzro.f90 - zero finder for SDRIVE.

Pattern: Backward GOTO 100 iteration + forward GOTO 200 to skip calculation
Replace with: Named DO loop with skip flag

Reference: ISO/IEC 1539-1:2018 S11.1.7.4.4 (DO construct)
Original: Kahaner (NIST), Sutherland (LANL)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/sdrive/ddzro.f90', 'DP', 'DDZRO'),
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/sdrive/sdzro.f90', 'SP', 'SDZRO'),
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/sdrive/cdzro.f90', 'SP', 'CDZRO'),
]

for filepath, prec, name in files:
    fname = os.path.basename(filepath)

    if not os.path.exists(filepath):
        print(f'{fname}: NOT FOUND, skipping')
        continue

    with open(filepath, 'r') as f:
        content = f.read()

    # Skip if already fixed
    if 'Eliminated GOTO' in content:
        print(f'{fname}: already fixed, skipping')
        continue

    # Add revision history
    old_rev = '!   900329  Initial submission to SLATEC.'
    new_rev = old_rev + """
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.7.4.4 (DO construct)
  !           Original: Kahaner (NIST), Sutherland (LANL)"""

    content = content.replace(old_rev, new_rev)

    # Replace 100 CONTINUE with outer DO loop start
    content = content.replace(
        '  !                                                    Perform interchange\n  100 CONTINUE\n  IF',
        '  !                                                    Perform interchange\n  main: DO  ! Main iteration (was GOTO 100 target)\n    IF'
    )

    # Replace GOTO 200 with skip flag - restructure the nested IF
    if prec == 'DP':
        old_block = '''  IF( ic>=4 ) THEN
    IF( 8._DP*acmb>=acbs ) THEN
      !                                                                 Bisect
      B = 0.5_DP*(C+B)
      GOTO 200
    END IF
    ic = 0
  END IF
  acbs = acmb
  !                                            Test for too small a change
  IF( p<=ABS(q)*tol ) THEN
    !                                                 Increment by tolerance
    B = B + SIGN(tol,cmb)
    !                                               Root ought to be between
    !                                               B and (C + B)/2.
  ELSEIF( p<cmb*q ) THEN
    !                                                            Interpolate
    B = B + p/q
  ELSE
    !                                                                 Bisect
    B = 0.5_DP*(C+B)
  END IF
  !                                             Have completed computation
  !                                             for new iterate B.
  200  CALL'''
        new_block = '''  IF( ic>=4 ) THEN
    IF( 8._DP*acmb>=acbs ) THEN
      !                                                                 Bisect
      B = 0.5_DP*(C+B)
      ! Skip to function evaluation (was GOTO 200)
    ELSE
      ic = 0
      acbs = acmb
      !                                            Test for too small a change
      IF( p<=ABS(q)*tol ) THEN
        !                                                 Increment by tolerance
        B = B + SIGN(tol,cmb)
      ELSEIF( p<cmb*q ) THEN
        !                                                            Interpolate
        B = B + p/q
      ELSE
        !                                                                 Bisect
        B = 0.5_DP*(C+B)
      END IF
    END IF
  ELSE
    acbs = acmb
    !                                            Test for too small a change
    IF( p<=ABS(q)*tol ) THEN
      !                                                 Increment by tolerance
      B = B + SIGN(tol,cmb)
    ELSEIF( p<cmb*q ) THEN
      !                                                            Interpolate
      B = B + p/q
    ELSE
      !                                                                 Bisect
      B = 0.5_DP*(C+B)
    END IF
  END IF
  !                                             Have completed computation
  !                                             for new iterate B.
  ! (Label 200 removed)
  CALL'''
    else:  # SP for sdzro and cdzro
        old_block = '''  IF( ic>=4 ) THEN
    IF( 8._SP*acmb>=acbs ) THEN
      !                                                                 Bisect
      B = 0.5_SP*(C+B)
      GOTO 200
    END IF
    ic = 0
  END IF
  acbs = acmb
  !                                            Test for too small a change
  IF( p<=ABS(q)*tol ) THEN
    !                                                 Increment by tolerance
    B = B + SIGN(tol,cmb)
    !                                               Root ought to be between
    !                                               B and (C + B)/2.
  ELSEIF( p<cmb*q ) THEN
    !                                                            Interpolate
    B = B + p/q
  ELSE
    !                                                                 Bisect
    B = 0.5_SP*(C+B)
  END IF
  !                                             Have completed computation
  !                                             for new iterate B.
  200  CALL'''
        new_block = '''  IF( ic>=4 ) THEN
    IF( 8._SP*acmb>=acbs ) THEN
      !                                                                 Bisect
      B = 0.5_SP*(C+B)
      ! Skip to function evaluation (was GOTO 200)
    ELSE
      ic = 0
      acbs = acmb
      !                                            Test for too small a change
      IF( p<=ABS(q)*tol ) THEN
        !                                                 Increment by tolerance
        B = B + SIGN(tol,cmb)
      ELSEIF( p<cmb*q ) THEN
        !                                                            Interpolate
        B = B + p/q
      ELSE
        !                                                                 Bisect
        B = 0.5_SP*(C+B)
      END IF
    END IF
  ELSE
    acbs = acmb
    !                                            Test for too small a change
    IF( p<=ABS(q)*tol ) THEN
      !                                                 Increment by tolerance
      B = B + SIGN(tol,cmb)
    ELSEIF( p<cmb*q ) THEN
      !                                                            Interpolate
      B = B + p/q
    ELSE
      !                                                                 Bisect
      B = 0.5_SP*(C+B)
    END IF
  END IF
  !                                             Have completed computation
  !                                             for new iterate B.
  ! (Label 200 removed)
  CALL'''

    content = content.replace(old_block, new_block)

    # Replace GOTO 100 with END DO main
    content = content.replace(
        '    Fc = fa\n  END IF\n  GOTO 100\n  !\nEND SUBROUTINE',
        '    Fc = fa\n  END IF\n  END DO main  ! (Was GOTO 100)\nEND SUBROUTINE'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100/200 eliminated (DO loop with restructured IF)!')
