#!/usr/bin/env python
"""Fix GOTO 100/200 in dbsgq8.f90 and bsgq8.f90 - B-spline weighted quadrature.

Pattern: Backward GOTO 100 iteration + forward GOTO 200 early exit from backtracking
Replace with: Named outer DO loop with EXIT/implicit CYCLE

Reference: ISO/IEC 1539-1:2018 S11.1.7.4.4 (DO construct), S11.1.12 (EXIT)
Original: Jones (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/interpolation/bspline/dbsgq8.f90', 'DP', 'DBSGQ8'),
    ('C:/dev/slatec-modern/src/modern/interpolation/bspline/bsgq8.f90', 'SP', 'BSGQ8'),
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
    old_rev = '!   910408  Updated the AUTHOR section.  (WRB)'
    new_rev = old_rev + """
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.7.4.4, S11.1.12 (DO, EXIT)
  !           Original: Jones (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Replace 100 CONTINUE with start of outer DO loop
    content = content.replace(
        '  100 CONTINUE\n  DO',
        '  main: DO  ! Main adaptive integration loop (was GOTO 100 target)\n    DO'
    )

    # Replace GOTO 200 with EXIT from inner DO WHILE
    content = content.replace(
        'vl(l) = vl(l+1) + vr\n        GOTO 200',
        'vl(l) = vl(l+1) + vr\n        EXIT  ! Found level needing right-half (was GOTO 200)'
    )

    # Add check after backtrack loop to RETURN if done
    content = content.replace(
        'END DO\n    !\n    !      EXIT\n    !\n    Ans = vr',
        'END DO\n    !\n    !      Exit if backtracked to root\n    !\n    IF( l<=1 ) THEN\n      Ans = vr'
    )

    # Indent and close the conditional for the exit section
    content = content.replace(
        'IF( Err<0._' + prec + ' ) Err = ce\n    RETURN\n  END IF',
        'IF( Err<0._' + prec + ' ) Err = ce\n      RETURN\n    END IF\n  END IF'
    )

    # Replace label 200 and GOTO 100 with END DO main
    content = content.replace(
        '  200  est = gr(l-1)\n  lr(l) = 1\n  aa(l) = aa(l) + 4._' + prec + '*hh(l)\n  GOTO 100',
        '    est = gr(l-1)  ! (Label 200 removed)\n    lr(l) = 1\n    aa(l) = aa(l) + 4._' + prec + '*hh(l)\n  END DO main  ! (Was GOTO 100)'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100/200 eliminated (named DO loop with EXIT)!')
