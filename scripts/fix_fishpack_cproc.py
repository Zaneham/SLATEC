#!/usr/bin/env python
"""Fix GOTO 100 in FISHPACK cproc/cprod family - tridiagonal system solvers.

Pattern: Backward GOTO 100 for iteration loop
Replace with: DO loop with EXIT

Reference: ISO/IEC 1539-1:2018 S11.1.4.2 (DO construct)
Original: Swarztrauber (NCAR)
"""

import os

files = [
    'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/cproc.f90',
    'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/cprocp.f90',
    'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/cprod.f90',
    'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/cprodp.f90',
]

for filepath in files:
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
    old_rev = '!   920308  Made changes as suggested by Amos.'
    new_rev = """!   920308  Made changes as suggested by Amos.
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.4.2 (DO construct)
  !           Original: Swarztrauber (NCAR)"""

    content = content.replace(old_rev, new_rev)

    # Replace label 100 CONTINUE with DO
    content = content.replace(
        '  100 CONTINUE\n  IFlg = 0',
        '  ! Main iteration loop (was GOTO 100 target)\n  DO\n    IFlg = 0'
    )

    # Replace conditional GOTO with EXIT/CYCLE
    content = content.replace(
        'IF( iflg>0 ) GOTO 100\n      RETURN',
        'IF( iflg<=0 ) EXIT  ! Done (was RETURN)\n      CYCLE  ! Continue iteration'
    )

    # Replace unconditional GOTO 100 with END DO
    content = content.replace(
        '  iflg = 1\n  GOTO 100\n  !\n  RETURN\nEND SUBROUTINE',
        '  iflg = 1\n  END DO  ! Main iteration loop\n  !\nEND SUBROUTINE'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (DO loop)!')
