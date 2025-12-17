#!/usr/bin/env python
"""Fix GOTO 100/200 in ckscl.f90 and zkscl.f90 - Bessel K scaling.

Pattern: Alternative value computation (GOTO 100/200)
Replace with: Early EXIT with flag

Reference: ISO/IEC 1539-1:2018 S11.1.1 (IF construct)
Original: Amos (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/ckscl.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/zkscl.f90', 'DP'),
]

for filepath, prec in files:
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
    old_rev = '!   910415  Prologue converted to Version 4.0 format.  (BAB)'
    new_rev = """!   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.1 (IF construct)
  !           Original: Amos (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Add flag declaration (find INTEGER :: line and add LOGICAL)
    content = content.replace(
        'INTEGER :: i, ic, k, kk, nn, nw\n  REAL',
        'INTEGER :: i, ic, k, kk, nn, nw\n  LOGICAL :: early_exit\n  REAL'
    )

    # Initialize flag before loop - find the DO i = 1, nn
    content = content.replace(
        '  DO i = 1, nn\n    Nz = Nz + 1',
        '  early_exit = .FALSE.\n  DO i = 1, nn\n    Nz = Nz + 1'
    )

    # Replace GOTO 100 with flag + EXIT
    content = content.replace(
        'IF( ic==(kk-1) ) GOTO 100',
        'IF( ic==(kk-1) ) THEN\n          early_exit = .TRUE.  ! (Was GOTO 100)\n          EXIT\n        END IF'
    )

    # Replace end-of-loop code with conditional
    content = content.replace(
        '  END DO\n  Nz = N\n  IF( ic==N ) Nz = N - 1\n  GOTO 200\n  100  Nz = kk - 2\n  200 CONTINUE\n  DO k = 1, Nz',
        '  END DO\n  ! Compute Nz (labels 100, 200 removed)\n  IF( early_exit ) THEN\n    Nz = kk - 2\n  ELSE\n    Nz = N\n    IF( ic==N ) Nz = N - 1\n  END IF\n  DO k = 1, Nz'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100/200 eliminated (flag + EXIT)!')
