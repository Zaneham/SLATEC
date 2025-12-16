#!/usr/bin/env python
"""Fix GOTO 20 in cunhj.f90 and zunhj.f90 - Bessel uniform H asymptotic.

Pattern: Convergence exit from DO loop
Replace with: EXIT + loop index check

Reference: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
Original: Amos (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/cunhj.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/zunhj.f90', 'DP'),
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
    old_rev = '!   930501  Made changes as suggested by Amos.'
    new_rev = """!   930501  Made changes as suggested by Amos.
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 20 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
  !           Original: Amos (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 20 with EXIT
    content = content.replace(
        'IF( ap(k)<Tol ) GOTO 20\n        END DO\n        k = 30\n      END IF\n      20  kmax = k',
        'IF( ap(k)<Tol ) EXIT  ! Converged (was GOTO 20)\n        END DO\n        ! Check if we exited early or exhausted iterations\n        IF( k > 30 ) k = 30\n      END IF\n      kmax = k  ! (Label 20 removed)'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 20 eliminated (EXIT)!')
