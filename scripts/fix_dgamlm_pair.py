#!/usr/bin/env python
"""Fix GOTO 100/200 in dgamlm.f90 and gamlim.f90 - Gamma function limits.

Pattern: Convergence exit from Newton iteration
Replace with: EXIT + loop index check

Reference: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
Original: Fullerton (LANL)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/dgamlm.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/gamlim.f90', 'SP'),
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
  !           Ref: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
  !           Original: Fullerton (LANL)"""

    content = content.replace(old_rev, new_rev)

    # Fix first convergence loop (GOTO 100)
    content = content.replace(
        ") GOTO 100\n  END DO\n  ERROR STOP 'D",
        ") EXIT  ! Converged\n  END DO\n  IF( i > 10 ) ERROR STOP 'D"
    )
    content = content.replace(
        ") GOTO 100\n  END DO\n  ERROR STOP 'G",
        ") EXIT  ! Converged\n  END DO\n  IF( i > 10 ) ERROR STOP 'G"
    )

    # Fix label 100
    content = content.replace(
        "FIND XMIN'\n  !\n  100  Xmin = -Xmin",
        "FIND XMIN'\n  Xmin = -Xmin  ! (Label 100 removed)"
    )

    # Fix second convergence loop (GOTO 200)
    content = content.replace(
        ") GOTO 200\n  END DO\n  ERROR STOP 'D",
        ") EXIT  ! Converged\n  END DO\n  IF( i > 10 ) ERROR STOP 'D"
    )
    content = content.replace(
        ") GOTO 200\n  END DO\n  ERROR STOP 'G",
        ") EXIT  ! Converged\n  END DO\n  IF( i > 10 ) ERROR STOP 'G"
    )

    # Fix label 200
    content = content.replace(
        "FIND XMAX'\n  !\n  200  Xmax = Xmax",
        "FIND XMAX'\n  Xmax = Xmax  ! (Label 200 removed)"
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100/200 eliminated (EXIT + check)!')
