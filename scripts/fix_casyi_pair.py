#!/usr/bin/env python
"""Fix GOTO 20/100 in casyi.f90 and zasyi.f90 - Bessel I asymptotic expansion.

Pattern: Convergence exit (GOTO 20) and exhaustion error (GOTO 100)
Replace with: EXIT + loop index check

Reference: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
Original: Amos (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/casyi.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/zasyi.f90', 'DP'),
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
  !   251217  Eliminated GOTO 20/100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
  !           Original: Amos (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 20 with EXIT
    content = content.replace(
        'IF( aa<=atol ) GOTO 20\n      END DO\n      GOTO 100\n      20  s2 = cs1',
        'IF( aa<=atol ) EXIT  ! Converged\n      END DO\n      IF( j > jl ) THEN  ! Loop exhausted - no convergence\n        Nz = -2\n        RETURN\n      END IF\n      s2 = cs1  ! (Label 20 removed)'
    )

    # Remove label 100 at end
    content = content.replace(
        '\n  100  Nz = -2\n  !\nEND SUBROUTINE',
        '\n  ! (Label 100 removed - handled inline)\nEND SUBROUTINE'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 20/100 eliminated (EXIT + check)!')
