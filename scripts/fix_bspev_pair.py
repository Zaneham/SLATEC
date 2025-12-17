#!/usr/bin/env python
"""Fix GOTO 20/100 in bspev.f90 and dbspev.f90 - B-spline evaluation.

Pattern: Search loop with early exit (GOTO 20) and error exit (GOTO 100)
Replace with: Named loop EXIT + inline error

Reference: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
Original: de Boor (NIST)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/interpolation/bspline/bspev.f90', 'SP', 'BSPEV'),
    ('C:/dev/slatec-modern/src/modern/interpolation/bspline/dbspev.f90', 'DP', 'DBSPEV'),
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
    old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
    new_rev = """!   920501  Reformatted the REFERENCES section.  (WRB)
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 20/100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
  !           Original: de Boor (NIST)"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 100 with inline error
    content = content.replace(
        f"IF( X>T(i) ) GOTO 100\n        DO WHILE",
        f"IF( X>T(i) ) ERROR STOP '{name} : X IS NOT IN T(K)<=X<=T(N+1)'  ! (Was GOTO 100)\n        search: DO WHILE"
    )

    # Replace GOTO 20 with EXIT search
    content = content.replace(
        "IF( X/=T(i) ) GOTO 20\n        END DO\n        ERROR STOP",
        "IF( X/=T(i) ) EXIT search  ! Found valid i (was GOTO 20)\n        END DO search\n        ERROR STOP"
    )

    # Remove label 20
    content = content.replace(
        "END IF\n      !\n      !- I* HAS BEEN FOUND",
        "END IF\n      !\n      !- I* HAS BEEN FOUND (label 20 removed)"
    )
    content = content.replace(
        f"      20  kp1mn",
        f"      kp1mn"
    )

    # Remove label 100 at end
    content = content.replace(
        f"\n  100  ERROR STOP '{name} : X IS NOT IN T(K)<=X<=T(N+1)'\nEND SUBROUTINE",
        f"\n  ! (Label 100 removed - handled inline)\nEND SUBROUTINE"
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 20/100 eliminated (named EXIT + inline error)!')
