#!/usr/bin/env python
"""Fix GOTO 100 in dxlegf.f90 and xlegf.f90 - Legendre functions.

Pattern: Error exit via GOTO 100
Replace with: Inline ERROR STOP

Reference: ISO/IEC 1539-1:2018 S11.2.3 (ERROR STOP)
Original: Smith, Amos (Legendre functions)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/dxlegf.f90', 'DP', 'DXLEGF'),
    ('C:/dev/slatec-modern/src/modern/special_functions/xlegf.f90', 'SP', 'XLEGF'),
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
    old_rev = '!   920618  Removed unnecessary variable (WRB).'
    new_rev = """!   920618  Removed unnecessary variable (WRB).
  !   211001  Converted to free-form, added INTENT.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.3 (ERROR STOP)
  !           Original: Smith, Amos (Legendre functions)"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 100 with inline ERROR STOP
    content = content.replace(
        ') GOTO 100',
        f") ERROR STOP '{name} : DNU1, NUDIFF, MU1, MU2, or ID not valid'"
    )

    # Remove the label 100 line
    content = content.replace(
        f"\n  100 ERROR STOP '{name} : DNU1, NUDIFF, MU1, MU2, or ID not valid'",
        '\n  ! (Label 100 removed - error handled inline)'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (inline ERROR STOP)!')
