#!/usr/bin/env python
"""Fix backward GOTO loops in SPLP pincw files.

Pattern: GOTO 100 at end looping back to 100 CONTINUE
Replace with: DO ... EXIT ... END DO

Reference: ISO/IEC 1539-1:2018 S11.1.4.2 (DO construct)
Original: Hanson & Hiebert (SPLP)
"""

import os

# Only dpincw/spincw have the backward GOTO 100 loop pattern
files = [
    ('C:/dev/slatec-modern/src/modern/optimization/splp/dpincw.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/optimization/splp/spincw.f90', 'SP'),
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
    old_rev = '!   900328  Added TYPE section.  (WRB)'
    new_rev = """!   900328  Added TYPE section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.4.2 (DO construct)
  !           Original: Hanson & Hiebert (SPLP)"""

    content = content.replace(old_rev, new_rev)

    # Replace 100 CONTINUE with DO
    content = content.replace(
        '  100 CONTINUE',
        '  ! Main search loop (was GOTO 100 loop)\n  DO'
    )

    # Replace GOTO 100 with EXIT condition
    content = content.replace(
        'IF( nnegrc<Npp .AND. j/=Jstrt ) GOTO 100',
        'IF( nnegrc>=Npp .OR. j==Jstrt ) EXIT\n  END DO'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (DO construct)!')
