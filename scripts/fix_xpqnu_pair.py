#!/usr/bin/env python
"""Fix GOTO 100 in dxpqnu.f90 and xpqnu.f90 - Legendre function recurrence.

Pattern: Backward GOTO for infinite loop with RETURN exits
Replace with: DO ... END DO

Reference: ISO/IEC 1539-1:2018 S11.1.4.2 (DO construct)
Original: Smith, NBS (Legendre functions)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/dxpqnu.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/xpqnu.f90', 'SP'),
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
    old_rev = '!   920127  Revised PURPOSE section of prologue.  (DWL)'
    new_rev = """!   920127  Revised PURPOSE section of prologue.  (DWL)
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.4.2 (DO construct)
  !           Original: Smith, NBS (Legendre functions)"""

    content = content.replace(old_rev, new_rev)

    # Replace label 100 with DO
    content = content.replace(
        '  100  pq1 = pq',
        '  ! Forward nu-wise recurrence loop (was GOTO 100)\n  DO\n    pq1 = pq'
    )

    # Replace GOTO 100 with END DO
    content = content.replace(
        '  GOTO 100\n  !\nEND SUBROUTINE',
        '  END DO\n  !\nEND SUBROUTINE'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (DO loop)!')
