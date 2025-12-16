#!/usr/bin/env python
"""Fix GOTO 100 in dxqnu.f90 and xqnu.f90 - Legendre Q function recurrence.

Pattern: GOTO 100 to restart outer loop with reinitialization
Replace with: Named outer loop + CYCLE

Reference: ISO/IEC 1539-1:2018 S11.2.1 (CYCLE statement)
Original: Smith, NBS (Legendre functions)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/dxqnu.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/xqnu.f90', 'SP'),
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
  !   251217  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (CYCLE statement)
  !           Original: Smith, NBS (Legendre functions)"""

    content = content.replace(old_rev, new_rev)

    # Wrap the label 100 line with named outer loop
    if prec == 'DP':
        content = content.replace(
            '  100  mu = 1\n  dmu = 1._DP\n  DO',
            '  ! Forward mu-recurrence loop (was GOTO 100 target)\n  mu_recurrence: DO\n    mu = 1\n    dmu = 1._DP\n    DO'
        )
    else:
        content = content.replace(
            '  100  mu = 1\n  dmu = 1._SP\n  DO',
            '  ! Forward mu-recurrence loop (was GOTO 100 target)\n  mu_recurrence: DO\n    mu = 1\n    dmu = 1._SP\n    DO'
        )

    # Replace GOTO 100 with CYCLE mu_recurrence
    content = content.replace(
        '        GOTO 100\n      END IF\n    END IF\n  END DO',
        '        CYCLE mu_recurrence  ! Restart mu recurrence (was GOTO 100)\n      END IF\n    END IF\n    END DO\n  END DO mu_recurrence'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (named CYCLE)!')
