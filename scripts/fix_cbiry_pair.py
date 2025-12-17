#!/usr/bin/env python
"""Fix GOTO 50/100 in cbiry.f90 and zbiry.f90 - Bessel Bi (Airy function).

Pattern: Error exits (overflow GOTO 50, other error GOTO 100)
Replace with: Inline error returns

Reference: ISO/IEC 1539-1:2018 S11.2.2 (RETURN statement)
Original: Amos (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/cbiry.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/zbiry.f90', 'DP'),
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
    old_rev = '!   920128  Category corrected.  (WRB)'
    new_rev = """!   920128  Category corrected.  (WRB)
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 50/100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.2 (RETURN statement)
  !           Original: Amos (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 50 with inline overflow error return
    content = content.replace(
        'IF( bb>elim ) GOTO 50\n        END IF',
        'IF( bb>elim ) THEN\n            nz = 0\n            Ierr = 2  ! Overflow (was GOTO 50)\n            RETURN\n          END IF\n        END IF'
    )

    # Replace GOTO 100 with inline error return
    content = content.replace(
        'GOTO 100\n      END IF\n    END IF\n    50  nz = 0\n    Ierr = 2\n    RETURN',
        'nz = 0\n        Ierr = 5  ! Error from CBINU (was GOTO 100)\n        RETURN\n      END IF\n    END IF\n    ! (Label 50 removed - handled inline)'
    )

    # Remove label 100 at end
    content = content.replace(
        '\n  100  nz = 0\n  Ierr = 5\n  !\n  RETURN\nEND SUBROUTINE',
        '\n  ! (Label 100 removed - handled inline)\nEND SUBROUTINE'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 50/100 eliminated (inline error returns)!')
