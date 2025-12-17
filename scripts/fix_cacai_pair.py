#!/usr/bin/env python
"""Fix GOTO 100 in cacai.f90 and zacai.f90 - Bessel K analytic continuation.

Pattern: Error exit GOTO 100 after subroutine calls
Replace with: Inline error return

Reference: ISO/IEC 1539-1:2018 S11.2.2 (RETURN statement)
Original: Amos (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/cacai.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/zacai.f90', 'DP'),
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
  !   251217  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.2 (RETURN statement)
  !           Original: Amos (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Replace first GOTO 100 with inline error return
    content = content.replace(
        'IF( nw<0 ) GOTO 100\n    ELSE',
        'IF( nw<0 ) THEN\n        Nz = -1\n        IF( nw==(-2) ) Nz = -2\n        RETURN\n      END IF\n    ELSE'
    )

    # Replace second GOTO 100 with inline error return
    content = content.replace(
        'IF( nw<0 ) GOTO 100\n    END IF',
        'IF( nw<0 ) THEN\n        Nz = -1\n        IF( nw==(-2) ) Nz = -2\n        RETURN\n      END IF\n    END IF'
    )

    # Remove label 100 section at end (but keep the closing)
    content = content.replace(
        '  END IF\n  100  Nz = -1\n  IF( nw==(-2) ) Nz = -2\n  !\nEND SUBROUTINE',
        '  END IF\n  ! (Label 100 removed - error handled inline)\nEND SUBROUTINE'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (inline error return)!')
