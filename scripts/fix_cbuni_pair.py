#!/usr/bin/env python
"""Fix GOTO 100 in cbuni.f90 and zbuni.f90 - Bessel uniform I/K asymptotic.

Pattern: Success skip forward GOTO
Replace with: Inline IF block with RETURN

Reference: ISO/IEC 1539-1:2018 S11.1.1 (IF construct)
Original: Amos (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/cbuni.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/zbuni.f90', 'DP'),
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
  !   251217  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.1 (IF construct)
  !           Original: Amos (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 100 with inline IF block
    content = content.replace(
        '    IF( nw>=0 ) GOTO 100\n  ELSE',
        '    IF( nw>=0 ) THEN\n      Nz = nw  ! Success (was GOTO 100)\n      RETURN\n    END IF\n  ELSE'
    )

    # Remove label 100 section
    content = content.replace(
        '\n  100  Nz = nw\n  !\n  RETURN\nEND SUBROUTINE',
        '\n  ! (Label 100 removed - handled inline)\nEND SUBROUTINE'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (inline IF)!')
