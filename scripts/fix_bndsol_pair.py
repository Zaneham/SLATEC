#!/usr/bin/env python
"""Fix GOTO 100 in bndsol.f90 and dbndsl.f90 - banded system solver.

Pattern: Error exit GOTO 100 for zero diagonal
Replace with: Inline error handling

Reference: ISO/IEC 1539-1:2018 S11.2.3 (ERROR STOP)
Original: Hanson (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/linear/bndsol.f90', 'SP', 'BNDSOL'),
    ('C:/dev/slatec-modern/src/modern/linear/dbndsl.f90', 'DP', 'DBNDSL'),
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
  !   251217  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.3 (ERROR STOP)
  !           Original: Hanson (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # The error message (for replacement)
    err_msg = f"'{name} : A ZERO DIAGONAL TERM IS IN THE N BY N UPPER TRIANGULAR MATRIX.'"

    # Replace first GOTO 100
    content = content.replace(
        'IF( G(j,l+1)==0 ) GOTO 100\n        X(j)',
        f'IF( G(j,l+1)==0 ) THEN\n          nerr = 1\n          iopt = 2\n          ERROR STOP {err_msg}\n        END IF\n        X(j)'
    )

    # Replace second GOTO 100
    content = content.replace(
        'IF( G(i,l+1)==0 ) GOTO 100\n    X(i)',
        f'IF( G(i,l+1)==0 ) THEN\n      nerr = 1\n      iopt = 2\n      ERROR STOP {err_msg}\n    END IF\n    X(i)'
    )

    # Remove label 100 section
    content = content.replace(
        f"\n  100  nerr = 1\n  iopt = 2\n  ERROR STOP {err_msg}\n\nEND SUBROUTINE",
        '\n  ! (Label 100 removed - error handled inline)\nEND SUBROUTINE'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (inline ERROR STOP)!')
