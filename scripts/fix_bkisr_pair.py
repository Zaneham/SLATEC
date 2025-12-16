#!/usr/bin/env python
"""Fix GOTO 100 in bkisr.f90 and dbkisr.f90 - Bessel K integral series.

Pattern: Convergence exit from DO loop via GOTO 100
Replace with: EXIT and loop index check

Reference: ISO/IEC 1539-1:2018 S11.2.1 (EXIT statement)
Original: D.E. Amos, SNLA (Bessel functions)
"""

import os
import re

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/bkisr.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/dbkisr.f90', 'DP'),
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
    old_rev = '!   910722  Updated AUTHOR section.  (ALS)'
    new_rev = """!   910722  Updated AUTHOR section.  (ALS)
  !   211001  Converted to free-form, added INTENT, ELEMENTAL.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT on convergence)
  !           Original: D.E. Amos, SNLA"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 100 with EXIT
    content = re.sub(
        r'\) GOTO 100',
        ') EXIT  ! Converged',
        content
    )

    # Replace "END DO / Ierr = 2 / RETURN" with check for non-convergence
    old_pattern = """    END DO
    Ierr = 2
    RETURN"""
    new_pattern = """    END DO
    IF( k > 20 ) THEN
      Ierr = 2
      RETURN
    END IF"""
    content = content.replace(old_pattern, new_pattern)

    # Remove label 100 and add comment
    content = re.sub(
        r'\n\s*100\s+Summ = ',
        '\n  ! (Label 100 removed - convergent path continues here)\n  Summ = ',
        content
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (convergence EXIT)!')
