#!/usr/bin/env python
"""Fix GOTO 50 in dplpfe.f90 and splpfe.f90 - SPLP find entering variable.

Pattern: GOTO 50 to skip to loop increment (CYCLE equivalent)
Replace with: CYCLE in restructured DO loop

Reference: ISO/IEC 1539-1:2018 S11.2.1 (CYCLE statement)
Original: Hanson & Hiebert (SPLP)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/optimization/splp/dplpfe.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/optimization/splp/splpfe.f90', 'SP'),
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
  !   251216  Eliminated GOTO 50 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (CYCLE statement)
  !           Original: Hanson & Hiebert (SPLP)"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 50 with CYCLE (the DO WHILE will handle the increment)
    # But first we need to restructure the loop slightly

    # Change the DO WHILE with manual increment to a proper DO loop
    content = content.replace(
        '  n20002 = Mrelas + Nvars\n  DO WHILE( (n20002-i)>=0 )',
        '  ! Search for variable to enter basis\n  DO i = Mrelas + 1, Mrelas + Nvars'
    )

    # Replace GOTO 50 with CYCLE
    content = content.replace(
        ') GOTO 50',
        ') CYCLE  ! Skip equation variables (was GOTO 50)'
    )

    # Remove the label and manual increment
    content = content.replace(
        '    50  i = i + 1\n  END DO',
        '    ! (Label 50 removed - handled by CYCLE)\n  END DO'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 50 eliminated (CYCLE in DO loop)!')
