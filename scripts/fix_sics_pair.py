#!/usr/bin/env python
"""Fix GOTO 50 in dsics.f90 and ssics.f90 - SLAP incomplete Cholesky.

Pattern: Nested loop success exit via GOTO 50
Replace with: Named loops with CYCLE

Reference: ISO/IEC 1539-1:2018 S11.2.1 (named CYCLE)
Original: Seager, LLNL (SLAP)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/linear/slap/dsics.f90', 'DP', 'DSICS'),
    ('C:/dev/slatec-modern/src/modern/linear/slap/ssics.f90', 'SP', 'SSICS'),
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
    old_rev = '!   930701  Updated CATEGORY section.  (FNF, WRB)'
    new_rev = """!   930701  Updated CATEGORY section.  (FNF, WRB)
  !   211001  Converted to free-form, added INTENT.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 50 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (named CYCLE)
  !           Original: Seager, LLNL (SLAP)"""

    content = content.replace(old_rev, new_rev)

    # Name the outer loop (DO irr = ...)
    content = content.replace(
        '    DO irr = irbgn, irend',
        '    irr_loop: DO irr = irbgn, irend'
    )

    # Replace GOTO 50 with CYCLE irr_loop
    content = content.replace(
        'GOTO 50',
        'CYCLE irr_loop  ! Success - continue to next irr (was GOTO 50)'
    )

    # Replace the 50 CONTINUE and END DO
    content = content.replace(
        '      50 CONTINUE\n    END DO',
        '      ! (Label 50 removed - success handled by named CYCLE)\n    END DO irr_loop'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 50 eliminated (named CYCLE)!')
