#!/usr/bin/env python
"""Fix remaining GOTO 300 in dpchid.f90 and pchid.f90 - PCHIP integral.

Pattern: Bounds check error exit GOTO 300
Replace with: Inline IF-THEN ERROR STOP

Reference: ISO/IEC 1539-1:2018 S11.2.3 (ERROR STOP)
Original: Fritsch, LLNL (PCHIP)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/interpolation/pchip/dpchid.f90', 'DP', 'DPCHID'),
    ('C:/dev/slatec-modern/src/modern/interpolation/pchip/pchid.f90', 'SP', 'PCHID'),
]

for filepath, prec, name in files:
    fname = os.path.basename(filepath)

    if not os.path.exists(filepath):
        print(f'{fname}: NOT FOUND, skipping')
        continue

    with open(filepath, 'r') as f:
        content = f.read()

    # Skip if GOTO 300 already fixed
    if 'GOTO 300' not in content:
        print(f'{fname}: GOTO 300 already fixed, skipping')
        continue

    # Update revision history
    old_rev = '!           Original: Fritsch, LLNL (PCHIP)'
    new_rev = '!           Original: Fritsch, LLNL (PCHIP)\n  !   251217  Eliminated GOTO 300 per MODERNISATION_GUIDE.md S1. (ZH)'

    content = content.replace(old_rev, new_rev)

    # Replace the two GOTO 300 statements with inline error
    content = content.replace(
        'IF( (Ia<1) .OR. (Ia>N) ) GOTO 300\n  IF( (Ib<1) .OR. (Ib>N) ) GOTO 300',
        f"IF( (Ia<1) .OR. (Ia>N) .OR. (Ib<1) .OR. (Ib>N) ) THEN\n    ERROR STOP '{name} : IA OR IB OUT OF RANGE'  ! (Was GOTO 300)\n  END IF"
    )

    # Remove label 300 at end
    content = content.replace(
        f"\n  !     IA OR IB OUT OF RANGE RETURN.\n  300 ERROR STOP '{name} : IA OR IB OUT OF RANGE'",
        '\n  ! (Label 300 removed - handled inline)'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 300 eliminated (inline error check)!')
