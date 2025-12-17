#!/usr/bin/env python
"""Fix GOTO 100 in dcov.f90 and scov.f90 - MINPACK covariance matrix.

Pattern: Error exit GOTO 100 to cleanup
Replace with: BLOCK with EXIT

Reference: ISO/IEC 1539-1:2018 S8.1.6 (BLOCK construct)
Original: More, Garbow, Hillstrom (ANL)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/approximation/minpack/dcov.f90', 'DP', 'DCOV'),
    ('C:/dev/slatec-modern/src/modern/approximation/minpack/scov.f90', 'SP', 'SCOV'),
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
  !           Ref: ISO/IEC 1539-1:2018 S8.1.6 (BLOCK construct)
  !           Original: More, Garbow, Hillstrom (ANL)"""

    content = content.replace(old_rev, new_rev)

    # Add BLOCK after the early returns
    content = content.replace(
        '  !\n  !     THE CALCULATION OF SIGMA',
        '  !\n  process: BLOCK  ! Main computation block (was GOTO 100 target)\n  !\n  !     THE CALCULATION OF SIGMA'
    )

    # Replace first GOTO 100 with EXIT process
    content = content.replace(
        'IF( iflag<0 ) GOTO 100\n          temp = Fvec(i)',
        'IF( iflag<0 ) EXIT process  ! Error (was GOTO 100)\n          temp = Fvec(i)'
    )

    # Replace second GOTO 100 with EXIT process
    content = content.replace(
        'IF( iflag<0 ) GOTO 100\n        !\n        !     COMPUTE THE QR',
        'IF( iflag<0 ) EXIT process  ! Error (was GOTO 100)\n        !\n        !     COMPUTE THE QR'
    )

    # Replace label 100 CONTINUE with END BLOCK
    content = content.replace(
        '  END IF\n  !\n  100 CONTINUE\n  IF( M<=0',
        '  END IF\n  END BLOCK process  ! (Label 100 removed)\n  !\n  IF( M<=0'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (BLOCK/EXIT)!')
