#!/usr/bin/env python
"""Fix GOTO 100 in dlssud.f90 and lssuds.f90 - underdetermined system solver.

Pattern: Success path GOTO 100 skipping error handling
Replace with: BLOCK construct with EXIT

Reference: ISO/IEC 1539-1:2018 S8.1.6 (BLOCK construct)
Original: Watts (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/bvsup/dlssud.f90', 'DP', 'DLSSUD'),
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/bvsup/lssuds.f90', 'SP', 'LSSUDS'),
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
  !           Original: Watts (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Add BLOCK after FIRST EXECUTABLE STATEMENT
    content = content.replace(
        '!* FIRST EXECUTABLE STATEMENT',
        '!* FIRST EXECUTABLE STATEMENT\n  success: BLOCK  ! Main processing block (was GOTO 100 target)'
    )

    # Replace first GOTO 100 (after storing divisors)
    content = content.replace(
        'END DO\n          !        .........EXIT\n          GOTO 100',
        'END DO\n          EXIT success  ! Success (was GOTO 100)'
    )
    # SP version has different comment structure
    content = content.replace(
        'END DO\n          GOTO 100\n        ELSE',
        'END DO\n          EXIT success  ! Success (was GOTO 100)\n        ELSE'
    )

    # Replace second GOTO 100 (Iflag==1 case)
    content = content.replace(
        'ELSEIF( Iflag==1 ) THEN\n        GOTO 100\n      END IF',
        'ELSEIF( Iflag==1 ) THEN\n        EXIT success  ! Already decomposed (was GOTO 100)\n      END IF'
    )
    content = content.replace(
        'ELSEIF( Iflag==1 ) THEN\n      GOTO 100\n    END IF',
        'ELSEIF( Iflag==1 ) THEN\n      EXIT success  ! Already decomposed (was GOTO 100)\n    END IF'
    )

    # Replace 100 CONTINUE with END BLOCK
    content = content.replace(
        '  100 CONTINUE\n  IF( Irank>0 )',
        '  END BLOCK success  ! (Label 100 removed)\n  IF( Irank>0 )'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (BLOCK/EXIT)!')
