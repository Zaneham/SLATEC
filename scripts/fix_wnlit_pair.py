#!/usr/bin/env python3
"""Fix GOTO in dwnlit.f90 and wnlit.f90 - WNNLS subsidiary"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/approximation/dwnlit.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/approximation/wnlit.f90', 'SP'),
]

for filepath, prec in files:
    with open(filepath, 'r') as f:
        content = f.read()

    # Add revision history
    old_rev = '!   910408  Updated AUTHOR section.  (WRB)'
    new_rev = '''!   910408  Updated AUTHOR section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT with flag)
  !           Original: Hanson & Haskell, SNLA'''

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 100 with EXIT and flag - the pattern is early exit from nested loops
    old_goto = '''      ELSE
        krank = i - 1
        GOTO 100
      END IF
      EXIT
    END DO
  END DO
  krank = l1
  !
  100 CONTINUE'''

    new_goto = '''      ELSE
        krank = i - 1
        EXIT  ! Early exit - krank set
      END IF
      EXIT
    END DO
    EXIT  ! Propagate early exit
  END DO
  IF( krank < 0 ) krank = l1  ! Only set if not early exit
  !
  ! (Label 100 removed)'''

    content = content.replace(old_goto, new_goto)

    # Initialize krank to -1 to detect early exit (find the declaration)
    if 'krank = -1' not in content:
        # Add initialization after declaration
        old_init = 'krank = 0'
        new_init = 'krank = -1  ! Sentinel for early exit detection'
        content = content.replace(old_init, new_init)

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 100 eliminated!')
