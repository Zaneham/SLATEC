#!/usr/bin/env python
"""Fix GOTO 100 in cscale.f90 and dcscal.f90 - BVP column scaling.

Pattern: Early exit GOTO 100 from condition check loop
Replace with: EXIT + loop index check

Reference: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
Original: Watts, Shampine (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/bvsup/cscale.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/bvsup/dcscal.f90', 'DP'),
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
    old_rev = '!   910722  Replaced DATA statements with I1MACH calls.'
    new_rev = """!   910722  Replaced DATA statements with I1MACH calls.
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
  !           Original: Watts, Shampine (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 100 with EXIT
    content = content.replace(
        ') GOTO 100\n      IF( (cs<1._' + prec + '/ten20) .OR. (cs>ten20) ) GOTO 100\n    END DO\n  END IF\n  !\n  DO k = 1, Ncol\n    Scales(k) = 1._' + prec + '\n  END DO\n  RETURN\n  !\n  100  alog2',
        ') EXIT  ! Needs scaling\n      IF( (cs<1._' + prec + '/ten20) .OR. (cs>ten20) ) EXIT  ! Needs scaling\n    END DO\n  END IF\n  !\n  IF( k > Ncol ) THEN  ! Loop completed - no scaling needed\n    DO k = 1, Ncol\n      Scales(k) = 1._' + prec + '\n    END DO\n    RETURN\n  END IF\n  !\n  ! Scaling required (was label 100)\n  alog2'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (EXIT + check)!')
