#!/usr/bin/env python
"""Fix GOTO 20 in bvalu.f90 and dbvalu.f90 - B-spline evaluation.

Pattern: Search loop with GOTO exit on success
Replace with: EXIT and loop index check

Reference: ISO/IEC 1539-1:2018 S11.2.1 (EXIT statement)
Original: de Boor, K. (A Practical Guide to Splines)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/interpolation/bspline/bvalu.f90', 'SP', 'BVALU'),
    ('C:/dev/slatec-modern/src/modern/interpolation/bspline/dbvalu.f90', 'DP', 'DBVALU'),
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
  !   211001  Converted to free-form, added INTENT.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 20 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT statement)
  !           Original: de Boor (A Practical Guide to Splines)"""

    content = content.replace(old_rev, new_rev)

    # Replace GOTO 20 with EXIT
    content = content.replace(
        'IF( X/=T(i) ) GOTO 20',
        'IF( X/=T(i) ) EXIT  ! Found suitable index'
    )

    # Add error check after the loop
    old_error = """          END DO
          ERROR STOP '""" + name + """ : A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)'
        END IF
      END IF
      !
      !- ** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
      !     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
      !
      20  imk = i - K"""

    new_error = """          END DO
          IF( i==K ) ERROR STOP '""" + name + """ : A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)'
        END IF
      END IF
      !
      !- ** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
      !     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
      !     (Label 20 removed - search exit handled by EXIT)
      !
      imk = i - K"""

    content = content.replace(old_error, new_error)

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 20 eliminated (EXIT + check)!')
