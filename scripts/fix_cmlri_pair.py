#!/usr/bin/env python
"""Fix GOTO 100/200 in cmlri.f90 and zmlri.f90 - Bessel I Miller algorithm.

Pattern: Convergence exit from loops (GOTO 100/200)
Replace with: EXIT + loop index check

Reference: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
Original: Amos (SNLA)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/special_functions/cmlri.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/special_functions/zmlri.f90', 'DP'),
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
    old_rev = '!   910415  Prologue converted to Version 4.0 format.  (BAB)'
    new_rev = """!   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
  !           Original: Amos (SNLA)"""

    content = content.replace(old_rev, new_rev)

    # Fix first GOTO 100 pattern
    content = content.replace(
        'IF( ap>tst*ak*ak ) GOTO 100\n    ak = ak + 1._' + prec + '\n  END DO\n  Nz = -2\n  RETURN\n  100  i = i + 1',
        'IF( ap>tst*ak*ak ) EXIT  ! Converged (was GOTO 100)\n    ak = ak + 1._' + prec + '\n  END DO\n  IF( i > 80 ) THEN  ! Loop exhausted\n    Nz = -2\n    RETURN\n  END IF\n  i = i + 1  ! (Label 100 removed)'
    )

    # Fix second GOTO 200 pattern
    content = content.replace(
        'IF( itime==2 ) GOTO 200\n        ack = ABS(ck)\n        flam = ack + SQRT(ack*ack-1._' + prec + ')\n        fkap = ap/ABS(p1)\n        rho = MIN(flam,fkap)\n        tst = tst*SQRT(rho/(rho*rho-1._' + prec + '))\n        itime = 2\n      END IF\n    END DO\n    Nz = -2\n    RETURN\n  END IF\n  !-----------------------------------------------------------------------\n  !     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION\n  !-----------------------------------------------------------------------\n  200  k = k + 1',
        'IF( itime==2 ) EXIT  ! Converged (was GOTO 200)\n        ack = ABS(ck)\n        flam = ack + SQRT(ack*ack-1._' + prec + ')\n        fkap = ap/ABS(p1)\n        rho = MIN(flam,fkap)\n        tst = tst*SQRT(rho/(rho*rho-1._' + prec + '))\n        itime = 2\n      END IF\n    END DO\n    IF( k > 80 ) THEN  ! Loop exhausted\n      Nz = -2\n      RETURN\n    END IF\n  END IF\n  !-----------------------------------------------------------------------\n  !     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION\n  !-----------------------------------------------------------------------\n  k = k + 1  ! (Label 200 removed)'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100/200 eliminated (EXIT + check)!')
