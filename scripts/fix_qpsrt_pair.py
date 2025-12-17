#!/usr/bin/env python
"""Fix GOTO 100/200/300 in dqpsrt.f90 and qpsrt.f90 - QUADPACK priority sorting.

Pattern: Early exit from search loops + skip to end
Replace with: Named loops with EXIT + loop index check

Reference: ISO/IEC 1539-1:2018 S11.1.7.4.4 (DO construct), S11.1.12 (EXIT)
Original: Piessens, de Doncker (K. U. Leuven)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/diff_integ/quadpack/dqpsrt.f90', 'DP', 'DQPSRT'),
    ('C:/dev/slatec-modern/src/modern/diff_integ/quadpack/qpsrt.f90', 'SP', 'QPSRT'),
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
    old_rev = '!   900328  Added TYPE section.  (WRB)'
    new_rev = old_rev + """
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200/300 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.7.4.4, S11.1.12 (DO, EXIT)
  !           Original: Piessens, de Doncker (K. U. Leuven)"""

    content = content.replace(old_rev, new_rev)

    # Name the first search loop
    content = content.replace(
        'DO i = ibeg, jbnd\n        isucc = Iord(i)\n        !- **JUMP OUT OF DO-LOOP\n        IF( errmax>=Elist(isucc) ) GOTO 100',
        'search_errmax: DO i = ibeg, jbnd\n        isucc = Iord(i)\n        IF( errmax>=Elist(isucc) ) EXIT search_errmax  ! Found insertion point (was GOTO 100)'
    )

    # Close the first search loop
    content = content.replace(
        'Iord(i-1) = isucc\n      END DO\n    END IF\n    Iord(jbnd) = Maxerr\n    Iord(jupbn) = Last\n  ELSE\n    Iord(1) = 1\n    Iord(2) = 2\n  END IF\n  GOTO 300',
        'Iord(i-1) = isucc\n      END DO search_errmax\n      IF( i<=jbnd ) THEN  ! Found early (was label 100)\n        Iord(i-1) = Maxerr\n        k = jbnd\n        search_errmin: DO j = i, jbnd\n          isucc = Iord(k)\n          IF( errmin<Elist(isucc) ) EXIT search_errmin  ! Found (was GOTO 200)\n          Iord(k+1) = isucc\n          k = k - 1\n        END DO search_errmin\n        IF( j<=jbnd ) THEN\n          Iord(k+1) = Last  ! (Was label 200)\n        ELSE\n          Iord(i) = Last  ! (Was after second loop)\n        END IF\n      ELSE  ! First loop completed\n        Iord(jbnd) = Maxerr\n        Iord(jupbn) = Last\n      END IF\n    END IF\n  ELSE\n    Iord(1) = 1\n    Iord(2) = 2\n  END IF'
    )

    # Remove the old errmin section (labels 100, 200, GOTO 300)
    content = content.replace(
        '\n  !\n  !           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.\n  !\n  100  Iord(i-1) = Maxerr\n  k = jbnd\n  DO j = i, jbnd\n    isucc = Iord(k)\n    !- **JUMP OUT OF DO-LOOP\n    IF( errmin<Elist(isucc) ) GOTO 200\n    Iord(k+1) = isucc\n    k = k - 1\n  END DO\n  Iord(i) = Last\n  GOTO 300\n  200  Iord(k+1) = Last',
        ''
    )

    # Remove label 300
    content = content.replace(
        '  300  Maxerr = Iord(Nrmax)',
        '  ! Set MAXERR and ERMAX (label 300 removed)\n  Maxerr = Iord(Nrmax)'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100/200/300 eliminated (named loops with EXIT)!')
