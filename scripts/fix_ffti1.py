#!/usr/bin/env python3
"""Fix GOTO in FFTPACK initialization files - prime factorization loop"""

import os

files = [
    'C:/dev/slatec-modern/src/modern/integ_trans/fftpack/cffti1.f90',
    'C:/dev/slatec-modern/src/modern/integ_trans/fftpack/rffti1.f90',
    'C:/dev/slatec-modern/src/modern/integ_trans/fftpack/ezfft1.f90',
]

for filepath in files:
    if not os.path.exists(filepath):
        print(f'{filepath}: NOT FOUND, skipping')
        continue

    with open(filepath, 'r') as f:
        content = f.read()

    # Skip if already fixed
    if 'Eliminated GOTO' in content:
        continue

    # Skip if no GOTO 100 pattern
    if 'GOTO 100' not in content:
        continue

    # Add revision history
    old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
    new_rev = """!   920501  Reformatted the REFERENCES section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S8.1.10.2 (Named DO construct, CYCLE)
  !           Original: Swarztrauber, NCAR (FFTPACK)"""

    content = content.replace(old_rev, new_rev)

    # The pattern: label 100 starts outer loop, GOTO 100 continues it
    # Convert to named DO construct with CYCLE
    old_pattern = """  100  j = j + 1
  IF( j<=4 ) THEN
    ntry = ntryh(j)
  ELSE
    ntry = ntry + 2
  END IF
  DO
    nq = nl/ntry
    nr = nl - ntry*nq
    IF( nr/=0 ) GOTO 100"""

    new_pattern = """  ! Prime factorization loop (was GOTO 100)
  factor_loop: DO
    j = j + 1
    IF( j<=4 ) THEN
      ntry = ntryh(j)
    ELSE
      ntry = ntry + 2
    END IF
    inner: DO
      nq = nl/ntry
      nr = nl - ntry*nq
      IF( nr/=0 ) CYCLE factor_loop  ! Try next trial factor"""

    content = content.replace(old_pattern, new_pattern)

    # Also need to fix the EXIT at the end of the loop
    # Change "EXIT" to "EXIT factor_loop" to exit the named outer loop
    content = content.replace('      EXIT\n    END IF\n  END DO', '      EXIT factor_loop  ! Factorization complete\n    END IF\n    END DO inner\n  END DO factor_loop')

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 100 eliminated (factor loop -> named CYCLE)!')
