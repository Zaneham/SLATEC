#!/usr/bin/env python3
"""Fix GOTO in dpchim.f90 and pchim.f90 - PCHIP monotone interpolation"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/interpolation/pchip/dpchim.f90', 'DP', 'DPCHIM'),
    ('C:/dev/slatec-modern/src/modern/interpolation/pchip/pchim.f90', 'SP', 'PCHIM'),
]

for filepath, prec, name in files:
    with open(filepath, 'r') as f:
        content = f.read()

    # Add revision history
    old_rev = '!   920429  Revised format and order of references.  (WRB,FNF)'
    new_rev = """!   920429  Revised format and order of references.  (WRB,FNF)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 50 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT with error flag)
  !           Original: Fritsch, LLNL (PCHIP); SIAM J. Sci. Stat. Comp. 5, 1984"""

    content = content.replace(old_rev, new_rev)

    # The GOTO 50 pattern: exit loop when X(i) <= X(i-1)
    # Convert to EXIT with flag check after loop
    if prec == 'DP':
        old_loop = """      DO i = 2, N
        IF( X(i)<=X(i-1) ) GOTO 50
      END DO"""
        new_loop = """      DO i = 2, N
        IF( X(i)<=X(i-1) ) EXIT
      END DO
      IF( i <= N ) THEN
        ! X-array not strictly increasing (was GOTO 50)
        Ierr = -3
        ERROR STOP 'DPCHIM : X-ARRAY NOT STRICTLY INCREASING'
      END IF"""
    else:
        old_loop = """      DO i = 2, N
        IF( X(i)<=X(i-1) ) GOTO 50
      END DO"""
        new_loop = """      DO i = 2, N
        IF( X(i)<=X(i-1) ) EXIT
      END DO
      IF( i <= N ) THEN
        ! X-array not strictly increasing (was GOTO 50)
        Ierr = -3
        ERROR STOP 'PCHIM : X-ARRAY NOT STRICTLY INCREASING'
      END IF"""

    content = content.replace(old_loop, new_loop)

    # Remove label 50 and its error handling
    if prec == 'DP':
        old_label = """    END IF
    !
    !     X-ARRAY NOT STRICTLY INCREASING.
    50  Ierr = -3
    ERROR STOP 'DPCHIM : X-ARRAY NOT STRICTLY INCREASING'
  END IF"""
        new_label = """    END IF
    ! (Label 50 removed - error handled inline after loop)
  END IF"""
    else:
        old_label = """    END IF
    !
    !     X-ARRAY NOT STRICTLY INCREASING.
    50 Ierr = -3
    ERROR STOP 'PCHIM : X-ARRAY NOT STRICTLY INCREASING'
  END IF"""
        new_label = """    END IF
    ! (Label 50 removed - error handled inline after loop)
  END IF"""

    content = content.replace(old_label, new_label)

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 50 eliminated!')
