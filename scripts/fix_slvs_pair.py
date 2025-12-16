#!/usr/bin/env python3
"""Fix GOTO in dslvs.f90 and slvs.f90 - DEBDF linear system solver"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/depac/dslvs.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/depac/slvs.f90', 'SP'),
]

for filepath, prec in files:
    with open(filepath, 'r') as f:
        content = f.read()

    # Add revision history
    old_rev = '!   920422  Changed DIMENSION statement.  (WRB)'
    new_rev = '''!   920422  Changed DIMENSION statement.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.3 (RETURN statement)
  !           Original: Watts, SNLA'''

    content = content.replace(old_rev, new_rev)

    # The GOTO 100 sets ier_com = -1 and returns (singular matrix)
    # Replace the GOTO with direct setting and RETURN
    if prec == 'DP':
        old_goto = '''          IF( ABS(di)==0._DP ) GOTO 100
          Wm(i+2) = 1._DP/di'''
        new_goto = '''          IF( ABS(di)==0._DP ) THEN
            ier_com = -1  ! Singular matrix
            RETURN
          END IF
          Wm(i+2) = 1._DP/di'''
    else:
        old_goto = '''          IF( ABS(di)==0._SP ) GOTO 100
          Wm(i+2) = 1._SP/di'''
        new_goto = '''          IF( ABS(di)==0._SP ) THEN
            ier_com = -1  ! Singular matrix
            RETURN
          END IF
          Wm(i+2) = 1._SP/di'''

    content = content.replace(old_goto, new_goto)

    # Remove label 100 and its error handling (now inline)
    old_label = '''  100  ier_com = -1
  !----------------------- END OF SUBROUTINE'''
    new_label = '''  ! (Label 100 removed - error exit now inline)
  !----------------------- END OF SUBROUTINE'''
    content = content.replace(old_label, new_label)

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 100 eliminated!')
