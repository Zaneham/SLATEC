#!/usr/bin/env python3
"""Fix GOTO in dplint.f90 and polint.f90 - polynomial interpolation"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/interpolation/dplint.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/interpolation/polint.f90', 'SP'),
]

for filepath, prec in files:
    with open(filepath, 'r') as f:
        content = f.read()

    # Add revision history
    old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
    new_rev = """!   920501  Reformatted the REFERENCES section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.3 (RETURN/ERROR STOP)
  !           Original: Huddleston, SNLL; Shampine & Davenport, SLA-74-0270"""

    content = content.replace(old_rev, new_rev)

    # The GOTO 100 is an error exit for non-distinct abscissas
    # Replace with direct ERROR STOP
    if prec == 'DP':
        old_goto = """        IF( dif==0._DP ) GOTO 100
        C(k) = (C(i)-C(k))/dif
      END DO
    END DO
    RETURN
  END IF
  100  ERROR STOP 'DPLINT : THE ABSCISSAS ARE NOT DISTINCT.'"""
        new_goto = """        IF( dif==0._DP ) ERROR STOP 'DPLINT : THE ABSCISSAS ARE NOT DISTINCT.'
        C(k) = (C(i)-C(k))/dif
      END DO
    END DO
    ! (Label 100 removed - error handled inline)
  END IF"""
    else:
        # SP version uses 0.0 not 0._SP
        old_goto = """        IF( dif==0.0 ) GOTO 100
        C(k) = (C(i)-C(k))/dif
      END DO
    END DO
    RETURN
  END IF
  100 ERROR STOP 'POLINT : THE ABSCISSAS ARE NOT DISTINCT.'"""
        new_goto = """        IF( dif==0.0 ) ERROR STOP 'POLINT : THE ABSCISSAS ARE NOT DISTINCT.'
        C(k) = (C(i)-C(k))/dif
      END DO
    END DO
    ! (Label 100 removed - error handled inline)
  END IF"""

    content = content.replace(old_goto, new_goto)

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 100 eliminated!')
