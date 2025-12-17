#!/usr/bin/env python
"""Fix GOTO 50 in dqawfe.f90 and qawfe.f90 - QUADPACK Fourier integral.

Pattern: Early exit GOTO 50 to partial sum result
Replace with: Inline result assignment + RETURN

Reference: ISO/IEC 1539-1:2018 S11.2.1 (RETURN)
Original: Piessens, de Doncker (K. U. Leuven)
"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/diff_integ/quadpack/dqawfe.f90', 'DP', 'DQAWFE'),
    ('C:/dev/slatec-modern/src/modern/diff_integ/quadpack/qawfe.f90', 'SP', 'QAWFE'),
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
    old_rev = '!   891214  Prologue converted to Version 4.0 format.  (BAB)'
    new_rev = old_rev + """
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 50 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (RETURN)
  !           Original: Piessens, de Doncker (K. U. Leuven)"""

    content = content.replace(old_rev, new_rev)

    # Define the inline exit code
    if prec == 'DP':
        exit_code = """Result = psum(numrl2)  ! (Was GOTO 50)
          Abserr = errsum + drl
          RETURN"""
        # First GOTO 50 (inside loop, accuracy test)
        content = content.replace(
            'IF( (errsum+drl)<=Epsabs .AND. Lst>=6 ) GOTO 50',
            'IF( (errsum+drl)<=Epsabs .AND. Lst>=6 ) THEN\n          ' + exit_code + '\n        END IF'
        )
        # Second GOTO 50 (inside loop, error recovery)
        content = content.replace(
            'IF( Ier==7 .AND. (errsum+drl)<=correc*0.1D+02 .AND. Lst>5 ) GOTO 50',
            'IF( Ier==7 .AND. (errsum+drl)<=correc*0.1D+02 .AND. Lst>5 ) THEN\n          ' + exit_code + '\n        END IF'
        )
        # Third GOTO 50 (after loop, comparison)
        content = content.replace(
            'IF( Abserr>errsum ) GOTO 50\n        IF( psum(numrl2)==0._DP ) RETURN',
            'IF( Abserr>errsum ) THEN\n          ' + exit_code + '\n        END IF\n        IF( psum(numrl2)==0._DP ) RETURN'
        )
    else:
        exit_code = """Result = psum(numrl2)  ! (Was GOTO 50)
          Abserr = errsum + drl
          RETURN"""
        # First GOTO 50
        content = content.replace(
            'IF( errsum+drl<=Epsabs .AND. Lst>=6 ) GOTO 50',
            'IF( errsum+drl<=Epsabs .AND. Lst>=6 ) THEN\n          ' + exit_code + '\n        END IF'
        )
        # Second GOTO 50
        content = content.replace(
            'IF( Ier==7 .AND. (errsum+drl)<=correc*10._SP .AND. Lst>5 ) GOTO 50',
            'IF( Ier==7 .AND. (errsum+drl)<=correc*10._SP .AND. Lst>5 ) THEN\n          ' + exit_code + '\n        END IF'
        )
        # Third GOTO 50
        content = content.replace(
            'IF( Abserr>errsum ) GOTO 50\n        IF( psum(numrl2)==0._SP ) RETURN',
            'IF( Abserr>errsum ) THEN\n          ' + exit_code + '\n        END IF\n        IF( psum(numrl2)==0._SP ) RETURN'
        )

    # Remove label 50 and convert to regular statements (for fall-through case)
    content = content.replace(
        '    50  Result = psum(numrl2)\n    Abserr = errsum + drl\n  END IF',
        '    ! (Label 50 removed - exits inlined)\n    Result = psum(numrl2)\n    Abserr = errsum + drl\n  END IF'
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 50 eliminated (inline exit)!')
