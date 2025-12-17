#!/usr/bin/env python
"""Fix GOTO 100/200/300 in pois3d.f90 - FISHPACK 3D Poisson solver.

Pattern: Validation loop with early exit error, skip validation
Replace with: Conditional block with named loop EXIT

Reference: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
Original: Adams (NCAR)
"""

import os

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/pois3d.f90'
fname = os.path.basename(filepath)

if not os.path.exists(filepath):
    print(f'{fname}: NOT FOUND')
    exit(1)

with open(filepath, 'r') as f:
    content = f.read()

# Skip if already fixed
if 'Eliminated GOTO' in content:
    print(f'{fname}: already fixed, skipping')
    exit(0)

# Add revision history
old_rev = '!   891214  Prologue converted to Version 4.0 format.  (BAB)'
new_rev = old_rev + """
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200/300 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.12 (EXIT statement)
  !           Original: Adams (NCAR)"""

content = content.replace(old_rev, new_rev)

# Replace the GOTO pattern with structured code
old_pattern = '''  IF( np/=1 ) GOTO 200
  DO k = 1, N
    IF( A(k)/=C(1) ) GOTO 100
    IF( C(k)/=C(1) ) GOTO 100
    IF( B(k)/=B(1) ) GOTO 100
  END DO
  GOTO 300
  100  Ierror = 9
  200 CONTINUE
  IF( Nperod==1 .AND. (A(1)/=0. .OR. C(N)/=0._SP) ) Ierror = 10
  300 CONTINUE'''

new_pattern = '''  IF( np==1 ) THEN  ! Check coefficient consistency (was GOTO 200 skip)
    check_coeff: DO k = 1, N
      IF( A(k)/=C(1) .OR. C(k)/=C(1) .OR. B(k)/=B(1) ) THEN
        Ierror = 9  ! (Was GOTO 100)
        EXIT check_coeff
      END IF
    END DO check_coeff
  END IF
  ! (Labels 100, 200, 300 removed)
  IF( Nperod==1 .AND. (A(1)/=0. .OR. C(N)/=0._SP) ) Ierror = 10'''

content = content.replace(old_pattern, new_pattern)

with open(filepath, 'w') as f:
    f.write(content)

print(f'{fname}: GOTO 100/200/300 eliminated (named loop with EXIT)!')
