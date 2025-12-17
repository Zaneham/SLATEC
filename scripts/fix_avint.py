#!/usr/bin/env python
"""Fix GOTO 100/200 in avint.f90 - tabulated function integration.

Pattern: Error exits GOTO 100 (< 3 values) and GOTO 200 (not increasing)
Replace with: Inline error handling

Reference: ISO/IEC 1539-1:2018 S11.2.3 (ERROR STOP)
Original: Jones (SNLA)
"""

import os

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ/avint.f90'
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
old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
new_rev = old_rev + """
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.3 (ERROR STOP)
  !           Original: Jones (SNLA)"""

content = content.replace(old_rev, new_rev)

# Fix GOTO 200 (abscissas not increasing) - make it a named loop with EXIT
content = content.replace(
    '''      DO i = 2, N
        IF( X(i)<=X(i-1) ) GOTO 200
        IF( X(i)>Xup ) EXIT
      END DO''',
    '''      check_order: DO i = 2, N
        IF( X(i)<=X(i-1) ) THEN
          Ierr = 4
          ERROR STOP 'AVINT : THE ABSCISSAS WERE NOT STRICTLY INCREASING. &
            &MUST HAVE X(I-1) < X(I) FOR ALL I.'  ! (Was GOTO 200)
          RETURN
        END IF
        IF( X(i)>Xup ) EXIT check_order
      END DO check_order'''
)

# Fix first two GOTO 100 (data range checks)
content = content.replace(
    '''      IF( N>=3 ) THEN
        IF( X(N-2)<Xlo ) GOTO 100
        IF( X(3)>Xup ) GOTO 100''',
    '''      IF( N>=3 ) THEN
        IF( X(N-2)<Xlo .OR. X(3)>Xup ) THEN
          Ierr = 3
          ERROR STOP 'AVINT : THERE WERE LESS THAN THREE FUNCTION VALUES BETWEEN THE &
            &LIMITS OF INTEGRATION.'  ! (Was GOTO 100)
          RETURN
        END IF'''
)

# Fix third GOTO 100 (too few points in range)
content = content.replace(
    '''        IF( (inrt-inlft)<2 ) GOTO 100
        istart = inlft''',
    '''        IF( (inrt-inlft)<2 ) THEN
          Ierr = 3
          ERROR STOP 'AVINT : THERE WERE LESS THAN THREE FUNCTION VALUES BETWEEN THE &
            &LIMITS OF INTEGRATION.'  ! (Was GOTO 100)
          RETURN
        END IF
        istart = inlft'''
)

# Remove label 100 and 200 at end
content = content.replace(
    '''  100  Ierr = 3
  ERROR STOP 'AVINT : THERE WERE LESS THAN THREE FUNCTION VALUES BETWEEN THE &
    &LIMITS OF INTEGRATION.'
  RETURN
  200  Ierr = 4
  ERROR STOP 'AVINT : THE ABSCISSAS WERE NOT STRICTLY INCREASING. &
    &MUST HAVE X(I-1) < X(I) FOR ALL I.'
  !
  RETURN''',
    '''  ! (Labels 100, 200 removed - errors handled inline)'''
)

with open(filepath, 'w') as f:
    f.write(content)

print(f'{fname}: GOTO 100/200 eliminated (inline errors)!')
