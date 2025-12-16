#!/usr/bin/env python
"""Fix GOTO 100 in mpadd3.f90 - MP extended precision addition.

Pattern: Early exit from loop skipping carry-off handling
Replace with: Flag and EXIT

Reference: ISO/IEC 1539-1:2018 S11.2.1 (EXIT with flag)
Original: Bailey & Brent (MP package)
"""

import os

filepath = 'C:/dev/slatec-modern/src/modern/linear/dqdot/mpadd3.f90'

with open(filepath, 'r') as f:
    content = f.read()

# Skip if already fixed
if 'Eliminated GOTO' in content:
    print('mpadd3.f90: already fixed')
    exit(0)

# Add revision history
old_rev = '!   930701  Updated CATEGORY section.  (FNF, WRB)'
new_rev = """!   930701  Updated CATEGORY section.  (FNF, WRB)
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT with flag)
  !           Original: Bailey & Brent (MP package)"""

content = content.replace(old_rev, new_rev)

# Add flag declaration
old_decl = '  INTEGER :: i, i2, i2p, j, Med'
new_decl = '  INTEGER :: i, i2, i2p, j, Med\n  LOGICAL :: early_exit'

content = content.replace(old_decl, new_decl)

# Replace the DO WHILE with EXIT and flag
old_loop = """    DO WHILE( i>0 )
      c = Y(i+2) + c
      IF( c<b_com ) GOTO 100
      r_com(i) = 0
      c = 1
      i = i - 1
    END DO
    IF( c==0 ) RETURN
    ! MUST SHIFT RIGHT HERE AS CARRY OFF END
    i2p = i2 + 1
    DO j = 2, i2
      i = i2p - j
      r_com(i+1) = r_com(i)
    END DO
    r_com(1) = 1
    Re = Re + 1
    RETURN
  END IF
  100  r_com(i) = c"""

new_loop = """    early_exit = .FALSE.
    DO WHILE( i>0 )
      c = Y(i+2) + c
      IF( c<b_com ) THEN
        early_exit = .TRUE.
        EXIT  ! No carry, exit early (was GOTO 100)
      END IF
      r_com(i) = 0
      c = 1
      i = i - 1
    END DO
    IF( .NOT. early_exit ) THEN
      IF( c==0 ) RETURN
      ! MUST SHIFT RIGHT HERE AS CARRY OFF END
      i2p = i2 + 1
      DO j = 2, i2
        i = i2p - j
        r_com(i+1) = r_com(i)
      END DO
      r_com(1) = 1
      Re = Re + 1
      RETURN
    END IF
  END IF
  ! (Label 100 removed - early exit handled by flag)
  r_com(i) = c"""

content = content.replace(old_loop, new_loop)

with open(filepath, 'w') as f:
    f.write(content)

print('mpadd3.f90: GOTO 100 eliminated (EXIT with flag)!')
