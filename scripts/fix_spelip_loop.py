#!/usr/bin/env python
"""Fix GOTO 100 in spelip.f90 - deferred correction loop.

Pattern: Backward GOTO for iteration (2nd order, then 4th order)
Replace with: DO loop with EXIT

Reference: ISO/IEC 1539-1:2018 S11.1.4.2 (DO construct)
Original: Swarztrauber (NCAR)
"""

import os

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/spelip.f90'
fname = os.path.basename(filepath)

if not os.path.exists(filepath):
    print(f'{fname}: NOT FOUND')
    exit(1)

with open(filepath, 'r') as f:
    content = f.read()

# Skip if already fixed properly
if 'Deferred correction loop' in content:
    print(f'{fname}: already fixed, skipping')
    exit(0)

# Update revision history to note second fix
old_comment = '! (Label 100 removed - np==1 case handled by empty CASE)'
new_comment = '! Deferred correction loop (was GOTO 100 target)'

content = content.replace(old_comment, new_comment)

# Replace CONTINUE with DO
content = content.replace(
    '! Deferred correction loop (was GOTO 100 target)\n  CONTINUE',
    '! Deferred correction loop (was GOTO 100 target)\n  DO'
)

# Replace IF( iord==2 ) RETURN with EXIT
content = content.replace(
    'IF( iord==2 ) RETURN\n  iord = 2',
    'IF( iord==2 ) EXIT  ! Done with correction loop\n    iord = 2'
)

# Replace GOTO 100 with END DO
content = content.replace(
    'CALL DEFER(COFX,COFY,Idmn,Usol,Grhs)\n  GOTO 100\n  !\nEND SUBROUTINE',
    'CALL DEFER(COFX,COFY,Idmn,Usol,Grhs)\n  END DO  ! Deferred correction loop\n  !\nEND SUBROUTINE'
)

with open(filepath, 'w') as f:
    f.write(content)

print(f'{fname}: GOTO 100 eliminated (DO loop)!')
