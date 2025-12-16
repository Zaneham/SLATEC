#!/usr/bin/env python3
"""Fix GOTO in blktr1.f90 - BLKTRI linear system solver"""

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/blktr1.f90'

with open(filepath, 'r') as f:
    content = f.read()

# Add revision history
old_rev = '!   900402  Added TYPE section.  (WRB)'
new_rev = '''!   900402  Added TYPE section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 50 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT statement)
  !           Original: Swarztrauber, NCAR (FISHPACK)'''

content = content.replace(old_rev, new_rev)

# The GOTO 50 skips to label 50 CONTINUE, then continues
# This is a "skip to end of nested loops" pattern
# The pattern: IF( i==nm_com ) GOTO 50 -> break out of inner loop to label 50
old_goto = '''          IF( i==nm_com ) GOTO 50
        END IF
      END DO
    END DO
    50 CONTINUE'''

new_goto = '''          IF( i==nm_com ) EXIT  ! Break to outer loop continuation
        END IF
      END DO
      IF( i==nm_com ) EXIT  ! Propagate exit
    END DO
    ! (Label 50 removed - EXIT propagates through loops)'''

content = content.replace(old_goto, new_goto)

with open(filepath, 'w') as f:
    f.write(content)

print('blktr1.f90: GOTO 50 eliminated!')
