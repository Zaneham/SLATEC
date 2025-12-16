#!/usr/bin/env python
"""Fix GOTO 100 in iswap.f90 - BLAS vector interchange.

Pattern: Skip forward to unrolled loop after cleanup
Replace with: Inline the unrolled loop in the appropriate branch

Reference: ISO/IEC 1539-1:2018 S11.1.8 (IF construct restructuring)
Original: Lawson, Hanson, Kincaid, Krogh (BLAS)
"""

import os

filepath = 'C:/dev/slatec-modern/src/modern/linear/iswap.f90'

with open(filepath, 'r') as f:
    content = f.read()

# Skip if already fixed
if 'Eliminated GOTO' in content:
    print('iswap.f90: already fixed')
    exit(0)

# Add revision history
old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
new_rev = """!   920501  Reformatted the REFERENCES section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.8 (IF construct)
  !           Original: Lawson, Hanson, Kincaid, Krogh (BLAS)"""

content = content.replace(old_rev, new_rev)

# The fix: Remove GOTO 100 and inline the loop
# Original structure has GOTO 100 at line 73, then ELSE branch, then END IF,
# then unequal increment code, then label 100 with main loop

# Replace the GOTO 100 and END IF with inline unrolled loop
old_pattern = """      END IF
      GOTO 100
    ELSE"""

new_pattern = """      END IF
      ! Main unrolled loop (was GOTO 100)
      mp1 = m + 1
      DO i = mp1, N, 3
        itemp1 = Ix(i)
        itemp2 = Ix(i+1)
        itemp3 = Ix(i+2)
        Ix(i) = Iy(i)
        Ix(i+1) = Iy(i+1)
        Ix(i+2) = Iy(i+2)
        Iy(i) = itemp1
        Iy(i+1) = itemp2
        Iy(i+2) = itemp3
      END DO
      RETURN
    ELSE"""

content = content.replace(old_pattern, new_pattern)

# Remove the label 100 and its duplicate loop (it's now inlined)
old_label = """  RETURN
  100  mp1 = m + 1
  DO i = mp1, N, 3
    itemp1 = Ix(i)
    itemp2 = Ix(i+1)
    itemp3 = Ix(i+2)
    Ix(i) = Iy(i)
    Ix(i+1) = Iy(i+1)
    Ix(i+2) = Iy(i+2)
    Iy(i) = itemp1
    Iy(i+1) = itemp2
    Iy(i+2) = itemp3
  END DO

  RETURN"""

new_label = """  ! (Label 100 removed - unrolled loop inlined above)
  RETURN"""

content = content.replace(old_label, new_label)

with open(filepath, 'w') as f:
    f.write(content)

print('iswap.f90: GOTO 100 eliminated (inlined unrolled loop)!')
