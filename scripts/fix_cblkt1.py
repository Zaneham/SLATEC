#!/usr/bin/env python
"""Fix GOTO 50 in cblkt1.f90 - FISHPACK complex block tridiagonal solver.

Pattern: Nested loop EXIT via GOTO 50
Replace with: Named loops with EXIT

Reference: ISO/IEC 1539-1:2018 S11.2.1 (EXIT with construct name)
Original: Swarztrauber & Sweet, NCAR (FISHPACK)
"""

import os

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/cblkt1.f90'

with open(filepath, 'r') as f:
    content = f.read()

# Skip if already fixed
if 'Eliminated GOTO' in content:
    print('cblkt1.f90: already fixed')
    exit(0)

# Add revision history
old_rev = '!   891214  Prologue converted to Version 4.0 format.  (BAB)'
new_rev = """!   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 50 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (named loop EXIT)
  !           Original: Swarztrauber & Sweet, NCAR (FISHPACK)"""

content = content.replace(old_rev, new_rev)

# Name the outer loop
old_outer = '    DO ll = 2, k_com'
new_outer = '    reduction_loop: DO ll = 2, k_com'
content = content.replace(old_outer, new_outer)

# Replace GOTO 50 with EXIT reduction_loop
old_goto = 'IF( i==nm_com ) GOTO 50'
new_goto = 'IF( i==nm_com ) EXIT reduction_loop  ! Exit both loops (was GOTO 50)'
content = content.replace(old_goto, new_goto)

# Update corresponding END DO
old_enddo = """      END DO
    END DO
    50 CONTINUE"""
new_enddo = """      END DO
    END DO reduction_loop
    ! (Label 50 removed - exit handled by named loop)"""
content = content.replace(old_enddo, new_enddo)

with open(filepath, 'w') as f:
    f.write(content)

print('cblkt1.f90: GOTO 50 eliminated (named loop EXIT)!')
