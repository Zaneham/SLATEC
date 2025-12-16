#!/usr/bin/env python
"""Fix convergence EXIT pattern in FNLIB special functions (Fullerton, LANL).

Pattern: GOTO 100 on convergence, ERROR STOP on non-convergence
Replace with: EXIT on convergence, check loop index after

Reference: ISO/IEC 1539-1:2018 S11.2.1 (EXIT statement)
Original: W. Fullerton, LANL (FNLIB), 1977
"""

import os
import re

# FNLIB files with single GOTO 100 (convergence exit pattern)
fnlib_files = [
    ('d9chu.f90', 'D9CHU', 'DP', 300),
    ('r9chu.f90', 'R9CHU', 'SP', 300),
    ('d9gmic.f90', 'D9GMIC', 'DP', 200),
    ('r9gmic.f90', 'R9GMIC', 'SP', 200),
    ('d9gmit.f90', 'D9GMIT', 'DP', 200),
    ('r9gmit.f90', 'R9GMIT', 'SP', 200),
    ('d9lgic.f90', 'D9LGIC', 'DP', 300),
    ('r9lgic.f90', 'R9LGIC', 'SP', 300),
    ('d9lgit.f90', 'D9LGIT', 'DP', 200),
    ('r9lgit.f90', 'R9LGIT', 'SP', 200),
]

base_dir = 'C:/dev/slatec-modern/src/modern/special_functions'

for fname, name, prec, max_iter in fnlib_files:
    filepath = os.path.join(base_dir, fname)

    if not os.path.exists(filepath):
        print(f'{fname}: NOT FOUND, skipping')
        continue

    with open(filepath, 'r') as f:
        content = f.read()

    # Skip if already fixed
    if 'Eliminated GOTO' in content:
        print(f'{fname}: already fixed, skipping')
        continue

    # Skip if no GOTO 100 pattern
    if 'GOTO 100' not in content:
        print(f'{fname}: no GOTO 100, skipping')
        continue

    # Add revision history - find last revision line
    rev_patterns = [
        '!   900720  Routine changed from user-callable to subsidiary.  (WRB)',
        '!   891214  Prologue converted to Version 4.0 format.  (BAB)',
    ]

    new_rev_added = False
    for old_rev in rev_patterns:
        if old_rev in content:
            new_rev = old_rev + """
  !   211001  Converted to free-form, added INTENT, ELEMENTAL.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT on convergence)
  !           Original: Fullerton, LANL (FNLIB)"""
            content = content.replace(old_rev, new_rev)
            new_rev_added = True
            break

    if not new_rev_added:
        print(f'{fname}: revision pattern not found, skipping')
        continue

    # Replace GOTO 100 with EXIT
    content = re.sub(
        r'\) GOTO 100',
        ') EXIT  ! Converged',
        content
    )

    # Find and fix the ERROR STOP / label 100 pattern
    # Pattern: ERROR STOP 'msg'\n  !\n  100  result = ...
    content = re.sub(
        r"(ERROR STOP '[^']+')([\s!]*\n\s*100\s+)",
        r'IF( i > ' + str(max_iter) + r' ) \1\n  ! (Label 100 removed - convergence handled by EXIT)\n  ',
        content
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (convergence EXIT)!')
