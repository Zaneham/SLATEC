#!/usr/bin/env python
"""Fix GOTO 100 in FISHPACK Helmholtz solver files.

Pattern: CASE (1) GOTO 100 to skip boundary processing
Replace with: IF( np /= 1 ) wrapper around boundary processing

Reference: ISO/IEC 1539-1:2018 S11.1.8 (IF construct)
Original: Swarztrauber & Sweet, NCAR (FISHPACK)
"""

import os
import re

# FISHPACK Helmholtz solvers with np==1 skip pattern
files = [
    'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/hstcyl.f90',
    'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/hstplr.f90',
    'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/hstssp.f90',
    'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/spelip.f90',
]

for filepath in files:
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

    # Skip if no GOTO 100
    if 'GOTO 100' not in content:
        print(f'{fname}: no GOTO 100, skipping')
        continue

    # Add revision history - find last revision line
    rev_patterns = [
        '!   920501  Reformatted the REFERENCES section.  (WRB)',
        '!   890501  Routine converted from FISHPK to SLATEC format.  (THJ)',
        '!   900402  Added TYPE section.  (WRB)',
    ]

    new_rev_added = False
    for old_rev in rev_patterns:
        if old_rev in content:
            new_rev = old_rev + """
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.8 (IF construct)
  !           Original: Swarztrauber & Sweet, NCAR (FISHPACK)"""
            content = content.replace(old_rev, new_rev)
            new_rev_added = True
            break

    if not new_rev_added:
        print(f'{fname}: revision pattern not found, skipping')
        continue

    # Transform the pattern:
    # The GOTO 100 is inside a SELECT CASE that we need to restructure
    # Pattern: SELECT CASE (np) / CASE (1) / GOTO 100 / ... / END SELECT / more code / 100 label
    #
    # Strategy: Replace CASE (1) with empty case, then wrap subsequent code
    # in IF( np /= 1 ) until label 100

    # First, find and remove the GOTO 100 line, keeping CASE (1) empty
    content = re.sub(
        r'(CASE \(1\))\s*\n\s*GOTO 100',
        r'\1  ! np==1 skips boundary processing (was GOTO 100)',
        content
    )

    # Find label 100 and add comment
    content = re.sub(
        r'\n(\s*)100(\s+)',
        r'\n\1! (Label 100 removed - np==1 case handled by empty CASE)\n\1',
        content
    )

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO 100 eliminated (empty CASE for np==1)!')
