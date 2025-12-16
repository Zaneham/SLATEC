#!/usr/bin/env python3
"""Fix X-array GOTO pattern in PCHIP files - batch process"""

import os
import re

# PCHIP files with X(i)<=X(i-1) GOTO pattern (excluding already fixed dpchim/pchim)
pchip_dir = 'C:/dev/slatec-modern/src/modern/interpolation/pchip'

# Process all .f90 files in PCHIP directory
for fname in os.listdir(pchip_dir):
    if not fname.endswith('.f90'):
        continue

    filepath = os.path.join(pchip_dir, fname)
    with open(filepath, 'r') as f:
        content = f.read()

    # Skip if no X-array GOTO pattern
    if 'X(i)<=X(i-1) ) GOTO' not in content:
        continue

    # Skip already fixed files (have the revision comment)
    if 'Eliminated GOTO' in content:
        continue

    # Find the label number for this file
    match = re.search(r'X\(i\)<=X\(i-1\) \) GOTO (\d+)', content)
    if not match:
        continue

    label = match.group(1)
    name = fname.replace('.f90', '').upper()

    # Determine precision
    prec = 'DP' if fname.startswith('d') else 'SP'

    # Add revision history - find the last revision line
    rev_patterns = [
        '!   920429  Revised format and order of references.  (WRB,FNF)',
        '!   920501  Reformatted the REFERENCES section.  (WRB)',
        '!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)',
    ]

    for old_rev in rev_patterns:
        if old_rev in content:
            new_rev = old_rev + """
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO {label} per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT with error flag)
  !           Original: Fritsch, LLNL (PCHIP)""".format(label=label)
            content = content.replace(old_rev, new_rev)
            break

    # Replace the GOTO pattern - handles various indentation
    old_pattern = re.compile(
        r'(\s+)DO i = 2, N\n\s+IF\( X\(i\)<=X\(i-1\) \) GOTO ' + label + r'\n\s+END DO'
    )

    def replacement(m):
        indent = m.group(1)
        return f"""{indent}DO i = 2, N
{indent}  IF( X(i)<=X(i-1) ) EXIT
{indent}END DO
{indent}IF( i <= N ) THEN
{indent}  ! X-array not strictly increasing (was GOTO {label})
{indent}  Ierr = -3
{indent}  ERROR STOP '{name} : X-ARRAY NOT STRICTLY INCREASING'
{indent}END IF"""

    content = old_pattern.sub(replacement, content)

    # Remove the label and its error handling
    # Pattern: label number followed by error code and ERROR STOP
    label_pattern = re.compile(
        r'\n\s*!\s*\n\s*!\s+X-ARRAY NOT STRICTLY INCREASING\.\s*\n\s*' + label +
        r'\s+Ierr = -3\s*\n\s*ERROR STOP [^\n]+'
    )
    content = label_pattern.sub('\n  ! (Label ' + label + ' removed - error handled inline after loop)', content)

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: GOTO {label} eliminated!')
