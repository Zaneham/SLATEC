#!/usr/bin/env python3
"""Fix GOTO in lsei.f90 - LSEI constrained least squares (single precision)"""

filepath = 'C:/dev/slatec-modern/src/modern/approximation/lsei.f90'

with open(filepath, 'r') as f:
    content = f.read()

# Add revision history
old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
new_rev = '''!   920501  Reformatted the REFERENCES section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.3 (RETURN statement)
  !           Original: Hanson & Haskell, SNLA'''

content = content.replace(old_rev, new_rev)

# Replace GOTO 100 with RETURN
old_goto = '''          Mode = Mode + 2
          GOTO 100'''
new_goto = '''          Mode = Mode + 2
          RETURN  ! Early exit on constraint violation'''

content = content.replace(old_goto, new_goto)

# Remove label 100 CONTINUE
old_label = '  100 CONTINUE\n  !'
new_label = '  ! (Label 100 removed - was just CONTINUE before END)\n  !'
content = content.replace(old_label, new_label)

with open(filepath, 'w') as f:
    f.write(content)

print('lsei.f90: GOTO 100 eliminated!')
