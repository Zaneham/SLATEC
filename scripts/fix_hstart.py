#!/usr/bin/env python3
"""Fix GOTO in hstart.f90 - DEPAC starting step size (single precision)"""

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ_eq/depac/hstart.f90'

with open(filepath, 'r') as f:
    content = f.read()

# Add revision history
old_rev = '!   910722  Updated AUTHOR section.  (ALS)'
new_rev = '''!   910722  Updated AUTHOR section.  (ALS)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT statement)
  !           Original: Watts, SNLA'''

content = content.replace(old_rev, new_rev)

# The GOTO 100 skips to end of loop - use EXIT instead
old_goto = '''    IF( k==lk ) GOTO 100'''
new_goto = '''    IF( k==lk ) EXIT  ! Last iteration complete'''

content = content.replace(old_goto, new_goto)

# Remove label 100
old_label = '''  100  ydpb = dfdxb + dfdub*fbnd'''
new_label = '''  ! (Label 100 removed - reached via EXIT or loop completion)
  ydpb = dfdxb + dfdub*fbnd'''

content = content.replace(old_label, new_label)

with open(filepath, 'w') as f:
    f.write(content)

print('hstart.f90: GOTO 100 eliminated!')
