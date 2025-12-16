#!/usr/bin/env python3
"""Fix GOTO in qelg.f90 - QUADPACK epsilon algorithm (single precision)"""

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ/quadpack/qelg.f90'

with open(filepath, 'r') as f:
    content = f.read()

# Add revision history - find the pattern
old_rev = '!   900328  Added TYPE section.  (WRB)'
new_rev = '''!   900328  Added TYPE section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT statement)
  !           Original: Piessens & de Doncker, K.U. Leuven (QUADPACK)'''

content = content.replace(old_rev, new_rev)

# Add converged flag - for SP version uses SP kind
if 'LOGICAL :: converged' not in content:
    old_vars = 'REAL(SP) :: delta1, delta2, delta3'
    new_vars = 'LOGICAL :: converged  ! Flag for early convergence exit\n  REAL(SP) :: delta1, delta2, delta3'
    content = content.replace(old_vars, new_vars)

# Initialize flag
if 'converged = .FALSE.' not in content:
    old_init = '!* FIRST EXECUTABLE STATEMENT  QELG\n  epmach'
    new_init = '!* FIRST EXECUTABLE STATEMENT  QELG\n  converged = .FALSE.\n  epmach'
    content = content.replace(old_init, new_init)

# Replace GOTO 100 with EXIT and flag
old_goto = '''        Result = res
        Abserr = err2 + err3
        !- **JUMP OUT OF DO-LOOP
        GOTO 100'''

new_goto = '''        Result = res
        Abserr = err2 + err3
        ! Early convergence - e0, e1, e2 equal within machine accuracy
        converged = .TRUE.
        EXIT'''

content = content.replace(old_goto, new_goto)

# Replace label 100
old_label = '  100  Abserr = MAX(Abserr,5._SP*epmach*ABS(Result))'
new_label = '  ! Final error bound (replaces label 100)\n  Abserr = MAX(Abserr,5._SP*epmach*ABS(Result))'
content = content.replace(old_label, new_label)

with open(filepath, 'w') as f:
    f.write(content)

print('qelg.f90: GOTO 100 eliminated!')
