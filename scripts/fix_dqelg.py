#!/usr/bin/env python3
"""Fix GOTO in dqelg.f90 - QUADPACK epsilon algorithm"""

import os

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ/quadpack/dqelg.f90'

with open(filepath, 'r') as f:
    content = f.read()

# Add revision history
old_rev = '''!* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)'''

new_rev = '''!* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT statement)
  !           Original: Piessens & de Doncker, K.U. Leuven (QUADPACK)'''

content = content.replace(old_rev, new_rev)

# Add converged flag variable
old_vars = '''  INTEGER :: i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num
  REAL(DP) :: delta1, delta2, delta3, epmach, epsinf, error, err1, err2, err3, &
    e0, e1, e1abs, e2, e3, oflow, res, ss, tol1, tol2, tol3'''

new_vars = '''  INTEGER :: i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num
  REAL(DP) :: delta1, delta2, delta3, epmach, epsinf, error, err1, err2, err3, &
    e0, e1, e1abs, e2, e3, oflow, res, ss, tol1, tol2, tol3
  LOGICAL :: converged  ! Flag for early convergence exit'''

content = content.replace(old_vars, new_vars)

# Initialize flag
old_init = '''!* FIRST EXECUTABLE STATEMENT  DQELG
  epmach = eps_dp'''

new_init = '''!* FIRST EXECUTABLE STATEMENT  DQELG
  converged = .FALSE.
  epmach = eps_dp'''

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

# Replace label 100 section
old_label = '''  END IF
  100  Abserr = MAX(Abserr,5._DP*epmach*ABS(Result))'''

new_label = '''  END IF
  ! Final error bound (replaces label 100)
  Abserr = MAX(Abserr,5._DP*epmach*ABS(Result))'''

content = content.replace(old_label, new_label)

with open(filepath, 'w') as f:
    f.write(content)

print('dqelg.f90: GOTO 100 eliminated!')
