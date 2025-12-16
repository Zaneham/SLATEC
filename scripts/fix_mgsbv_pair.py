#!/usr/bin/env python3
"""Fix GOTO in dmgsbv.f90 and mgsbv.f90 - BVSUP orthogonalization"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/bvsup/dmgsbv.f90', 'DP'),
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/bvsup/mgsbv.f90', 'SP'),
]

for filepath, prec in files:
    with open(filepath, 'r') as f:
        content = f.read()

    # Add revision history
    old_rev = '!   910722  Updated AUTHOR section.  (ALS)'
    new_rev = '''!   910722  Updated AUTHOR section.  (ALS)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 50 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.2.1 (EXIT with flag)
  !           Original: Watts, SNLA'''

    content = content.replace(old_rev, new_rev)

    # The GOTO 50 pattern: early exit when y <= eps*S(ix) -> sets Iflag=2, Wcnd=eps
    # We'll use a logical flag to handle this
    if prec == 'DP':
        # Add flag declaration (after existing INTEGER declarations)
        old_decl = '''INTEGER :: i, ip1, ix, iz, j, jk, jp, jq, jy, jz, k, kd, kj, kp, l, lix, lr, &
    m2, nivn, nmnr, nn, np1, nr, nrm1'''
        new_decl = '''INTEGER :: i, ip1, ix, iz, j, jk, jp, jq, jy, jz, k, kd, kj, kp, l, lix, lr, &
    m2, nivn, nmnr, nn, np1, nr, nrm1
  LOGICAL :: lin_dep  ! Flag for linear dependence detection'''
        content = content.replace(old_decl, new_decl)

        # Initialize flag after Wcnd = 1._DP
        old_init = '''Wcnd = 1._DP
    nivn = Niv'''
        new_init = '''Wcnd = 1._DP
    lin_dep = .FALSE.  ! Initialize linear dependence flag
    nivn = Niv'''
        content = content.replace(old_init, new_init)

        # Replace GOTO 50 with flag set and EXIT
        old_goto = '''          IF( y<=eps_com*S(ix) ) GOTO 50'''
        new_goto = '''          IF( y<=eps_com*S(ix) ) THEN
            lin_dep = .TRUE.  ! Linear dependence detected
            EXIT
          END IF'''
        content = content.replace(old_goto, new_goto)

        # After the main DO loop and before the Inhomo check, need to check flag
        # and handle the label 50 case
        old_label = '''    END DO
      !              *********************************************************
      !
      !                  TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
      !
      !        ......EXIT
      IF( Inhomo/=1 ) RETURN
      IF( (N>1) .AND. (S(np1)<1._SP) ) RETURN
      vnorm = NORM2(V(1:M))**2
      IF( S(np1)/=0._DP ) Wcnd = MIN(Wcnd,vnorm/S(np1))
      !        ......EXIT
      IF( vnorm>=eps_com*S(np1) ) RETURN
    END IF
    50  Iflag = 2
    Wcnd = eps_com'''
        new_label = '''    END DO
      !              *********************************************************
      !
      !              Check for linear dependence detected during loop
      IF( .NOT. lin_dep ) THEN
        !
        !                  TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
        !
        IF( Inhomo/=1 ) RETURN
        IF( (N>1) .AND. (S(np1)<1._SP) ) RETURN
        vnorm = NORM2(V(1:M))**2
        IF( S(np1)/=0._DP ) Wcnd = MIN(Wcnd,vnorm/S(np1))
        IF( vnorm>=eps_com*S(np1) ) RETURN
      END IF
    END IF
    ! (Label 50 removed - linear dependence sets flag)
    Iflag = 2
    Wcnd = eps_com'''
        content = content.replace(old_label, new_label)
    else:  # SP version
        # Add flag declaration
        old_decl = '''INTEGER :: i, ip1, ix, iz, j, jk, jp, jq, jy, jz, k, kd, kj, kp, l, lix, lr, &
    m2, nivn, nmnr, nn, np1, nr, nrm1'''
        new_decl = '''INTEGER :: i, ip1, ix, iz, j, jk, jp, jq, jy, jz, k, kd, kj, kp, l, lix, lr, &
    m2, nivn, nmnr, nn, np1, nr, nrm1
  LOGICAL :: lin_dep  ! Flag for linear dependence detection'''
        content = content.replace(old_decl, new_decl)

        # Initialize flag
        old_init = '''Wcnd = 1._SP
    nivn = Niv'''
        new_init = '''Wcnd = 1._SP
    lin_dep = .FALSE.  ! Initialize linear dependence flag
    nivn = Niv'''
        content = content.replace(old_init, new_init)

        # Replace GOTO 50
        old_goto = '''          IF( y<=eps_com*S(ix) ) GOTO 50'''
        new_goto = '''          IF( y<=eps_com*S(ix) ) THEN
            lin_dep = .TRUE.  ! Linear dependence detected
            EXIT
          END IF'''
        content = content.replace(old_goto, new_goto)

        # Replace label section
        old_label = '''    END DO
      !- *********************************************************************
      !
      !     TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
      !
      IF( Inhomo/=1 ) RETURN
      IF( (N>1) .AND. (S(np1)<1._SP) ) RETURN
      vnorm = NORM2(V(1:M))**2
      IF( S(np1)/=0. ) Wcnd = MIN(Wcnd,vnorm/S(np1))
      IF( vnorm>=eps_com*S(np1) ) RETURN
    END IF
    50  Iflag = 2
    Wcnd = eps_com'''
        new_label = '''    END DO
      !- *********************************************************************
      !
      !     Check for linear dependence detected during loop
      IF( .NOT. lin_dep ) THEN
        !     TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
        IF( Inhomo/=1 ) RETURN
        IF( (N>1) .AND. (S(np1)<1._SP) ) RETURN
        vnorm = NORM2(V(1:M))**2
        IF( S(np1)/=0. ) Wcnd = MIN(Wcnd,vnorm/S(np1))
        IF( vnorm>=eps_com*S(np1) ) RETURN
      END IF
    END IF
    ! (Label 50 removed - linear dependence sets flag)
    Iflag = 2
    Wcnd = eps_com'''
        content = content.replace(old_label, new_label)

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 50 eliminated (error exit -> flag)!')
