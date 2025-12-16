#!/usr/bin/env python3
"""Fix GOTO in qmomo.f90 - QUADPACK modified Chebyshev moments (single precision)"""

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ/quadpack/qmomo.f90'

with open(filepath, 'r') as f:
    content = f.read()

# Add revision history
old_rev = '!   900328  Added TYPE section.  (WRB)'
new_rev = '''!   900328  Added TYPE section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.8 (IF construct)
  !           Original: Piessens & de Doncker, K.U. Leuven (QUADPACK)'''

content = content.replace(old_rev, new_rev)

# Replace GOTO 100 with IF construct (SP version uses _SP suffix)
old_pattern = '''      IF( Integr==2 ) GOTO 100
    END IF
    !
    !           COMPUTE RH USING A FORWARD RECURRENCE RELATION.
    !
    Rh(1) = -Rj(1)/betp1
    Rh(2) = -(rbet+rbet)/(betp2*betp2) - Rh(1)
    an = 2._SP
    anm1 = 1._SP
    im1 = 2
    DO i = 3, 25
      Rh(i) = -(an*(an-betp2)*Rh(im1)-an*Rj(im1)+anm1*Rj(i))/(anm1*(an+betp1))
      anm1 = an
      an = an + 1._SP
      im1 = i
    END DO
    DO i = 2, 25, 2
      Rh(i) = -Rh(i)
    END DO
  END IF
  100 CONTINUE'''

new_pattern = '''      IF( Integr/=2 ) THEN
        !
        !           COMPUTE RH USING A FORWARD RECURRENCE RELATION.
        !           (Skipped when Integr==2)
        !
        Rh(1) = -Rj(1)/betp1
        Rh(2) = -(rbet+rbet)/(betp2*betp2) - Rh(1)
        an = 2._SP
        anm1 = 1._SP
        im1 = 2
        DO i = 3, 25
          Rh(i) = -(an*(an-betp2)*Rh(im1)-an*Rj(im1)+anm1*Rj(i))/(anm1*(an+betp1))
          anm1 = an
          an = an + 1._SP
          im1 = i
        END DO
        DO i = 2, 25, 2
          Rh(i) = -Rh(i)
        END DO
      END IF
    END IF
  END IF
  ! (Label 100 removed - was just CONTINUE)'''

content = content.replace(old_pattern, new_pattern)

with open(filepath, 'w') as f:
    f.write(content)

print('qmomo.f90: GOTO 100 eliminated!')
