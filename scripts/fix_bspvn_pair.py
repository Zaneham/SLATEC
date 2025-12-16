#!/usr/bin/env python3
"""Fix GOTO in bspvn.f90 and dbspvn.f90 - B-spline evaluation"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/interpolation/bspline/bspvn.f90', 'SP'),
    ('C:/dev/slatec-modern/src/modern/interpolation/bspline/dbspvn.f90', 'DP'),
]

for filepath, prec in files:
    if not os.path.exists(filepath):
        print(f'{filepath}: NOT FOUND, skipping')
        continue

    with open(filepath, 'r') as f:
        content = f.read()

    # Add revision history
    old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
    new_rev = '''!   920501  Reformatted the REFERENCES section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.8 (IF construct)
  !           Original: de Boor (SIAM J. Numer. Anal. 14, 1977); Amos, SNLA'''

    content = content.replace(old_rev, new_rev)

    # The GOTO 100 skips to RETURN when Iwork>=Jhigh
    # The pattern: IF( Iwork>=Jhigh ) GOTO 100 followed by a DO loop
    # We can restructure this as an IF-ELSE
    if prec == 'SP':
        old_pattern = '''    IF( Indexx/=2 ) THEN
      Iwork = 1
      Vnikx(1) = 1._SP
      IF( Iwork>=Jhigh ) GOTO 100
    END IF
    DO
      !
      ipj = Ileft + Iwork
      Work(Iwork) = T(ipj) - X
      imjp1 = Ileft - Iwork + 1
      Work(K+Iwork) = X - T(imjp1)
      vmprev = 0._SP
      jp1 = Iwork + 1
      DO l = 1, Iwork
        jp1ml = jp1 - l
        vm = Vnikx(l)/(Work(l)+Work(K+jp1ml))
        Vnikx(l) = vm*Work(l) + vmprev
        vmprev = vm*Work(K+jp1ml)
      END DO
      Vnikx(jp1) = vmprev
      Iwork = jp1
      IF( Iwork>=Jhigh ) EXIT
    END DO
  END IF
  !
  100  RETURN'''
        new_pattern = '''    IF( Indexx/=2 ) THEN
      Iwork = 1
      Vnikx(1) = 1._SP
    END IF
    ! Only iterate if Iwork < Jhigh (was GOTO 100)
    DO WHILE( Iwork < Jhigh )
      ipj = Ileft + Iwork
      Work(Iwork) = T(ipj) - X
      imjp1 = Ileft - Iwork + 1
      Work(K+Iwork) = X - T(imjp1)
      vmprev = 0._SP
      jp1 = Iwork + 1
      DO l = 1, Iwork
        jp1ml = jp1 - l
        vm = Vnikx(l)/(Work(l)+Work(K+jp1ml))
        Vnikx(l) = vm*Work(l) + vmprev
        vmprev = vm*Work(K+jp1ml)
      END DO
      Vnikx(jp1) = vmprev
      Iwork = jp1
    END DO
  END IF
  ! (Label 100 removed - restructured with DO WHILE)'''
    else:  # DP
        old_pattern = '''    IF( Indexx/=2 ) THEN
      Iwork = 1
      Vnikx(1) = 1._DP
      IF( Iwork>=Jhigh ) GOTO 100
    END IF
    DO
      !
      ipj = Ileft + Iwork
      Work(Iwork) = T(ipj) - X
      imjp1 = Ileft - Iwork + 1
      Work(K+Iwork) = X - T(imjp1)
      vmprev = 0._DP
      jp1 = Iwork + 1
      DO l = 1, Iwork
        jp1ml = jp1 - l
        vm = Vnikx(l)/(Work(l)+Work(K+jp1ml))
        Vnikx(l) = vm*Work(l) + vmprev
        vmprev = vm*Work(K+jp1ml)
      END DO
      Vnikx(jp1) = vmprev
      Iwork = jp1
      IF( Iwork>=Jhigh ) EXIT
    END DO
  END IF
  !
  100  RETURN'''
        new_pattern = '''    IF( Indexx/=2 ) THEN
      Iwork = 1
      Vnikx(1) = 1._DP
    END IF
    ! Only iterate if Iwork < Jhigh (was GOTO 100)
    DO WHILE( Iwork < Jhigh )
      ipj = Ileft + Iwork
      Work(Iwork) = T(ipj) - X
      imjp1 = Ileft - Iwork + 1
      Work(K+Iwork) = X - T(imjp1)
      vmprev = 0._DP
      jp1 = Iwork + 1
      DO l = 1, Iwork
        jp1ml = jp1 - l
        vm = Vnikx(l)/(Work(l)+Work(K+jp1ml))
        Vnikx(l) = vm*Work(l) + vmprev
        vmprev = vm*Work(K+jp1ml)
      END DO
      Vnikx(jp1) = vmprev
      Iwork = jp1
    END DO
  END IF
  ! (Label 100 removed - restructured with DO WHILE)'''

    content = content.replace(old_pattern, new_pattern)

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 100 eliminated!')
