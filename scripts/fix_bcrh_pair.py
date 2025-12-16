#!/usr/bin/env python3
"""Fix GOTO in bcrh.f90 and bsrh.f90 - FISHPACK bisection root finder"""

import os

files = [
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/bcrh.f90', 'BCRH'),
    ('C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/bsrh.f90', 'BSRH'),
]

for filepath, fname_upper in files:
    with open(filepath, 'r') as f:
        content = f.read()

    # Add revision history
    old_rev = '!   900402  Added TYPE section.  (WRB)'
    new_rev = '''!   900402  Added TYPE section.  (WRB)
  !   211001  Converted to free-form, added INTENT, PURE.  (Mehdi Chinoune)
  !   251216  Eliminated GOTO 100 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S8.1.6.2 (DO WHILE construct)
  !           Original: Swarztrauber, NCAR (FISHPACK)'''

    content = content.replace(old_rev, new_rev)

    # Replace the backward GOTO loop with DO WHILE
    # The pattern is: label 100, then IF at end with GOTO 100
    old_loop = f'''  100  x = 0.5_SP*(xl+xr)
  IF( Sgn*F(x,Iz,C,A,Bh)<0 ) THEN
    xl = x
  ELSEIF( Sgn*F(x,Iz,C,A,Bh)==0 ) THEN
    {fname_upper} = 0.5_SP*(xl+xr)
    RETURN
  ELSE
    xr = x
  END IF
  dx = 0.5_SP*dx
  IF( dx>cnv_com ) GOTO 100
  {fname_upper} = 0.5_SP*(xl+xr)'''

    new_loop = f'''  ! Bisection iteration (was backward GOTO 100)
  DO WHILE( dx > cnv_com )
    x = 0.5_SP*(xl+xr)
    IF( Sgn*F(x,Iz,C,A,Bh) < 0 ) THEN
      xl = x
    ELSEIF( Sgn*F(x,Iz,C,A,Bh) == 0 ) THEN
      {fname_upper} = 0.5_SP*(xl+xr)
      RETURN
    ELSE
      xr = x
    END IF
    dx = 0.5_SP*dx
  END DO
  {fname_upper} = 0.5_SP*(xl+xr)'''

    content = content.replace(old_loop, new_loop)

    with open(filepath, 'w') as f:
        f.write(content)

    fname = os.path.basename(filepath)
    print(f'{fname}: GOTO 100 eliminated (bisection loop -> DO WHILE)!')
