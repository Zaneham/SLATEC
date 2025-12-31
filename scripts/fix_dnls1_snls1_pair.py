#!/usr/bin/env python
"""Fix GOTOs in dnls1.inc and snls1.inc - MINPACK nonlinear least squares.

Pattern:
  - Early validation GOTO 200 -> inline validation with early exit
  - Label 100 CONTINUE -> outer_loop: DO with CYCLE
  - GOTO 100 (back to outer loop) -> CYCLE outer_loop
  - GOTO 200 (error/termination) -> EXIT outer_loop

Reference: ISO/IEC 1539-1:2018 S11.1.7.4.3 (EXIT/CYCLE)
Original: Hiebert, K. L. (SNLA) - MINPACK
"""

import os
import re

def fix_file(filepath, precision_suffix):
    """Fix GOTOs in a single file."""
    fname = os.path.basename(filepath)

    if not os.path.exists(filepath):
        print(f'{fname}: NOT FOUND')
        return False

    with open(filepath, 'r') as f:
        content = f.read()

    # Skip if already fixed
    if 'outer_loop: DO' in content or 'Eliminated GOTO' in content:
        print(f'{fname}: already fixed, skipping')
        return True

    # Count original GOTOs
    original_gotos = len(re.findall(r'GOTO\s+\d+', content))
    print(f'{fname}: found {original_gotos} GOTOs')

    # Add revision history
    old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
    new_rev = old_rev + """
  !   260101  Eliminated GOTO using outer_loop DO/EXIT/CYCLE.  (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.7.4.3 (EXIT/CYCLE)"""
    content = content.replace(old_rev, new_rev)

    # 1. Replace early validation GOTOs with validation block
    # These are before the outer loop starts

    # Pattern: multi-line validation check
    old_validation = f'''IF( Iopt<1 .OR. Iopt>3 .OR. N<=0 .OR. M<N .OR. Ldfjac<N .OR. Ftol<0._{precision_suffix} .OR. &
    Xtol<0._{precision_suffix} .OR. Gtol<0._{precision_suffix} .OR. Maxfev<=0 .OR. Factor<=0._{precision_suffix} ) GOTO 200
  IF( Iopt<3 .AND. Ldfjac<M ) GOTO 200'''

    new_validation = f'''validation: BLOCK
    LOGICAL :: invalid_input
    invalid_input = Iopt<1 .OR. Iopt>3 .OR. N<=0 .OR. M<N .OR. Ldfjac<N .OR. Ftol<0._{precision_suffix} .OR. &
      Xtol<0._{precision_suffix} .OR. Gtol<0._{precision_suffix} .OR. Maxfev<=0 .OR. Factor<=0._{precision_suffix}
    IF( Iopt<3 .AND. Ldfjac<M ) invalid_input = .TRUE.
    IF( invalid_input ) EXIT validation'''

    content = content.replace(old_validation, new_validation)

    # Pattern: Diag validation inside Mode==2 block
    old_diag_check = f'''IF( Mode==2 ) THEN
    DO j = 1, N
      IF( Diag(j)<=0._{precision_suffix} ) GOTO 200
    END DO
  END IF'''

    new_diag_check = f'''IF( Mode==2 ) THEN
      DO j = 1, N
        IF( Diag(j)<=0._{precision_suffix} ) EXIT validation
      END DO
    END IF'''

    content = content.replace(old_diag_check, new_diag_check)

    # Pattern: iflag check after function evaluation
    content = content.replace(
        'IF( iflag<0 ) GOTO 200\n  fnorm = NORM2(Fvec)',
        'IF( iflag<0 ) EXIT validation\n    fnorm = NORM2(Fvec)'
    )

    # Close the validation block and add outer_loop
    old_outer_start = '''!
  !     BEGINNING OF THE OUTER LOOP.
  !
  !
  !        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
  !
  100 CONTINUE'''

    new_outer_start = '''  END BLOCK validation
  !
  !     BEGINNING OF THE OUTER LOOP.
  !
  outer_loop: DO
    !
    !        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
    !'''

    content = content.replace(old_outer_start, new_outer_start)

    # 2. Replace GOTO 200 inside the loop with EXIT outer_loop
    content = re.sub(r'IF\( iflag<0 \) GOTO 200', 'IF( iflag<0 ) EXIT outer_loop', content)

    # 3. Replace GOTO 100 with CYCLE outer_loop
    content = content.replace('IF( ratio>=p0001 ) GOTO 100', 'IF( ratio>=p0001 ) CYCLE outer_loop')

    # 4. Close the outer loop and fix termination section
    old_termination = '''    END DO
  END IF
  !
  !     TERMINATION, EITHER NORMAL OR USER IMPOSED.
  !
  200 CONTINUE'''

    new_termination = '''      EXIT outer_loop
    END DO
  END IF
  EXIT outer_loop
  END DO outer_loop
  !
  !     TERMINATION, EITHER NORMAL OR USER IMPOSED.
  !'''

    content = content.replace(old_termination, new_termination)

    # Count remaining GOTOs
    remaining_gotos = len(re.findall(r'GOTO\s+\d+', content))

    with open(filepath, 'w') as f:
        f.write(content)

    print(f'{fname}: {original_gotos} -> {remaining_gotos} GOTOs (eliminated {original_gotos - remaining_gotos})')
    return remaining_gotos == 0


if __name__ == '__main__':
    base_path = 'C:/dev/slatec-modern/src/modern/approximation/minpack'

    success_d = fix_file(f'{base_path}/dnls1.inc', 'DP')
    success_s = fix_file(f'{base_path}/snls1.inc', 'SP')

    if success_d and success_s:
        print('\nBoth files fixed successfully!')
    else:
        print('\nSome files may need manual review.')
