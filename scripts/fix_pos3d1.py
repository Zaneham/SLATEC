#!/usr/bin/env python
"""Fix GOTO 100/200 in pos3d1.f90 - FISHPACK 3D Poisson solver helper.

Pattern: SELECT CASE with GOTO to skip alternate initialization
Replace with: IF-ELSE block for alternate paths

Reference: ISO/IEC 1539-1:2018 S11.1.9 (SELECT CASE)
Original: Adams, Swarztrauber, Sweet (NCAR)
"""

import os

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/pos3d1.f90'
fname = os.path.basename(filepath)

if not os.path.exists(filepath):
    print(f'{fname}: NOT FOUND')
    exit(1)

with open(filepath, 'r') as f:
    content = f.read()

# Skip if already fixed
if 'Eliminated GOTO' in content:
    print(f'{fname}: already fixed, skipping')
    exit(0)

# Add revision history
old_rev = '!   900402  Added TYPE section.  (WRB)'
new_rev = old_rev + """
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.9 (SELECT CASE)
  !           Original: Adams, Swarztrauber, Sweet (NCAR)"""

content = content.replace(old_rev, new_rev)

# Fix X initialization section (GOTO 100)
old_xinit = '''  SELECT CASE (Lp)
    CASE (1)
      Xrt(1) = 0._SP
      Xrt(lr) = -4._SP*C1
      DO i = 3, lr, 2
        Xrt(i-1) = -4._SP*C1*(SIN((i-1)*dx))**2
        Xrt(i) = Xrt(i-1)
      END DO
      CALL RFFTI(lr,Wx)
      GOTO 100
    CASE (2)
      di = 0._SP
    CASE (4)
      di = 1._SP
    CASE DEFAULT
      di = 0.5
      scalx = 2._SP*scalx
  END SELECT
  DO i = 1, lr
    Xrt(i) = -4._SP*C1*(SIN((i-di)*dx))**2
  END DO
  scalx = 2._SP*scalx
  SELECT CASE (Lp)
    CASE (1)
    CASE (3)
      CALL SINQI(lr,Wx)
    CASE (4)
      CALL COSTI(lr,Wx)
    CASE (5)
      CALL COSQI(lr,Wx)
    CASE DEFAULT
      CALL SINTI(lr,Wx)
  END SELECT
  100  mrdel'''

new_xinit = '''  IF( Lp==1 ) THEN  ! Special X initialization (was GOTO 100 skip)
    Xrt(1) = 0._SP
    Xrt(lr) = -4._SP*C1
    DO i = 3, lr, 2
      Xrt(i-1) = -4._SP*C1*(SIN((i-1)*dx))**2
      Xrt(i) = Xrt(i-1)
    END DO
    CALL RFFTI(lr,Wx)
  ELSE
    SELECT CASE (Lp)
      CASE (2)
        di = 0._SP
      CASE (4)
        di = 1._SP
      CASE DEFAULT
        di = 0.5
        scalx = 2._SP*scalx
    END SELECT
    DO i = 1, lr
      Xrt(i) = -4._SP*C1*(SIN((i-di)*dx))**2
    END DO
    scalx = 2._SP*scalx
    SELECT CASE (Lp)
      CASE (3)
        CALL SINQI(lr,Wx)
      CASE (4)
        CALL COSTI(lr,Wx)
      CASE (5)
        CALL COSQI(lr,Wx)
      CASE DEFAULT
        CALL SINTI(lr,Wx)
    END SELECT
  END IF
  ! (Label 100 removed)
  mrdel'''

content = content.replace(old_xinit, new_xinit)

# Fix Y initialization section (GOTO 200)
old_yinit = '''  SELECT CASE (Mp)
    CASE (1)
      Yrt(1) = 0._SP
      Yrt(mr) = -4._SP*C2
      DO j = 3, mr, 2
        Yrt(j-1) = -4._SP*C2*(SIN((j-1)*dy))**2
        Yrt(j) = Yrt(j-1)
      END DO
      CALL RFFTI(mr,Wy)
      GOTO 200
    CASE (2)
      dj = 0._SP
    CASE (4)
      dj = 1._SP
    CASE DEFAULT
      dj = 0.5
      scaly = 2._SP*scaly
  END SELECT
  DO j = 1, mr
    Yrt(j) = -4._SP*C2*(SIN((j-dj)*dy))**2
  END DO
  scaly = 2._SP*scaly
  SELECT CASE (Mp)
    CASE (1)
    CASE (3)
      CALL SINQI(mr,Wy)
    CASE (4)
      CALL COSTI(mr,Wy)
    CASE (5)
      CALL COSQI(mr,Wy)
    CASE DEFAULT
      CALL SINTI(mr,Wy)
  END SELECT
  200  nrdel'''

new_yinit = '''  IF( Mp==1 ) THEN  ! Special Y initialization (was GOTO 200 skip)
    Yrt(1) = 0._SP
    Yrt(mr) = -4._SP*C2
    DO j = 3, mr, 2
      Yrt(j-1) = -4._SP*C2*(SIN((j-1)*dy))**2
      Yrt(j) = Yrt(j-1)
    END DO
    CALL RFFTI(mr,Wy)
  ELSE
    SELECT CASE (Mp)
      CASE (2)
        dj = 0._SP
      CASE (4)
        dj = 1._SP
      CASE DEFAULT
        dj = 0.5
        scaly = 2._SP*scaly
    END SELECT
    DO j = 1, mr
      Yrt(j) = -4._SP*C2*(SIN((j-dj)*dy))**2
    END DO
    scaly = 2._SP*scaly
    SELECT CASE (Mp)
      CASE (3)
        CALL SINQI(mr,Wy)
      CASE (4)
        CALL COSTI(mr,Wy)
      CASE (5)
        CALL COSQI(mr,Wy)
      CASE DEFAULT
        CALL SINTI(mr,Wy)
    END SELECT
  END IF
  ! (Label 200 removed)
  nrdel'''

content = content.replace(old_yinit, new_yinit)

with open(filepath, 'w') as f:
    f.write(content)

print(f'{fname}: GOTO 100/200 eliminated (IF-ELSE restructure)!')
