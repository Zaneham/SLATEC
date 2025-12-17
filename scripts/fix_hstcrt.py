#!/usr/bin/env python
"""Fix GOTO 100/200 in hstcrt.f90 - FISHPACK Helmholtz solver.

Pattern: SELECT CASE with GOTO to skip boundary processing
Replace with: IF block wrapping the processing

Reference: ISO/IEC 1539-1:2018 S11.1.9 (SELECT CASE)
Original: Adams, Swarztrauber, Sweet (NCAR)
"""

import os

filepath = 'C:/dev/slatec-modern/src/modern/diff_integ_eq/fishpack/hstcrt.f90'
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
old_rev = '!   920501  Reformatted the REFERENCES section.  (WRB)'
new_rev = old_rev + """
  !   211001  Converted to free-form.  (Mehdi Chinoune)
  !   251217  Eliminated GOTO 100/200 per MODERNISATION_GUIDE.md S1. (ZH)
  !           Ref: ISO/IEC 1539-1:2018 S11.1.9 (SELECT CASE)
  !           Original: Adams, Swarztrauber, Sweet (NCAR)"""

content = content.replace(old_rev, new_rev)

# Fix X boundary section (GOTO 100)
old_xbound = '''  SELECT CASE (mp)
    CASE (1)
      GOTO 100
    CASE (4,5)
      DO j = 1, N
        F(1,j) = F(1,j) + Bda(j)*twdelx
      END DO
      W(id2+1) = W(id2+1) + W(1)
    CASE DEFAULT
      DO j = 1, N
        F(1,j) = F(1,j) - Bda(j)*delxsq
      END DO
      W(id2+1) = W(id2+1) - W(1)
  END SELECT
  SELECT CASE (mp)
    CASE (1)
    CASE (3,4)
      DO j = 1, N
        F(M,j) = F(M,j) - Bdb(j)*twdelx
      END DO
      W(id3) = W(id3) + W(1)
    CASE DEFAULT
      DO j = 1, N
        F(M,j) = F(M,j) - Bdb(j)*delxsq
      END DO
      W(id3) = W(id3) - W(1)
  END SELECT
  !
  !     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
  !
  100 CONTINUE'''

new_xbound = '''  IF( mp/=1 ) THEN  ! Process X boundaries (was GOTO 100 skip)
    SELECT CASE (mp)
      CASE (4,5)
        DO j = 1, N
          F(1,j) = F(1,j) + Bda(j)*twdelx
        END DO
        W(id2+1) = W(id2+1) + W(1)
      CASE DEFAULT
        DO j = 1, N
          F(1,j) = F(1,j) - Bda(j)*delxsq
        END DO
        W(id2+1) = W(id2+1) - W(1)
    END SELECT
    SELECT CASE (mp)
      CASE (3,4)
        DO j = 1, N
          F(M,j) = F(M,j) - Bdb(j)*twdelx
        END DO
        W(id3) = W(id3) + W(1)
      CASE DEFAULT
        DO j = 1, N
          F(M,j) = F(M,j) - Bdb(j)*delxsq
        END DO
        W(id3) = W(id3) - W(1)
    END SELECT
  END IF
  !
  !     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
  !
  ! (Label 100 removed)'''

content = content.replace(old_xbound, new_xbound)

# Fix Y boundary section (GOTO 200)
old_ybound = '''  SELECT CASE (np)
    CASE (1)
      GOTO 200
    CASE (4,5)
      DO i = 1, M
        F(i,1) = F(i,1) + Bdc(i)*twdely
      END DO
    CASE DEFAULT
      DO i = 1, M
        F(i,1) = F(i,1) - Bdc(i)*twdysq
      END DO
  END SELECT
  SELECT CASE (np)
    CASE (1)
    CASE (3,4)
      DO i = 1, M
        F(i,N) = F(i,N) - Bdd(i)*twdely
      END DO
    CASE DEFAULT
      DO i = 1, M
        F(i,N) = F(i,N) - Bdd(i)*twdysq
      END DO
  END SELECT
  200 CONTINUE'''

new_ybound = '''  IF( np/=1 ) THEN  ! Process Y boundaries (was GOTO 200 skip)
    SELECT CASE (np)
      CASE (4,5)
        DO i = 1, M
          F(i,1) = F(i,1) + Bdc(i)*twdely
        END DO
      CASE DEFAULT
        DO i = 1, M
          F(i,1) = F(i,1) - Bdc(i)*twdysq
        END DO
    END SELECT
    SELECT CASE (np)
      CASE (3,4)
        DO i = 1, M
          F(i,N) = F(i,N) - Bdd(i)*twdely
        END DO
      CASE DEFAULT
        DO i = 1, M
          F(i,N) = F(i,N) - Bdd(i)*twdysq
        END DO
    END SELECT
  END IF
  ! (Label 200 removed)'''

content = content.replace(old_ybound, new_ybound)

with open(filepath, 'w') as f:
    f.write(content)

print(f'{fname}: GOTO 100/200 eliminated (IF block restructure)!')
