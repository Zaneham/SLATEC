# Bugs Found During SLATEC Modernisation

This document records bugs discovered and fixed during the GOTO elimination and modernisation effort.

## 1. DBESK/BESK Numerical Drift (K Bessel Function)

**Discovery Date:** 2025-12-20
**Severity:** High (numerical accuracy)
**Affected Routines:** `dbsknu.f90`, `besknu.f90`
**Reference:** Temme, N.M., J. Comp. Physics 19, 1975, pp. 324-337

### Symptoms
- K Bessel function values showed relative error of approximately 4×10⁻⁵
- Example: K₄(1) computed as 44.2324 instead of correct value 44.2341

### Root Cause
During GOTO elimination, the control flow for the power series computation (X ≤ 2) was incorrectly structured. The original SLATEC code used:

```fortran
! After series computation
GOTO 20  ! Skip to forward recursion
...
! Coefficient/asymptotic/Miller code
...
20  ck = (dnu+dnu+2._DP)/X  ! Forward recursion
```

The modernised code removed the GOTO but failed to skip the coefficient/asymptotic/Miller section. This caused the already-computed series values (s1, s2) to be overwritten by the coefficient calculation code.

### Fix
Added a `series_done` logical flag:

```fortran
LOGICAL :: series_done
...
series_done = .FALSE.
...
! After series computation
series_done = .TRUE.
...
! Coefficient/asymptotic/Miller section
IF( .NOT. series_done ) THEN
  ! ... coefficient/asymptotic/Miller code ...
END IF
```

### Verification
After the fix, the relative error dropped from ~4×10⁻⁵ to ~1×10⁻¹⁴ (machine precision).

---

## 2. DCOV/SCOV Missing BLOCK Statement

**Discovery Date:** 2025-12-20
**Severity:** Medium (compilation failure)
**Affected Routines:** `dcov.inc`, `scov.inc` (MINPACK covariance routines)

### Symptoms
Compilation error: `Name 'process' in EXIT statement at (1) is unknown`

### Root Cause
During GOTO elimination, a named BLOCK construct was introduced to replace a labelled error exit, but the opening `process: BLOCK` statement was missing whilst `END BLOCK process` remained.

### Fix
Added the missing `process: BLOCK` statement before the main IF block.

---

## 3. DPSORT/HPSORT/IPSORT/SPSORT Unclosed DO Loop

**Discovery Date:** 2025-12-20
**Severity:** Medium (compilation failure)
**Affected Routines:** Data handling sort routines

### Symptoms
Compilation error: `Syntax error in END DO statement`

### Root Cause
A nested DO loop was missing its `END DO` statement. The code had:

```fortran
main_sort: DO
  ...
  DO            ! Inner loop
    ...
  END DO main_sort  ! Wrong! This closes outer loop, not inner
```

### Fix
Added the missing `END DO` for the inner loop:

```fortran
main_sort: DO
  ...
  DO            ! Inner loop
    ...
  END DO        ! Closes inner loop
END DO main_sort  ! Now correctly closes outer loop
```

---

## 4. PCHID/PCHIA Undeclared Variable

**Discovery Date:** 2025-12-20
**Severity:** Low (compilation failure)
**Affected Routines:** PCHIP interpolation routines

### Symptoms
Compilation error: `Symbol 'ierr' at (1) has no IMPLICIT type`

### Root Cause
During GOTO elimination, error return code assignments (`Ierr = -3`) were left in place but the subsequent line was changed to `ERROR STOP`. The variable assignment became redundant and the variable was never declared.

### Fix
Removed the redundant `Ierr` assignments since the `ERROR STOP` statements make them unreachable.

---

## Pre-existing Issues (Not Introduced by Modernisation)

The following issues were discovered but predate the modernisation effort:

1. **PCHSP.inc**: Contains GOTOs (20, 50) that were never converted
2. **PCHFE.inc**: Contains multiple unconverted GOTOs (100, 200, 300, 320, 350, 400, 600)
3. **Various PCHIP routines**: Unused labels (100, 200) generating warnings

These require additional GOTO elimination work in the interpolation module.

---

## Build System Issues

The following structural issues were addressed to enable compilation:

1. **Include files as .f90**: Files intended for inclusion via Fortran `include` statements were named `.f90`, causing fpm to attempt standalone compilation. Renamed to `.inc`.

2. **Module dependencies**: Several modules (blas, lapack, linpack) needed to be extracted to separate files to ensure correct compilation order.

3. **Missing USE statements**: Some module files lacked explicit `USE` statements for their dependencies, preventing fpm from determining the correct build order.
