# SLATEC Modernisation Guide

## Document Control
- **Version**: 0.1
- **Date**: 2025-12-16
- **Author**: Zane Hambly
- **Based On**: Mehdi Chinoune's SLATEC modernisation (2021)

## Primary References

### Original SLATEC Documentation
1. **SLATEC Guide** (July 1993) - Fong, Jefferson, Suyehiro, Walton
   - Location: `docs/guide`
   - Defines FORTRAN 77 coding standards used in original library
   - Section 6: "Coding Guidelines - General Requirements for SLATEC"
   - Section 7: "Source Code Format"

2. **SLATEC Table of Contents**
   - Location: `docs/toc`
   - GAMS classification of all routines

3. **Individual Routine Documentation**
   - Each routine contains prologue comments citing original authors and references
   - Format defined in SLATEC Guide Section 8: "Prologue Format for Subprograms"

### Academic References (from SLATEC Guide)
- **[SL82]** Vandevender & Haskell, SIGNUM Newsletter 17(3), Sept 1982, pp. 16-21.
- **[JK83]** Jones & Kahaner, Software P&E 13(3), Mar 1983, pp. 251-257.
- **[BHK85]** Boisvert et al, ACM TOMS 11(4), Dec 1985, pp. 313-355 (GAMS).
- **[QUADPACK]** Piessens et al, Springer-Verlag, 1983.
- **[PCHIP84]** Fritsch & Butland, SIAM J. Sci. Stat. Comp. 5(2), 1984.
- **[FFT82]** Swarztrauber, in Parallel Computations, Academic Press, 1982.

### Modern Fortran Standards
4. **ISO/IEC 1539-1:2018** - Fortran 2018 Standard
   - Primary reference for modern language features

5. **ISO/IEC 1539-1:2010** - Fortran 2008 Standard
   - Introduced BLOCK construct, DO CONCURRENT, etc.

6. **ISO/IEC 1539-1:2004** - Fortran 2003 Standard
   - Introduced: allocatable components, INTENT for pointers, ASSOCIATE

7. **Metcalf, Reid, Cohen** - "Modern Fortran Explained" (2018)
   - Authoritative guide to Fortran 2008/2018 features

### Prior Modernisation Work
8. **Mehdi Chinoune** - https://github.com/MehdiChinoune/SLATEC (2021)
   - Converted to free-form source
   - Added INTENT attributes
   - Added KIND parameters (SP, DP for single/double precision)
   - Marked PURE/ELEMENTAL where applicable
   - Organised into modules by category

---

## Modernisation Transformations

### 1. GOTO Elimination

**Rationale**: GOTO statements create "spaghetti code" that is difficult to read, maintain, and verify. Modern Fortran provides structured alternatives.

**Reference**:
- SLATEC Guide Section 6 states routines should be "modular in structure" and "composed of reasonably small subprograms which in turn are made up of easily understandable blocks."
- ISO/IEC 1539-1:2018 §11.2 provides EXIT, CYCLE for loop control
- ISO/IEC 1539-1:2008 §8.1.6 provides BLOCK construct for scoped exits

**Transformation Patterns**:

| Original Pattern | Modern Replacement | Reference |
|------------------|-------------------|-----------|
| `GO TO label` at end of IF | Remove, restructure with IF-THEN-ELSE | F2018 §11.1.8 |
| `GO TO label` for error exit | `RETURN` with error code set | F2018 §11.2.3 |
| Computed `GO TO (l1,l2,l3), I` | `SELECT CASE (I)` | F2018 §11.1.9 |
| `GO TO` inside DO loop | `EXIT` or `CYCLE` with named loops | F2018 §11.2.1 |
| Arithmetic IF | `IF-THEN-ELSE IF-ELSE` | F2018 §11.1.8 |

**Example** (from `diff_integ/davint.f90`):
```fortran
! BEFORE (F77 style):
DO i = 2, N
  IF( X(i)<=X(i-1) ) GOTO 50  ! Jump to error handler
END DO
...
50 Ierr = 4
   ERROR STOP 'DAVINT: Abscissas not increasing'

! AFTER (F2018 style):
DO i = 2, N
  IF( X(i)<=X(i-1) ) THEN
    Ierr = 4
    ERROR STOP 'DAVINT: Abscissas not increasing'
    RETURN
  END IF
END DO
```

---

### 2. Assumed-Size to Assumed-Shape Arrays

**Rationale**: Assumed-size arrays `(*)` don't carry bounds information, preventing runtime bounds checking and requiring explicit size passing.

**Reference**:
- ISO/IEC 1539-1:2018 §8.5.8.5 - Assumed-shape arrays
- SLATEC Guide Section 6.4: "Array dimensions should be passed as arguments"

**Transformation**:
```fortran
! BEFORE:
SUBROUTINE FOO(X, N)
  INTEGER, INTENT(IN) :: N
  REAL, INTENT(INOUT) :: X(*)  ! Assumed-size

! AFTER:
SUBROUTINE FOO(X)
  REAL, INTENT(INOUT) :: X(:)  ! Assumed-shape
  INTEGER :: N
  N = SIZE(X)  ! Get size from array itself
```

---

### 3. Generic Interfaces for Precision Variants

**Rationale**: SLATEC provides separate routines for single (S), double (D), and complex (C) precision. Modern Fortran allows generic interfaces.

**Reference**:
- ISO/IEC 1539-1:2018 §15.4.3.4 - Generic interfaces
- SLATEC naming convention (Guide Section 7): "S indicates single precision, D double precision, C complex"

**Transformation**:
```fortran
! Modern module with generic interface
MODULE gamma_functions
  INTERFACE gamma
    MODULE PROCEDURE sgamma   ! Single precision
    MODULE PROCEDURE dgamma   ! Double precision
    MODULE PROCEDURE cgamma   ! Complex
  END INTERFACE gamma
END MODULE
```

---

## Documentation Requirements

### For Each Modernised Routine

1. **Preserve Original Prologue** - Keep all original author credits and references
2. **Add Modernisation Note** - Document changes made
3. **Cite Transformations** - Reference this guide section

**Template**:
```fortran
!* MODERNISATION HISTORY
!   2025-12-16  Eliminated GOTO per MODERNISATION_GUIDE.md §1
!               Replaced assumed-size arrays per §2
!   Original: SLATEC v4.1 (1993), modernised by Mehdi Chinoune (2021)
```

---

## Testing Requirements

### Golden Regression Tests

For each modernised routine, we must verify numerical equivalence with the original.

**Reference**:
- SLATEC Guide Section 10: "SLATEC Quick Check Philosophy"
- Original quick checks in `test/` directory

**Procedure**:
1. Capture original routine output for test cases
2. Store as "golden" reference data
3. Compare modernised routine output
4. Tolerance: Machine epsilon (from D1MACH/R1MACH)

---

## Change Log
| Date | File | Change | Reference |
|------|------|--------|-----------|

| 2025-12-16 | Project setup | Copied Mehdi base | Prior Work |
| 2025-12-17 | special/dxlegf,xlegf | GOTO 100 → ERROR STOP | §1, Smith |
| 2025-12-17 | special/dxpqnu,xpqnu | GOTO 100 → DO loop | §1, Smith |
| 2025-12-17 | special/dxqnu,xqnu | GOTO 100 → named CYCLE | §1, Smith |
| 2025-12-17 | special/cbuni,zbuni | GOTO 100 → inline IF | §1, Amos |
| 2025-12-17 | special/cunhj,zunhj,cunik,zunik | GOTO 20 → EXIT | §1, Amos |
| 2025-12-17 | splp/dplpfe,splpfe | GOTO 50 → CYCLE | §1, Hanson |
| 2025-12-17 | fishpack/spelip | GOTO 100 → DO loop | §1, Swarztrauber |
| 2025-12-17 | special/cacai,zacai | GOTO 100 → error return | §1, Amos |
| 2025-12-17 | special/casyi,zasyi | GOTO 20/100 → EXIT | §1, Amos |
| 2025-12-17 | special/cbiry,zbiry | GOTO 50/100 → error returns | §1, Amos |
| 2025-12-17 | special/ckscl,zkscl | GOTO 100/200 → flag+EXIT | §1, Amos |
| 2025-12-17 | special/cmlri,zmlri | GOTO 100/200 → EXIT | §1, Amos |
| 2025-12-18 | Build architecture | Renamed ~1047 .f90 → .inc for fpm | Build fix |
| 2025-12-18 | special/drc3jj,rc3jj | Fixed END DO/END IF nesting | S1, Stone |
| 2025-12-18 | special/drc6j,rc6j | Fixed END DO/END IF nesting | S1, Stone |
| 2025-12-18 | special/dasyjy,asyjy | Removed orphaned END IF | S1, Amos |
| 2025-12-18 | special/zkscl,ckscl | Added early_exit declaration | S1, Amos |
| 2025-12-18 | special/r9lg* d9lg* r9gm* d9gm* | Added missing i declaration | S1 |
| 2025-12-18 | special/dexint,exint | 10 GOTOs each -> 0, internal subroutines | S1, Amos |
| 2025-12-19 | special/dbesi,besi | 26 GOTOs each -> 0, state machine + internal subs | S1, Amos/Daniel |
| 2025-12-19 | special/dbskin,bskin | Removed PURE (calls non-PURE dexint/exint) | Build fix |
| 2025-12-19 | special/dbesj,besj | 23 GOTOs each -> 0, state machine + internal subs | S1, Amos/Daniel/Weston |
