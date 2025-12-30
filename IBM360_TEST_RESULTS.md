# SLATEC Library Test Results

## The Great SLATEC Verification Expedition

*Wherein we discover that mathematical libraries from the Reagan era still produce correct numbers when run on compilers from the Johnson administration.*

---

## Test Environment

| Component | Version | Origin Story |
|-----------|---------|--------------|
| **Compiler** | IBM FORTRAN G Level 21 (IEYFORT) | 1966 - older than most of the code it's compiling |
| **System** | MVS 3.8j (TK4-) | The mainframe that refuses to die |
| **Emulator** | Hercules 3.07 | Making virtual iron from actual silicon |
| **Library** | SLATEC 4.1 | Sandia/Los Alamos/Air Force - the nuclear option for numerical computing |

---

## Machine Constants Verification

IBM System/360 uses hexadecimal floating-point arithmetic because *of course it does*.

### I1MACH Results (Integer Machine Constants)
```
I1MACH(1)  = 5       Standard input unit
I1MACH(2)  = 6       Standard output unit
I1MACH(10) = 16      Base for floating-point (hex!)
I1MACH(11) = 6       Hex digits in single precision mantissa
I1MACH(14) = 14      Hex digits in double precision mantissa
```

### R1MACH Results (Single Precision Machine Constants)
```
R1MACH(1) = 0.5398E-78   Smallest positive magnitude
R1MACH(2) = 0.7237E+76   Largest magnitude
R1MACH(4) = 9.54E-07     Machine epsilon (unit roundoff)
```

### D1MACH Results (Double Precision Machine Constants)
```
D1MACH(1) = 5.4E-79      Smallest positive magnitude
D1MACH(2) = 7.2E+75      Largest magnitude
D1MACH(4) = 2.22E-16     Machine epsilon (unit roundoff)
```

---

## Test Results Summary

| Routine | Category | Precision | Status | Max Error |
|---------|----------|-----------|--------|-----------|
| **GAMLN** | Special Functions | Single | ✅ PASS | ~10⁻⁶ |
| **DGAMLN** | Special Functions | Double | ✅ PASS | ~10⁻¹⁵ |
| **PYTHAG** | Numerical Utilities | Single | ✅ PASS | ~10⁻⁶ |
| **ENORM** | Vector Operations | Single | ✅ PASS | 0 |
| **DENORM** | Vector Operations | Double | ✅ PASS | ~10⁻¹⁵ |
| **CDIV** | Complex Arithmetic | Single | ✅ PASS | ~10⁻⁷ |
| **CSROOT** | Complex Arithmetic | Single | ✅ PASS | 0 |
| **VNWRMS** | Vector Operations | Single | ✅ PASS | ~10⁻⁶ |
| **HVNRM** | Vector Operations | Single | ✅ PASS | 0 |
| **INTRV** | Interpolation | Integer | ✅ PASS | exact |
| **RC** | Elliptic Integrals | Single | ✅ PASS | ~10⁻⁶ |
| **DHVNRM** | Vector Operations | Double | ✅ PASS | 0 |
| **DVNRMS** | Vector Operations | Double | ✅ PASS | 0 |
| **DRC** | Elliptic Integrals | Double | ✅ PASS | ~10⁻¹⁶ |

**14 routines tested. 14 routines passed. 0 surprises.**

---

## GAMLN (Single Precision Log-Gamma) Test Results

### Test Date: 30 December 2025

The GAMLN function computes ln(Γ(z)) where Γ is the gamma function.
For positive integers: GAMLN(n) = ln((n-1)!)

### Integer Argument Tests (Table Lookup Path)

| Input | Computed | Expected | Error | Status |
|-------|----------|----------|-------|--------|
| 1.0 | 0.000000 | 0.000000 | 0.0 | PASS |
| 2.0 | 0.000000 | 0.000000 | 0.0 | PASS |
| 3.0 | 0.6931471 | 0.6931472 | 6.0E-08 | PASS |
| 5.0 | 3.178053 | 3.178053 | 0.0 | PASS |
| 10.0 | 12.80183 | 12.80183 | 9.5E-07 | PASS |

### Non-Integer Argument Tests (Asymptotic Expansion Path)

| Input | Computed | Expected | Error | Status |
|-------|----------|----------|-------|--------|
| 0.5 | 0.5723624 | 0.5723649 | 2.4E-06 | PASS |
| 1.5 | -0.1207842 | -0.1207822 | 2.0E-06 | PASS |

---

## DGAMLN (Double Precision Log-Gamma) Test Results

### Test Date: 30 December 2025

Double precision version using 14 hexadecimal digits (~16 decimal digits).

### Integer Argument Tests (Table Lookup Path)

| Input | Computed | Expected | Error | Status |
|-------|----------|----------|-------|--------|
| 1.0 | 0.0000000000 | 0.0000000000 | 0.0 | PASS |
| 3.0 | 0.6931471806 | 0.6931471806 | 1.0E-12 | PASS |
| 5.0 | 3.1780538303 | 3.1780538303 | 0.0 | PASS |
| 10.0 | 12.8018274801 | 12.8018274801 | 1.1E-15 | PASS |

### Non-Integer Argument Tests (Asymptotic Expansion Path)

| Input | Computed | Expected | Error | Status |
|-------|----------|----------|-------|--------|
| 0.5 | 0.5723649429 | 0.5723649429 | 2.2E-16 | PASS |
| 1.5 | -0.1207822376 | -0.1207822376 | 2.2E-16 | PASS |
| 100.0 | 359.1342054 | 359.1342054 | 0.0 | PASS |

---

## PYTHAG (Overflow-Safe Hypotenuse) Test Results

### Test Date: 30 December 2025

Computes √(a² + b²) without overflow or destructive underflow using Moler-Morrison algorithm.

| Input (a, b) | Computed | Expected | Error | Status |
|--------------|----------|----------|-------|--------|
| (3.0, 4.0) | 5.0000 | 5.0 | 9.5E-07 | PASS |
| (5.0, 12.0) | 13.0000 | 13.0 | 0.0 | PASS |
| (1.0, 1.0) | 1.4142 | √2 | 9.5E-07 | PASS |
| (0.0, 5.0) | 5.0000 | 5.0 | 0.0 | PASS |
| **(1E30, 1E30)** | **1.414E30** | **√2 × 10³⁰** | **1.7E-06** | **PASS** |

The last test is the important one - naive computation of √(a² + b²) with a=b=10³⁰
would overflow, but PYTHAG handles it correctly.

---

## ENORM (Single Precision Euclidean Norm) Test Results

### Test Date: 30 December 2025

Computes ‖x‖₂ = √(Σxᵢ²) with overflow protection using three-accumulator algorithm.

| Input Vector | Computed | Expected | Error | Status |
|--------------|----------|----------|-------|--------|
| [3, 4] | 5.0000 | 5.0 | 0.0 | PASS |
| [1, 2, 2] | 3.0000 | 3.0 | 0.0 | PASS |
| [1, 1, 1, 1, 1] | 2.2361 | √5 | 9.5E-07 | PASS |

---

## DENORM (Double Precision Euclidean Norm) Test Results

### Test Date: 30 December 2025

Double precision version of ENORM.

| Input Vector | Computed | Expected | Error | Status |
|--------------|----------|----------|-------|--------|
| [3, 4] | 5.0000000000 | 5.0 | 0.0 | PASS |
| [1, 2, 2] | 3.0000000000 | 3.0 | 0.0 | PASS |
| [1, 1, 1, 1, 1] | 2.2360679775 | √5 | 2.2E-16 | PASS |
| **(1E30, 1E30)** | **1.414213562E30** | **√2 × 10³⁰** | **0.0** | **PASS** |

---

## CDIV (Complex Division) Test Results

### Test Date: 30 December 2025

Computes (a + bi) / (c + di) with scaling to avoid overflow.

| Numerator | Denominator | Computed | Expected | Error | Status |
|-----------|-------------|----------|----------|-------|--------|
| 2 + 3i | 1 + 1i | 2.5 + 0.5i | 2.5 + 0.5i | 0.0 | PASS |
| 4 + 0i | 2 + 0i | 2.0 + 0.0i | 2.0 + 0.0i | 0.0 | PASS |
| 1 + 1i | 0 + 1i | 1.0 - 1.0i | 1.0 - 1.0i | 0.0 | PASS |
| 3 + 4i | 3 - 4i | -0.28 + 0.96i | -0.28 + 0.96i | 6.0E-08 | PASS |

---

## CSROOT (Complex Square Root) Test Results

### Test Date: 30 December 2025

Computes √(x + yi) using PYTHAG for numerical stability.

| Input | Computed | Expected | Error | Status |
|-------|----------|----------|-------|--------|
| 4 + 0i | 2.0 + 0.0i | 2 + 0i | 0.0 | PASS |
| -1 + 0i | 0.0 + 1.0i | 0 + i | 0.0 | PASS |
| 0 + 4i | 1.414 + 1.414i | √2 + √2i | 0.0 | PASS |
| 3 + 4i | 2.0 + 1.0i | 2 + i | 0.0 | PASS |
| 0 - 4i | 1.414 - 1.414i | √2 - √2i | 0.0 | PASS |

---

## VNWRMS (Weighted Root-Mean-Square Norm) Test Results

### Test Date: 30 December 2025

Computes √((1/n) × Σ(vᵢ/wᵢ)²) for ODE integrator error control.

| V Vector | W Vector | Computed | Expected | Error | Status |
|----------|----------|----------|----------|-------|--------|
| [1,1,1,1,1] | [1,1,1,1,1] | 1.0000 | 1.0 | 0.0 | PASS |
| [2,2,2,2,2] | [2,2,2,2,2] | 1.0000 | 1.0 | 0.0 | PASS |
| [3, 4] | [1, 1] | 3.5355 | √12.5 | 9.5E-07 | PASS |
| [6, 8] | [2, 2] | 3.5355 | √12.5 | 9.5E-07 | PASS |

---

## HVNRM (Maximum/Infinity Norm) Test Results

### Test Date: 30 December 2025

Computes ‖v‖∞ = max|vᵢ| for ODE integrator step control.

| Input Vector | Computed | Expected | Error | Status |
|--------------|----------|----------|-------|--------|
| [1, 2, 3, 4, 5] | 5.0 | 5.0 | 0.0 | PASS |
| [-7, 3, 5, -2, 1] | 7.0 | 7.0 | 0.0 | PASS |
| [0.1, 0.2, 0.3] | 0.3 | 0.3 | 0.0 | PASS |

---

## INTRV (Interval Search for B-Splines) Test Results

### Test Date: 30 December 2025

Binary search to find interval containing a value in a knot vector.
Written by Carl de Boor for the classic B-spline package.

Knot vector: XT = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

| X | ILEFT | MFLAG | Meaning | Status |
|---|-------|-------|---------|--------|
| 3.5 | 3 | 0 | In interval [3,4) | PASS |
| 7.0 | 7 | 0 | In interval [7,8) | PASS |
| 0.5 | 1 | -1 | Below range | PASS |
| 10.5 | 10 | 1 | Above range | PASS |

Sequential search with ILO caching also verified.

---

## RC (Carlson Elliptic Integral) Test Results

### Test Date: 30 December 2025

Carlson's RC function computes the degenerate elliptic integral:

```
RC(x,y) = (1/2) ∫₀^∞ (t+x)^(-1/2) (t+y)^(-1) dt
```

This elegant function can express π, logarithms, and inverse trigonometric
functions through various identities. The duplication theorem is iterated
until convergence, then a Taylor series is evaluated.

| Test | Computed | Identity | Expected | Error | Status |
|------|----------|----------|----------|-------|--------|
| RC(0, 1/4) | **3.1415930** | **π** | 3.1415927 | 9.5E-07 | PASS |
| RC(1, 2) | 0.7853984 | arctan(1) = π/4 | 0.7853982 | 1.8E-07 | PASS |
| RC(1/16, 1/8) | **3.1415930** | **π** (alternate) | 3.1415927 | 9.5E-07 | PASS |
| RC(9/4, 2) | 0.6931474 | ln(2) | 0.6931472 | 1.8E-07 | PASS |

*Yes, we're computing π using the Carlson duplication theorem on 1966 hardware.*

---

## DHVNRM (Double Precision Maximum Norm) Test Results

### Test Date: 30 December 2025

Double precision version of HVNRM.

| Input Vector | Computed | Expected | Error | Status |
|--------------|----------|----------|-------|--------|
| [1, 2, 3, 4, 5] | 5.0000000000 | 5.0 | 0.0 | PASS |
| [-7, 3, 5, -2, 1] | 7.0000000000 | 7.0 | 0.0 | PASS |

---

## DVNRMS (Double Precision Weighted RMS Norm) Test Results

### Test Date: 30 December 2025

Double precision weighted root-mean-square norm for ODE integrators.

| V Vector | W Vector | Computed | Expected | Error | Status |
|----------|----------|----------|----------|-------|--------|
| [1,1,1,1,1] | [1,1,1,1,1] | 1.0000000000 | 1.0 | 0.0 | PASS |
| [3, 4] | [1, 1] | 3.5355339059 | √12.5 | 0.0 | PASS |

---

## DRC (Double Precision Carlson Elliptic Integral) Test Results

### Test Date: 30 December 2025

Double precision version of RC. Computing transcendental constants to 15-16 significant digits!

| Test | Computed | Identity | Expected | Error | Status |
|------|----------|----------|----------|-------|--------|
| DRC(0, 1/4) | **3.14159265358979** | **π** | 3.14159265358979324 | 2.4E-15 | PASS |
| DRC(1, 2) | 0.785398163397448 | arctan(1) = π/4 | 0.78539816339744831 | 6.9E-16 | PASS |
| DRC(1/16, 1/8) | **3.14159265358979** | **π** | 3.14159265358979324 | 2.4E-15 | PASS |
| DRC(9/4, 2) | 0.693147180559945 | ln(2) | 0.69314718055994531 | **5.6E-17** | PASS |

The ln(2) computation achieves an error of 5.6×10⁻¹⁷ - essentially at machine epsilon.
*Computing transcendentals to 16 significant figures on vintage iron. Not bad for 1966.*

---

## Analysis

### Accuracy Notes

**Single Precision (6 hex digits = 24 bits ≈ 7.2 decimal digits)**
- Observed errors: 10⁻⁵ to 10⁻⁷
- Within expected tolerance for IBM 360 hexadecimal floating-point

**Double Precision (14 hex digits = 56 bits ≈ 16.8 decimal digits)**
- Observed errors: 10⁻¹⁵ to 10⁻¹⁶
- Approaching machine epsilon, as expected

### FORTRAN IV Compatibility Fixes Applied

| Issue | Original | Fixed For FORTRAN G |
|-------|----------|---------------------|
| Deck markers | `*DECK NAME` | `C     DECK NAME` |
| Generic intrinsics | `MAX`, `MIN`, `ABS` | `AMAX1`, `AMIN1`, `ABS` (SP) |
| Generic intrinsics | `MAX`, `MIN`, `ABS`, `LOG` | `DMAX1`, `DMIN1`, `DABS`, `DLOG` (DP) |
| Assumed-size arrays | `X(*)` | `X(N)` |
| SAVE statement | `SAVE VAR` | (removed - DATA provides static init) |
| Hex constants | `Z'00100000'` | Decimal via EQUIVALENCE |
| Block IF | `IF (...) THEN` | `IF (...) GO TO label` |

---

## Files Used

```
tests/slatec/
├── Machine Constants
│   ├── i1mach.f         # Integer machine constants for IBM 360
│   ├── r1mach.f         # Single precision constants (hex FP)
│   └── d1mach.f         # Double precision constants (hex FP)
│
├── Special Functions
│   ├── gamln.f          # Single precision log-gamma
│   ├── dgamln.f         # Double precision log-gamma
│   ├── test_gamln.f     # GAMLN test driver
│   └── test_dgamln.f    # DGAMLN test driver
│
├── Numerical Utilities
│   ├── pythag.f         # Overflow-safe √(a²+b²)
│   ├── enorm.f          # Single precision Euclidean norm
│   ├── denorm.f         # Double precision Euclidean norm
│   ├── vnwrms.f         # Weighted RMS norm (SP)
│   ├── dvnrms.f         # Weighted RMS norm (DP)
│   ├── hvnrm.f          # Maximum norm (SP)
│   ├── dhvnrm.f         # Maximum norm (DP)
│   └── test_numutil.f   # PYTHAG/ENORM test driver
│
├── Complex Arithmetic
│   ├── cdiv.f           # Complex division
│   ├── csroot.f         # Complex square root
│   └── test_complex_all.f  # Combined complex tests
│
├── Interpolation
│   ├── intrv.f          # Interval search (Carl de Boor)
│   └── test_intrv_hvnrm.f  # INTRV/HVNRM/RC test driver
│
├── Elliptic Integrals
│   ├── rc.f             # Carlson RC single precision
│   ├── drc.f            # Carlson RC double precision (computes π to 15 digits!)
│   └── test_dp_final.f  # Final DP tests (DRC, DHVNRM, DVNRMS)
│
└── Combined Tests
    └── test_gamln_vnwrms.f  # GAMLN + VNWRMS combined
```

---

## Verdict

**ALL 14 ROUTINES TESTED. ALL 14 PASSED.**

The 1966 IBM FORTRAN G compiler successfully compiled and executed SLATEC
mathematical library code, producing numerically correct results within
the expected precision limits of IBM 360 hexadecimal floating-point arithmetic.

Key achievements:
- ✅ Special functions (log-gamma) in both single and double precision
- ✅ Overflow-protected vector norms (Euclidean, weighted RMS, infinity)
- ✅ Complex arithmetic with proper branch cuts
- ✅ Numerical utilities for ODE integrators
- ✅ B-spline interval search algorithm
- ✅ Carlson elliptic integrals computing π and ln(2) via duplication theorem
- ✅ Double precision routines achieving 10⁻¹⁶ accuracy (machine epsilon!)

*Take that, "modern" numerics libraries.*

---

## Historical Context

- **SLATEC**: 1982 joint effort by Sandia, Los Alamos, and Air Force labs
- **FORTRAN G**: 1966 IBM compiler, predating SLATEC by 16 years
- **IBM 360**: 1964 architecture, still running (virtually) in 2025

The fact that 1982 library code runs perfectly on a 1966 compiler is a
testament to FORTRAN's legendary backward compatibility. They literally
don't make 'em like they used to - because they don't have to.

---

*"In FORTRAN we trust, but we verify on vintage iron."*
