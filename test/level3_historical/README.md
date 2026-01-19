# Level 3: Historical Baseline Tests

*"Does the output match what users would have seen on the original hardware?"*

## Purpose

Level 3 tests establish the **last known good** reference. SLATEC was developed and validated on IBM mainframes running FORTRAN G/H compilers. If our modernised code produces different output than the original would have on IBM hardware, we need to understand *why*.

## Platform Specification

- **Hardware**: IBM System/360/370 (emulated via Hercules)
- **Operating System**: MVT (Multiple Variable Tasks) via TK4-
- **Compilers**:
  - IBM FORTRAN G (IEYFORT) - 1966
  - IBM FORTRAN H (IEKQFORT) - 1969
- **Floating-Point**: IBM Hexadecimal (base-16, not IEEE 754)

## Why IBM Hexadecimal Matters

IBM hexadecimal floating-point differs fundamentally from IEEE 754:

| Property | IBM Hex | IEEE 754 |
|----------|---------|----------|
| Base | 16 | 2 |
| Double precision mantissa | 56 bits (14 hex digits) | 52 bits |
| Gradual underflow | No | Yes |
| NaN/Infinity | No | Yes |
| "Wobbling precision" | 0-3 bits lost | None |

A deviation at Level 3 is not necessarily a bug — it may be an unavoidable consequence of the floating-point representation change. But it must be **documented and explained** in [DEVIATIONS.md](../../DEVIATIONS.md).

## Test Location

Actual IBM 360 tests are maintained in:

```
/c/dev/fortran360/tests/slatec/minpack/
```

This is the fortran360 project which provides:
- Hercules emulator integration
- JCL job submission
- FORTRAN G/H compilation
- Authentic IBM 360 execution environment

## Coverage

### MINPACK

| Test File | Routines | Status |
|-----------|----------|--------|
| test_enorm.f | DENORM (Euclidean norm) | **7/7 PASS** |
| test_qrfac.f | DQRFAC (QR factorisation) | Pending |
| test_dnls1.f | DNLS1 (Nonlinear least squares) | Pending |
| test_dnsq.f | DNSQ (Nonlinear equations) | Pending |

### BLAS

| Test File | Routines | Status |
|-----------|----------|--------|
| test_daxpy.f | DAXPY (Y = alpha*X + Y) | **5/5 PASS** |
| test_drotg.f | DROTG (Givens rotation) | **4/4 PASS** |

BLAS tests use Pythagorean triples (3-4-5, 5-12-13, 8-15-17, 7-24-25) for exact integer verification.

### Interpolation

| Test File | Routines | Status |
|-----------|----------|--------|
| test_l3_interpolation.f90 | DPLINT, DPOLVL | ✓ **4/4 PASS** (polynomial) |
| test_l3_interpolation.f90 | DPCHIM, DPCHFE | N/A (post-1980 algorithm) |
| test_l3_interpolation.f90 | DBINT4, DBVALU | N/A (post-1978 + L2 issues) |

**Polynomial Interpolation**: Golden values captured from IBM 360 via Hercules/TK4- on 18 January 2026.
Test case: y = x³ interpolated through x = 0,1,2,3,4, evaluated at midpoints.

**PCHIP**: Fritsch & Carlson (1980) - algorithm post-dates IBM 360 era. No L3 golden values possible.

**B-spline**: de Boor (1978) - algorithm post-dates IBM 360 era. Additionally has L2 mathematical issues.

### Diff_Integ

| Test File | Routines | Status |
|-----------|----------|--------|
| test_l3_diff_integ.f90 | Trapezoidal rule (1960s method) | ✓ **2/2 PASS** |
| test_l3_diff_integ.f90 | QUADPACK (DQAGS, DQAGI, DQNG) | N/A (post-1983) |
| test_l3_diff_integ.f90 | DGAUS8 | N/A (1980s implementation) |

**Trapezoidal Rule**: Basic numerical integration method available in 1960s FORTRAN. Golden values captured from IBM 360.

**QUADPACK**: Piessens, de Doncker, et al. (1983) - algorithms post-date IBM 360 era. No L3 golden values possible.

## Running Tests

### IBM 360 (Hercules)
```bash
cd /c/dev/fortran360
python -m src.fortran360.cli run tests/slatec/minpack/test_enorm.f
python -m src.fortran360.cli run tests/slatec/blas/test_daxpy.f
python -m src.fortran360.cli run tests/slatec/blas/test_drotg.f
```

### Modern Fortran Comparison Tests
```bash
cd /c/dev/slatec-modern
gfortran -o test_l3_blas test/level3_historical/test_l3_linear_blas.f90
./test_l3_blas

gfortran -o test_l3_minpack test/level3_historical/test_l3_minpack.f90
./test_l3_minpack
```

## FORTRAN IV Compatibility Notes

Level 3 tests must be written in FORTRAN IV (FORTRAN 66) compatible syntax:

1. **No IF-THEN-ELSE** — Use arithmetic IF and logical IF with GOTO
2. **No IMPLICIT NONE** — I-N integer convention applies
3. **Fixed-form source** — Columns 1-72 only
4. **No modules** — Subroutines and COMMON blocks only
5. **No DO-WHILE** — Use labelled DO loops with EXIT via GOTO

## Interpreting Failures

If Level 3 fails but Level 2 passes:
- The mathematics is correct
- The platform differs from IBM 360
- Document the deviation with root cause analysis

---

*Last verified on IBM System/360: 1 January 2026*
