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

## Current MINPACK Coverage

| Test File | Routines | Status |
|-----------|----------|--------|
| test_enorm.f | DENORM (Euclidean norm) | **7/7 PASS** |
| test_qrfac.f | DQRFAC (QR factorisation) | Pending |
| test_dnls1.f | DNLS1 (Nonlinear least squares) | Pending |
| test_dnsq.f | DNSQ (Nonlinear equations) | Pending |

## Running Tests

```bash
cd /c/dev/fortran360
python -m src.fortran360.cli run tests/slatec/minpack/test_enorm.f
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
