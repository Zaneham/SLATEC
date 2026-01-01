# SLATEC-Modern Test Results

*Comprehensive test execution log — 1 January 2026*

## Executive Summary

| Suite | L1 | L2 | L3 | L4 | Total |
|-------|-----|-----|-----|------|-------|
| **BLAS** | 18/18 | 16/16 | 9/9 | 65/65 | **108/108** |
| **MINPACK** | 9/9 | 17/17 | 7/7 | 45/45 | **78/78** |
| **LINPACK** | 6/6 | 13/13 | 9/9 | 18/18 | **46/46** |
| **Combined** | 33/33 | 46/46 | 25/25 | 128/128 | **232/232** |

**All tests pass with safe compiler flags (`-O2` or `-O3`).**

---

## Compiler Flag Matrix

### gfortran 13.2.0 on Windows (MSYS2)

| Flag Combination | BLAS L4 | MINPACK L4 | LINPACK L4 | Safe? |
|------------------|---------|------------|------------|-------|
| `-O2` | 65/65 ✓ | 45/45 ✓ | 18/18 ✓ | **YES** |
| `-O3` | 65/65 ✓ | 45/45 ✓ | 18/18 ✓ | **YES** |
| `-O3 -march=native` | 65/65 ✓ | 45/45 ✓ | 18/18 ✓ | **YES** |
| `-O2 -flto` | 65/65 ✓ | 45/45 ✓ | 18/18 ✓ | **YES** |
| `-O2 -ffast-math` | 59/65 ✗ | 44/45 ✗ | 18/18 ⚠ | **NO** |
| `-Ofast` | 59/65 ✗ | 44/45 ✗ | 18/18 ⚠ | **NO** |
| `-O2 -ffinite-math-only` | 61/65 ✗ | — | 18/18 ⚠ | **NO** |
| `-O2 -funsafe-math-optimizations` | 63/65 ✗ | — | 18/18 ✓ | **NO** |

**Note:** LINPACK L4 tests pass with all flags but show warnings (⚠) for Inf/NaN handling when `-ffinite-math-only` is active. The tests detect but don't fail on altered IEEE behavior.

---

## Detailed Failure Analysis

### `-ffast-math` Failures (BLAS: 6 failures, MINPACK: 1 failure)

#### BLAS Failures

| Category | Test | Reason |
|----------|------|--------|
| Subnormal Handling | DAXPY flushed subnormals to zero | FTZ (Flush-to-Zero) mode enabled |
| Subnormal Handling | DAXPY with subnormal alpha flushed | FTZ mode enabled |
| Inf/NaN Propagation | DAXPY infinity propagation unexpected | `-ffinite-math-only` assumes no Inf |
| Inf/NaN Propagation | DAXPY did not propagate NaN | `-ffinite-math-only` assumes no NaN |
| NaN Variants | Quiet NaN not detected | NaN checks optimized away |
| NaN Variants | NaN propagation incorrect | NaN handling undefined |

#### MINPACK Failures

| Category | Test | Reason |
|----------|------|--------|
| Edge Cases | Subnormal inputs: non-finite result | FTZ flushes input to zero, causing division by zero |

### `-ffinite-math-only` Failures (BLAS: 4 failures)

| Category | Test | Reason |
|----------|------|--------|
| Inf/NaN Propagation | DAXPY infinity propagation | Assumes no infinities exist |
| Inf/NaN Propagation | DAXPY NaN propagation | Assumes no NaNs exist |
| NaN Variants | Quiet NaN not detected | `ieee_is_nan()` returns false |
| NaN Variants | NaN propagation incorrect | NaN comparisons undefined |

### `-funsafe-math-optimizations` Failures (BLAS: 2 failures)

| Category | Test | Reason |
|----------|------|--------|
| Subnormal Handling | DAXPY flushed subnormals | Includes FTZ behavior |
| Subnormal Handling | DAXPY with subnormal alpha | Includes FTZ behavior |

---

## What `-ffast-math` Actually Does

The `-ffast-math` flag is equivalent to:

```
-fno-math-errno
-funsafe-math-optimizations
-ffinite-math-only
-fno-rounding-math
-fno-signaling-nans
-fcx-limited-range
-fexcess-precision=fast
```

### Breaking Down the Damage

| Sub-flag | Effect | Tests Affected |
|----------|--------|----------------|
| `-ffinite-math-only` | Assumes no NaN or Inf | Inf/NaN propagation, NaN variants |
| `-funsafe-math-optimizations` | Enables FTZ, DAZ | Subnormal handling |
| `-fno-signed-zeros` | Treats -0 as +0 | Signed zero tests (warnings) |
| `-fno-trapping-math` | No FP exceptions | Exception detection |
| `-fassociative-math` | Reorders operations | Associativity tests |

---

## Test Categories Explained

### Level 1: Regression (33 tests)
- **Purpose**: Does the code work?
- **Stability**: May be modified when code changes
- **Pass Rate**: 100%
- **Breakdown**: BLAS (18), MINPACK (9), LINPACK (6)

### Level 2: Mathematical (46 tests)
- **Purpose**: Does output match mathematical truth?
- **Stability**: Read-only (reference values)
- **Pass Rate**: 100%
- **Breakdown**: BLAS (16), MINPACK (17), LINPACK (13)

### Level 3: Historical (25 tests)
- **Purpose**: Does output match IBM 360?
- **Stability**: Read-only (golden values from Hercules/TK4-)
- **Pass Rate**: 100%
- **Breakdown**: BLAS (9), MINPACK (7), LINPACK (9)

### Level 4: Hostile (128 tests)
- **Purpose**: What breaks under stress?
- **Stability**: Read-only (platform detection)
- **Pass Rate**: 100% with safe flags
- **Breakdown**: BLAS (65), MINPACK (45), LINPACK (18)

---

## Level 4 Test Category Details

### BLAS Level 4 (65 tests)

| Category | Tests | What It Detects |
|----------|-------|-----------------|
| Subnormal Handling | 4 | FTZ/DAZ, `-ffast-math` flush |
| Signed Zero | 4 | IEEE 754 ±0 compliance |
| Inf/NaN Propagation | 4 | `-ffinite-math-only` effects |
| Extreme Values | 5 | Near overflow/underflow |
| Accumulation Precision | 3 | Precision loss in summations |
| SIMD Edge Cases | 5 | Vectorization boundaries |
| Rounding Modes | 5 | All four IEEE modes |
| FMA Detection | 4 | Fused multiply-add |
| Catastrophic Cancellation | 4 | Precision loss |
| Associativity | 4 | Operation reordering |
| Extended Precision (x87) | 3 | 80-bit register effects |
| Reproducibility | 3 | Determinism across runs |
| Compiler Flag Detection | 5 | Dangerous flag warnings |
| NaN Variants | 4 | qNaN, sNaN, propagation |
| ULP Accuracy | 4 | Mathematical precision |
| Memory Alignment | 4 | Various access patterns |

### MINPACK Level 4 (45 tests)

| Category | Tests | What It Detects |
|----------|-------|-----------------|
| ULP Precision | 4 | Accuracy vs constants |
| Edge Cases | 5 | Subnormals, overflow |
| Associativity | 4 | Operation reordering |
| Rounding Modes | 4 | IEEE rounding sensitivity |
| FMA Effects | 3 | Fused multiply-add |
| Catastrophic Cancellation | 3 | Norm computation precision |
| Condition Number | 3 | Ill-conditioned inputs |
| Convergence Reproducibility | 3 | Determinism |
| Compiler Flag Detection | 4 | Dangerous flag warnings |
| Extended Precision (x87) | 3 | 80-bit precision leakage |
| Scaling Sensitivity | 3 | Homogeneity preservation |
| Jacobian Computation | 3 | Finite difference accuracy |
| LM Parameter Sensitivity | 3 | Trust region effects |

### LINPACK Level 4 (18 tests)

| Category | Tests | What It Detects |
|----------|-------|-----------------|
| Ill-Conditioned Matrices | 2 | Hilbert, Vandermonde stability |
| Subnormal Elements | 2 | FTZ/DAZ effects on matrix elements |
| Extreme Scaling | 3 | Near overflow/underflow |
| Pivoting Stress | 2 | Partial pivoting correctness |
| Near-Singular | 2 | Detection of rank deficiency |
| Inf/NaN Handling | 2 | IEEE special value propagation |
| Cholesky Edge Cases | 3 | SPD detection, conditioning |
| Reproducibility | 2 | Deterministic results |

---

## Recommendations

### Safe Compiler Flags

```bash
# Production builds
gfortran -O2 -o myprogram myprogram.f90

# Maximum safe optimization
gfortran -O3 -march=native -o myprogram myprogram.f90

# With LTO (safe)
gfortran -O2 -flto -o myprogram myprogram.f90
```

### Dangerous Flags — DO NOT USE

```bash
# NEVER use these with SLATEC
-ffast-math        # Breaks subnormals, Inf, NaN
-Ofast             # Includes -ffast-math
-ffinite-math-only # Breaks Inf/NaN handling
-funsafe-math-optimizations  # Breaks subnormals
```

### If You Must Use Fast Math

If performance is critical and you understand the risks:

1. **Test thoroughly** with Level 4 tests
2. **Avoid subnormal inputs** — scale your data
3. **Never rely on NaN/Inf detection** — check inputs beforehand
4. **Document the limitation** in your code

---

## Platform Information

Tests run on:
- **OS**: Windows (MSYS2/MinGW64)
- **Compiler**: GCC/gfortran 13.2.0
- **Architecture**: x86_64
- **IEEE 754**: Full support
- **Date**: 1 January 2026

### Machine Constants

| Constant | Value |
|----------|-------|
| `EPSILON(1.0d0)` | 2.220E-16 |
| `TINY(1.0d0)` | 2.225E-308 |
| `HUGE(1.0d0)` | 1.798E+308 |
| Infinity | Supported |
| NaN | Supported |

---

## Historical Comparison

### IBM System/360 vs Modern IEEE 754

| Property | IBM 360 | IEEE 754 | Test Impact |
|----------|---------|----------|-------------|
| Base | 16 (hexadecimal) | 2 (binary) | Different rounding |
| Mantissa (double) | 56 bits | 52 bits | Slight precision difference |
| Subnormals | No | Yes | L4 subnormal tests N/A on 360 |
| NaN | No | Yes | L4 NaN tests N/A on 360 |
| Infinity | No | Yes | L4 Inf tests N/A on 360 |
| Wobbling precision | 0-3 bits | None | May see 1-3 ULP difference |

**Level 3 tests verify our modern output matches IBM 360 golden values.**

---

## Running All Tests

```bash
cd /c/dev/slatec-modern

# Level 1: Regression
gfortran -O2 -o test_l1_blas test/level1_regression/test_l1_linear_blas.f90 && ./test_l1_blas
gfortran -O2 -o test_l1_minpack test/level1_regression/test_l1_minpack.f90 && ./test_l1_minpack
gfortran -O2 -o test_l1_linpack test/level1_regression/test_l1_linear_linpack.f90 && ./test_l1_linpack

# Level 2: Mathematical
gfortran -O2 -o test_l2_blas test/level2_mathematical/test_l2_linear_blas.f90 && ./test_l2_blas
gfortran -O2 -o test_l2_minpack test/level2_mathematical/test_l2_minpack_mgh.f90 && ./test_l2_minpack
gfortran -O2 -o test_l2_linpack test/level2_mathematical/test_l2_linear_linpack.f90 && ./test_l2_linpack

# Level 3: Historical
gfortran -O2 -o test_l3_blas test/level3_historical/test_l3_linear_blas.f90 && ./test_l3_blas
gfortran -O2 -o test_l3_minpack test/level3_historical/test_l3_minpack.f90 && ./test_l3_minpack
gfortran -O2 -o test_l3_linpack test/level3_historical/test_l3_linear_linpack.f90 && ./test_l3_linpack

# Level 4: Hostile (safe flags)
gfortran -O2 -o test_l4_blas test/level4_hostile/test_l4_linear_blas.f90 && ./test_l4_blas
gfortran -O2 -o test_l4_minpack test/level4_hostile/test_l4_minpack.f90 && ./test_l4_minpack
gfortran -O2 -o test_l4_linpack test/level4_hostile/test_l4_linear_linpack.f90 && ./test_l4_linpack

# Level 4: Hostile (with hostile flags — expected failures)
gfortran -Ofast -o test_l4_hostile test/level4_hostile/test_l4_linear_blas.f90 && ./test_l4_hostile
```

---

*"The first principle is that you must not fool yourself — and you are the easiest person to fool."* — Richard Feynman
