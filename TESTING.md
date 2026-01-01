# SLATEC-Modern Test Coverage

*A Four-Level Testing Framework for Numerical Software, or: How to Find Where Things Went Wrong*

## The Four Levels

This testing framework is designed not merely to detect failures, but to *locate* them. When a deviation occurs, the level at which it manifests tells us where to look.

| Level | Name | Question Answered | What We Change |
|-------|------|-------------------|----------------|
| **1** | Regression | Does the code work? Does it still work after changes? | ✓ Unilateral changes permitted |
| **2** | Mathematical | Does the algorithm match the mathematics? | ✗ Read-only |
| **3** | Historical | Does output match what IBM 360/370 users saw? | ✗ Read-only |
| **4** | Hostile | Do compilers, processors, or OSes change anything? | ✗ Read-only |

### Interpreting Failures

- **Level 1 fails** → Code bug. Fix the code.
- **Level 2 fails** → Algorithm doesn't match mathematics. Either the original was wrong, or the modernisation broke something. Investigate.
- **Level 3 fails** → Deviation from historical baseline. Document in [DEVIATIONS.md](DEVIATIONS.md) with explanation.
- **Level 4 fails** → Portability issue. Compiler/platform-specific behaviour. Document and potentially ifdef.

The beauty of this structure is that a failure at Level N but success at Level N-1 immediately narrows the search space. If Level 2 passes but Level 3 fails, we know the mathematics is correct but something about the platform differs from the original IBM hardware.

---

## Level 1: Regression Tests

*"Does the code work? If we change things, does it still work?"*

These are the tests we run continuously. They verify basic functionality and catch regressions when code is modified. **This is the only level where we make unilateral changes.**

| Module | Status | Location |
|--------|--------|----------|
| service | ⏳ Pending | `test/level1_regression/` |
| special_functions | Partial (Bessel) | `test/test_bessel_regression.f90` |
| linear (BLAS) | ✓ **18/18 PASS** | `test/level1_regression/test_l1_linear_blas.f90` |
| linear (LINPACK) | ⏳ In Progress | `test/level1_regression/test_l1_linear_linpack.f90` |
| diff_integ | Partial (QUADPACK) | `test/test_diff_integ_quadpack.f90` |
| diff_integ_eq | Partial | `test/test_diff_integ_eq.f90` |
| interpolation | ⏳ Pending | — |
| integ_trans | ⏳ Pending | — |
| approximation (MINPACK) | ✓ **9/9 PASS** | `test/level1_regression/test_l1_minpack.f90` |
| nonlin_eq | ⏳ Pending | — |
| optimisation | ⏳ Pending | — |
| data_handling | ⏳ Pending | — |

---

## Level 2: Mathematical Verification

*"Does the algorithm match the original mathematics?"*

These tests compare computed values against authoritative mathematical references. If Level 2 fails, either:
1. The original SLATEC algorithm was incorrect (rare but documented historically)
2. The modernisation introduced a mathematical error
3. The reference value is wrong (check multiple sources)

**References:**
- Abramowitz & Stegun (1964) — *Handbook of Mathematical Functions*
- NIST Digital Library of Mathematical Functions — https://dlmf.nist.gov/
- Closed-form solutions where available (e.g. the integral of x squared from 0 to 1 equals one-third)

| Module | Status | Reference Source |
|--------|--------|------------------|
| special_functions | Partial | A&S Tables 9.1, 9.8 (Bessel) |
| linear (BLAS) | ✓ **16/16 PASS** | Lawson et al., ACM TOMS 5(3), 1979 |
| diff_integ | Partial | Closed-form integrals |
| approximation (MINPACK) | ✓ **17/17 PASS** | Moré, Garbow, Hillstrom (1981) |

**Test Files:**
- `test/test_bessel_golden.f90` — Bessel J, I, K vs A&S tables
- `test/test_bessel_reference.f90` — Extended validation
- `test/test_specfun_bessel.f90` — Edge cases
- `test/level2_mathematical/test_l2_linear_blas.f90` — BLAS linearity, orthogonality, Pythagorean properties
- `test/level2_mathematical/test_l2_minpack_mgh.f90` — MGH test functions (Rosenbrock, Powell, Helical Valley, Freudenstein-Roth)

---

## Level 3: Historical Baseline (IBM 360/370)

*"Does the output match what users would have seen on the original hardware?"*

This is the **"last known good"** reference. SLATEC was developed and validated on IBM mainframes. If our modernised code produces different output than the original would have on IBM hardware, we need to understand *why*.

**Platform:** IBM System/360 via Hercules emulation (TK4-/MVT)
**Compiler:** IBM FORTRAN G (1966) and FORTRAN H (1969)
**Floating-Point:** IBM Hexadecimal (base-16, not IEEE 754)

### Why This Matters

IBM hexadecimal floating-point differs from IEEE 754:
- Base 16 vs base 2
- Different mantissa lengths (56 bits vs 52 bits for double)
- No gradual underflow
- No NaN or Infinity
- "Wobbling precision" (0-3 bits lost depending on normalisation)

A deviation at Level 3 is not necessarily a bug, it may be an unavoidable consequence of the floating-point representation change. But it must be **documented and explained**.

### Test Location

`/c/dev/fortran360/tests/slatec/`

### Coverage

| Module | Routines Tested | Status |
|--------|-----------------|--------|
| service | D1MACH, R1MACH, I1MACH | ✓ PASS |
| special_functions | GAMLN, DGAMLN, CDIV, CSROOT | ✓ PASS |
| linear (BLAS) | DAXPY, DROTG | ✓ **9/9 PASS** |
| linear | ENORM, PYTHAG | ✓ PASS |
| approximation (MINPACK) | DENORM | ✓ **7/7 PASS** |
| approximation (MINPACK) | DQRFAC | ⏳ Pending |
| approximation (MINPACK) | DNLS1 | ⏳ Pending |
| approximation (MINPACK) | DNSQ | ⏳ Pending |

*BLAS tested with Pythagorean triple golden values. FORTRAN IV source in `/c/dev/fortran360/tests/slatec/blas/`.*
*DENORM tested on IBM System/360 (Hercules/TK4-) using FORTRAN G, 1 January 2026.*

---

## Level 4: Hostile Tests

*"Do compilers, processors, or operating systems change anything?"*

These are the portability stress tests. The same code compiled with different compilers, on different processors, under different operating systems, should produce the same results (within floating-point tolerance). When it doesn't, we need to know.

### Compiler Matrix

| Compiler | Platform | Status |
|----------|----------|--------|
| gfortran 13.x | Linux (x86_64) | CI tested |
| gfortran 13.x | macOS (ARM64) | CI tested |
| gfortran 13.x | Windows (MSYS2) | CI tested |
| ifort/ifx | Linux | ⏳ Untested |
| flang | Linux | ⏳ Untested |
| nvfortran | Linux | ⏳ Untested |

### BLAS Level 4 Results

| Test Category | Tests | Default | `-ffast-math` |
|---------------|-------|---------|---------------|
| Subnormals | 3 | ✓ PASS | May FAIL (DAZ/FTZ) |
| Signed Zero | 3 | ✓ PASS | ✓ PASS |
| Inf/NaN Propagation | 3 | ✓ PASS | ✓ PASS |
| Extreme Values | 4 | ✓ PASS | ✓ PASS |
| Accumulation Precision | 3 | ✓ PASS | ✓ PASS |
| SIMD Edge Cases | 4 | ✓ PASS | ✓ PASS |
| Rounding Modes | 4 | ✓ PASS | ✓ PASS |
| FMA Detection | 3 | ✓ PASS | ✓ PASS |
| Catastrophic Cancellation | 3 | ✓ PASS | May vary |
| Associativity | 4 | ✓ PASS | May FAIL |
| Extended Precision (x87) | 4 | ✓ PASS | ✓ PASS |
| Reproducibility | 3 | ✓ PASS | ✓ PASS |
| Compiler Flag Detection | 4 | ✓ PASS | FAIL (expected) |
| NaN Variants | 4 | ✓ PASS | May FAIL |
| ULP Accuracy | 4 | ✓ PASS | ✓ PASS |
| Memory Alignment | 12 | ✓ PASS | ✓ PASS |
| **Total** | **65** | **65/65 PASS** | May vary |

### MINPACK Level 4 Results

| Test Category | Tests | Default | `-ffast-math` |
|---------------|-------|---------|---------------|
| ULP Precision | 4 | ✓ PASS | ✓ PASS |
| Edge Cases | 5 | ✓ PASS | May FAIL |
| Associativity | 4 | ✓ PASS | May FAIL |
| Rounding Modes | 4 | ✓ PASS | ✓ PASS |
| FMA Effects | 3 | ✓ PASS | ✓ PASS |
| Catastrophic Cancellation | 3 | ✓ PASS | May vary |
| Condition Number | 3 | ✓ PASS | ✓ PASS |
| Convergence Reproducibility | 3 | ✓ PASS | ✓ PASS |
| Compiler Flag Detection | 4 | ✓ PASS | FAIL (expected) |
| Extended Precision (x87) | 3 | ✓ PASS | ✓ PASS |
| Scaling Sensitivity | 3 | ✓ PASS | ✓ PASS |
| Jacobian Computation | 3 | ✓ PASS | May vary |
| LM Parameter Sensitivity | 3 | ✓ PASS | ✓ PASS |
| **Total** | **45** | **45/45 PASS** | May FAIL |

### Known Hostile Behaviours

| Issue | Affected | Mitigation |
|-------|----------|------------|
| `-ffast-math` FTZ | **DENORM subnormals** | Subnormals flushed to zero. DENORM returns 0 instead of correct norm. **Do not use.** |
| `-ffast-math` | Bessel K functions | 2-4 ULP deviation. **Do not use.** |
| Aggressive optimisation | Some QUADPACK routines | May reorder operations. Test carefully. |
| ARM vs x86 FPU | Potentially all | Under investigation |

See [DEVIATIONS.md](DEVIATIONS.md#critical-deviation-subnormal-flush-with--ffast-math) for full analysis.

### Test Files

- `test/test_bessel_portability.f90` — Cross-platform ULP measurement
- `test/level4_hostile/test_l4_linear_blas.f90` — BLAS hostile tests (subnormals, Inf/NaN, SIMD edges)
- `test/level4_hostile/test_l4_minpack.f90` — MINPACK portability (detects FTZ)

---

## Module Coverage Summary

*Last updated: 1 January 2026*

| Module | Routines | L1 | L2 | L3 | L4 | Overall |
|--------|----------|----|----|----|----|---------|
| **service** | 3 | — | — | ✓ | — | ~33% |
| **special_functions** | 270 | ~20 | ~20 | 4 | ~20 | ~9% |
| **linear (BLAS)** | 217 | 18 | 16 | 9 | 65 | ~50% |
| **diff_integ** | 81 | ~10 | ~10 | — | — | ~12% |
| **diff_integ_eq** | 225 | ~5 | — | — | — | ~2% |
| **interpolation** | 80 | — | — | — | — | 0% |
| **integ_trans** | 48 | — | — | — | — | 0% |
| **approximation** | 78 | 9 | 17 | 7 | 45 | ~100% |
| **nonlin_eq** | 15 | — | — | — | — | 0% |
| **optimisation** | 46 | — | — | — | — | 0% |
| **data_handling** | 16 | — | — | — | — | 0% |
| **TOTAL** | **1,079** | ~62 | ~63 | ~22 | ~130 | **~18%** |

---

## MINPACK Test Details

*Testing the routines we just modernised (eliminated 27 GOTOs)*

### Complete 4-Level Results (1 January 2026)

| Level | Tests | Result | Notes |
|-------|-------|--------|-------|
| **L1 Regression** | 9 | ✓ **9/9 PASS** | DENORM (7) + ENORM (2) |
| **L2 Mathematical** | 17 | ✓ **17/17 PASS** | Rosenbrock, Powell, Helical, F-R, DENORM |
| **L3 Historical** | 7 | ✓ **7/7 PASS** | IBM 360 golden values |
| **L4 Hostile** | 45 | ✓ **45/45 PASS** | 13 categories, comprehensive portability |
| **L4 with -ffast-math** | 45 | May FAIL | Subnormal flush, associativity detected |
| **TOTAL** | **78** | **78/78 PASS** | (without hostile flags) |

### Test Files

| Level | File | Description |
|-------|------|-------------|
| L1 | `test/level1_regression/test_l1_minpack.f90` | Modern Fortran regression |
| L2 | `test/level2_mathematical/test_l2_minpack_mgh.f90` | MGH test functions |
| L3 | `test/level3_historical/test_l3_minpack.f90` | IBM 360 golden comparisons |
| L3 | `/c/dev/fortran360/tests/slatec/minpack/test_enorm.f` | FORTRAN IV for Hercules |
| L4 | `test/level4_hostile/test_l4_minpack.f90` | Portability stress tests |

### IBM 360 Historical Tests (Level 3)

| Test File | Routines | Problems | Status |
|-----------|----------|----------|--------|
| test_enorm.f | DENORM | 7 norm tests | ✓ **7/7 PASS** |
| test_qrfac.f | DQRFAC | 4 QR tests | ⏳ Pending |
| test_dnls1.f | DNLS1 | Rosenbrock, Powell, Freudenstein-Roth | ⏳ Pending |
| test_dnsq.f | DNSQ | Helical, Trigonometric, Broyden | ⏳ Pending |

---

## BLAS Test Details

*Testing core linear algebra routines (DAXPY, DSCAL, DCOPY, DSWAP, DROT, DROTG)*

### Complete 4-Level Results (1 January 2026)

| Level | Tests | Result | Notes |
|-------|-------|--------|-------|
| **L1 Regression** | 18 | ✓ **18/18 PASS** | DAXPY, DSCAL, DCOPY, DSWAP, DROT, DROTG |
| **L2 Mathematical** | 16 | ✓ **16/16 PASS** | Linearity, orthogonality, Pythagorean |
| **L3 Historical** | 9 | ✓ **9/9 PASS** | IBM 360 golden values (Pythagorean triples) |
| **L4 Hostile** | 65 | ✓ **65/65 PASS** | 16 categories, comprehensive portability |
| **TOTAL** | **108** | **108/108 PASS** | |

### Test Files

| Level | File | Description |
|-------|------|-------------|
| L1 | `test/level1_regression/test_l1_linear_blas.f90` | Modern Fortran regression |
| L2 | `test/level2_mathematical/test_l2_linear_blas.f90` | Mathematical property verification |
| L3 | `test/level3_historical/test_l3_linear_blas.f90` | IBM 360 golden comparisons |
| L3 | `/c/dev/fortran360/tests/slatec/blas/test_daxpy.f` | FORTRAN IV for Hercules |
| L3 | `/c/dev/fortran360/tests/slatec/blas/test_drotg.f` | FORTRAN IV for Hercules |
| L4 | `test/level4_hostile/test_l4_linear_blas.f90` | Hostile environment tests |

### Mathematical Properties Tested (Level 2)

| Property | Test |
|----------|------|
| Distributivity | (a+b)*X = a*X + b*X |
| Distributivity | a*(X+Y) = a*X + a*Y |
| Orthogonality | c² + s² = 1 for Givens rotations |
| Pythagorean | r = √(a² + b²) for Givens |
| Norm preservation | \|\|Rx\|\| = \|\|x\|\| for rotations |
| Triangle inequality | \|\|x+y\|\| ≤ \|\|x\|\| + \|\|y\|\| |

### Hostile Tests (Level 4) — 65 Tests in 16 Categories

| Category | Tests | What It Catches |
|----------|-------|-----------------|
| Subnormals | 3 | DAZ/FTZ from `-ffast-math` |
| Signed Zero | 3 | IEEE 754 ±0 handling |
| Inf/NaN | 3 | Special value propagation |
| Extreme Values | 4 | Near overflow/underflow |
| Accumulation | 3 | Precision loss in summations |
| SIMD Edges | 4 | Vectorization boundary bugs |
| Rounding Modes | 4 | All four IEEE rounding modes |
| FMA Detection | 3 | Fused multiply-add vs separate |
| Catastrophic Cancellation | 3 | Nearly-equal value subtraction |
| Associativity | 4 | Summation reordering, `-ffast-math` |
| Extended Precision (x87) | 4 | 80-bit register leakage |
| Reproducibility | 3 | Identical results across runs |
| Compiler Flag Detection | 4 | Automatic dangerous flag detection |
| NaN Variants | 4 | qNaN vs sNaN, payload preservation |
| ULP Accuracy | 4 | Mathematical precision deviation |
| Memory Alignment | 12 | Non-aligned, strided, reverse access |

---

## Contributing Tests

When adding tests, consider which level they belong to:

1. **Level 1**: Basic functionality. Does `DGEMM(A, B)` return something reasonable?
2. **Level 2**: Mathematical correctness. Does `DGEMM(I, A)` return `A`?
3. **Level 3**: Historical baseline. Does our `DGEMM` match IBM 360 output?
4. **Level 4**: Portability. Does `DGEMM` give the same answer on ARM vs x86?

Only Level 1 tests should be modified when code changes. Levels 2-4 are reference tests — if they fail, investigate before changing them.

---

## Deviations

All deviations from expected values are documented in [DEVIATIONS.md](DEVIATIONS.md), including:
- Which level detected the deviation
- Root cause analysis
- Whether it's a bug, a floating-point artefact, or expected behaviour
- Remediation (if any)

---

*"The purpose of testing is not to prove the code works. It's to find where it doesn't."*

— Someone who has clearly spent too long debugging numerical software
