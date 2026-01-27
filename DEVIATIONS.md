# SLATEC-Modern Floating-Point Deviation Report

## Purpose

This document records **exact** floating-point values produced by SLATEC-Modern
under different compiler configurations. This is critical for:

- **Reproducibility**: Scientific results must be reproducible across platforms
- **Validation**: Comparing against published reference values (A&S, NIST)
- **Audit trails**: Financial and regulatory compliance
- **Error bounds**: Understanding worst-case deviations for error propagation analysis

## Test Environment

- **Compiler**: gfortran 13.2.0
- **Platform**: Windows (MSYS2/MinGW-w64)
- **CPU**: x86-64 (specific model not recorded)
- **Date**: 2025-12-23

## Machine Constants

| Precision | Digits | Epsilon |
|-----------|--------|---------|
| SP (float32) | 6 | 1.192E-07 |
| DP (float64) | 15 | 2.220E-16 |

## Configurations Tested

1. **-O0**: No optimisation (baseline)
2. **-O2**: Standard optimisation
3. **-O3 -ffast-math**: Aggressive optimisation with relaxed IEEE compliance

## Deviation Summary

### Identical Across -O0, -O2, -O3 -ffast-math

The following functions produced **bit-identical** results across all configurations:

| Function | Precision | Value | Hex |
|----------|-----------|-------|-----|
| J_0(1.0) | SP | 0.7651976943016052 | 3F43E3FF |
| J_0(1.0) | DP | 0.7651976865579666 | 3FE87C7FDBD7B8F0 |
| J_0(2.0) | SP | 0.2238907665014267 | 3E65439F |
| J_0(2.0) | DP | 0.2238907791412356 | 3FCCA873FB24CEF6 |
| J_0(5.0) | SP | -0.1775967627763748 | BE35DBED |
| J_0(5.0) | DP | -0.1775967713143385 | BFC6BB7DB255CB92 |
| J_0(50.0) | SP | 0.0558124817907810 | 3D649BA1 |
| J_0(50.0) | DP | 0.0558123276692521 | 3FAC936EF41B2C72 |
| J_0(100.0) | SP | 0.0199856813997030 | 3CA3B903 |
| J_0(100.0) | DP | 0.0199858503042233 | 3F94772BB5C1EFA9 |
| J_1(1.0) | SP | 0.4400506019592285 | 3EE14E50 |
| J_1(1.0) | DP | 0.4400505857449336 | 3FDC29C9EE970C6D |
| I_0(1.0) | SP | 1.2660659551620483 | 3FA20E73 |
| I_0(1.0) | DP | 1.2660658777520082 | 3FF441CE4B386C2C |
| I_0(2.0) | SP | 2.2795851230621338 | 4011E4B9 |
| I_0(2.0) | DP | 2.2795853023360668 | 40023C97380FCE31 |
| I_1(1.0) | SP | 0.5651590824127197 | 3F10AE44 |
| I_1(1.0) | DP | 0.5651591039924850 | 3FE215C88B95E67E |
| K_1(1.0) | SP | 0.6019072532653809 | 3F1A1698 |
| K_1(1.0) | DP | 0.6019072301972346 | 3FE342D2F39D89C2 |
| Y_0(1.0) | DP | 0.0882569642156769 | 3FB6980226F358DE |
| Y_0(2.0) | DP | 0.5103756726497451 | 3FE054FF5CD68C8D |

### DEVIATIONS DETECTED

The following functions produced **different** results with -O3 -ffast-math:

#### K_0(1.0) - Single Precision

| Config | Decimal Value | Hex | Deviation from -O0 |
|--------|--------------|-----|-------------------|
| -O0 | 0.4210244417190552 | 3ED79084 | (baseline) |
| -O2 | 0.4210244417190552 | 3ED79084 | 0 |
| -O3 -ffast-math | 0.4210243821144104 | **3ED79082** | **5.96E-08** |

**Analysis**: 2 ULPs difference in SP. Hex differs by 2 in last digit.

#### K_0(1.0) - Double Precision

| Config | Decimal Value | Hex | Deviation from -O0 |
|--------|--------------|-----|-------------------|
| -O0 | 0.4210244382407085 | 3FDAF2107C43E11C | (baseline) |
| -O2 | 0.4210244382407085 | 3FDAF2107C43E11C | 0 |
| -O3 -ffast-math | 0.4210244382407082 | **3FDAF2107C43E118** | **3E-16** |

**Analysis**: 4 ULPs difference in DP. Hex differs by 4 in last digit (0x1C vs 0x18).

#### K_0(2.0) - Single Precision

| Config | Decimal Value | Hex | Deviation from -O0 |
|--------|--------------|-----|-------------------|
| -O0 | 0.1138938963413239 | 3DE94134 | (baseline) |
| -O2 | 0.1138938963413239 | 3DE94134 | 0 |
| -O3 -ffast-math | 0.1138938665390015 | **3DE94130** | **2.98E-08** |

**Analysis**: 4 ULPs difference in SP. Hex differs by 4 in last digit.

#### K_0(2.0) - Double Precision

| Config | Decimal Value | Hex | Deviation from -O0 |
|--------|--------------|-----|-------------------|
| -O0 | 0.1138938727495334 | 3FBD28261AAC8D54 | (baseline) |
| -O2 | 0.1138938727495334 | 3FBD28261AAC8D54 | 0 |
| -O3 -ffast-math | 0.1138938727495334 | **3FBD28261AAC8D58** | **~1E-16** |

**Analysis**: 4 ULPs difference in DP. Hex differs by 4 in last digit (0x54 vs 0x58).

## Large Argument Bessel Functions - Verified Correct

SLATEC's Bessel function implementations have been verified against mpmath (arbitrary
precision arithmetic, 100+ digits) for large arguments:

| Function | SLATEC Computed | mpmath (100 digit) | Status |
|----------|-----------------|---------------------|--------|
| J_0(500) | -0.0341005568807320 | -0.0341005568807320 | ✓ Correct |
| J_0(1000) | +0.0247866861524202 | +0.0247866861524202 | ✓ Correct |
| J_9(20) | 0.1251262546479942 | 0.1251262546479942 | ✓ Correct |
| Y_0(50) | -0.0980649954700771 | -0.0980649954700771 | ✓ Correct |
| Y_0(100) | -0.0772443133650832 | -0.0772443133650832 | ✓ Correct |

**Note:** An earlier version of this document contained incorrect reference values
that were mistakenly attributed to NIST. The SLATEC implementations are accurate.

---

## Comparison to Reference Values

### A&S Table 9.1 (Bessel J)

| Function | A&S Reference | SLATEC DP | Abs Error | Rel Error |
|----------|---------------|-----------|-----------|-----------|
| J_0(1.0) | 0.7651976865579666 | 0.7651976865579666 | 0 | 0 |
| J_0(2.0) | 0.2238907791412357 | 0.2238907791412356 | 1E-16 | 4E-16 |
| J_0(5.0) | -0.1775967713143383 | -0.1775967713143385 | 2E-16 | 1E-15 |
| J_1(1.0) | 0.4400505857449335 | 0.4400505857449336 | 1E-16 | 2E-16 |

### A&S Table 9.8 (Bessel I, K)

| Function | A&S Reference | SLATEC DP | Abs Error | Rel Error |
|----------|---------------|-----------|-----------|-----------|
| I_0(1.0) | 1.2660658777520084 | 1.2660658777520082 | 2E-16 | 2E-16 |
| I_0(2.0) | 2.2795853023360673 | 2.2795853023360668 | 5E-16 | 2E-16 |
| I_1(1.0) | 0.5651591039924851 | 0.5651591039924850 | 1E-16 | 2E-16 |
| K_0(1.0) | 0.4210244382407084 | 0.4210244382407085 | 1E-16 | 2E-16 |
| K_0(2.0) | 0.1138938727495334 | 0.1138938727495334 | 0 | 0 |
| K_1(1.0) | 0.6019072301972346 | 0.6019072301972346 | 0 | 0 |

## Impact Assessment

-- Please note: I am not an expert in this field. Please exercise your own judgement and caution. 

### For Physics/Mathematics

The maximum DP deviation observed (4 ULPs in K_0) is:
- **Acceptable** for most numerical work
- **May accumulate** over 10^6+ iterations
- **Recommend** using -O2 without -ffast-math for reproducibility

### For Finance

With -O3 -ffast-math:
- K_0(1.0) SP deviation: 5.96E-08 = **$59.60 per $1B**
- K_0(2.0) SP deviation: 2.98E-08 = **$29.80 per $1B**

**Recommendation**: Never use -ffast-math for financial calculations.

### For Reproducibility

- -O0 and -O2 produce **bit-identical** results
- -O3 -ffast-math introduces measurable deviations in Bessel K functions
- Cross-platform reproducibility requires:
  1. Same compiler version
  2. Same optimisation flags
  3. Verification of hex values against this document

## Recommendations

1. **Default build**: Use -O2 (bit-identical to -O0, faster)
2. **Never use -ffast-math** for numerical libraries
3. **Document compiler/flags** in any publication using these functions
4. **Verify critical results** by checking hex values against reference

---

## QUADPACK Integration Routines

QUADPACK integration routines (QAGS, DQAGS, QAGI, DQAGI, QNG, DQNG) were tested
against integrals with known closed-form solutions.

### Reference Integrals Tested

| Integral | Exact Value | Description |
|----------|-------------|-------------|
| ∫₀¹ x² dx | 1/3 | Polynomial |
| ∫₀^π sin(x) dx | 2 | Trigonometric |
| ∫₀¹ exp(x) dx | e - 1 | Exponential |
| ∫₀¹ 1/(1+x²) dx | π/4 | Arctan |
| ∫₀^∞ exp(-x²) dx | √π/2 | Gaussian (infinite) |

### QUADPACK Results: Bit-Identical Across All Configurations

**QUADPACK produces bit-identical results across -O2 and -O3 -ffast-math.**

This is in contrast to Bessel K_0 which shows 2-4 ULP deviations with -ffast-math.

| Integral | Prec | Computed | Hex | Rel Error |
|----------|------|----------|-----|-----------|
| x² | DP | 0.3333333333333334 | 3FD5555555555556 | 1.67E-16 |
| x² | SP | 0.3333333433 | 3EAAAAAB | 0 |
| sin(x) | DP | 2.0000000000000000 | 4000000000000000 | 0 |
| sin(x) | SP | 2.0000002384 | 40000001 | 1.19E-07 |
| exp(x) | DP | 1.7182818284590453 | 3FFB7E151628AED3 | 1.29E-16 |
| exp(x) | SP | 1.7182817459 | 3FDBF0A8 | 0 |
| 1/(1+x²) | DP | 0.7853981633974484 | 3FE921FB54442D19 | 1.41E-16 |
| 1/(1+x²) | SP | 0.7853981256 | 3F490FDA | 7.59E-08 |
| exp(-x²) | DP | 0.8862269254527579 | 3FEC5BF891B4EF6A | 0 |
| exp(-x²) | SP | 0.8862268925 | 3F62DFC4 | 6.73E-08 |

### QUADPACK Accuracy Assessment

- **Double precision**: All integrals computed to within 2 ULPs of exact value
- **Single precision**: All integrals computed to within 1-2 ULPs of exact value
- **No optimisation-dependent deviations** detected
- **Infinite interval integration** (DQAGI) works correctly for Gaussian integral

### QUADPACK Recommendations

1. QUADPACK is **safe to use** with any optimisation level, including -ffast-math
2. Results are reproducible across optimisation configurations
3. For critical applications, verify hex values match expected results

## Raw Data Files

- `outputs/deviation_O0.txt` - Full report with -O0
- `outputs/deviation_O2.txt` - Full report with -O2
- `outputs/deviation_O3_fastmath.txt` - Full report with -O3 -ffast-math
- `outputs/quadpack_O2.txt` - QUADPACK deviation report with -O2
- `outputs/quadpack_O3_fastmath.txt` - QUADPACK deviation report with -O3 -ffast-math

## References

1. Abramowitz, M. & Stegun, I.A. (1964). Handbook of Mathematical Functions. NBS.
2. NIST Digital Library of Mathematical Functions. https://dlmf.nist.gov/
3. IEEE 754-2019. Standard for Floating-Point Arithmetic.

---

## MINPACK (Approximation Module)

*Tested: 1 January 2026*

### DENORM (Double Precision Euclidean Norm)

DENORM computes the Euclidean norm using a three-accumulator algorithm designed to prevent overflow and underflow. It partitions input values into three ranges:
- **Large components** (> rgiant/n): Scaled accumulation to prevent overflow
- **Intermediate components**: Direct accumulation
- **Small components** (< rdwarf): Scaled accumulation to prevent underflow

#### Results Across Configurations

| Test Case | -O0 | -O2 | -O3 -ffast-math |
|-----------|-----|-----|-----------------|
| Unit vector | ✓ 1.0 | ✓ 1.0 | ✓ 1.0 |
| Ones vector (√10) | ✓ 3.162277660168... | ✓ identical | ✓ identical |
| 3-4-5 triangle | ✓ 5.0 | ✓ 5.0 | ✓ 5.0 |
| Large (10¹⁵) | ✓ 1.732...×10¹⁵ | ✓ identical | ✓ identical |
| Small (10⁻¹⁵) | ✓ 1.732...×10⁻¹⁵ | ✓ identical | ✓ identical |
| Mixed (10⁻¹⁸ to 10¹⁸) | ✓ 1.0×10¹⁸ | ✓ identical | ✓ identical |
| Zero vector | ✓ 0.0 | ✓ 0.0 | ✓ 0.0 |
| **Subnormal (10⁻³⁰⁸)** | ✓ 2.49×10⁻³⁰⁸ | ✓ identical | **✗ 0.0** |

### CRITICAL DEVIATION: Subnormal Flush with -ffast-math

**Detected by**: Level 4 Hostile Test  
**Severity**: Catastrophic (silent data corruption)  
**Affected**: Any computation involving subnormal (denormalised) floating-point numbers

#### The Problem

When compiled with `-ffast-math`, GCC enables **Flush-To-Zero (FTZ)** mode. FTZ causes the CPU to silently replace subnormal numbers with zero.

#### Observed Behaviour

```
Input vector: x = (1.11×10⁻³⁰⁸, 1.11×10⁻³⁰⁸, 1.11×10⁻³⁰⁸, 1.11×10⁻³⁰⁸, 1.11×10⁻³⁰⁸)

Expected:     ||x|| = √5 × 1.11×10⁻³⁰⁸ ≈ 2.49×10⁻³⁰⁸
Without -ffast-math:  DENORM = 2.48770820×10⁻³⁰⁸  ✓
With -ffast-math:     DENORM = 0.0                 ✗
```

#### Root Cause Analysis

1. The subnormal value `tiny(1.0_dp) / 2.0_dp = 1.11×10⁻³⁰⁸` is created
2. With FTZ enabled, the CPU flushes this to zero before any arithmetic
3. DENORM receives a vector of zeros and correctly returns zero
4. The algorithm's underflow protection is **bypassed at the hardware level**

DENORM's three-accumulator design handles values down to `rdwarf = 3.834×10⁻²⁰`, but subnormals (below `tiny ≈ 2.225×10⁻³⁰⁸`) are flushed before DENORM ever sees them.

#### Hex Analysis

| Compiler Flags | Subnormal Value | Hex Representation |
|----------------|-----------------|-------------------|
| -O2 | 1.112536929253601×10⁻³⁰⁸ | 0x000FFFFFFFFFFFFF |
| -O3 -ffast-math | 0.0 (flushed) | 0x0000000000000000 |

#### Impact Assessment

**Scientific Computing**: Simulations involving very small quantities (e.g., quantum mechanics, particle physics) may silently produce zeros instead of small values.

**Numerical Linear Algebra**: Condition number estimation, null space computation, and iterative refinement may fail silently.

**Optimisation Algorithms**: MINPACK's Levenberg-Marquardt algorithm uses DENORM for step size computation. Subnormal Jacobian entries would cause premature termination.

#### Mitigation

1. **Never use `-ffast-math`** for numerical code
2. If aggressive optimisation is required, use `-O3` without `-ffast-math`
3. Test with Level 4 hostile tests before deployment

### Comparison to IBM 360 Historical Values (Level 3)

DENORM produces **bit-identical** results to IBM System/360 (Hercules/TK4-/FORTRAN G) for all test cases within IEEE 754 representable range:

| Test Case | IBM 360 (Hex FP) | Modern (IEEE 754) | Match |
|-----------|------------------|-------------------|-------|
| Unit vector | 0.1000000000000D+01 | 1.0000000000000000 | ✓ |
| √10 | 0.3162277660168D+01 | 3.1622776601683795 | ✓ |
| 3-4-5 | 0.5000000000000D+01 | 5.0000000000000000 | ✓ |
| Large | 0.1732050807569D+16 | 1.7320508075688772×10¹⁵ | ✓ |
| Small | 0.1732050807569D-14 | 1.7320508075688772×10⁻¹⁵ | ✓ |
| Mixed | 0.1000000000000D+19 | 1.0000000000000000×10¹⁸ | ✓ |
| Zero | 0 | 0.0 | ✓ |

Note: IBM 360 uses hexadecimal floating-point (base 16) with different mantissa length. Values match within representation limits of both formats.

---

## BLAS (Linear Algebra Primitives)

*Tested: 1 January 2026*

### Level 4 Hostile Test Results

BLAS routines (DAXPY, DSCAL, DCOPY, DSWAP, DROT, DROTG) were subjected to 65 hostile environment tests across 16 categories.

#### Results by Compiler Flag

| Flag | Passed | Failed | Failure Categories |
|------|--------|--------|-------------------|
| `-O2` | 65/65 | 0 | None |
| `-O3` | 65/65 | 0 | None |
| `-O3 -march=native` | 65/65 | 0 | None |
| `-O2 -flto` | 65/65 | 0 | None |
| `-O2 -ffast-math` | 59/65 | 6 | Subnormals (2), Inf/NaN (2), NaN variants (2) |
| `-Ofast` | 59/65 | 6 | Subnormals (2), Inf/NaN (2), NaN variants (2) |
| `-O2 -ffinite-math-only` | 61/65 | 4 | Inf/NaN (2), NaN variants (2) |
| `-O2 -funsafe-math-optimizations` | 63/65 | 2 | Subnormals (2) |

### CRITICAL DEVIATION: DAXPY Subnormal Flush

**Detected by**: Level 4 Subnormal Handling Test
**Severity**: Catastrophic (silent data corruption)
**Affected**: Any DAXPY computation with subnormal inputs or results

#### Observed Behaviour

```fortran
! With -O2 (correct)
call daxpy(n, 1.0_dp, subnormal_x, 1, y, 1)
! y = subnormal values preserved

! With -ffast-math (incorrect)
call daxpy(n, 1.0_dp, subnormal_x, 1, y, 1)
! y = 0.0 (subnormals flushed to zero)
```

#### Impact

DAXPY is the fundamental operation `Y = alpha*X + Y` used throughout numerical linear algebra. Subnormal flush affects:

1. **Iterative refinement** — Small corrections lost
2. **Null space computation** — Near-zero vectors corrupted
3. **Gradient descent** — Small gradients vanish
4. **Condition estimation** — Small singular values lost

### CRITICAL DEVIATION: DAXPY NaN Non-Propagation

**Detected by**: Level 4 Inf/NaN Propagation Test
**Severity**: High (incorrect results go undetected)
**Affected**: Error detection in iterative algorithms

#### Observed Behaviour

```fortran
! Input with NaN
x = [1.0, NaN, 3.0]
call daxpy(3, 1.0_dp, x, 1, y, 1)

! With -O2: y contains NaN (correct - propagates error)
! With -ffast-math: y = finite values (NaN silently dropped)
```

#### Impact

NaN values typically indicate:
- Division by zero
- 0 * Inf
- Inf - Inf
- sqrt of negative number

When NaN doesn't propagate, errors go undetected and produce plausible-looking but wrong results.

### CRITICAL DEVIATION: Infinity Arithmetic

**Detected by**: Level 4 Inf/NaN Propagation Test
**Severity**: High
**Affected**: Overflow detection

#### Observed Behaviour

```fortran
! With -O2 (IEEE 754 compliant)
Inf - Inf = NaN  ✓
0 * Inf = NaN    ✓

! With -ffast-math (non-compliant)
Inf - Inf = undefined (not NaN)
0 * Inf = undefined (not NaN)
```

### Additional Warnings (Not Failures)

The following were detected but classified as warnings, not failures:

| Test | -O2 Result | -ffast-math Result | Status |
|------|------------|-------------------|--------|
| Signed zero distinction | +0 ≠ -0 | +0 = -0 | WARN |
| 1/(+0) vs 1/(-0) | +Inf vs -Inf | Same | WARN |
| (1+1e-15) - 1 | 1.11e-16 | 1.11e-16 | OK (same) |
| Associativity | Violated (IEEE correct) | Violated | OK |

### BLAS Safe Usage Summary

| Routine | Safe Flags | Dangerous Flags | Notes |
|---------|------------|-----------------|-------|
| DAXPY | `-O2`, `-O3` | `-ffast-math` | Subnormal, NaN issues |
| DSCAL | `-O2`, `-O3` | `-ffast-math` | Subnormal issues |
| DCOPY | Any | None | Pure memory copy |
| DSWAP | Any | None | Pure memory swap |
| DROT | `-O2`, `-O3` | `-ffast-math` | May affect precision |
| DROTG | `-O2`, `-O3` | `-ffast-math` | Extreme value handling |

---

## Comprehensive Flag Impact Matrix

*What exactly breaks with each dangerous flag?*

| Flag | Subnormals | Signed Zero | Inf | NaN | Associativity | Extended Precision |
|------|------------|-------------|-----|-----|---------------|-------------------|
| `-ffast-math` | FLUSH | Merged | Broken | Broken | May reorder | OK |
| `-ffinite-math-only` | OK | OK | Broken | Broken | OK | OK |
| `-funsafe-math-optimizations` | FLUSH | May merge | OK | OK | May reorder | May leak |
| `-fno-signed-zeros` | OK | Merged | OK | OK | OK | OK |
| `-fassociative-math` | OK | OK | OK | OK | Reorders | OK |
| `-freciprocal-math` | OK | OK | OK | OK | OK | OK |
| `-fno-trapping-math` | OK | OK | OK | OK | OK | OK |

### Legend

- **FLUSH**: Values flushed to zero (silent data corruption)
- **Merged**: +0 and -0 treated as identical
- **Broken**: IEEE 754 semantics violated
- **Reorders**: (a+b)+c may not equal a+(b+c)
- **May leak**: x87 80-bit precision may affect results
- **OK**: Behaviour unchanged from strict IEEE 754
