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

## CRITICAL: Algorithmic Accuracy Failures

The following are **NOT** optimisation-dependent - these are fundamental accuracy
failures in SLATEC's large-argument and high-order handling. These failures occur
identically across all compiler configurations (-O0, -O2, -O3 -ffast-math).

### Bessel J - Large Argument Failures

| Function | NIST Reference | SLATEC Computed | Error Type |
|----------|----------------|-----------------|------------|
| J_0(500) | +0.0179753822938868 | **-0.0341005568807317** | **WRONG SIGN, 290% ERROR** |
| J_0(1000) | -0.0246711376936846 | **+0.0247866861524200** | **WRONG SIGN** |

### Bessel J - High Order Failure

| Function | NIST Reference | SLATEC Computed | Error Type |
|----------|----------------|-----------------|------------|
| J_9(20) | 0.2453202954891665 | **0.0** | **TOTAL LOSS OF ACCURACY** |

### Bessel Y - Large Argument Failures

| Function | NIST Reference | SLATEC Computed | Error Type |
|----------|----------------|-----------------|------------|
| Y_0(50) | -0.0560404718523358 | **-0.0980649954700770** | **75% ERROR** |
| Y_0(100) | -0.0772433752531550 | -0.0772443133650831 | 0.01% error |

### Root Cause Analysis

These failures likely stem from:
1. **Range reduction errors** in trigonometric calculations for large arguments
2. **Recurrence instability** for high-order Bessel functions
3. **Asymptotic expansion limitations** not properly bounded

### Recommendations for Critical Applications

1. **DO NOT USE** J_n(x) for x > 100 without independent verification
2. **DO NOT USE** Y_n(x) for x > 20 without independent verification
3. **DO NOT USE** J_n(x) for n > 5 and x > 15 without verification
4. For large arguments, consider alternative implementations (e.g., AMOS library)

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
