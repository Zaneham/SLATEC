# Level 4: Hostile/Portability Tests

*"Do compilers, processors, or operating systems change anything?"*

## Purpose

Level 4 tests detect platform-specific behaviour. The same code compiled with different compilers, on different processors, under different operating systems, should produce the same results (within floating-point tolerance). When it doesn't, we need to know.

**Do not change Level 4 tests — they detect platform bugs, not code bugs.**

## Compiler Matrix

| Compiler | Platform | Status |
|----------|----------|--------|
| gfortran 13.x | Linux (x86_64) | CI tested |
| gfortran 13.x | macOS (ARM64) | CI tested |
| gfortran 13.x | Windows (MSYS2) | CI tested |
| ifort/ifx | Linux | Untested |
| flang | Linux | Untested |
| nvfortran | Linux | Untested |

## Known Hostile Behaviours

| Issue | Affected | Mitigation |
|-------|----------|------------|
| `-ffast-math` | Bessel K, associativity, subnormals | **Do not use** |
| `-Ofast` | All numerical code | Use `-O3` instead |
| `/fp:fast` (MSVC) | All numerical code | Use `/fp:precise` |
| FTZ/DAZ modes | Subnormal handling | Disable in MXCSR |
| Aggressive FMA | Summation order | Usually acceptable |
| ARM vs x86 FPU | Potentially all | Under investigation |
| SIMD vectorisation | Reduction order | Usually 1-2 ULP |
| x87 extended precision | Intermediate results | Use SSE/AVX |

## What We Test

1. **ULP Deviation** — How many floating-point values apart are results?
2. **Edge Cases** — Subnormals, near-overflow, extreme ranges
3. **Associativity** — Does (a+b)+c equal a+(b+c)? (Detects `-ffast-math`)
4. **Signed Zero** — IEEE 754 ±0 handling
5. **Inf/NaN Propagation** — Special value handling
6. **SIMD Boundaries** — Vectorization edge cases
7. **Rounding Modes** — Sensitivity to ieee_nearest, ieee_up, ieee_down, ieee_to_zero
8. **FMA Detection** — Fused multiply-add vs separate operations
9. **Catastrophic Cancellation** — Precision loss in subtraction
10. **Extended Precision** — x87 80-bit register leakage
11. **Reproducibility** — Identical results across multiple runs
12. **Compiler Flag Detection** — Automatic detection of dangerous flags

## Coverage

| File | Module | Tests | Status |
|------|--------|-------|--------|
| test_minpack_portability.f90 | approximation | 45 | ✓ **45/45 PASS** |
| test_linear_blas_hostile.f90 | linear | 65 | ✓ **65/65 PASS** |

**Total: 110 Level 4 tests**

### MINPACK Hostile Test Categories

| Category | Tests | What It Detects |
|----------|-------|-----------------|
| ULP Precision | 4 | Floating-point accuracy vs mathematical constants |
| Edge Cases | 5 | Subnormals, near-overflow, extreme ranges |
| Associativity | 4 | `-ffast-math`, reordering optimizations |
| Rounding Modes | 4 | Sensitivity to rounding mode changes |
| FMA Effects | 3 | Fused multiply-add vs separate ops |
| Catastrophic Cancellation | 3 | Precision loss in norm computation |
| Condition Number | 3 | Sensitivity to ill-conditioned inputs |
| Convergence Reproducibility | 3 | Identical results across runs |
| Compiler Flag Detection | 4 | Automatic `-ffast-math`, FTZ detection |
| Extended Precision | 3 | x87 80-bit precision leakage |
| Scaling Sensitivity | 3 | Homogeneity preservation |
| Jacobian Computation | 3 | Finite difference accuracy |
| LM Parameter Sensitivity | 3 | Trust region regularization effects |

### BLAS Hostile Test Categories

| Category | Tests | What It Detects |
|----------|-------|-----------------|
| Subnormals | 3 | DAZ/FTZ modes, `-ffast-math` flush |
| Signed Zero | 3 | IEEE 754 ±0 compliance |
| Inf/NaN | 3 | Special value propagation |
| Extreme Values | 4 | Near overflow/underflow behaviour |
| Accumulation | 3 | Precision loss in long operations |
| SIMD Edges | 4 | Odd lengths, non-unit strides, negative strides |
| Rounding Modes | 4 | All four IEEE rounding modes |
| FMA Detection | 3 | FMA vs separate multiply-add |
| Catastrophic Cancellation | 3 | Nearly-equal value subtraction |
| Associativity | 4 | Summation order, -ffast-math |
| Extended Precision | 4 | x87 80-bit register effects |
| Reproducibility | 3 | Identical results across runs |
| Compiler Flags | 4 | Automatic dangerous flag detection |
| NaN Variants | 4 | qNaN vs sNaN, payload preservation |
| ULP Accuracy | 4 | Deviation from mathematical truth |
| Memory Alignment | 12 | Non-aligned, strided, reverse access |

## Running

```bash
cd /c/dev/slatec-modern

# MINPACK - Safe flags
gfortran -O2 -o test_l4_minpack test/level4_hostile/test_minpack_portability.f90
./test_l4_minpack

# BLAS - Safe flags
gfortran -O2 -o test_l4_blas test/level4_hostile/test_linear_blas_hostile.f90
./test_l4_blas

# With hostile flags (EXPECTED TO FAIL - validates detection)
gfortran -Ofast -ffast-math -o test_l4_hostile test/level4_hostile/test_linear_blas_hostile.f90
./test_l4_hostile
```

## Interpreting Failures

If Level 4 fails but Level 3 passes:
- Historical baseline from IBM 360 is correct
- Current platform differs from expectation
- Options:
  1. Document as platform-specific behaviour
  2. Add `#ifdef` for platform-specific code
  3. File compiler bug report

### Common Failure Patterns

| Failure | Likely Cause | Fix |
|---------|--------------|-----|
| Subnormal flush | `-ffast-math` or FTZ mode | Remove flag or set MXCSR |
| Associativity | `-ffast-math` reordering | Use `-O3` instead of `-Ofast` |
| NaN not propagating | `/fp:fast` or `-ffinite-math-only` | Use strict FP mode |
| Signed zero wrong | Non-IEEE optimization | Use `/fp:precise` |
| Reproducibility | Vectorization variance | Add `-ffp-contract=off` |
| Extended precision | x87 register spill | Use `-mfpmath=sse` |

---

*"It worked on my machine." — Famous last words*
