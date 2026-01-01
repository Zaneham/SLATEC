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
| `-ffast-math` | Bessel K, associativity | **Do not use** |
| `-Ofast` | All numerical code | Use `-O3` instead |
| `/fp:fast` (MSVC) | All numerical code | Use `/fp:precise` |
| Aggressive FMA | Summation order | Usually acceptable |
| ARM vs x86 FPU | Potentially all | Under investigation |
| SIMD vectorisation | Reduction order | Usually 1-2 ULP |

## What We Test

1. **ULP Deviation** — How many floating-point values apart are results?
2. **Edge Cases** — Subnormals, near-overflow, extreme ranges
3. **Associativity** — Does (a+b)+c equal a+(b+c)? (Detects `-ffast-math`)
4. **Signed Zero** — IEEE 754 ±0 handling
5. **Inf/NaN Propagation** — Special value handling
6. **SIMD Boundaries** — Vectorization edge cases

## Coverage

| File | Module | Tests | Status |
|------|--------|-------|--------|
| test_minpack_portability.f90 | approximation | 9 | ✓ **9/9 PASS** |
| test_linear_blas_hostile.f90 | linear | 20 | ✓ **20/20 PASS** |

### BLAS Hostile Test Categories

| Category | Tests | What It Detects |
|----------|-------|-----------------|
| Subnormals | 3 | DAZ/FTZ modes, `-ffast-math` flush |
| Signed Zero | 3 | IEEE 754 ±0 compliance |
| Inf/NaN | 3 | Special value propagation |
| Extreme Values | 4 | Near overflow/underflow behaviour |
| Accumulation | 3 | Precision loss in long operations |
| SIMD Edges | 4 | Odd lengths, non-unit strides, negative strides |

## Running

```bash
cd /c/dev/slatec-modern

# MINPACK - Safe flags
gfortran -O2 -o test_l4_minpack test/level4_hostile/test_minpack_portability.f90
./test_l4_minpack

# BLAS - Safe flags
gfortran -O2 -o test_l4_blas test/level4_hostile/test_linear_blas_hostile.f90
./test_l4_blas

# With hostile flags (may fail subnormal tests)
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

---

*"It worked on my machine." — Famous last words*
