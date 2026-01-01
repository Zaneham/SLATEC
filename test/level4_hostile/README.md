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

## Test File

| File | Purpose |
|------|---------|
| test_minpack_portability.f90 | Platform detection, ULP tests, edge cases, associativity |

## Running

```bash
# With safe flags
gfortran -O2 -o test_l4 test/level4_hostile/test_minpack_portability.f90
./test_l4

# With hostile flags (should fail associativity test)
gfortran -Ofast -ffast-math -o test_l4_hostile test/level4_hostile/test_minpack_portability.f90
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
