# Level 1: Regression Tests

*"Does the code work? If we change things, does it still work?"*

## Purpose

Level 1 tests verify basic functionality and catch regressions when code is modified. **This is the only level where we make unilateral changes** — if Level 1 fails, fix the code.

## What We Test

- Basic function calls work
- Expected outputs match
- Edge cases handled properly
- No crashes or exceptions

## Modifying These Tests

When code changes, Level 1 tests may need updating. This is expected and acceptable:

1. Run the tests
2. If they fail due to intentional code changes, update the tests
3. If they fail unexpectedly, the code has regressed — investigate

## Coverage

| Test File | Module | Routines | Tests | Status |
|-----------|--------|----------|-------|--------|
| test_l1_minpack.f90 | approximation | DENORM/ENORM | 9 | ✓ **9/9 PASS** |
| test_l1_linear_blas.f90 | linear | DAXPY, DSCAL, DCOPY, DSWAP, DROT, DROTG | 18 | ✓ **18/18 PASS** |
| test_l1_linear_linpack.f90 | linear | DGEFA/DGESL, DPOFA/DPOSL | — | ⏳ In Progress |
| test_l1_interpolation.f90 | interpolation | DBINT4, DBVALU, DPCHIM, DPCHFE, DPCHIA, DPLINT, DPOLCF, DPOLVL | 9 | ✓ **9/9 PASS** |
| test_l1_diff_integ.f90 | diff_integ | DGAUS8, DQAGS, DQAGI, DQNG | 12 | ✓ **12/12 PASS** |

## Running

```bash
cd /c/dev/slatec-modern

# Build library first (if needed)
cmake -B build && cmake --build build

# MINPACK
gfortran -O2 -I build/modules -o test_l1_minpack test/level1_regression/test_l1_minpack.f90 -L build/lib -lslatec
./test_l1_minpack

# BLAS
gfortran -O2 -I build/modules -o test_l1_blas test/level1_regression/test_l1_linear_blas.f90 -L build/lib -lslatec
./test_l1_blas

# Interpolation
gfortran -O2 -I build/modules -o test_l1_interpolation test/level1_regression/test_l1_interpolation.f90 -L build/lib -lslatec
./test_l1_interpolation

# Diff_Integ
gfortran -O2 -I build/modules -o test_l1_diff_integ test/level1_regression/test_l1_diff_integ.f90 -L build/lib -lslatec
./test_l1_diff_integ
```

Or via fpm:

```bash
fpm test
```

---

*"If it compiles and the tests pass, ship it." — Every developer ever, immediately before a production incident*
