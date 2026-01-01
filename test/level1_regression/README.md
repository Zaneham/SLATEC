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

## MINPACK Coverage

| Test File | Routines | Tests |
|-----------|----------|-------|
| test_minpack.f90 | DENORM/ENORM | 9 tests |

## Running

```bash
cd /c/dev/slatec-modern
gfortran -o test_l1 test/level1_regression/test_minpack.f90
./test_l1
```

Or via fpm:

```bash
fpm test
```

---

*"If it compiles and the tests pass, ship it." — Every developer ever, immediately before a production incident*
