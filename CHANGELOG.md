# CHANGELOG

All notable changes to SLATEC-Modern are documented in this file.

This project adheres to [Semantic Versioning](https://semver.org/).

---

## [1.0.0] - 2025-12-27

### Summary

First stable release of SLATEC-Modern: a complete modernisation of the SLATEC Common Mathematical Library (Version 4.1, July 1993) from FORTRAN 77 to Fortran 2018.

**Original Library:** 735 files, 168,216 lines of FORTRAN 77
**Modernised Library:** 34 modules, 1,079 routines, ~85,000 lines of Fortran 2018

---

### Transformations from Original SLATEC

The following systematic transformations were applied to the entire codebase:

| Transformation | Description |
|----------------|-------------|
| **Free-form source** | Converted from fixed-form (columns 1-72) to free-form Fortran 2018 |
| **GOTO elimination** | Replaced ~1,400 GOTO statements with structured control flow |
| **Module structure** | Organised routines into logical modules with explicit interfaces |
| **Arithmetic IF removal** | Replaced three-way arithmetic IFs with IF-THEN-ELSE constructs |
| **INTENT attributes** | Added INTENT(IN/OUT/INOUT) to all procedure arguments |
| **DATA → PARAMETER** | Replaced DATA statements with named PARAMETER constants |
| **IMPLICIT NONE** | Enforced explicit typing throughout |
| **Modern error handling** | Replaced XERMSG calls with ERROR STOP |

#### GOTO Elimination Patterns

The following patterns were used to replace GOTO statements whilst preserving the original algorithm semantics:

**State Machine (for complex multi-path control flow):**
```fortran
INTEGER, PARAMETER :: ST_INIT = 1, ST_COMPUTE = 2, ST_DONE = 99
state_loop: DO WHILE (istate /= ST_DONE)
  SELECT CASE (istate)
    CASE (ST_INIT)
      ! initialisation
      istate = ST_COMPUTE
    CASE (ST_COMPUTE)
      ! computation
      istate = ST_DONE
  END SELECT
END DO state_loop
```

**Logical Flags (for forward jumps):**
```fortran
LOGICAL :: skip_section
IF (condition) skip_section = .TRUE.
IF (.NOT. skip_section) THEN
  ! code that was previously skipped via GOTO
END IF
```

**Named DO Loops (for backward jumps/retry logic):**
```fortran
retry_loop: DO
  ! computation
  IF (converged) EXIT retry_loop
END DO retry_loop
```

**BLOCK/EXIT (for common exit points):**
```fortran
cleanup: BLOCK
  IF (error) EXIT cleanup
  ! normal processing
END BLOCK cleanup
! common cleanup code here
```

---

### Module Organisation

The modernised library is organised into the following modules:

| Module | Contents | Original Source |
|--------|----------|-----------------|
| `service` | Machine constants, utility routines | D1MACH, R1MACH, I1MACH |
| `special_functions` | Bessel, Gamma, Airy, Error functions, etc. | FNLIB, Amos algorithms |
| `linear` | BLAS, LINPACK, SLAP sparse solvers | BLAS, LINPACK, SLAP |
| `diff_integ` | Numerical integration (QUADPACK) | QUADPACK |
| `diff_integ_eq` | ODE/DAE solvers (DASSL, DEPAC, SDRIVE) | DASSL, DEPAC, SDRIVE |
| `interpolation` | Splines, PCHIP, polynomial interpolation | BSPLINE, PCHIP |
| `integ_trans` | Fast Fourier Transforms | FFTPACK |
| `approximation` | Curve fitting, least squares (MINPACK) | MINPACK |
| `nonlin_eq` | Nonlinear equation solvers | FZERO, SOSEQ |
| `optimization` | Linear/quadratic programming | SPLP |
| `data_handling` | Sorting, permutation | SLATEC utilities |

---

### Build System

#### CMake (Recommended)

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
ctest --output-on-failure
```

#### Fortran Package Manager (fpm)

```bash
fpm build --profile release
fpm test
```

#### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `SLATEC_BUILD_TESTS` | ON | Build test suite |
| `SLATEC_BUILD_SHARED` | OFF | Build shared library (experimental) |
| `SLATEC_USE_EXTERNAL_LAPACK` | OFF | Link against system LAPACK |

---

### Test Suite

The test suite validates correctness against published reference values:

| Test File | Coverage |
|-----------|----------|
| `test_bessel_golden.f90` | J, I, K Bessel functions against Abramowitz & Stegun Tables 9.1, 9.8 |
| `test_bessel_reference.f90` | Extended Bessel function validation |
| `test_bessel_portability.f90` | Cross-platform ULP error measurement |
| `test_bessel_regression.f90` | Regression tests against stored values |
| `test_diff_integ_quadpack.f90` | QUADPACK integration against closed-form solutions |
| `test_diff_integ_eq.f90` | ODE solver validation |
| `test_specfun_bessel.f90` | Special function edge cases |

**Golden Test Methodology:**

Tests compare computed values against:
1. Abramowitz & Stegun (1964) handbook tables
2. NIST Digital Library of Mathematical Functions (dlmf.nist.gov)
3. Known closed-form solutions (e.g., ∫₀¹ x² dx = 1/3)

All comparisons include relative error bounds and ULP (Units in Last Place) analysis.

---

### Numerical Accuracy

See `DEVIATIONS.md` for detailed floating-point analysis.

**Summary:**
- Double precision routines match A&S reference values to within 2 ULPs
- QUADPACK produces bit-identical results across optimisation levels
- Bessel K functions show 2-4 ULP deviation with `-ffast-math` (not recommended)

**Known Limitations:**
- J Bessel functions lose accuracy for arguments > 100
- High-order Bessel functions (n > 5) may fail for large arguments
- See `DEVIATIONS.md` for full analysis and recommendations

---

### Bugs Fixed During Modernisation

The following bugs were discovered and fixed during the modernisation process:

| Bug | Severity | Affected Routines | Resolution |
|-----|----------|-------------------|------------|
| DBESK numerical drift | High | dbsknu, besknu | Missing flag causing value overwrite |
| Missing BLOCK statement | Medium | dcov, scov | Added missing END BLOCK |
| Unclosed DO loop | Medium | sort routines | Fixed loop termination |
| Undeclared variable | Low | pchid, pchia | Added declaration |
| Mixed kind specifiers | Low | dpsifn, dasyjy | Corrected literal kinds |

---

### What Has NOT Changed

The following aspects of SLATEC remain unchanged:

1. **Algorithm implementations** - The mathematical algorithms are identical to the original
2. **Routine names** - All original routine names are preserved
3. **Calling conventions** - Argument order and semantics match the original
4. **Numerical results** - Output values match the original within floating-point precision

---

### Compatibility Notes

**Fortran Compilers Tested:**
- gfortran 13.x (Linux, macOS, Windows via MSYS2)
- Intel Fortran (ifort/ifx) - untested but should work

**Platforms:**
- Linux (Ubuntu, tested in CI)
- macOS (tested in CI)
- Windows (MSYS2/MinGW-w64, tested in CI)

---

### Known Issues

1. **Shared library builds** may produce duplicate symbol warnings due to module structure
2. **C interoperability** not yet implemented (planned for v1.1)
3. **OpenMP parallelisation** not yet implemented (planned for v1.2)

---

### Future Plans

**v1.1 (Planned):**
- C API via ISO_C_BINDING for key routines
- C header files for language bindings
- Modular shared libraries (libslatec_quadpack.so, etc.)

**v1.2 (Planned):**
- OpenMP parallelisation for suitable routines
- Extended test coverage
- FORD documentation generation

---

### Acknowledgements

SLATEC was developed by the US National Laboratories:
- Sandia National Laboratories
- Los Alamos National Laboratory
- Lawrence Livermore National Laboratory
- National Institute of Standards and Technology (NIST)
- Oak Ridge National Laboratory

The original library is in the public domain as a work of the US Government.

### References

- Abramowitz, M. & Stegun, I.A. (1964). *Handbook of Mathematical Functions*. NBS.
- NIST Digital Library of Mathematical Functions. https://dlmf.nist.gov/
- Amos, D.E. (1986). Algorithm 644: A Portable Package for Bessel Functions. ACM TOMS 12(3).
- Piessens, R. et al. (1983). QUADPACK: A Subroutine Package for Automatic Integration.
- Petzold, L. (1982). A Description of DASSL. SAND82-8637.

---

### Licence

Public Domain. The original SLATEC library is a work of the United States Government and is not subject to copyright protection in the United States.

