# Level 2: Mathematical Verification

*"Does the algorithm match the original mathematics?"*

## Purpose

Level 2 tests compare computed values against authoritative mathematical references. If Level 2 fails, either:

1. The original SLATEC algorithm was incorrect (rare but documented historically)
2. The modernisation introduced a mathematical error
3. The reference value is wrong (check multiple sources)

**Do not change Level 2 tests without mathematical justification.**

## Authoritative References

### MINPACK (Optimisation)
- Moré, J.J., Garbow, B.S., Hillstrom, K.E. (1981). "Testing Unconstrained Optimization Software." *ACM TOMS* 7(1), 17-41.
- Moré, J.J. et al. (1980). "User Guide for MINPACK-1." ANL-80-74. Argonne National Laboratory.

### Test Functions Used

| Function | Reference | Dimension | Solution |
|----------|-----------|-----------|----------|
| Rosenbrock | Rosenbrock (1960) | 2 | (1, 1) |
| Powell Singular | Powell (1962) | 4 | (0, 0, 0, 0) |
| Helical Valley | Fletcher & Powell (1963) | 3 | (1, 0, 0) |
| Freudenstein-Roth | Freudenstein & Roth (1963) | 2 | (5, 4) |
| Broyden Tridiagonal | Broyden (1965) | N | varies |

### Other References
- Abramowitz & Stegun (1964) — *Handbook of Mathematical Functions*
- NIST DLMF — https://dlmf.nist.gov/
- Closed-form solutions (e.g., ∫₀¹ x² dx = 1/3)

## Coverage

| File | Module | Properties Tested | Status |
|------|--------|-------------------|--------|
| test_l2_minpack_mgh.f90 | approximation | Rosenbrock, Powell, Helical Valley, Freudenstein-Roth, DENORM | ✓ **17/17 PASS** |
| test_l2_linear_blas.f90 | linear | Linearity, distributivity, orthogonality, Pythagorean, norms | ✓ **16/16 PASS** |
| test_l2_interpolation.f90 | interpolation | B-spline, PCHIP, polynomial interpolation | ⚠ **14/17 PASS** |
| test_l2_diff_integ.f90 | diff_integ | Gauss-Kronrod exactness, classical integrals, infinite intervals | ✓ **17/17 PASS** |

### BLAS Mathematical Properties

| Property | Mathematical Statement |
|----------|------------------------|
| Distributivity | (a+b)·X = a·X + b·X |
| Distributivity | a·(X+Y) = a·X + a·Y |
| Orthogonality | c² + s² = 1 (Givens rotation) |
| Pythagorean | r = √(a² + b²) (Givens) |
| Norm Preservation | \|\|Rx\|\| = \|\|x\|\| |
| Triangle Inequality | \|\|x+y\|\| ≤ \|\|x\|\| + \|\|y\|\| |

### Interpolation Mathematical Properties

| Property | Reference | Status |
|----------|-----------|--------|
| B-spline interpolation: S(x_i) = y_i | de Boor (1978) Ch. IX | **FAIL** |
| Clamped spline: S'(a) = f'(a), S'(b) = f'(b) | de Boor (1978) Thm IX.1 | ✓ PASS |
| Natural spline: S''(a) = S''(b) = 0 | de Boor (1978) | **FAIL** (S''(b) ≠ 0) |
| PCHIP interpolation: p(x_i) = y_i | Fritsch & Carlson (1980) | ✓ PASS |
| PCHIP monotonicity preservation | Fritsch & Carlson (1980) Thm 1 | ✓ PASS |
| Polynomial uniqueness theorem | Burden & Faires Thm 3.2 | ✓ PASS |
| Polynomial exactness for degree ≤ n-1 | Burden & Faires Cor 3.3 | ✓ PASS |

#### Interpolation Test Results (2025-01-18)

**B-spline Issues Identified:**
- DBINT4 natural spline does not interpolate exactly (max error: 3.56e-03)
- DBINT4 clamped spline does not interpolate exactly (max error: 0.36)
- Natural spline boundary condition S''(1) = 60.1 instead of 0

**Investigation Notes:**
The B-spline interpolation failure appears to be in DBINT4's handling of the
coefficient system. The boundary conditions ARE being applied correctly to the
linear system, but the spline does not pass through the data points as expected.
This may indicate an issue with the basis function evaluation or system solve.

PCHIP and polynomial interpolation work correctly.

### Diff_Integ Mathematical Properties

| Property | Reference | Status |
|----------|-----------|--------|
| Gauss-Kronrod exactness: QK15 exact for deg ≤ 29 | QUADPACK (1983) | ✓ PASS |
| Gauss-Kronrod exactness: QK21 exact for deg ≤ 41 | QUADPACK (1983) | ✓ PASS |
| Gauss-Kronrod exactness: QK31 exact for deg ≤ 61 | QUADPACK (1983) | ✓ PASS |
| Classical integrals: ∫sin(x)dx, ∫exp(-x)dx | A&S | ✓ PASS |
| Endpoint singularities: ∫1/√x dx, ∫log(x)dx | A&S | ✓ PASS |
| Infinite intervals: ∫₀^∞ exp(-x)dx = 1 | A&S | ✓ PASS |
| Gaussian integral: ∫_{-∞}^{∞} exp(-x²)dx = √π | A&S 7.1.1 | ✓ PASS |

## Running

```bash
cd /c/dev/slatec-modern

# MINPACK
gfortran -o test_l2_minpack test/level2_mathematical/test_l2_minpack_mgh.f90
./test_l2_minpack

# BLAS
gfortran -o test_l2_blas test/level2_mathematical/test_l2_linear_blas.f90
./test_l2_blas

# Interpolation
gfortran -O2 -I build/modules -o test_l2_interpolation test/level2_mathematical/test_l2_interpolation.f90 -L build/lib -lslatec
./test_l2_interpolation

# Diff_Integ
gfortran -O2 -I build/modules -o test_l2_diff_integ test/level2_mathematical/test_l2_diff_integ.f90 -L build/lib -lslatec
./test_l2_diff_integ
```

---

*"In God we trust. All others must bring data." — W. Edwards Deming*
