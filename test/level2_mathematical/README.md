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

### BLAS Mathematical Properties

| Property | Mathematical Statement |
|----------|------------------------|
| Distributivity | (a+b)·X = a·X + b·X |
| Distributivity | a·(X+Y) = a·X + a·Y |
| Orthogonality | c² + s² = 1 (Givens rotation) |
| Pythagorean | r = √(a² + b²) (Givens) |
| Norm Preservation | \|\|Rx\|\| = \|\|x\|\| |
| Triangle Inequality | \|\|x+y\|\| ≤ \|\|x\|\| + \|\|y\|\| |

## Running

```bash
cd /c/dev/slatec-modern

# MINPACK
gfortran -o test_l2_minpack test/level2_mathematical/test_l2_minpack_mgh.f90
./test_l2_minpack

# BLAS
gfortran -o test_l2_blas test/level2_mathematical/test_l2_linear_blas.f90
./test_l2_blas
```

---

*"In God we trust. All others must bring data." — W. Edwards Deming*
