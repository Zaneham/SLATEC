# SLATEC-Modern

Kia Ora! A complete modernisation of the SLATEC Common Mathematical Library from FORTRAN 77 to modern Fortran (2018+).

## Original Library

- **Version**: 4.1 (July 1993)
- **Source**: [Netlib SLATEC](https://netlib.org/slatec/)
- **Licence**: Public Domain (US Government work)
- **Statistics**: 735 files, 168,216 lines of FORTRAN 77

## What is SLATEC?

SLATEC (Sandia, Los Alamos, Air Force Weapons Laboratory Technical Exchange Committee) is a comprehensive mathematical library developed by US national laboratories:

- Sandia National Laboratories
- Los Alamos National Laboratory
- Lawrence Livermore National Laboratory
- Air Force Weapons Laboratory
- National Institute of Standards and Technology
- Oak Ridge National Laboratory

## Included Sub-Libraries

SLATEC incorporates several well-known numerical libraries:

| Library | Purpose |
|---------|---------|
| BLAS | Basic Linear Algebra Subprograms |
| LINPACK | Linear equation solving |
| EISPACK | Eigenvalue computation |
| FFTPACK | Fast Fourier Transforms |
| QUADPACK | Numerical integration |
| PCHIP | Piecewise Cubic Hermite Interpolation |
| SLAP | Sparse Linear Algebra |
| FNLIB | Special functions |
| BSPLINE | B-spline interpolation |

## Mathematical Categories

| Cat | Description | Routines |
|-----|-------------|----------|
| A | Arithmetic, error analysis | |
| C | Elementary & special functions | Bessel, Gamma, Airy, Error functions |
| D | Linear algebra | Vectors, matrices, eigenvalues, SVD |
| E | Interpolation | Splines, PCHIP, polynomial |
| F | Nonlinear equations | Root finding, systems |
| G | Optimisation | Linear/quadratic programming |
| H | Differentiation & integration | QUADPACK routines |
| I | Differential equations | ODEs, BVPs, PDEs |
| J | Integral transforms | FFT, Laplace |
| K | Approximation | Least squares, fitting |
| L | Statistics & probability | RNG, distributions |

## Modernisation Goals

- [ ] Convert to free-form Fortran 2018
- [ ] Create proper modules with explicit interfaces
- [ ] Replace GOTOs with structured control flow
- [ ] Replace arithmetic IFs
- [ ] Replace COMMON blocks with module variables
- [ ] Replace DATA statements with parameter constants
- [ ] Add `intent` attributes to all arguments
- [ ] Replace `EXTERNAL` with procedure interfaces
- [ ] Modern error handling (replacing XERMSG)
- [ ] Add optional OpenMP parallelisation where applicable
- [ ] Comprehensive test suite
- [ ] Documentation with FORD

## Project Structure

```
slatec-modern/
├── src/
│   ├── original/     # Original FORTRAN 77 source
│   └── modern/       # Modernised Fortran 2018+ source
│       ├── core/     # Machine constants, error handling
│       ├── blas/     # Basic Linear Algebra
│       ├── special/  # Special functions (Category C)
│       ├── linalg/   # Linear algebra (Category D)
│       ├── interp/   # Interpolation (Category E)
│       ├── nonlin/   # Nonlinear equations (Category F)
│       ├── optim/    # Optimisation (Category G)
│       ├── integ/    # Integration/differentiation (Category H)
│       ├── diffeq/   # Differential equations (Category I)
│       ├── fft/      # FFT routines (Category J)
│       ├── approx/   # Approximation (Category K)
│       └── stats/    # Statistics (Category L)
├── tests/
├── examples/
└── docs/
```

## Building

TODO: Add fpm.toml and CMake support

## Authors & Credits

**Original SLATEC**: US National Laboratories (Public Domain)

**Prior Modernisation Work**: [Mehdi Chinoune](https://github.com/MehdiChinoune/SLATEC) - Started SLATEC modernisation in 2021, converting to free-form, adding intents/KIND, marking pure/elemental procedures, and organising into modules. This project builds upon and continues his excellent foundational work.

**Current Modernisation**: Zane Hambly

## References

- [SLATEC Guide](https://www.netlib.org/slatec/guide)
- [SLATEC Table of Contents](https://www.netlib.org/slatec/toc)
- [John Burkardt's SLATEC page](https://people.math.sc.edu/Burkardt/f_src/slatec/slatec.html)
