# SLATEC-Modern

Kia ora! Welcome to a complete modernisation of the SLATEC Common Mathematical Library, taking this 1990s FORTRAN 77 codebase into the modern Fortran era (2018+).

## What's the Story?

Saw this massive old-school maths library sitting there with thousands of GOTOs and thought "screw it, why not give it a go?" Turns out SLATEC gets used in quite a bit of scientific code out there, built by US national laboratories (Sandia, Los Alamos, Lawrence Livermore, NIST, Oak Ridge). The maths is mint, but the code style is properly vintage and its kind of giving my old 1998 Kia Sportage.

This library is still used and is pulled into many scientific applications, it's embedded in scipy and referenced pretty much everywhere. But no one wants to touch it.

Now I know what you must be thinking: Why're the All Blacks losing their steam? Nah jokes! "Why not GSL?"

Well the issue with GSL is that it sits under the GPL licence, now I am not a lawyer but if you link against it then any code you write must also be released under the same licence. So its absolutely rubbish for people who're using it for commercial or scientific applications or even just regular people messing about. GSL is also just C with a FORTRAN wig on. SLATEC is FORTRAN native so theres far less friction!

So here we are, keeping the mathematical goodness while making the code actually readable.

## Original Library

- **Version**: 4.1 (July 1993)
- **Source**: [Netlib SLATEC](https://netlib.org/slatec/)
- **Licence**: Public Domain (US Government work)
- **Size**: 735 files, 168,216 lines of FORTRAN 77

## What's Inside?

SLATEC packs in some absolute classics:

| Library | What it Does |
|---------|--------------|
| BLAS | Basic Linear Algebra, the foundation of everything |
| LINPACK | Linear equation solving |
| EISPACK | Eigenvalue computation |
| FFTPACK | Fast Fourier Transforms |
| QUADPACK | Numerical integration (proper good stuff) |
| PCHIP | Piecewise Cubic Hermite Interpolation |
| SLAP | Sparse Linear Algebra |
| FNLIB | Special functions galore |
| BSPLINE | B-spline interpolation |

## Mathematical Categories

| Cat | Description | Examples |
|-----|-------------|----------|
| A | Arithmetic, error analysis | Machine constants |
| C | Elementary & special functions | Bessel, Gamma, Airy, Error functions |
| D | Linear algebra | Vectors, matrices, eigenvalues, SVD |
| E | Interpolation | Splines, PCHIP, polynomial |
| F | Nonlinear equations | Root finding, systems |
| G | Optimisation | Linear/quadratic programming |
| H | Differentiation & integration | QUADPACK routines |
| I | Differential equations | ODEs, BVPs, PDEs |
| J | Integral transforms | FFT, Laplace |
| K | Approximation | Least squares, curve fitting |
| L | Statistics & probability | RNG, distributions |

## Modernisation Progress

The big cleanup jobs:

- [x] Convert to free-form Fortran 2018
- [x] Create proper modules with explicit interfaces
- [x] **Replace GOTOs with structured control flow** *(~1,400 eliminated!)*
- [ ] Replace arithmetic IFs
- [ ] Replace COMMON blocks with module variables *(note: some .f90 files in subdirectories kept for reference)*
- [ ] Replace DATA statements with parameter constants
- [x] Add `intent` attributes to all arguments
- [ ] Replace `EXTERNAL` with procedure interfaces
- [x] Modern error handling (ERROR STOP replacing XERMSG calls)
- [x] CMake build system (static library)
- [ ] Add optional OpenMP parallelisation where applicable
- [ ] Comprehensive test suite
- [ ] Documentation with FORD

## Project Structure

```
slatec-modern/
├── src/
│   ├── original/          # Original FORTRAN 77 source (for reference)
│   └── modern/            # Modernised Fortran 2018+ source
│       ├── approximation/ # Curve fitting, least squares
│       ├── diff_integ/    # Differentiation & integration
│       ├── diff_integ_eq/ # Differential equations
│       ├── interpolation/ # Splines, PCHIP
│       ├── linear/        # Linear algebra
│       ├── service/       # Utility routines
│       └── special_functions/
├── scripts/               # Python fix scripts for batch GOTO elimination
├── tests/
├── examples/
└── docs/
```

## Building

Uses the Fortran Package Manager (fpm):

```bash
fpm build
fpm test
```

## Licence

The original SLATEC library is **public domain** software, developed by the United States Government. Works created by US Government employees within the scope of their employment are not subject to domestic copyright protection under 17 U.S.C. § 105.


**PCHIP Reference**: [Jacob Williams](https://github.com/jacobwilliams/PCHIP) created a clean modern Fortran PCHIP implementation in 2019. His control flow patterns (using `do`/`exit` with logical flags instead of GOTOs) informed our approach to modernising the PCHIP routines whilst preserving algorithmic correctness.

**Current Work**: Zane Hambly, continuing the GOTO elimination crusade and general cleanup.
=======
Modernised code in this repository is released under the same public domain terms, except where components incorporate code from other open-source projects (see Acknowledgements below for specific licence terms).

## Acknowledgements

This project builds upon the work of many contributors. I gratefully acknowledge:

### Original SLATEC Library

The SLATEC Common Mathematical Library was developed collaboratively by:

- **Sandia National Laboratories** (Albuquerque, NM)
- **Los Alamos National Laboratory** (Los Alamos, NM)
- **Air Force Weapons Laboratory** (Kirtland AFB, NM)
- **Lawrence Livermore National Laboratory** (Livermore, CA)
- **National Institute of Standards and Technology** (formerly NBS)
- **Oak Ridge National Laboratory** (Oak Ridge, TN)

"SLATEC" is an acronym for the Sandia, Los Alamos, Air Force Weapons Laboratory Technical Exchange Committee, established in 1974.

### Prior Modernisation Work

**Mehdi Chinoune** ([MehdiChinoune/SLATEC](https://github.com/MehdiChinoune/SLATEC))

Initiated the modernisation effort in 2021, providing the foundation for this project:
- Conversion from fixed-form to free-form Fortran
- Addition of `INTENT` attributes to procedure arguments
- Implementation of `KIND` parameters for portable precision
- Identification of `PURE` and `ELEMENTAL` procedures
- Organisation into cohesive module structures

### Jacob Williams' Modern Fortran Libraries

Several components incorporate modernisation patterns and code from Jacob Williams' excellent suite of modern Fortran libraries, released under the BSD-3-Clause licence:

**QUADPACK** ([jacobwilliams/quadpack](https://github.com/jacobwilliams/quadpack))
- Complete GOTO elimination in numerical integration routines
- Modern control flow patterns using `DO`/`EXIT` constructs
- BSD-3-Clause Licence, Copyright © 2021-2022 Jacob Williams

**PCHIP** ([jacobwilliams/PCHIP](https://github.com/jacobwilliams/PCHIP))
- Piecewise Cubic Hermite Interpolation Package
- Structured control flow patterns that informed our approach

**ddeabm** ([jacobwilliams/ddeabm](https://github.com/jacobwilliams/ddeabm))
- Adams-Bashforth-Moulton ODE solver modernisation
- Object-oriented Fortran design patterns

**Additional libraries**: [bspline-fortran](https://github.com/jacobwilliams/bspline-fortran), [polyroots-fortran](https://github.com/jacobwilliams/polyroots-fortran), [carlson-elliptic-integrals](https://github.com/jacobwilliams/carlson-elliptic-integrals)

### Current Development

**Zane Hambly** — Ongoing GOTO elimination, control flow modernisation, and test suite development.


## The use of generative AI

Generative LLM models such as Anthropic's Claude (Opus 4.5 and Sonnet) and my own offline models (Qwen coder with augmented epistemic detection) have been used for running tests and summarising their results. Code and documentation for better or worse is authored by me.

Apologies for the NZ English throughout :-)

## References

### Documentation
- [SLATEC Guide](https://www.netlib.org/slatec/guide)
- [SLATEC Table of Contents](https://www.netlib.org/slatec/toc)
- [John Burkardt's SLATEC page](https://people.math.sc.edu/Burkardt/f_src/slatec/slatec.html)

### Key Publications
- Piessens, R., de Doncker-Kapenga, E., Überhuber, C. W., and Kahaner, D. K. *QUADPACK: A Subroutine Package for Automatic Integration*. Springer-Verlag, 1983.
- Fritsch, F. N. and Carlson, R. E. "Monotone Piecewise Cubic Interpolation". *SIAM Journal on Numerical Analysis* 17(2), 1980, pp. 238–246.
- Amos, D. E. "A Subroutine Package for Bessel Functions of a Complex Argument and Nonnegative Order". Sandia National Laboratories Report SAND85-1018, 1985.

## Contributing

Got a GOTO that's doing your head in? Found a gnarly bit of code that needs sorting? PRs welcome.
