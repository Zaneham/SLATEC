# SLATEC-Modern

Kia ora! Welcome to a complete modernisation of the SLATEC Common Mathematical Library, taking this 1990s FORTRAN 77 codebase into the modern Fortran era (2018+).

## What's the Story?

Saw this massive old-school maths library sitting there with thousands of GOTOs and thought "screw it, why not give it a go?" Turns out SLATEC gets used in quite a bit of scientific code out there, built by US national laboratories (Sandia, Los Alamos, Lawrence Livermore, NIST, Oak Ridge). The maths is mint, but the code style is properly vintage and its kind of giving my old 1998 Kia Sportage.

This library is still used and is pulled into many scientific applications, it's embedded in scipy and referenced pretty much everywhere. But no one wants to touch it. 

Now I know what you must be thinking: Why're the All Blacks losing their steam? Nah jokes! "Why not GSL?"

Well the issue with GSL is that it sits under the GPL license, now I am not a lawyer but if you link against it then any code you write must also be released under the same license. So its absolutely rubbish for people who're using it for commercial or scientific applications or even just regular people messing about. GSL is also just C with a FORTRAN wig on. SLATEC is FORTRAN native so theres far less friction!

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
- [ ] **Replace GOTOs with structured control flow** *(in progress, started with ~2000, chipping away)*
- [ ] Replace arithmetic IFs
- [ ] Replace COMMON blocks with module variables
- [ ] Replace DATA statements with parameter constants
- [x] Add `intent` attributes to all arguments
- [ ] Replace `EXTERNAL` with procedure interfaces
- [x] Modern error handling (ERROR STOP replacing XERMSG calls)
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

TODO: Add fpm.toml and CMake support

## Credits

**Original SLATEC**: US National Laboratories (Public Domain)

**Prior Modernisation**: [Mehdi Chinoune](https://github.com/MehdiChinoune/SLATEC) kicked off the modernisation back in 2021, converting to free-form, adding intents and KIND parameters, marking pure/elemental procedures, and organising into modules. Massive props for the groundwork.

**PCHIP Reference**: [Jacob Williams](https://github.com/jacobwilliams/PCHIP) created a clean modern Fortran PCHIP implementation in 2019. His control flow patterns (using `do`/`exit` with logical flags instead of GOTOs) informed our approach to modernising the PCHIP routines whilst preserving algorithmic correctness.

**Current Work**: Zane Hambly, continuing the GOTO elimination crusade and general cleanup.

## The use of generative AI

Generative LLM models such as Anthropic's Claude (Opus 4.5 and Sonnet) and my own offline models (Qwen coder with augmented epistemic detection) have been used for running tests and summarising their results. Code and documentation for better or worse is authored by me.

Apologies for the NZ English throughout :-)

## References

- [SLATEC Guide](https://www.netlib.org/slatec/guide)
- [SLATEC Table of Contents](https://www.netlib.org/slatec/toc)
- [John Burkardt's SLATEC page](https://people.math.sc.edu/Burkardt/f_src/slatec/slatec.html)
- [Jacob Williams Modern PCHIP](https://github.com/jacobwilliams/PCHIP)
- Fritsch, F.N. and Carlson, R.E., "Monotone Piecewise Cubic Interpolation", SIAM J. Numer. Anal. 17(2), 1980, pp. 238-246
- Fritsch, F.N., "PCHIP Package", UCRL-87285, Lawrence Livermore National Laboratory, 1982

## Contributing

Got a GOTO that's doing your head in? Found a gnarly bit of code that needs sorting? PRs welcome.
