# SLATEC-Modern

> **Verified on Period Hardware**: This modernisation is tested against authentic 1966 IBM FORTRAN G compilers running on emulated System/360 mainframes. The original SLATEC targeted IBM hexadecimal floating-point (base-16), which has fundamentally different precision characteristics than modern IEEE 754. Our test suite validates numerical correctness across both vintage and contemporary architectures. See testing documentation for details.

Kia ora (hello)! Welcome to a complete modernisation of the SLATEC Common Mathematical Library, dragging this 1993 FORTRAN 77 codebase—kicking and screaming, in some cases—into the modern Fortran era (2018+).

---

## What's All This Then?

SLATEC is one of those libraries that quietly underpins half of scientific computing whilst receiving approximately zero recognition for it. Built by the combined might of Sandia, Los Alamos, Lawrence Livermore, NIST, and Oak Ridge National Laboratories, it contains some of the most battle-tested numerical algorithms ever committed to punch cards.

The mathematics? Absolutely sublime. The code style? Let's just say it has *character*. Over 8,000 GOTOs worth of character, to be precise, written in an era when structured programming was considered a bit avant-garde and "readable code" meant "has comments".

This library is embedded in SciPy, referenced in countless papers, and used in applications ranging from weather prediction to spacecraft trajectories. Yet nobody wants to touch it, presumably because doing so requires a working knowledge of FORTRAN 77, a high tolerance for spaghetti code, and the sort of patience normally reserved for bomb disposal technicians.

So here I am. Someone had to do it.

### But Why Not Just Use GSL?

A fair question. The GNU Scientific Library is comprehensive and well-documented. It also sits under the GPL licence, which means if you link against it, your code must also be released under GPL. This is, to use the technical term, *a bit of a faff* for commercial applications, proprietary research code, or anyone who simply wants to compute a Bessel function without entering into a licensing philosophy debate.

SLATEC, being a product of the US Government, is **public domain**. Use it for whatever you like. Print it out and wallpaper your bathroom with it. I won't judge.

Also, GSL is fundamentally C with a FORTRAN wig on. SLATEC is native Fortran, which means less friction if you're working in a Fortran codebase—and despite rumours to the contrary, rather a lot of scientific computing still is.

---

## The Original Library

| | |
|---|---|
| **Version** | 4.1 (July 1993) |
| **Source** | [Netlib SLATEC](https://netlib.org/slatec/) |
| **Licence** | Public Domain (US Government work) |
| **Size** | 735 files, 168,216 lines of FORTRAN 77 |
| **GOTOs** | 8,296 (I counted, eventually) |
| **Vibes** | Retro |

---

## What's Inside?

SLATEC is essentially a greatest-hits compilation of 1970s–80s numerical libraries:

| Library | What It Does |
|---------|--------------|
| **BLAS** | Basic Linear Algebra Subprograms—the foundation upon which all else rests |
| **LINPACK** | Linear equation solving, from the era before LAPACK existed |
| **EISPACK** | Eigenvalue computation, for when you absolutely must diagonalise a matrix |
| **FFTPACK** | Fast Fourier Transforms, making signals comprehensible since 1985 |
| **QUADPACK** | Numerical integration done properly |
| **PCHIP** | Piecewise Cubic Hermite Interpolation—smooth curves without the oscillations |
| **SLAP** | Sparse Linear Algebra Package, for matrices that are mostly zeros |
| **FNLIB** | Special functions: Bessel, Gamma, Airy, and friends |
| **BSPLINE** | B-spline interpolation and approximation |

### Mathematical Categories

SLATEC organises its routines using a hierarchical classification scheme:

| Category | Description | Examples |
|----------|-------------|----------|
| **A** | Arithmetic & error analysis | Machine constants, precision utilities |
| **C** | Elementary & special functions | Bessel, Gamma, Airy, error functions |
| **D** | Linear algebra | Vectors, matrices, eigenvalues, SVD |
| **E** | Interpolation | Splines, PCHIP, polynomial methods |
| **F** | Nonlinear equations | Root finding, systems of equations |
| **G** | Optimisation | Linear programming, quadratic programming |
| **H** | Differentiation & integration | QUADPACK numerical integration |
| **I** | Differential equations | ODEs, boundary value problems, PDEs |
| **J** | Integral transforms | FFT, convolution, Laplace |
| **K** | Approximation | Least squares, curve fitting |
| **L** | Statistics & probability | Random number generation, distributions |

---

## Modernisation Progress

The grand cleanup:

- [x] Convert to free-form Fortran 2018 (Cheers Mehdi!)
- [x] Create proper modules with explicit interfaces
- [x] **Eliminate GOTOs with structured control flow** (8,296 of 8,296 removed — 100%. Every last one.)
- [x] Replace arithmetic IF statements (remember those? No? Good.)
- [x] Replace DATA statements with parameter constants
- [x] Add `intent` attributes to all procedure arguments
- [x] Modern error handling (ERROR STOP replacing XERMSG)
- [x] CMake build system
- [x] Replace COMMON blocks with module variables (25 blocks across 112+ files — all converted)
- [ ] Replace `EXTERNAL` declarations with procedure interfaces
- [ ] Add optional OpenMP parallelisation
- [ ] Comprehensive test suite
- [ ] FORD documentation

---

## Project Structure

```
slatec-modern/
├── CMakeLists.txt           # CMake build configuration
├── fpm.toml                 # Fortran Package Manager config
├── CHANGELOG.md             # Detailed modernisation history
├── DEVIATIONS.md            # Floating-point deviation analysis
├── IBM360_TEST_RESULTS.md   # Period hardware verification
├── DISCLAIMER               # Public domain notice
├── src/
│   ├── original/            # Original FORTRAN 77 (for reference)
│   └── modern/              # Modernised Fortran 2018+ source
│       ├── approximation/       # Curve fitting, least squares
│       ├── data_handling/       # Sorting, permutation
│       ├── diff_integ/          # Differentiation & integration
│       ├── diff_integ_eq/       # Differential equations
│       ├── integ_trans/         # Integral transforms
│       ├── interpolation/       # Splines, PCHIP
│       ├── linear/              # Linear algebra (BLAS, LINPACK, SLAP)
│       ├── nonlin_eq/           # Nonlinear equation solvers
│       ├── optimization/        # Linear/quadratic programming
│       ├── service/             # Utilities, machine constants
│       └── special_functions/   # Bessel, Gamma, Airy, etc.
├── test/                    # Test suite
├── outputs/                 # Deviation test outputs
├── scripts/                 # Batch processing scripts
└── docs/                    # SLATEC guide and documentation
```

---

## Building

### Fortran Package Manager (fpm)

```bash
fpm build
fpm test
```

### CMake

```bash
mkdir build && cd build
cmake ..
make
ctest
```

Both approaches produce a static library. The fpm route is generally less fuss if you just want to get on with things.

---

## Licence

**TL;DR:** Do what you like with it.

The original SLATEC library is **public domain** software, developed by the United States Government. Works created by US Government employees within the scope of their employment are not subject to domestic copyright protection under 17 U.S.C. § 105.

Modernisation work (GOTO elimination, module restructuring, control flow rewriting) is released under the **Apache License 2.0**. See [LICENSE](LICENSE) for the full text.

Some components incorporate code from Jacob Williams' excellent modern Fortran libraries, released under the **BSD-3-Clause** licence. See Acknowledgements for details.

In practical terms: it's all permissively licensed. Use it commercially, modify it, distribute it, put it on a satellite. Just don't blame me if your Bessel functions come out wonky.

---

## Acknowledgements

This project stands on the shoulders of rather a lot of giants.

### The Original SLATEC Library

The SLATEC Common Mathematical Library was developed collaboratively by:

- **Sandia National Laboratories** (Albuquerque, NM)
- **Los Alamos National Laboratory** (Los Alamos, NM)
- **Air Force Weapons Laboratory** (Kirtland AFB, NM)
- **Lawrence Livermore National Laboratory** (Livermore, CA)
- **National Institute of Standards and Technology** (formerly NBS)
- **Oak Ridge National Laboratory** (Oak Ridge, TN)

"SLATEC" is an acronym for the **S**andia, **L**os **A**lamos, Air Force Weapons Laboratory **T**echnical **E**xchange **C**ommittee, established in 1974 to coordinate mathematical software development across the US national laboratories.

The original authors—too numerous to list individually—created algorithms that have stood the test of time for over four decades. Their work continues to underpin scientific computing worldwide.

### Prior Modernisation Work

**Mehdi Chinoune** ([MehdiChinoune/SLATEC](https://github.com/MehdiChinoune/SLATEC)) initiated the modernisation effort in 2021, providing the foundation for this project:

- Conversion from fixed-form to free-form Fortran
- Addition of `INTENT` attributes to procedure arguments
- Implementation of `KIND` parameters for portable precision
- Identification of `PURE` and `ELEMENTAL` procedures
- Organisation into cohesive module structures

### Jacob Williams' Modern Fortran Libraries

Several components incorporate modernisation patterns and code from Jacob Williams' excellent suite of modern Fortran libraries, released under the BSD-3-Clause licence:

| Library | Contribution |
|---------|-------------|
| [**quadpack**](https://github.com/jacobwilliams/quadpack) | Complete GOTO elimination in numerical integration. BSD-3-Clause, Copyright 2021-2022 Jacob Williams |
| [**PCHIP**](https://github.com/jacobwilliams/PCHIP) | Modern piecewise cubic Hermite interpolation patterns |
| [**ddeabm**](https://github.com/jacobwilliams/ddeabm) | Adams-Bashforth-Moulton ODE solver modernisation |
| [**bspline-fortran**](https://github.com/jacobwilliams/bspline-fortran) | B-spline interpolation |
| [**polyroots-fortran**](https://github.com/jacobwilliams/polyroots-fortran) | Polynomial root finding |
| [**carlson-elliptic-integrals**](https://github.com/jacobwilliams/carlson-elliptic-integrals) | Elliptic integral computation |

Jacob's work on eliminating GOTOs through structured `DO`/`EXIT` constructs with logical flags provided the template for much of our control flow modernisation.

### Current Development

**Zane Hambly** — Complete GOTO elimination (8,296 of 8,296!), control flow modernisation, IBM 360 verification testing, and general acts of code archaeology.

---

## Reproducing the IBM 360 Verification

The period hardware verification is fully reproducible. Here's how:

### Requirements

1. **Hercules** (v3.07+) — Open-source IBM mainframe emulator ([hercules-390.eu](http://www.hercules-390.eu/))
2. **TK4-** — Pre-built MVS 3.8j distribution with compilers ([tk4-.org](http://wotho.ethz.ch/tk4-/))
3. **FORTRAN 360** — Python toolchain for the vintage compilers ([github.com/Zaneham/Fortran360](https://github.com/Zaneham/Fortran360))

### Setup

```bash
# Clone the FORTRAN 360 toolchain
git clone https://github.com/Zaneham/Fortran360.git
cd Fortran360

# Start Hercules with TK4- (see TK4- documentation)
# The system exposes a card reader on port 3505

# Run SLATEC tests against IBM FORTRAN G (1966)
python -m src.fortran360.cli run tests/slatec/slatec_test.f --era=1966
```

### What Gets Tested

| Component | Verification |
|-----------|--------------|
| **I1MACH** | Integer machine constants (FP base=16, mantissa digits=6) |
| **R1MACH** | Real machine constants via EQUIVALENCE bit patterns |
| **GAMLN** | Log-gamma function across integer and non-integer arguments |

### Why This Matters

IBM 360 uses **hexadecimal floating-point** (base-16), which differs fundamentally from IEEE 754:

- Mantissa precision "wobbles" between 21–24 bits depending on leading digit
- No denormalised numbers, NaN, or infinity
- Truncation rounding (not round-to-nearest)

Numerical edge cases that pass on modern IEEE hardware may behave differently on the original architecture. Testing both ensures the modernised code is faithful to the original mathematical intent, not just the original hardware quirks.

See [IBM360_TEST_RESULTS.md](IBM360_TEST_RESULTS.md) for full test output and analysis.

---

## A Note on AI Assistance

Generative LLM models (Anthropic's Claude Opus and Sonnet, plus local models including Qwen Coder) have been used primarily as:

- **Test result summarisation** — Collating output from verification runs
- **Rubber duck debugging** — Talking through control flow transformations. I speak to the rubber duck, it speaks back to me.
- **Documentation drafting** — Under human review and editing

Code, architectural decisions, and the GOTO elimination work itself remain human-authored. The AI doesn't write the Fortran; it helps make sense of what the Fortran is doing.

Apologies for the NZ English throughout. Colour has a 'u' in it and that's simply how it is.

---

## References

### Documentation

- [SLATEC Guide](https://www.netlib.org/slatec/guide) — The official user guide
- [SLATEC Table of Contents](https://www.netlib.org/slatec/toc) — Complete routine listing
- [John Burkardt's SLATEC page](https://people.math.sc.edu/Burkardt/f_src/slatec/slatec.html) — Comprehensive reference

### Key Publications

- Piessens, R., de Doncker-Kapenga, E., Überhuber, C. W., and Kahaner, D. K. *QUADPACK: A Subroutine Package for Automatic Integration*. Springer-Verlag, 1983.
- Fritsch, F. N. and Carlson, R. E. "Monotone Piecewise Cubic Interpolation". *SIAM Journal on Numerical Analysis* 17(2), 1980, pp. 238–246.
- Amos, D. E. "A Subroutine Package for Bessel Functions of a Complex Argument and Nonnegative Order". Sandia National Laboratories Report SAND85-1018, 1985.

---

## Related Projects

[SLATEC.jl](https://github.com/Zaneham/SLATEC.jl) ports validated routines to Julia for those who prefer their numerical archaeology without the semicolons. Only routines passing all four test levels get ported. The rest remain here, documented but not shipped.

---

## Contributing

Found something that needs sorting? Spotted an issue or want to help with the remaining modernisation tasks? PRs welcome.

If you're feeling particularly brave, the differential equation solvers in `diff_integ_eq/` still have some rather creative control flow that could use attention. I recommend a stiff drink beforehand. Or several. I'm not here to judge.

The COMMON blocks are done (all 25 of them, converted to modules with `_com` suffix variables), so you'll need to find a different excuse to question your life choices. The `EXTERNAL` declarations are still waiting, if you're feeling masochistic.
