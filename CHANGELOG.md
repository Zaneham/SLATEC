# CHANGELOG - SLATEC Modernisation Project

## Overview

This project modernises the SLATEC Fortran 77 library to Fortran 2018 standards, with a primary focus on eliminating GOTO statements and replacing them with structured control flow constructs.

**Total GOTOs Eliminated:** ~1,400+ across 53 batches

---

## Recent Changes (Unpushed)

### Batch 53 - K Bessel Function Drivers (34 GOTOs)
- `dbesk.inc`: 17 GOTOs → 0
- `besk.inc`: 17 GOTOs → 0
- Pattern: Logical flags + named DO loops for underflow handling

### Batch 52 - I Bessel Functions (52 GOTOs)
- `dbesi.inc`: 26 GOTOs → 0 (state machine with 12 states)
- `besi.inc`: 26 GOTOs → 0
- Pattern: SELECT CASE state machine for complex multi-path control flow

### Bug Fix: Pre-existing Build Errors
Fixed compilation issues that crept in during earlier modernization:
- `dpsifn.inc`: Mixed 'd' exponent with explicit kind (0.5D-18_DP → 0.5E-18_DP)
- `dasyjy.inc` / `asyjy.inc`: Missing LOGICAL declarations
- `dxnrmp.inc`: Wrong SP/DP types for D-prefixed routine
- `dbesy.inc`: Typo DDBESYNU → DBSYNU

### Batch 51 - QUADPACK Integration (12 GOTOs)
- `dqagpe.inc`, `qagpe.inc`, `qagie.inc`: 4 GOTOs each → 0
- Pattern: `final_result: BLOCK` with EXIT

### Batch 50 - Nonlinear Equation Solvers (18 GOTOs)
- `dsoseq.inc`: 9 GOTOs → 0
- `soseqs.inc`: 9 GOTOs → 0
- Pattern: Nested DO WHILE loops + logical flags

### Batch 49 - LQ Factorisation (12 GOTOs)
- `du12ls.inc`, `u12ls.inc`: 6 GOTOs each → 0

### Batch 48 - Singleton Quicksort (39 GOTOs)
- `spsort.inc`, `dpsort.inc`, `ipsort.inc`, `hpsort.inc`

### Batch 47 - SPLP Sparse Matrix (4 GOTOs)
- `dpnnzr.inc`, `pnnzrs.inc`: Sparse nonzero retrieval routines

### Bug Fix: Restore Lost Code (210 GOTOs recovered)
- Commit `9df00fa`: Fixed regression from .f90 → .inc refactor
- Restored GOTO-free code that was accidentally overwritten

### Batch 46 - Oscillatory Integration (12 GOTOs)
- `dqawoe.inc`, `qawoe.inc`: Following Jacob Williams' quadpack pattern

### Batch 45 - QUADPACK Adaptive (26 GOTOs)
- `dqagse.inc`, `qagse.inc` and related routines

### Batch 44 - Wigner/Bessel Asymptotic (16 GOTOs)
- `d9b0mp.inc`, `d9b1mp.inc`, `drc3jm.inc`, `drc3jj.inc`

### Infrastructure: .f90 → .inc Refactor
- Renamed include files from `.f90` to `.inc` extension
- Prevents fpm from attempting standalone compilation

### Bug Fix: DBESK Numerical Drift
- **Severity:** High (numerical accuracy)
- K₄(1) computed as 44.2324 vs correct 44.2341 (~4×10⁻⁵ error)
- Root cause: series_done flag missing, causing value overwrite
- Fixed in `dbsknu.inc`, `besknu.inc`

### Batch 43 - Final Special Functions (41 GOTOs)
- `cuni2.inc`, `psifn.inc`, `dpsifn.inc`, `xnrmp.inc`, `dxnrmp.inc`
- Completed GOTO elimination in special_functions directory

### Batch 42 - I/Y Bessel Series (10 GOTOs)
- `cseri.inc`, `besy.inc`, `dbesy.inc`

### Batch 41 - Extended-Range & Bickley (30 GOTOs)
- `xadd.inc`, `dxadd.inc`, `bskin.inc`, `dbskin.inc`

### Batch 40 - K Bessel Kernel (32 GOTOs)
- `zbknu.inc`, `cbknu.inc`: Complex K Bessel main kernel

### Batch 39 - Single-Precision Bessel (50 GOTOs)
- `cunk1.inc`, `cunk2.inc`, `cbinu.inc`

### Batch 38 - zbinu Dispatcher (16 GOTOs)
- I Bessel computation method selection

### Batch 37 - Uniform Asymptotic (34 GOTOs)
- `zunk1.inc`, `zunk2.inc`: K Bessel analytic continuation

### Test Infrastructure
- `test_bessel_regression.f90`: Regression tests against stored values
- `test_bessel_reference.f90`: Tests against A&S tables
- `test_bessel_portability.f90`: Cross-platform ULP error measurement
- `test_utils.f90`: Drift analysis utilities

### Batch 36 - J Bessel Functions (46 GOTOs)
- `dbesj.inc`: 23 GOTOs → 0
- `besj.inc`: 23 GOTOs → 0

### Batch 35 - DEPAC RK Routines (20 GOTOs)
- `derkfs.inc`, `drkfs.inc`: Fehlberg RK4-5 steppers

### Batch 34 - DEPAC Steps/LSOD (36 GOTOs)
- `dsteps.inc`, `steps.inc`, `dlsod.inc`, `lsod.inc`

### Batch 33 - DEPAC Step Routines (47 GOTOs)
- `dstod.inc`, `stod.inc`: Adams-Bashforth/BDF steppers

### Batch 32 - SDRIVE ODE Package (162 GOTOs)
- `ddriv3.inc`, `sdriv3.inc`, `cdriv3.inc`
- `ddstp.inc`, `sdstp.inc`, `cdstp.inc`

### Batch 31 - DASSL DAE Solver (118 GOTOs)
- `ddassl.inc`, `sdassl.inc`: Main driver
- `ddastp.inc`, `sdastp.inc`: Step routines

---

## Bugs Found & Fixed

See `bugs_found.md` for detailed documentation. Summary:

| Bug | Severity | Routines | Status |
|-----|----------|----------|--------|
| DBESK numerical drift | High | dbsknu, besknu | Fixed |
| Missing BLOCK statement | Medium | dcov, scov | Fixed |
| Unclosed DO loop | Medium | sort routines | Fixed |
| Undeclared variable | Low | pchid, pchia | Fixed |
| Build system issues | Medium | various | Fixed |

---

## Patterns Used

### State Machine (Complex Multi-Path)
```fortran
INTEGER, PARAMETER :: ST_INIT = 1, ST_COMPUTE = 2, ST_DONE = 99
state_loop: DO WHILE( istate /= ST_DONE )
  SELECT CASE( istate )
    CASE( ST_INIT )
      ! ...
      istate = ST_COMPUTE
    CASE( ST_COMPUTE )
      ! ...
      istate = ST_DONE
  END SELECT
END DO state_loop
```

### Logical Flags (Forward Jumps)
```fortran
LOGICAL :: skip_section, done
IF( condition ) skip_section = .TRUE.
IF( .NOT. skip_section ) THEN
  ! ... code that was skipped via GOTO ...
END IF
```

### Named DO Loops (Backward Jumps)
```fortran
retry_loop: DO
  ! ... computation ...
  IF( converged ) EXIT retry_loop
  IF( max_iter ) EXIT retry_loop
END DO retry_loop
```

### BLOCK/EXIT (Forward to Common Exit)
```fortran
final_result: BLOCK
  IF( error ) EXIT final_result
  ! ... normal processing ...
END BLOCK final_result
! Common cleanup code here
```

---

## References

- Amos, D.E., Daniel, S.L., Weston, M.K. - ACM TOMS 3, 1977 (Bessel functions)
- Temme, N.M. - J. Comp. Physics 19, 1975 (K Bessel)
- Petzold, L. - SAND82-8637 (DASSL)
- Piessens, R. & de Doncker, E. - QUADPACK (Integration)
- Hanson, R.J. & Hiebert, K.L. - SAND81-0297 (SPLP)
