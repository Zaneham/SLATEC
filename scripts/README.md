# GOTO Elimination Scripts

This folder contains Python scripts used to automate GOTO elimination in SLATEC Fortran routines.

## Why Python?

These scripts are **not** replacing Fortran code with Python. They're used to:

1. **Write Fortran files** - Shell escaping of Fortran syntax (quotes, parentheses, special characters) is extremely error-prone. Python's triple-quoted strings handle this cleanly.

2. **Batch transformations** - Some patterns (like `_SP` to `_DP` conversions for single/double precision pairs) are easier to script than edit manually.

3. **Find GOTO candidates** - `find_goto_files.py` scans the codebase to identify files with the most GOTOs remaining.

## Script Naming Convention

- `fix_<routine>_pair.py` - Fixes both single and double precision versions (e.g., `besy.inc` and `dbesy.inc`)
- `fix_<routine>_trio.py` - Fixes single, double, and complex versions
- `fix_<routine>.py` - Fixes a single routine

## How They Work

Each script typically:
1. Contains the modernised Fortran code as a Python string
2. Writes it to the appropriate `.inc` file in `src/modern/`

Example usage:
```bash
python scripts/fix_besy_pair.py
```

## Output

All output is **Fortran code** (`.inc` files) that gets included into the module structure via Fortran's `INCLUDE` statement.

The `.inc` extension prevents fpm from trying to compile these files as standalone units.
