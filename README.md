# libfint

This is a Fortran port of Qiming Sun's [libcint](https://github.com/sunqm/libcint).
The intention of the port is exploring the Fortran language for performance
oriented, friendly API design for general electronic integrals over Gaussian
functions used in Quantum Chemistry programs. The code was validated against
the C version and a CI workflow is performed to evaluate drift from the
original package.

Additional intentions are meant to provide more mature libraries to be tested
with the [LFortran](https://github.com/lfortran/lfortran) compiler, simpler
integration into the [Fortran Package Manager](https://github.com/fortran-lang/fpm)
and just me being stubborn and curious.

If using libfint please also cite libcint since it is a derivative work.
Libfint will become the main integral engine for
[mqc](https://github.com/JorgeG94/metalquicha) and if you are developing
versatile quantum chemistry codes in Fortran and wish to stick to an all
Fortran environment, you can use libfint.

## Validation

Bit-identical, not "equivalent to within tolerance".  The suite walks every
entry point of both libraries over a shared basis set and compares every
double:

```
574 entry points, 2,751,264 values compared against libcint, 0 over tolerance
```

There is no tolerance at any Rys-root count. One difference is documented
rather than hidden, and it is visible in the check's own output: 228 values in
the r-polynomial families differ in the last bit, because the C folds a
coefficient into the product where the generator materialises the product and
scales afterwards. And
`int3c1e_ip1_r6_origk` differs by rather more against an unpatched libcint,
because the C reads a `g` block it never writes -- that one is a bug in the C,
and the patch is offered upstream.

## Status

620 of libcint's 635 entry points are ported -- the 15 that are not are ones
libcint itself refuses -- along with all 222 optimizers, plus `int2c2e_ipip2`,
which libcint does not have.

Performance is 1.03x to 1.18x libcint's, measured per entry point on an idle
machine. `doc/PERFORMANCE.md` has the table, the method, and a long list of
things that were tried and did not work, which is the more useful half.

## Building

With fpm:

```
fpm build
fpm test
```

With CMake:

```
cmake -S . -B build && cmake --build build -j
```

Requires a Fortran 2008 compiler (gfortran 11+ or ifx). Nothing else -- no
LAPACK, no BLAS, no C. The undefined symbols in `libfint.a` are `erf`, `erfc`,
`exp`, `lgamma`, `log` and `pow`, and that is the whole list.

The Wheeler quadrature needs a tridiagonal eigensolver, and libfint carries a
translation of the one libcint carries -- LAPACK's DSTEMR, which libcint
vendored and which its build always uses. `-DWITH_EXTERNAL_LAPACK=ON` swaps in
a system `dstemr` instead, for callers who already link LAPACK and would
rather have a vendor implementation. It is faster and it is not
bit-identical.

## Verifying the claim

libfint contains no C. Checking it against libcint therefore means fetching
libcint, which the build does for you -- from `sunqm/libcint` master, not a
fork and not a pin:

```
scripts/fetch_libcint.sh
cmake -S . -B build -DLIBFINT_BUILD_TESTS=ON && cmake --build build -j
cd build && ctest
```

Floating rather than pinned is deliberate. A pin would make this a regression
test on libfint; upstream's master makes it a drift detector, and a red run
says which kind of drift:

| failure | meaning |
|---|---|
| `catalogue_check` | upstream changed a kernel and the Fortran generator has not been told |
| `integral_coverage` | upstream added an integral; copy its description into `scripts/gen_f90_all.cl` |
| build | upstream changed its build or its headers |

Pass a ref to reproduce an older run: `scripts/fetch_libcint.sh v6.1.3`.

The catalogue is generated against whichever libcint was fetched rather than
committed, so it always compares the intersection of the two libraries -- and
`integral_coverage` is what tells you about the difference.

## Relationship to libcint

This is a derivative work under the Apache 2.0 license, which libcint grants
explicitly. It tracks a specific upstream version rather than diverging: when
libcint adds an integral, the description is copied into the Fortran
generator and the suite re-run.

One deliberate difference from stock libcint 6.1.3: a missing
`G1E_D_I(g76, ...)` in `src/cint3c1e_a.c`, which is a real bug in the C.
libfint has it fixed, so `int3c1e_a` differs from an unpatched libcint. The
patch has been offered upstream.

## License

Apache 2.0, the same as libcint. See `LICENSE` and `NOTICE`.
