# Obara–Saika with the Head-Gordon–Pople split

A third route to the four-centre Coulomb integrals, beside Rys quadrature
and the rotated-axis McMurchie–Davidson kernels. It exists because the
rotated-axis path stops paying above p: measured inside a threaded Fock
build, its d-touching quartets ran 2.5x slower than Rys while its s, p and
L quartets ran 1.75x faster.

Like the rotated-axis kernels, none of this is written by hand.
`scripts/hgp/` carries the recurrences, emits one Fortran module per class
into `src/hgp/`, and — before any Fortran exists —
`scripts/hgp/check_numeric.py` evaluates the recurrence graph numerically
and compares it against an independent McMurchie–Davidson reference. That
check agrees to 2e-14 over 4148 values across every class through (dd|dd).

## 1. The two recurrences

Write a primitive quartet with exponents a, b on the bra and c, d on the
ket, p = a + b, q = c + d, P and Q the Gaussian product centres,
rho = pq/(p+q), and W = (pP + qQ)/(p+q).

**Vertical (Obara–Saika).** Build integrals with every unit of angular
momentum on A and C, denoted [e0|f0], from the auxiliary series
[00|00]^(m) = pref · F_m(T), T = rho|P−Q|²:

    [(e+1_i)0|f0]^(m) = (P−A)_i [e0|f0]^(m) + (W−P)_i [e0|f0]^(m+1)
                      + (e_i/2p) ( [(e−1_i)0|f0]^(m) − (rho/p)[(e−1_i)0|f0]^(m+1) )
                      + (f_i/2(p+q)) [e0|(f−1_i)0]^(m+1)

and the mirror form incrementing f with (Q−C), (W−Q), 1/2q and rho/q.

**Horizontal (Head-Gordon–Pople).** Move momentum onto B and D:

    (a(b+1_i)|cd) = ((a+1_i)b|cd) + (A−B)_i (ab|cd)
    (ab|c(d+1_i)) = (ab|(c+1_i)d) + (C−D)_i (ab|cd)

## 2. Why the split is the point

The horizontal recurrence contains no exponent. It therefore commutes with
the contraction, and HGP runs it *after*: the vertical recurrence and the
contraction happen per primitive quartet, the transfers once per contracted
quartet. What crosses the contraction is only the [e0|f0] set, which is far
smaller than the final component count once momentum is spread over four
centres.

That number is the one that decided this. The rotated-axis attempt at d
carried a ket accumulator per distinct polynomial and reached 4351 of them
for (dd|dd), expanded into 96k terms; the table lookups that produced were
what lost to Rys.

| class | rotated-axis carried | HGP carries | components |
|---|---|---|---|
| pppp | 143 | 81 | 81 |
| LdLd | 1384 | 400 | 576 |
| ppdd | 1050 | 279 | 324 |
| dddd | 4351 | 961 | 1296 |

## 3. Nothing is expanded

The generator emits the recurrence, not its solution. Each intermediate is
a node computed once from earlier nodes, so common subexpressions are
shared by construction and the emitted code is straight-line arithmetic
whose length is the node count: 2330 vertical and 6254 transfer nodes for
(dd|dd), against the 96k expanded terms the other approach produced for
the same class. This is the lesson from that attempt applied directly.

## 4. Shell kinds

A slot is a kind — s, p, L or d — exactly as in the rotated-axis
generator, an L shell being one four-component slot carrying s and p on
shared exponents with a coefficient column for each. The vertical
recurrence is shared across an L shell's two angular momenta, and only the
transfers and the contraction weights differ per coefficient type, so an
(LL|LL) class runs one vertical recurrence and sixteen small transfer
blocks.

Classes are canonical with the **higher** kind first in each pair, since
that is the direction the transfers run. This is the reverse of the
rotated-axis convention, where the kinds ascend.

## 5. Contraction, and why it is staged

The contraction coefficients multiply in at the contraction rather than
riding in the primitive weight, which is what lets one vertical recurrence
serve every contraction column and every coefficient type. But *how* they
multiply in decides whether the path is usable at all on a generally
contracted basis.

Done flat — touching every (bra column, ket column) pair at every primitive
quartet — the accumulation costs the product of the two contraction counts
per quartet. On cc-pVDZ carbon, whose s shell is 9 primitives into 3
contractions and whose p shell is 4 into 2, an (ss|ss) means 81
accumulations of the whole target vector at each of 6561 primitive
quartets. It is invisible on a segmented basis, where both counts are one.

Staged, it is a sum instead of a product: sum the bra columns into a buffer
inside the primitive loop, then fold the ket columns in once per ket
primitive.

    inside the bra loop:   cb(ef, bcol)       += c_i c_j (bcol) · [e0|f0]
    once per ket primitive: c(ef, bcol, kcol) += c_k c_l (kcol) · cb(ef, bcol)

This is what libcint's `cint_prim_to_ctr` staging does, and it is the
single change that took this path from losing to winning on a generally
contracted basis. The geometric part — 2π^{5/2} exp(−μR_AB²) exp(−νR_CD²)
/ (p q √(p+q)) — rides in the pair tables and is shared throughout.

## 6. Numerics

The Boys function is the same guarded routine the rotated-axis kernels use,
and it now has three branches rather than two.

* Below the turnover point, the ascending series. libcint's closed form
  cancels catastrophically for F_1 at a tiny argument, which the Rys path
  never asks for.
* Above t = 50, a closed form with no library calls at all: exp(−t) is
  under 2e-22 there and erf(√t) is 1 to better than that, so
  F_0 = √π/(2√t) and the rest is (2i−1)/(2t) upward. This matters because
  a screened Fock build spends much of its time on distant pairs, and
  profiled inside one, libm's erf and exp were 12 s of 101 s on this path
  against a Rys path that calls neither.
* In between, the closed form as before.

## 7. The scatter

Written the obvious way -- one flat loop per slot, splitting the index
back into a component and a contraction with `mod` and `/` -- the scatter
costs eight integer divisions per output element. Profiled inside a Fock
build, `one_block` was then the largest single item in the run, 29% of CPU
and larger than any integral kernel. Walking component-then-contraction
directly and carrying the strides and normalisation products down the loop
nest removes every division, and it is what took the d classes from behind
Rys to ahead of it. The rotated-axis driver had the same scatter and the
same fix.

## 8. Status

* Coulomb (ab|cd) only, Cartesian and spherical, contracted and generally
  contracted, L shells. No range separation, no derivatives, no 3-centre.
* All 55 canonical classes through d are generated and pass `hgp_check`
  against the Rys path: 1,324,167 values, worst scaled difference 4.6e-14,
  tolerance 1e-12. It is a different algorithm, so it cannot be
  bit-identical.
* Public ABI: `libcint_2e_hgp_cart/sph` with the `libcint_2e_*` argument
  list (the optimizer argument is accepted and ignored) and
  `libcint_hgp_supported(shls, bas, nbas)`. No state between calls, so it
  is safe under an OpenMP loop over quartets. Which path a quartet takes is
  the caller's decision, not libfint's.

## 9. Timing, so far

Per quartet against `int2e_cart`, Cartesian, on a **generally contracted**
basis: carbon with 9 primitives into 3 s contractions and 4 into 2 p,
cc-pVDZ's shape, which is the case that separates the algorithms.

| total l | Rys µs | rotated-axis µs | HGP µs | HGP/Rys |
|---|---|---|---|---|
| 0 | 105.6 | 63.7 | 62.9 | 1.68 |
| 1 | 65.2 | 40.4 | 40.4 | 1.61 |
| 2 | 40.0 | 24.7 | 25.4 | 1.57 |
| 3 | 24.3 | 20.7 | 21.9 | 1.11 |
| 4 | 20.3 | 25.8 | 19.8 | 1.02 |

HGP is ahead of Rys at every angular momentum here and tracks the
rotated-axis path closely, overtaking it at total l 4 where that path
falls behind. Getting the d classes over the line took the scatter fix in
§7 rather than anything in the recurrences.

**Take any per-quartet harness with salt.** This one rated the rotated-axis
path at parity where a threaded Fock build measured 1.75x. It also missed
the staged-contraction defect in §5 entirely, because every basis it had
was segmented — the case only appears once a shell has more than one
contraction. A Fock build over a contracted basis is the verdict; nothing
should be routed by default until one says so.
