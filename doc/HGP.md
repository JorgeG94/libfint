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

## 5. Contraction

The contraction coefficients multiply in at the contraction rather than
riding in the primitive weight, which is what lets one vertical recurrence
serve every contraction column and every coefficient type:

    c(ef, bcol, kcol) += c_i c_j (bcol) · c_k c_l (kcol) · [e0|f0]

The geometric part — 2π^{5/2} exp(−μR_AB²) exp(−νR_CD²) / (p q √(p+q)) —
rides in the pair tables and is shared.

## 6. Numerics

The Boys function is the same guarded routine the rotated-axis kernels use,
for the same reason: libcint's closed form cancels catastrophically for F_1
at a tiny argument, which the Rys path never asks for.

## 7. Status

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

## 8. Timing, so far

Per quartet on C2H6 with 6-31G and 6-311G**, Cartesian, against
`int2e_cart`:

| total l | Rys µs | rotated-axis µs | HGP µs |
|---|---|---|---|
| 0 | 1.76 | 1.51 | 1.50 |
| 1 | 1.67 | 1.25 | 1.23 |
| 2 | 1.54 | 1.28 | 1.22 |
| 3 | 1.48 | 1.90 | 1.52 |
| 4 | 2.51 | 3.65 | 3.20 |

HGP matches the rotated-axis path at low angular momentum and is well
ahead of it where that path lost, but does not beat Rys above total l 2
here.

**This harness is not the verdict.** It rated the rotated-axis path at
parity where a threaded Fock build over a real basis measured 1.75x, and
the reason applies again: these quartets are small and the per-call cost
dominates. Two things also work against HGP specifically here. The d
shells in both bases are single primitives, so the contraction-before-
transfer split — the entire point — buys nothing on exactly the classes
being measured. And the transfer blocks have not been tuned at all.

The measurement that counts is a Fock build over a properly contracted
basis. Until that exists, nothing here should be routed by default.
