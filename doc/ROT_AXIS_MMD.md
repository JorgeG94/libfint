# Rotated-axis McMurchie–Davidson for (ab|cd)

A second route to the two-electron Coulomb integrals libfint otherwise
computes by Rys quadrature.  It exists for the low angular momenta — s, p,
the L shells of the Pople bases, and d — where a rotated frame collapses
most of the Hermite expansion and the primitive loop becomes a Boys
function and a few dozen multiply-adds.

The formulas are not written by hand.  `scripts/rotaxis_mmd/derive.py`
carries the derivation below as executable algebra, `factor.py` splits the
result into the three loop levels a kernel has, and `emit.py` writes one
module per class into `src/rotaxis/` plus the dispatcher
`src/cint_rotaxis_kernels.f90` and the CMake list beside the classes.  This note is the same derivation in prose,
so that the generated code can be read against something.

The hand-written GAMESS-lineage kernels in `gamess-libERI` (the `int0001`,
`int0101`, … files) compute exactly these quantities; the derivation
reproduces their accumulators term for term once the sign convention of the
Boys terms is accounted for (they carry `+2ρ F₁`, this derivation
`(−2α)ⁿ Fₙ`).  Those kernels were the check on the algebra before the Rys
path was.

## 1. The frame

Four centres A, B, C, D.  Put the origin at A, the z axis along A→B, and
choose x so that C→D lies in the xz plane:

    ẑ = (B − A)/|B − A|
    x̂ = w/|w|,   w = (D − C) − ((D − C)·ẑ) ẑ
    ŷ = ẑ × x̂

When A = B, ẑ is the lab z; when C→D is parallel to ẑ (or C = D), x̂ is the
lab axis least aligned with ẑ, projected.  The frame is always a proper
rotation, so p and d components transform as vectors and rank-2 tensors.

In this frame

    A = (0, 0, 0)          B = (0, 0, R_AB)
    C = (cx, cy, cz)       D = C + (Rs, 0, Rc)

with `Rs = R_CD sin γ`, `Rc = R_CD cos γ`, γ the angle between AB and CD.
In the code Rs and Rc are taken as `D_loc − C_loc` directly, so the only
approximation is dropping the y component of that difference, which is
rounding (≈ 1e-16 · R_CD).

## 2. What the frame buys

For a primitive bra pair with exponents a, b, p = a + b, the product centre
is

    P = (0, 0, za),   za = b R_AB / p,   zb = za − R_AB   (= P − B along z)

so `X_PA = X_PB = Y_PA = Y_PB = 0`.  For a ket pair c, d, q = c + d,
yc = c/q, yd = d/q,

    Q = C + yd (D − C) = (ax, cy, az),   ax = cx + yd Rs,   az = cz + yd Rc
    Q − C = yd (Rs, 0, Rc),   Q − D = −yc (Rs, 0, Rc)

so `Y_QC = Y_QD = 0`.  Every geometric quantity of a bra pair is one scalar
(za); of a ket pair, two (ax, az); the pair-independent geometry is three
constants (Rs, Rc, cy).  And P − Q, which the Coulomb integrals depend on,
is

    X_PQ = −ax,   Y_PQ = −cy,   Z_PQ = zq = za − az

where only zq varies over the bra primitives once a ket primitive is fixed.

## 3. The expansion

Helgaker, Jørgensen and Olsen, *Molecular Electronic-Structure Theory*,
eq. 9.9.33:

    (ab|cd) = (2π^{5/2} / (p q √(p+q))) K_ab K_cd
              × Σ_{tuv} E^{ab}_{tuv} Σ_{τνφ} (−1)^{τ+ν+φ} E^{cd}_{τνφ} R_{t+τ, u+ν, v+φ}(α, R_PQ)

with α = pq/(p+q), K_ab = exp(−(ab/p) R_AB²), the Hermite expansion
coefficients from the one-dimensional recursion (9.5.6–7)

    E^{00}_0 = 1
    E^{i+1,j}_t = (1/2p) E^{ij}_{t−1} + X_PA E^{ij}_t + (t+1) E^{ij}_{t+1}
    E^{i,j+1}_t = (1/2p) E^{ij}_{t−1} + X_PB E^{ij}_t + (t+1) E^{ij}_{t+1}

and the Hermite Coulomb integrals from (9.9.18–20)

    R^n_{000} = (−2α)^n F_n(α R_PQ²)
    R^n_{t+1,u,v} = t R^{n+1}_{t−1,u,v} + X_PQ R^{n+1}_{tuv}     (same for u, v)

All of it is polynomial in the frame quantities of §2 with integer
coefficients, and F_n only ever appears linearly.  The derivation script
builds each Cartesian component of a class as exactly this polynomial over
the symbols

    bra primitive:  ip = 1/(2p), za, zb, zq, B₀…B_L   (B_n stands for (−2α)^n F_n(T))
    ket primitive:  iq = 1/(2q), yc, yd, ax
    quartet:        Rs, Rc, cy

with L = la + lb + lc + ld the highest Boys order.

## 4. Factorisation into loops

Every term of every component is a product of a bra-primitive monomial, a
ket-primitive monomial and a quartet-constant monomial.  The kernel runs
ket primitives outside and bra primitives inside, so:

**Level 1, inside the bra loop.**  For each distinct bra monomial
`g(ip, za, zb) · B_n · zq^k` an accumulator

    S[g,n,k] += (c_a c_b K_ab / p) · g · B_n · zq^k

where B_n = (−2α)^n F_n(T) / √(p+q) and T = α (ax² + cy² + zq²).  The set of
(g, n, k) is what the class needs and nothing else; for (pp|pp) it is 43
accumulators over 5 distinct g.  Nothing here depends on which component
is being built, which is the whole saving: the inner loop is the Boys
function plus S.

**Level 2, once per ket primitive.**  Each component's polynomial, with the
S symbols substituted, is a polynomial in the ket symbols and the constants.
Group by constant monomial; the remaining factor P_j(iq, yc, yd, ax; S) is
accumulated with the ket weight:

    r[j] += (c_c c_d K_cd / q) · P_j

Components share r's whenever the factor is the same polynomial up to a
rational scale (a ket p_x and a bra p_x pull the same R_{100}, for instance).
(pp|pp) has 143 of them.

**Level 3, once per quartet.**

    (ab|cd)_local = 2π^{5/2} · Σ_j scale_j · (Rs^a Rc^b cy^c) · r[j]

The per-l normalisation libcint carries in `common_factor`
(`CINTcommon_fac_sp`) is applied by the driver, per component, so that an L
shell's s and p get theirs separately.

## 5. Shell kinds and L shells

A slot of the quartet is a *kind*, not just an l: s, p, d, or L.  An L shell
(libcint's `KAPPA_SP_SHELL`) is one slot of four components — s, then
p_x p_y p_z — whose s takes the s contraction coefficient and whose p take
the p one.  In the factorisation the coefficient type simply joins the S
and r keys, so an (L L|L L) kernel accumulates the s and p parts side by
side inside one primitive loop and one Boys evaluation.  The weights
`kab`, `kcd` carry one column per (type pair, contraction pair).

Classes are canonical: kinds ascend within each pair and the pairs ascend,
in the order s < p < L < d.  The driver permutes any quartet into that
order and scatters the result back through the permutation.

## 6. Back to the lab frame

Each index is transformed by the matrix that takes local Cartesian
components to lab ones: rot for p, the 6×6 built from products of rot's
entries for d (a lab x^i y^j z^k is the product of i factors of rot(1,:)·x,
j of rot(2,:)·x, k of rot(3,:)·x, expanded), 1 ⊕ rot for L.  For spherical
output the Cartesian-to-spherical matrix of `cint_cart2sph` is composed onto
it, so the four index transforms do both at once.

## 7. Numerics

Two things bit the first version and are worth knowing:

* `gamma_inc_like` (the Boys function libcint uses) switches to its closed
  form above `turnover_point(m)`, and that table is 0 for m = 0 and 1.  At
  T ≈ 1e-33 — which a quartet with two identical primitive pairs produces
  by rounding — the closed form's `(F₀ − e^{−T})/(2T)` cancels to 1e16.  The
  Rys path never asks it for F₁ there.  The kernels use the ascending
  series below T = 1 as well as below the table's turnover.

* A quartet that vanishes by symmetry gives the Rys path 1e-18 of rounding
  and this path an exact 0.  `rotaxis_check` scales its tolerance to the
  largest value in the quartet and floors that scale, for that reason.

`rotaxis_check` holds the path to the Rys one at 1e-12 scaled over every
quartet of the five reference bases up to the generated l, the adversarial
systems, and the packed L-shell basis of `l_shell_check`; measured, the two
agree to 5e-14.

## 8. Unrolled and table-driven kernels

A class below the generator's `--unroll-limit` (3000 ket-level terms, which
is every s/p/L class and the light d ones) is emitted fully unrolled: each
accumulator update is its own statement.  Above it the same three levels
are emitted as loops over DATA tables of indices and rational coefficients
-- (dd|dd) has 96k ket-level terms, and unrolled it is 32k lines gfortran
takes minutes over, while the tables compile in seconds.  GAMESS's CPU
kernels are loop-driven too; full unrolling was a GPU choice.  The limit is
a flag, so the benchmark can decide where the crossover really is.

## 9. Status and what is not here

* Coulomb (ab|cd) only, Cartesian and spherical, contracted and generally
  contracted, L shells.  No range separation, no derivatives, no 3-centre.
* Public ABI in `libfint.f90`: `libcint_2e_rotaxis_cart/sph` with the
  `libcint_2e_*` argument list (the optimizer argument is accepted and
  ignored) and `libcint_rotaxis_supported(shls, bas, nbas)`.  The path
  keeps no state between calls, so it is safe under an OpenMP loop over
  quartets exactly as the Rys path is.
* The committed kernels are every canonical class through d (55).  The
  per-quartet harness on C2H6 rated the d classes at 0.3-0.5x of Rys, but
  that harness also rated the s/p/L classes at parity where a threaded
  Fock build over a real basis (mqc, msn slice, 6-31G) measures 1.6x, so
  the d classes go in to be measured the same way.  The driver errors out
  beyond the generated classes rather than falling back;
  `rotaxis_supported` is the caller's question to ask first.
* Timed against `int2e_cart` over all quartets of 6-31G and 6-311G** on
  C2H6 (20 repeats each), Cartesian: total l 0 is 1.07x Rys, l 1 is 0.98x,
  l 2 0.81x, l 3 0.56x, l 4 0.57x.  Both paths sit near 1.7 us per quartet,
  which says the driver's fixed cost -- allocatable pair tables and an
  allocatable transform matrix per index, a generic scatter -- is what is
  being measured, not the kernel.  Fitting the driver into the workspace
  scheme the Rys path uses and specialising the s/p/L rotations is the next
  step, and the kernels themselves have not been profiled at all.
