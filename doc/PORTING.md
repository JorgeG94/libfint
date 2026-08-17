# Porting libcint to Fortran — a scoping study

Exploration only. Nothing here has been built. Numbers are line counts and
call-graph facts measured against this tree (`develop`, `3c78069`) and against
`../metalquicha` as it stands today.

The short version: the job is roughly **5× smaller than the repo's line count
suggests**, the hard part is concentrated in about 3,000 lines of Rys
quadrature, and the thing that decides whether this is clever or foolish is not
the porting cost at all — it is that you fork away from upstream forever.

---

## 1. Anatomy: what is actually in here

59,042 lines of C across `src/` and `include/`. That number is misleading in
three separate ways.

| Category | Lines | What it really is |
|---|---:|---|
| Machine-generated kernels (`src/autocode/`) | 23,743 | Emitted by `scripts/*.cl` (2,834 lines of Common Lisp) |
| Numeric data tables | ~9,200 | Emitted by `scripts/*.py`; Chebyshev/Wheeler/c2s coefficients |
| Hand-written logic | ~24,000 | The part a human actually has to think about |

The generated fraction matters enormously. `src/autocode/intor1.c` is 3,730
lines that nobody wrote — a Lisp program wrote them, from a symbolic
description of the integral. You do not port those. You retarget the emitter.

Same for the tables. `rys_wheeler.c` is 6,347 lines of which **5,392 are
numeric literals** — 85% of that file is a data table with a `.py` script
sitting in `scripts/` that produced it. Re-emitting it as a Fortran module with
`real(dp), parameter ::` arrays is a formatting change to a Python script, not a
port.

### The hand-written core, honestly measured

| File | Total | Table | Real logic |
|---|---:|---:|---:|
| `g2e.c` | 4,918 | 0 | **4,918** |
| `cart2sph.c` | 6,978 | 2,810 | 4,168 |
| `rys_roots.c` | 2,037 | 775 | 1,262 |
| `eigh.c` | 1,477 | 0 | 1,477 |
| `cint2e.c` | 1,229 | 0 | 1,229 |
| `rys_wheeler.c` | 6,347 | 5,392 | 955 |
| `cint3c2e.c` | 769 | 0 | 769 |
| `g1e.c`, `g1e_grids.c` | 1,054 | 0 | 1,054 |
| `fmt.c` | 580 | 40 | 540 |
| everything else in the closure | ~2,600 | — | ~2,600 |

`g2e.c` is the heart: the 2D vertical/horizontal recursion plus 70 hand-unrolled
`_g0_2d4d_XXXX` kernels (lines 695–2020) selected by a bitmask switch. Those 70
kernels are themselves mechanical — they are unrolled index arithmetic, exactly
the kind of thing that should be generated rather than written.

---

## 2. Scope collapse: what metalquicha actually needs

libcint declares **760 integral entry points**. Cross-referencing every
`libcint_*` and `LIBCINT_*` symbol referenced anywhere in `../metalquicha/{src,app,backends,test}`:

**Integral kernels actually called — 15 of 760:**

- 1e: `ovlp`, `nuc`, `kin` (cart + sph)
- 1e multipoles: `r`, `rr`, `rrr`
- 1e grids: `int1e_grids` (+ optimizer) — the ESP path
- 2e: `cint2e` (cart + sph) + optimizer
- 3c2e, 2c2e (cart + sph) + optimizers

**Helpers:** `cgto_cart/sph`, `tot_cgto_cart/sph`, `gto_norm`,
`del_optimizer`, and the `ATM_SLOTS`/`BAS_SLOTS` offset constants.

**Not used *today*:** every spinor and spin-orbit path, Breit, Gaunt, DKB,
LRESC, F12, 4c1e, 3c1e, and the entire gradient/Hessian autocode. The bridge is
explicit about the last one — `mqc_libcint_bridge.f90:110` errors out on
`want_gradient` with "the CPU backend has no gradients yet".

⚠️ **"Today" is doing real work in that sentence.** Relativistic capability is a
stated goal for metalquicha, which puts the spinor path back in scope. See §4.5
— it turns out to cost much less than its line count suggests, but it is not
zero, and it changes what the emitter has to handle from day one.

The gradient/Hessian autocode (13,453 lines) is genuinely deferrable, because it
is pure generator output and comes back for free once §4 is done.

**Dependency closure for what metalquicha uses: 32,380 of 57,165 `.c` lines** —
and of those, ~9,200 are tables and 3,730 are generated. So the genuine
hand-port surface is **roughly 17,000 lines**, and three of those files should
be deleted rather than ported:

- `eigh.c` (1,477) — a hand-rolled MRRR/dqds symmetric tridiagonal eigensolver,
  present only so libcint has no LAPACK dependency. metalquicha already has
  `pic_lapack_interfaces`/`pic_syev`. **Delete, call LAPACK.**
- `fblas.c` (220) — same story. metalquicha has `pic_blas_interfaces`. **Delete.**
- `c2f.c` (150) — the C↔Fortran shim. In a native Fortran libcint it does not
  exist. **Delete.**

Net: **~15,000 lines of real hand-porting** for everything metalquicha does today.

---

## 3. The genuine blockers, ranked

### 3.1 Extended precision in the Rys roots — the one real constraint

`CINTrys_roots` escalates by root count:

| nroots | Path |
|---|---|
| ≤ 5 | closed-form `rys_root1..5`, double |
| 6–8 | Jacobi/Schmidt, double → long double |
| 9–12 | long double (`CINTlrys_*`) |
| > 12 | `__float128` if available |

Two things make this less scary than it looks.

First, **the quad path is already optional upstream.** `rys_roots.c:35-39`:

```c
#ifndef HAVE_QUADMATH_H
#define CINTqrys_schmidt        CINTlrys_schmidt
#define CINTqrys_laguerre       CINTlrys_laguerre
#define CINTqrys_jacobi         CINTlrys_jacobi
#endif
```

Without libquadmath, everything falls back to `long double`. So quad is a
precision refinement, not a correctness requirement.

Second, **for metalquicha's basis sets you barely reach it.** For a 4-centre ERI,
`nroots = (li+lj+lk+ll)/2 + 1`. All-f is 7 roots; all-g is 9. You touch the
long-double path at g functions and the quad path essentially never.

In Fortran: use `real(real128)` from `iso_fortran_env` for both the `l` and `q`
ladders. That is *more* accurate than the C's 80-bit `long double`, and it
collapses two code paths into one. Available on gfortran and ifx — and
metalquicha's CMake already constrains itself to `GNU|Intel|IntelLLVM`
(`CMakeLists.txt:176`), so this costs you nothing you have.

⚠️ **nvfortran has no `real128`.** If a GPU path via nvfortran is ever the plan,
the high-nroots ladder needs a different answer there. Worth deciding early
rather than late.

### 3.2 Anonymous unions in `CINTEnvVars` — trivial

```c
union {FINT nfk; FINT grids_offset;};
union {FINT nfl; FINT ngrids;};
union {double *rl; double *grids;};
```

Fortran has no unions. These are pure readability aliases — the same slot
called two names depending on whether you are doing a 4-centre integral or a
grid integral. Pick one name per slot, add a comment. Non-issue.

### 3.3 Function pointers — clean in F2003, with one ordering trap

Three dispatch slots, declared K&R-style (`FINT (*f_g0_2e)();`) so C never
checks them. The actual callees *are* uniform per slot:

| Slot | Signature | Callees |
|---|---|---|
| `f_gout` | `(gout, g, idx, envs, empty)` | the `CINTgout1e_*`/`CINTgout2e_*` family |
| `f_g0_2e` | `(g, rij, rkl, cutoff, envs) → FINT` | `CINTg0_2e`, `_stg`, `_yp` |
| `f_g0_2d4d` | `(g, bc, envs)` | 12 variants |

So: three `abstract interface` blocks + `procedure(...), pointer` components.
Standard F2003.

The trap is module circularity — the interfaces take `type(CINTEnvVars)`, and
the gout kernels live in modules that need the envs type. Put the derived type
*and* the abstract interfaces in one base module that everything else uses. If
you discover this late it is a painful untangling; if you plan for it, it costs
nothing.

### 3.4 `CINTOpt` ragged pointer arrays — moderate

```c
FINT **index_xyz_array;   // LMAX1**4 pointers
PairData **pairdata;      // NULL = not-initialized, NO_VALUE = skippable
```

Fortran wants a wrapper type with an allocatable component. The genuine work is
that libcint overloads `NULL` to mean two different things ("never built" vs.
"built and screened out"). That needs an explicit `logical :: initialized` flag
rather than a sentinel. The result is *better* than the C, but it is a
restructure, not a transliteration.

### 3.5 The cache protocol — an opportunity, not a blocker

Every driver implements the two-call idiom: `if (out == NULL) return cache_size;`
then the caller allocates and calls again. It exists so libcint never mallocs
inside a hot loop and stays thread-safe.

Fortran maps this fine, and modern Fortran can do better — but **do not
naively replace it with automatic arrays.** `cint2e.c:805` computes cache sizes
that get checked against `INT32_MAX`; those buffers are large enough that
stack-allocating them per shell quartet will blow the stack on real molecules.
The right answer is a persistent workspace derived type allocated once per
thread. Design this deliberately.

### 3.6 Indexing — the actual source of bugs

The whole `g`-array machinery is flat index arithmetic (`g_stride_i`,
`g_stride_k`, `g_stride_l`, `g_stride_j`, `g_size`), and the `atm`/`bas`/`env`
slot constants are 0-based throughout.

**Strong recommendation: keep 0-based in the port.** Declare `env(0:)`,
`atm(0:,0:)`, `bas(0:,0:)` so every line of ported code has a
literal one-to-one correspondence with the C you are porting from. Across ~5,000
lines of index arithmetic, off-by-one is the failure mode that will cost you
weeks, and it produces *plausible* wrong numbers rather than crashes. Flip to
1-based later if you want, as a separate mechanical pass with the test suite green.

### 3.7 Things that are non-issues

- **`complex.h`** — used only in spinor paths metalquicha does not touch, and
  Fortran has native `complex` anyway.
- **`gout2e_simd.c`** (44K of intrinsics) — not referenced in `CMakeLists.txt`
  at all. Ignore it.
- **No threading in libcint** — no `#pragma omp`, no pthreads anywhere in
  `src/`. Parallelism is entirely the caller's problem. One less thing.
- **`static inline`** → module procedure; the compiler inlines within a module.
- **`double ai[1]`** — size-1 arrays used as addressable scalars. Just scalars.

---

## 4. The leverage point: retarget the Lisp generator

This is the part that turns "insane" into "actually tractable".

```
scripts/parser.cl      684 lines   symbolic integral expression parsing
scripts/derivator.cl   409 lines   symbolic differentiation
scripts/gen-code.cl  1,447 lines   emits C
scripts/auto_intor.cl  247 lines   the list of integrals to generate
scripts/utility.cl      47 lines
```

`parser.cl` and `derivator.cl` — 1,093 lines — are **language-agnostic
symbolic math**. They do not care what you emit.

`gen-code.cl` is the emitter, and it is C-flavoured `format` strings:

```lisp
(format fout "gout[n*~a+~a] +=~a;~%" comp gid s)
```

Retargeting that to Fortran is maybe a 30–40% rewrite of one 1,447-line file.
What you buy: **all 760 entry points regenerate as Fortran**, instead of hand-porting
23,743 lines of generated C. And the 70 unrolled `_g0_2d4d_*` kernels in `g2e.c`
should move under the generator too, rather than being ported by hand.

Do this and the port stops being "translate libcint" and becomes "port the
runtime, regenerate the kernels".

---

## 4.5 The relativistic path — the generator already knows the physics

This is the part that most rewards owning the generator, and it is worth being
precise about, because the line counts mislead badly in *both* directions.

### All of it is generated

| File | Lines | Origin |
|---|---:|---|
| `autocode/breit1.c` | 1,616 | `gen-code.cl` |
| `autocode/dkb.c` | 1,036 | `gen-code.cl` |
| `autocode/lresc.c` | 1,056 | `gen-code.cl` |
| `autocode/gaunt1.c` | 965 | `gen-code.cl` |
| `cint3c1e_a.c` | 675 | `gen-code.cl` |
| `cint1e_a.c` | 317 | `gen-code.cl` |

*(I initially filed the two `_a` files as hand-written relativistic runtime.
They are not — both carry the `generated by gen-code.cl` header.)*

**5,665 lines of relativistic kernels, none of it hand-written.** Nobody ports
those. They regenerate.

### The physics lives in the Lisp, not the C

`scripts/parser.cl` is a symbolic operator algebra that already understands the
relativistic vocabulary:

```lisp
(defparameter *two-electron-operator* '(r12 ccc2e nabla-r12 gaunt breit-r1 breit-r2))
(defparameter *nabla-not-comutable*   '(rinv nuc grids nabla-rinv r12 nabla-r12
                                        gaunt breit-r1 breit-r2))
(defparameter *act-left-right*        '(nabla-rinv nabla-r12 breit-r1 breit-r2))
;; p = -i∇, carried symbolically as a complex coefficient:
((p) (make-vec (make-cell #C(0 -1) '() '(nabla x)) ...))
((sigma) ...)   ; the Pauli vector
((g) ...)       ; GIAO factor, i/2 (R_m - R_n) × r0
```

And `auto_intor.cl` declares integrals in that DSL rather than in code:

```lisp
'("int1e_spsp"               (sigma dot p \| sigma dot p))
'("int2e_ssp1ssp2"           ( \, sigma dot p \| gaunt \| \, sigma dot p))
'("int2e_gauge_r1_ssp1ssp2"  ( \, sigma dot p \| breit-r1 \| \, sigma dot p))
```

Commutation rules, the σ·p algebra, Gaunt and Breit gauge terms, GIAO factors —
all of it is in 1,093 lines of language-agnostic Lisp (`parser.cl` +
`derivator.cl`) that a Fortran port does not touch. **Your instinct is right:
port the emitter well and the relativistic integrals are not a separate
project.**

### What spinor support actually costs at runtime

The generator emits `_cart`, `_sph`, and `_spinor` from a single integral
description, and they differ only in which cart-to-spherical function pointer
gets handed to a shared driver (`gen-code.cl:319-329, 498-539`):

```lisp
(format fout "return CINT1e_drv(out, dims, &envs, cache, &c2s_cart_1e, ~d);
(format fout "return CINT1e_drv(out, dims, &envs, cache, &c2s_sph_1e,  ~d);
(format fout "return CINT1e_spinor_drv(out, dims, &envs, cache, ~a, ~d);
```

So spinor is nearly free *at the emitter*, and the cost is a fixed, bounded
chunk of hand-written runtime:

| Piece | Lines | What |
|---|---:|---|
| `cart2sph.c` spinor section (3920–6978) | ~3,050 | `cart2spinor` transforms; `_sf`/`_si`/`_zf`/`_zi` variants; 59 `double complex` entry points |
| `breit.c` | 317 | Gauge assembly — combines the `breit-r1`/`breit-r2` terms; hand-written |
| `CINT{1e,2e,3c2e,2c2e}_spinor_drv` | ~250 | Complex-output drivers |
| `kappa`/spinor counting in `cint_bas.c` | ~60 | `CINTlen_spinor`, `CINTcgto_spinor`, offsets |

**~3,700 lines of runtime, once, and it unlocks all 190 spinor entry points.**
That is a far better ratio than any other part of this project.

### Two notes, one good and one annoying

**Fortran is genuinely better here.** The C hand-splits complex values into
parallel real/imaginary arrays (`gspR`/`gspI`, coefficient tables
`cart2j_gt_lR`/`cart2j_gt_lI`) because C's `double complex` optimises poorly and
the transform tables are stored split. Fortran has native `complex(dp)` and lets
you pick either representation — split for vectorisation in the hot transform,
native for the API. This code should come out *shorter and clearer* in Fortran,
which is not something I can say about `g2e.c`.

**`CINT2c2e_spinor_drv` is not implemented upstream** (`cint2c2e.c:294-298` —
it prints `"not implemented"` and bails). If the relativistic density-fitting
path matters to you, that is a gap you would be writing yourself either way. Mildly
in favour of owning the code.

### Consequence for the plan

If relativistic is the destination, **do not defer the complex/spinor output
path to a later phase.** Retrofitting complex returns and the `_sf`/`_si`
dispatch into an emitter and a driver layer designed only for real output is the
kind of rework that costs more than doing it up front. Design the driver
signatures and the emitter's c2s selection to carry the spinor case from the
start, even if you only wire up `_cart`/`_sph` initially.

---

## 5. Proposed phasing

**Phase 0 — the oracle. Mostly already written; two gaps to close.**

`examples/fortran/fortran_time_c2h6.F90` already has the hard part: the correct
unique shell-quartet traversal (`i≥j`, `ij≥kl`), OpenMP structure, buffer
sizing via `CINTcgto_spheric`, and **five progressively larger basis setups —
6-31G, 6-311G\*\*, cc-pVDZ, cc-pVTZ, cc-pVQZ** (`setup_*_basis`, lines
310–861). That is l = 0..4, s through g. Better coverage than a from-scratch
harness would likely have got, and the g functions matter specifically: an
all-g quartet needs 9 Rys roots, which reaches the long-double ladder in §3.1.
`fortran_time_c2h6_pure.F90` is the same sweep through the high-level API.
Reuse this rather than writing a new one.

Two gaps, though:

**(a) It discards every buffer.** Lines 174–176, 225–227, 276–278 are all
`allocate → call → deallocate` with `buf` never read. It is a pure timing
benchmark, not a validating one. Retaining results and comparing is a small
edit, but nothing is being checked today.

**(b) ⚠️ The Python suite validates by a permutation-invariant fingerprint.**
`test_cint.py:149-166` accumulates `abs(buf).sum()` across all quartets into a
scalar and compares to a stored `vref`:

```python
def close(v1, vref, count, place):
    return round(abs(v1-vref)/count**.5, place) == 0
```

For regression-testing the C that is fine. **For validating a port it is the
wrong instrument**, because `abs().sum()` is invariant under reordering *within*
a buffer — so a transposed or mis-strided result (§3.6: the single most likely
porting bug, and the one that produces plausible wrong numbers rather than
crashes) passes silently. The suite also needs PySCF and only ever exercises the
C side.

The pattern you want is already in that file — `test_comp1e_spinor`
(`test_cint.py:203-233`) compares two implementations **elementwise**
(`dd = abs(op - op_ref)`) and reports the offending shell pair and index. That
is the check to generalise.

**So Phase 0 is: take the C2H6 sweep, retain the buffers, call C libcint and
Fortran libcint side by side, and diff elementwise to ~1e-12 reporting
`(shell quartet, flat index)` on failure.** Days rather than weeks, and mostly
editing code that exists. Still do it first — but it is cheaper than I first
estimated.

The one thing the existing setups do *not* stress is numerically nasty geometry:
near-coincident centres, very diffuse and very tight exponents in the same
shell pair. Worth adding a handful of adversarial cases, since that is where the
Rys escalation ladder actually gets exercised.

**Phase 1 — tables (~1 week).** Re-emit `cart2sph`, `rys_wheeler`,
`polyfits`, `roots_for_x0.dat` as Fortran `parameter` modules by editing the
existing Python. Verify by loading both and diffing. Low risk, gets ~9,200
lines done almost for free, and builds confidence.

**Phase 2 — the Rys core (the risk concentration).** `rys_roots.c`,
`rys_wheeler.c`, `fmt.c`, `find_roots.c`, `polyfits.c` — ~3,300 lines of real
logic. Numerically delicate, has dedicated tests (`test_rys_roots.py`,
`test_rys_wheeler.py`), and is where the `real128` decision lands. **Budget
~40% of total project risk here.** It is also fully self-contained — you can
port and validate it against the C with no other Fortran in place.

**Phase 3 — support layer.** `cart2sph` transform loops, `cint_bas`, `misc`,
`optimizer`. Delete `eigh`/`fblas`/`c2f` in favour of `pic_lapack`/`pic_blas`.

**Phase 4 — the engine.** `g1e`, `g1e_grids`, `g2e`, and the `cint1e`/`cint2e`/
`cint3c2e`/`cint2c2e` drivers. The core, and where the indexing discipline from
§3.6 earns its keep. **Design the driver signatures for complex output now**
(§4.5), even though nothing uses it yet.

**Phase 5 — retarget `gen-code.cl`, early and on a small target.** Do not save
this for the end. As soon as one driver works, retarget the emitter and generate
*one* integral family (say `intor1.c`'s `int1e_kin`) rather than hand-porting
it. That tells you whether the emitter strategy holds while the blast radius is
still small. Everything after this point gets dramatically cheaper — this is the
phase that decides whether the project is 6 months or 18.

**Phase 6 — the spinor runtime** (~3,700 lines, §4.5): the `cart2spinor`
transforms, the complex drivers, `breit.c`, kappa counting. Unlocks all 190
spinor entry points via the generator.

**Phase 7 — regenerate everything else.** Gradients, Hessians, Gaunt, Breit,
DKB, LRESC, F12. This is now a matter of running the generator and validating,
not of writing code.

**Phase 8 — optimizer/screening** and performance work.

### Honest cost

For a person fluent in both the physics and both languages, working full time:
**4–7 months** to bitwise-comparable results across metalquicha's current
non-relativistic subset. Phases 0–2 are ~6 weeks of that and de-risk most of the
rest.

Add roughly **6–10 weeks** for Phase 6 to reach relativistic parity — and note
that this is where the estimate is most favourable, because 5,665 lines of
relativistic kernels arrive from the generator rather than from a keyboard.

Call it **6–9 months to a Fortran libcint with relativistic capability.** The
marginal cost of relativistic *on top of* the non-relativistic port is small,
which is unusual and worth exploiting.

The estimate assumes Phase 0 exists. Without a differential oracle, multiply by
something unpleasant.

---

## 6. What you gain, what you lose

**Gain**

- Single-language build. No C toolchain, no `iso_c_binding` marshalling at the
  hot boundary, one set of flags, one debugger story.
- **Bounds checking on the `g` array.** Given what §3.6 says about index
  arithmetic, `-fcheck=bounds` on a debug build of the integral engine is worth
  a genuinely large amount.
- **A real shot at `do concurrent` / OpenMP-target over the shell-quartet
  loop.** libcint has no threading at all, and the kernels are pure flat-array
  arithmetic with no aliasing. This is the most interesting prize on the list,
  and it is exactly the pattern you have GPU-ported before.
- Fortran's dummy-argument non-aliasing rules are stronger than C's. The
  unrolled 2d4d kernels plausibly get *better* codegen, not worse.

**Lose**

- **Upstream. This is the real cost and it dwarfs the porting effort.** libcint
  is actively developed and PySCF-adjacent; every numerical fix, new integral
  type, and performance improvement stops flowing to you the day you fork. You
  are signing up to own an integrals library forever.
- 745 entry points you would not have on day one.
- The C is battle-tested by an enormous user base. Yours would not be.

**Performance read:** no reason Fortran should be slower. The one intrinsics
file is not even in the build. The hot loop is flat indexing and FMA chains,
which both compilers handle identically. The risk is not the language — it is
losing the hand-tuned cache/workspace discipline (§3.5) by being careless.

---

## 7. Recommendation

The engineering is tractable and the scope is far smaller than the repo looks.
Phases 0–2 would tell you within six weeks whether you want the rest.

One thing worth stating plainly, because it sharpens the argument rather than
weakening it: **wanting relativistic capability is not by itself a reason to
port.** libcint already has all 190 spinor entry points, Gaunt, Breit, DKB and
LRESC. If the goal is "metalquicha can run DHF or X2C", the cheapest path by an
enormous margin is to expose the existing C entry points through the interface
you have been extending — days of work, not months.

The real argument is the one you made, and it is about the generator, not the
integrals:

> *if we port the generator well enough we should be able to do anything*

That is the correct read of this codebase. What libcint actually is, under the
23,743 lines of emitted C, is a **symbolic integral compiler** — 1,093 lines of
operator algebra that knows σ·p, Gaunt, Breit gauge terms, GIAO factors and
their commutation rules, plus a 1,447-line backend that happens to emit C.
Owning that in your own toolchain means you can emit integrals that *do not
exist upstream*, in the language and for the target you choose. Nothing you can
get by linking the C gives you that.

So the three framings, honestly ranked:

- **Language purity alone** — a permanent maintenance liability bought for an
  aesthetic property. Doesn't pay for itself.
- **Access to relativistic integrals** — real goal, wrong mechanism. Bind to the
  C; it is already there and already correct.
- **Owning a symbolic integral generator that emits Fortran for the backend of
  your choice** — this is the one that justifies the cost, and relativistic
  capability is strong *evidence* for it rather than the reason itself. It is
  also where the GPU story lives: libcint has no threading and no device path,
  so an emitter that can target `do concurrent` is going somewhere upstream is
  not. You would not be duplicating libcint. You would be building the thing
  libcint's generator could have been.

The fork cost from §6 is unchanged and still real — but under this third framing
it is a price for a capability rather than a tax on a rewrite.

**Suggested next steps, in order:**

1. **Close the two gaps in the existing harness** (§5, Phase 0) — retain the
   buffers in the C2H6 sweep and diff elementwise rather than by fingerprint.
   The traversal and the s-through-g basis coverage are already written. Days,
   not weeks. Useful immediately for validating the existing Fortran bridge, and
   a prerequisite for everything else. Do this even if the port never happens.
2. **Spend a week on `gen-code.cl` before committing to anything.** Read the
   emitter, and try to make it produce Fortran for one trivial integral. The
   entire strategic case rests on that file being retargetable, and that
   assumption is cheap to test and expensive to be wrong about. If the emitter
   turns out to be too entangled with C semantics to retarget cleanly, the whole
   project reverts to "hand-port 15,000 lines" and the answer changes.
3. Decide after Phase 2, when the Rys core has told you what the real numerical
   difficulty is.

Steps 1 and 2 together are about two weeks and retire most of the risk in the
estimate. That is a good trade for a decision of this size.

---

## Appendix: measured facts

- `develop` @ `3c78069`; metalquicha as of 2026-08-14.
- 59,042 lines C in `src/` + `include/`; 57,165 in `.c` alone.
- 760 integral entry points in `include/cint_funcs.h`
  (110 × 1e, 61 × 2e, 12 × 3c2e, 4 × 2c2e, 3 × 3c1e families, × cart/sph/spinor).
- Dependency closure for metalquicha's usage: 32,380 lines; out-of-scope: 24,785.
- `src/autocode/`: 23,743 generated lines from 2,834 lines of Common Lisp.
- `g2e.c`: 70 unrolled `_g0_2d4d_*` kernels, lines 695–2020.
- No `#pragma omp`, no pthreads, in `src/`.
- `gout2e_simd.c` is not referenced by `CMakeLists.txt`.
- `MXRYSROOTS 32`, `ANG_MAX 15` (`include/cint.h.in:123-125`).
- 190 `_spinor` entry points in `include/cint_funcs.h`.
- Relativistic kernels are 100% generated: `breit1.c` (1,616), `dkb.c` (1,036),
  `lresc.c` (1,056), `gaunt1.c` (965), `cint3c1e_a.c` (675), `cint1e_a.c` (317)
  — 5,665 lines, all carrying the `generated by gen-code.cl` header.
- Hand-written relativistic runtime: `cart2sph.c:3920-6978` (~3,050),
  `breit.c` (317), the four `*_spinor_drv` functions (~250), kappa counting in
  `cint_bas.c` (~60).
- `CINT2c2e_spinor_drv` is a stub that prints "not implemented"
  (`cint2c2e.c:294-298`).
- `parser.cl` operator vocabulary: `ovlp rinv nuc nabla-rinv ccc1e grids` (1e);
  `r12 ccc2e nabla-r12 gaunt breit-r1 breit-r2` (2e); plus `sigma`, `p = -i∇`,
  GIAO `g`, and commutation/act-left-right rules.
- `fortran_time_c2h6.F90` basis setups: 6-31G, 6-311G\*\*, cc-pVDZ, cc-pVTZ,
  cc-pVQZ (lines 310–861) → l = 0..4. Buffers allocated and freed unread at
  lines 174–176, 225–227, 276–278.
- ctest registers only `cinttest` (`test_cint.py`) and `cint3c2etest`
  (`test_3c2e.py`), both gated on finding a Python interpreter
  (`CMakeLists.txt:238-245`).
- `test_cint.py:149` validates via `abs(v1-vref)/sqrt(count)` on an
  `abs(buf).sum()` fingerprint — permutation-invariant. `test_comp1e_spinor`
  (line 203) is the one elementwise comparison in the suite.
