# Closing the gap to the C — method, findings, and what is left

This is a working blueprint, not a report. It exists so the next person
attacking the port's performance does not repeat the two weeks of guessing
that produced the first half of it.

Everything here is measured on one target — `int2e_sph`, the four-centre ERI,
over the five reference bases of `cint_test_systems` — because a single class
driven to the bottom teaches more than six classes sampled. The method
generalises; the specific numbers do not.

**State: 1.128× like-for-like at `-O3`; 1.10× with PGO against a stock C build.**
The port was 1.64× at the true Release baseline.

---

## 1. Do not profile with a clock

The single most useful change was to stop timing and start counting.

```
valgrind --tool=callgrind --callgrind-out-file=cg.f.out ./bench port 1
valgrind --tool=callgrind --callgrind-out-file=cg.c.out ./bench c 1
callgrind_annotate --threshold=90 cg.f.out
callgrind_annotate --auto=yes cg.f.out    # per source line
```

Instructions retired, first measurement:

| | instructions | ratio |
|---|---:|---|
| C | 2,770,553,559 | — |
| port | 3,584,701,617 | 1.294× |

The wall-clock ratio at the time was 1.262×. **They agree.** That single fact
reframed the problem: the port was not stalling, not missing cache, not
mis-scheduling. It was executing 29% more instructions. Instructions can be
found and removed; "it feels slow" cannot.

Three practical notes:

- Callgrind is ~50× slower than native. One rep of the benchmark is enough —
  the counts are exact, not sampled, so there is no averaging to do.
- **Function totals are trustworthy; line attribution is not.** Twice a line
  that looked expensive turned out to be the loop overhead folded onto it, and
  a "fix" changed the instruction count by zero. Always confirm at the
  function level.
- VTune's software sampling was what found the first two wins, but it cannot
  do this. On this Broadwell its hardware sampling does not run at all.

---

## 2. The pattern that produced most of the wins

**libcint hand-specialises small angular momenta, and a faithful transcription
of the general arm is not a faithful transcription.**

The C is full of `switch` statements that write out `l = 0, 1, 2` (and `mi =
1, 3, 5, 6, 7`) as straight-line code with immediate offsets, falling through
to the general loop only above that. Those are the cases that dominate every
real basis set. Port the general arm alone and you are slower on 90% of the
work and level on the rest.

The tell is unmistakable once you know to look for it: **isolate the kernel and
sweep `l`.** If the port is level with the C at high `l` and behind at low `l`,
a specialisation is missing. `CINTg2e_index_xyz`:

| | l=0 | l=1 | l=2 | l=3 | l=4 |
|---|---|---|---|---|---|
| before | 0.72× | 1.54× | 1.75× | 1.06× | 1.09× |
| after | 0.67× | 1.11× | 1.18× | 1.01× | 0.96× |

Level at f and g both times — that is where the C stops specialising.

Find them all with:

```
grep -n "switch (" src/*.c
```

At the time of writing, the ones that matter are in `CINTgout2e` (on
`nrys_roots`), `CINTg2e_index_xyz` (on `i_l`), `dcopy_iklj` (on `mi`),
`sph2e_inner` (on `l`), and the `_2d4d_unrolled` pair (on `type_ijkl`). All
are now ported. **Check this list again for any arity you work on next** — the
one-, three- and two-centre paths have not been swept this way.

---

## 3. The other pattern: reaching through a pointer that the C hoisted

The C caches base pointers outside the loop nest:

```c
double *ai = env + bas(PTR_EXP, i_sh);   /* once */
... envs->ai[0] = ai[ip];                /* one load */
```

The port reached through the `envs%env` pointer component every time:

```fortran
envs%ai = envs%env(ai_o + ip)   ! descriptor, offset, load
```

One line, measured: **30.4M instructions against the C's 5.0M** for the same
assignment. Binding the rows once — `aip(0:) => envs%env(ai_o:)` — closed it.

Related, same cause: dummy arguments that are always small and fixed should be
**explicit-shape**, not assumed-shape. `rij(0:)` made `rij(2)-rkl(2)` cost 10
instructions to the C's 3. `rij(0:2)` costs what the C costs.

### 3a. The same pattern inside the `g` array, and in the generator

`g` is one flat array holding several blocks. The C opens every kernel that
touches more than one with

```c
double *g0 = g;
double *g1 = g0 + envs->g_size * 3;
... g1[ix+0] ...
```

The port wrote `g(g1+ix+0)` — an extra add on every access, and a gout term
is three of them. Declare the dummy `target`, bind one
`real(dp), pointer, contiguous :: g1p(:)` per block **after** the operation
sequence that assigns the offsets, and read `g1p(ix+0)`.

This is wired into the backend rather than applied by hand: `gref`,
`emit-gptr-decls` and `emit-gptr-binds` in `scripts/f90-backend.cl` do it for
all four gout emitters, above no more than `*f90-gptr-max*` blocks. 174 of
the generated kernels take it. `int2e_ip1_sph` 1.346x → 1.252x.

The same shape appears by hand in the 2D recursions — `cint_g1e_ovlp`,
`cint_g1e_nuc`, `cint_g0_2e_2d` and the four `_2d_4d` transfers. Bind
`gxp/gyp/gzp` there. `int1e_ipnuc_sph` 1.254x → 1.189x, almost all of it
from `cint_g1e_nuc` alone: 137.7M instructions → 124.0M.

### 3c. Whether the hoist pays depends on the arity, so measure per loop

All six primitive loops now hoist their exponent and coefficient rows, bind
the gout kernel's three arguments, and expand `prim2ctr` inline. What that
was worth is not uniform, and the reason is the shape of the nest:

| loop | before | after |
|---|---|---|
| `cint_3c1e_loop` | 1.494x | 1.273x |
| `cint_3c2e_loop_nopt` | 1.264x | 1.172x |
| `cint_1e_loop` | 1.254x | 1.189x |

`cint_3c1e_loop` gained most because its innermost body names the three
exponents six times. `cint_1e_loop` gained almost nothing from the exponent
hoist specifically — **25,040 instructions out of 506 million**, which is
noise — and that part is deliberately not there; the win in the 1e path came
from binding the gout arguments and from §3a. A two-deep nest wrapped around
a whole 2D recursion has nothing to save. Do not assume the pattern
transfers; each loop is a separate measurement.

### 3d. Assumed-size at a procedure-pointer boundary

The largest single win in the whole exercise, and the least obvious.

The gout kernel is reached through `envs%f_gout`, matching the C's
`(*envs->f_gout)(...)`. Neither side can inline it; both pay a real
indirect call. So the per-call overhead §7b measured — loop bodies
instruction-identical, all the difference at the boundary — could not be a
macro-versus-procedure tax. It was **what crosses the boundary**.

The C passes three bare pointers. `gout(0:)` and `idx(0:)` are
assumed-shape: the caller builds an array descriptor, the callee unpacks
one, on a call that happens once per primitive quartet — and for
`int1e_nuc` once per *nucleus* per primitive pair.

```fortran
real(dp), intent(inout) :: gout(0:*)   ! address, nothing else
integer,  intent(in)    :: idx(0:*)
```

`int1e_ipnuc_sph` went 1.189x → **1.048x** on that change alone. On that
path the port's gout and its nabla together are now *ahead* of the C's,
221.6M instructions to 226.4M.

Two things to know before reaching for it. An assumed-size array cannot be
passed as a whole array, so a forwarding routine passes `gout(0)` and lets
sequence association carry the address — the F77 idiom, and here the
correct one. And the extent is not lost, only unenforced: `gout` is
`nf*ncomp`, `idx` is `nf*3`, both derivable from the `envs` that travels
with them. Say so at the interface, which is where a reader looks.

**Not taken: the same change for `g`.** After the above, the gout boundary
costs nothing — the port is ahead of the C there. `g` would cascade into
sixteen operation routines, need an explicit upper bound at every pointer
binding, and weaken the interface everywhere for a boundary that is no
longer costing anything. Measured first, then declined.

A preprocessor — `cpp` via `-cpp`, or fypp — was considered for this and
is not needed. The cost was never macro expansion; you cannot expand an
indirect call in either language. It was the descriptors.

### 3e. Costs that only matter on cheap integrals

Everything above was found on `int2e_sph`, which is an expensive integral.
The cheap ones — `int2c2e`, `int3c1e`, `int1e_ipovlp` — were the worst in
the port for a while, and for a different reason: fixed per-call overhead
that an expensive integral amortises and a cheap one does not.

**Declaring the envs.** `type(cint_env_vars) :: envs` default-initialises
52 components at every declaration; `CINTEnvVars envs;` costs nothing.
Measured at 17.6M on the two-electron path, which is 0.6% and was declined
there. On `int2c2e` the same line is **4,352,000 instructions, 3.8% of the
run and 27% of the whole gap**. Now removed: no component carries a
default except `g0_kind`, which has no C counterpart and is assigned only
by the F12 setup.

The reason this was declined the first time is worth keeping. The risk is
that a component gets read before it is written and the sweep passes
anyway, because undefined memory usually reads as zero. Do not argue
about that — test it. Build a second time with

```
-finit-derived -finit-integer=-99999 -finit-real=snan -finit-logical=true
```

which poisons exactly the components that lost their initialisers, and run
the whole suite against it. And check the test has teeth: dropping
`g0_kind`'s default as well makes `catalogue_check` crash outright.

**Zeroing a buffer the C does not zero.** None of `CINT1e_drv`,
`CINT3c2e_drv`, `CINT2c2e_drv`, `CINT3c1e_drv` or `CINT1e_grids_drv`
memsets `gctr`; each threads an `empty` flag so the first contraction sets
where later ones accumulate. The port threaded the same flag *and* memset,
at five sites. cachegrind put 24% of the port's writes in the driver on
`int2c2e`; this was most of them. Same discipline: a memset removal passes
by luck when the buffer happens to hold zeros, so patch `ws_alloc_d` to
fill every allocation with a quiet NaN and re-run before believing it.

**An unused allocatable component is not free.** `integer, allocatable ::
idx(:)` in `cint_env_vars` was referenced by nothing, and cost a
nullify-on-entry and a deallocate on every return path of every routine
that declared an envs.

---

## 3b. Fortran-only costs with no C counterpart

These are the highest-yield findings in the whole exercise, because the C
never pays them at all and so nothing in a line-by-line reading of the C
suggests looking for them. Two of the three were found by noticing library
calls in the profile that had no business being there.

**Automatic arrays with runtime bounds go on the heap.** `integer ::
zeroidx(0:ictr-1)` in `cint_non0coeff_byshell` was calling malloc and free
**414,720 times** — four per shell quartet, about 75 million instructions.
The C writes `FINT zeroidx[ictr]`, a C99 VLA, and pays a stack adjustment.
Removing the scratch entirely (two passes over the coefficients instead of
one plus a buffer) was worth **2.2% on its own**. Sweep for these: any local
array whose bound is a dummy or a local is one.

**A section assignment with the same array on both sides allocates a
temporary.** `ws%d(o2:o2+n-1) = ws%d(o1:o1+n-1)` cost a malloc and a memcpy
per call, 14,766 of them, because gfortran cannot prove the ranges disjoint.
An explicit loop says what the programmer knows.

**Copying where the C repoints.** `c2s_sph_2e1` threads a `tmp1` pointer
through four transform stages: each returns either its output buffer or its
input untouched, and the next stage reads whichever came back. The port
copied to normalise the buffer instead. Thread an offset.

**A macro is not a subroutine.** `PRIM2CTR` is a macro in the C. As a
Fortran subroutine its wrapper — prologue, empty-flag test, epilogue — cost
32 million instructions over 350,616 calls that the C does not pay. Expand
it at the call sites. Note `!GCC$ ATTRIBUTES always_inline` is rejected by
gfortran 13 for a module subroutine, so this is by hand. gfortran will not
inline it on its own either: it takes `ws` by reference and writes `flag`.
All four `prim2ctr` wrappers in the port are now gone; on `int3c2e_sph` the
one in `cint_3c2e` was 18.3M instructions that vanished entirely.

**Declaring the envs costs 170 instructions.** `type(cint_env_vars) :: envs`
default-initialises 52 components that the init routine immediately
overwrites; `CINTEnvVars envs;` costs nothing. Measured at 17.6M. **Not
taken**: only 25 of the 52 are set by every one of the base init routines,
and `g0_kind` — read by the two-electron loop — is set only by the F12
setup. An undefined read there would very likely pass every sweep, because
undefined usually reads as zero. Half a percent is not worth that.

## 4. Things that are not the problem (measured, do not retry)

Each of these looked obviously right and each was worth nothing or worse.

| idea | result |
|---|---|
| `-march=native` | **slower**, and breaks bit-identity in 371,752 values — gfortran contracts FMAs where gcc does not |
| `-funroll-loops` | 3.213 s against 2.825 s |
| `-mno-gather` | 3.365 s — worse than the blunt `-fno-tree-vectorize` |
| assumed-size dummies on the gout kernels | **exactly 1 instruction per call**, twice measured, once with inlining possible and once through a procedure pointer. Would have meant changing 197 generated kernels |
| `contiguous` on **dummy arguments** | 30% worse — gfortran emits copy-in/copy-out |
| `contiguous` on **pointer declarations** | **worth 0.6%** — not the same thing as the row above, and the distinction matters: on a pointer it is a promise the compiler can use, on a dummy it forces a repack |
| striding `kij` alongside the loop instead of `jp*i_prim+ip` | zero |
| writing the `rkl` three-vector componentwise | zero |
| replacing the envs setup's four coordinate section assignments, each a memcpy call | 5.7M, below noise |
| binding ws%d / ws%i / envs%env once and passing the pointers | 0.5% **worse** — the extra pointers cost more register pressure in the loop than the descriptor builds they save |
| `contiguous` on the **gout kernels' dummies** | **produces wrong answers.** 366,672 values over tolerance, worst 6.9e56. At least one caller's actual cannot be proven contiguous, gfortran inserts a repack, and the kernels' writes into `g` — the derivative forms build their intermediates there — do not survive it. The 0.5% it appeared to buy was measured on a broken build |
| `contiguous` **dummies** with contiguous actuals | still 1 instruction per call. The disassembly does show the stride test disappear, and assumed-size removes the base loads too — but the caller pays back what the callee saves, three separate measurements agree |
| one-array-plus-offsets through the Rys layer | the `g(w:)` section costs 0.27 ns a call |
| hoisting `nf*j` out of the contraction inner loop | zero instructions — gfortran had already done it |
| `-ftree-vectorize` on the contraction kernels | zero — the one-array idiom blocks vectorisation regardless |

The lesson in the middle rows: **a microbenchmark that lets the compiler inline
the callee measures nothing.** The assumed-size test had to be redone through a
procedure pointer before it meant anything — and then gave the same answer, which
is why it is recorded as settled.

---

## 5. Flags

`-fno-tree-vectorize` and LTO are worth 12% and 15% respectively, and the
reasoning for both is in `fortran/CMakeLists.txt` beside the flags. Do not
treat `-fno-tree-vectorize` as a blanket good: it is right for the gout
kernels, which index `g` by three values out of `idx` and get AVX gathers that
lose to scalar loads on this Broadwell. It is wrong for anything shaped like an
axpy. `cint_g1e.f90` already overrides it.

Intel builds **require** `-fp-model=precise`. Without it `icx`/`ifx` disagree
with the C in 75,593 of `int2e_check`'s 209,254 values, worst 4.2e-10. The
guard is in `fortran/CMakeLists.txt`.

---

## 5b. Profile-guided optimisation

Worth knowing about, and worth not overselling. Two-phase build, and the
profile data path embeds the build directory, so generate and use in the
*same* one or almost every file silently misses its data:

```
cmake -S . -B b -DCMAKE_Fortran_FLAGS="-fprofile-generate=/tmp/pgo"
cmake --build b && ./bench port 2
cmake -S . -B b -DCMAKE_Fortran_FLAGS="-fprofile-use=/tmp/pgo -fprofile-correction"
cmake --build b
```

| build | time |
|---|---:|
| port, `-O3` tuned | 2.194 s |
| port, same + PGO | **2.144 s** |
| C, `-O3` | 1.945 s |
| C, same + PGO | **1.841 s** |

It changes no results — the catalogue is bit-identical — and it is worth
2.3% to the port. But it is worth 5.3% to the C, so like-for-like it makes
the ratio *worse*, 1.165× against 1.128×. PGO is an orthogonal lever both
sides can pull, not a way to close the gap. Against a stock libcint build,
a PGO'd port is 1.10×, which is the honest number for someone deploying it
next to a distribution package.

Not wired into CMakeLists: it needs a representative workload, and which
workload is the caller's business.

## 6. The gap is not the language

Worth settling early, because it is the first thing anyone will suspect. All
four of these builds are bit-identical against the C:

| toolchain | flags | C | Fortran | ratio |
|---|---|---:|---:|---|
| gcc / gfortran | `-O3` | 1.945 s | 3.199 s | 1.64× |
| icx / ifx | `-O3 -fp-model=precise` | 1.987 s | 3.093 s | 1.56× |
| gcc / gfortran | tuned | 1.945 s | 2.480 s | 1.27× |
| icx / ifx | `-fp-model=precise -ipo` | 2.025 s | 2.754 s | 1.36× |

Two unrelated compiler families land within 5% of each other, and both respond
to the same lever (cross-module inlining: LTO, `-ipo`). Whatever is left is in
the Fortran *formulation*, not in gfortran.

---

## 7. Where every arity stands

`int2e_sph` was the one taken furthest, and §7b breaks it down. The rest were
swept afterwards with §3a and §3c, and are now within a few points of it:

c2h6 over the reference basis set, gfortran `-O3` tuned against gcc
`-O3`, every one bit-identical.  Wall clock is min-of-5 on an idle
40-core box (load average 0.16), alternating the two implementations so
any drift hits both equally.

| entry point | wall | instructions |
|---|---:|---:|
| `int1e_ipnuc_sph` | **1.038x** | 1.035x |
| `int2e_sph` | **1.109x** | 1.103x |
| `int3c2e_sph` | **1.123x** | 1.114x |
| `int1e_ipovlp_sph` | **1.126x** | 1.122x |
| `int2c2e_sph` | **1.159x** | 1.116x |
| `int2e_ip1_sph` | **1.171x** | 1.181x |
| `int3c1e_sph` | **1.175x** | 1.126x |

`int2e_ip1_sph` is the one entry that is not a Fortran problem — see the
vectoriser finding in §7b.

Five of the seven agree within a point, which is the useful methodological
result: on those paths instruction count is a valid stand-in for time, and
can be taken on a machine somebody else is using.

`int2c2e` and `int3c1e` do not agree — 4 to 5 points slower than their
instruction counts predict, and it reproduces exactly. Three explanations
are ruled out. **Cache**: total misses are in the thousands against 100M+
references, and the port has *fewer* than the C on both. **Branches**: the
port has fewer conditional mispredicts and far fewer indirect ones —
17,857 against 121,123 on `int2c2e`. **Expensive instruction mix**:
divisions and square roots are identical in count, and `CINTg3c1e_ovlp`
and `cint_g3c1e_ovlp` divide by `aijk` in the same three places. What is
left is instruction-level parallelism, and it is not decomposed. It is a
4-point residual on the two cheapest integrals in the set.
The two things left undone are named where they belong: the small-`l`
specialisation sweep (§2) has still only been done on the four-centre path,
and the leads in §7b. Note that `int2e_sph` is the only one whose *gout* is
hand-written; every other number above runs a generated kernel, so §3a's
generator wiring is what moves them.

## 7b. What is left on `int2e_sph`

Instructions: 3,113,640,910 against the C's 2,770,553,559. Self costs,
grouped so that differing inlining decisions do not distort the comparison
— the port inlines `prim_to_ctr` and the 2D→4D transfers where the C does
not, and callgrind charges those to the caller.

| group | delta | note |
|---|---:|---|
| `cint_2e_loop_nopt` | +83M | largest remaining; its screening is already at parity |
| `cint_gout2e` | +47M | loop bodies were **instruction-identical**, all per-call overhead over 2.34M calls — **answered by §3d**, it was the argument descriptors |
| driver + c2s | +40M | after threading the offset |
| envs declaration | +18M | measured, deliberately not taken — see §3b |
| `cint_g0_2e_2d` | +18M | after binding the three g rows |
| Rys roots | +11M | |
| `prim_to_ctr` | +24M | the C aliases `gc` and `gp` deliberately when `i_ctr == 1`, so the dummies cannot be split |
| `cint_g2e_index_xyz` | −1M | **done** — at parity with the C plus `CINTcart_comp` |
| `cint_g0_2e` + unrolled + transfers | −35M | **done** — the port is ahead |

**`int2e_ip1_sph` is the outlier, and the answer is gcc's vectoriser.**
Its generated gout plus the nabla it inlines cost 1,834M instructions to
the C's 1,487M over an identical 1,219,372 calls — 285 per call, on
source that is line-for-line the same shape. Line profiling could not
explain it. A disassembly diff could, immediately:

```
C     unpcklpd %xmm2,%xmm3     gfortran   mulsd  %xmm3,%xmm4
      mulpd    %xmm2,%xmm3                mulsd  %xmm5,%xmm3
      addpd    %xmm2,%xmm3                addsd  %xmm3,%xmm0
```

gcc vectorises the unrolled `nroots` arms **across the three s
components** — two of the three multiply-accumulates in one packed
instruction, from loads it already has. gfortran emits scalar code for
the identical Fortran.

This is not fixable with the flag, and the reason is worth recording
because the flag's own comment had it wrong. `-fno-tree-vectorize` was
justified in this file by bit-identity: vectorising made 371,752 values
differ. That was misattributed. The culprit was FMA contraction, not
vectorisation, and `-ffp-contract=off` separates them — with
`-ftree-vectorize -ffp-contract=off` the whole suite is bit-identical.
But measured that way, gfortran's vectoriser is a large regression
everywhere, including on the very kernel gcc vectorises well:

| path | scalar | vectorised |
|---|---:|---:|
| `int2e_sph` | 3,055M | 3,964M |
| `int3c2e_sph` | 245M | 291M |
| `int2e_ip1_sph` gout alone | 1,386M | 1,419M |

So the flag stays, now for the right reason: not that vectorising breaks
the contract — it does not — but that gfortran 13 vectorises this shape
badly. This is the one place found so far where §6 does not hold and the
gap really is the toolchain rather than the formulation.

**Ruled out along the way**, so nobody repeats them:

- *The inlined nabla is causing register pressure.* Plausible — the port
  inlines `cint_nabla1i_2e` into the gout where the C keeps it separate,
  and the port spills at 112 distinct sites against the C's 18, 244M
  instructions touching the stack against 59M. But `!GCC$ ATTRIBUTES
  noinline` (which gfortran does accept on a module procedure, unlike
  `always_inline`) made it **worse**: 3,750M → 3,791M, and the gout body
  did not improve at all. The spills are not from the nabla.

One other named lead, still not chased:

- `cint_g0_2e` on that same path is 479M to the C's 302M, where on
  `int2e_sph` the port is *ahead* of the C for the same function. The
  difference is angular momentum: the derivative raises `l` and a
  different branch runs. That is a §2 sweep, on the branch the four-centre
  sweep did not reach.
- The port dispatches the 2D→4D transfer with `select case`; the C uses a
  function pointer, and callgrind charges ~5 instructions a call to the
  comparison chain. The one place found so far where the C's indirection is
  cheaper.

## 8. Recipe

1. Pick one entry point. Build a two-mode benchmark that runs the port or the C
   from the same driver, so the comparison shares the harness.
2. Count instructions with callgrind. Confirm the ratio matches wall clock; if
   it does not, the problem is memory or scheduling and none of this applies.
3. Rank the functions by **delta**, not by share. A function at 2× that is 1%
   of runtime is not worth what a function at 1.2× that is 20% is.
4. For the top item, isolate it in a head-to-head microbenchmark and **sweep
   `l`**. Level at high `l` means a missing specialisation (§2).
5. Otherwise annotate by line and compare against the C's same lines. Look for
   pointer reaches the C hoisted (§3).
6. Verify with `catalogue_check` after every change, and read the timing only
   after it passes. Bit-identity is the acceptance criterion and the fastest
   way to know a "harmless" rewrite was not — `contiguous` on the gout dummies
   looked like a 0.5% win and was corrupting 366,672 values at the same time.
7. Record what did not work, with the number. Half of this document is that.
8. If the fix is in a gout kernel, put it in `scripts/f90-backend.cl`, not in
   the 197 generated files. §3a is the worked example.
9. Re-measure per arity before claiming the fix generalised. §3c is the worked
   counter-example: the same hoist was worth 15% on one loop and 0.005% on
   another.
