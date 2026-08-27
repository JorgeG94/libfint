#!/usr/bin/env python3
"""Re-emit the extended-precision tables as three-double quad rebuilds.

Runs on the original real128 output of c_tables_to_fortran.py (recover it
with `git show <rev>:src/cint_tab_jacobi_ext.f90` from before the dd
conversion), because that stage has already carried the C literal's exact
67 digits across. This stage reads those digits and splits each into THREE
doubles: a double holds 53 bits and binary128 needs 113, so two doubles are
seven bits short -- the very shortfall that made double-double insufficient
at n=12, x=200 -- and three are exact with room over.

**Each part is the remainder against the exact binary value of the previous
parts, never against their printed form.** `repr(x)` is the shortest decimal
that round-trips, which is a different number; splitting against it leaves
the low word wrong by ~2e-17 and nothing notices until the n=12, x=200 tail
of rys_check. Decimal(float) is the exact value; that is what is used.

Every triple is verified here, in the generator, before it is emitted, and
the check is the real criterion rather than a tolerance: the exact rational
hi+mid+lo must round to the SAME binary128 (round-to-nearest-even, 113-bit
significand, computed with Fraction arithmetic) as the source decimal does.
If any triple fails, nothing is written.

A quad cannot be a `parameter` -- rebuilding one is a call into the C shim --
so unlike the dd stage this cannot emit initialised arrays. It emits module
variables and an idempotent, OpenMP-guarded init that fills them element by
element through quad_from3.
"""
import re
import sys
from decimal import Decimal as D, getcontext
from fractions import Fraction as F

getcontext().prec = 90

CHUNK = 416          # assignments per fill subroutine, so none is huge
MAXCOL = 130         # wrap an assignment that would pass gfortran's 132

DECL = re.compile(
    r"real\(qp\), parameter :: (\w+)\(([^)]*)\) = \[(.*?)\]", re.S)

EXPECTED = {
    "lJACOBI_ALPHA": 48, "lJACOBI_BETA": 48,
    "lJACOBI_RN_PART2": 88, "lJACOBI_SN": 88, "lJACOBI_COEF": 1176,
    "qJACOBI_ALPHA": 64, "qJACOBI_BETA": 64,
    "qJACOBI_RN_PART2": 64, "qJACOBI_SN": 64, "qJACOBI_COEF": 2080,
}


def split3(dec):
    hi = float(dec)
    mid = float(dec - D(hi))            # exact binary value of hi, not repr
    lo = float(dec - D(hi) - D(mid))
    return hi, mid, lo


def round_binary128(x):
    """Round an exact rational to binary128, round-to-nearest-even.

    Returns (sign, integer significand in [2^112, 2^113), exponent); the
    tables are far from the subnormal range, so no clamping is needed.
    """
    if x == 0:
        return (0, 0, 0)
    s = 1 if x > 0 else -1
    a = abs(x)
    e = a.numerator.bit_length() - a.denominator.bit_length()
    if F(2) ** e > a:
        e -= 1
    while F(2) ** (e + 1) <= a:
        e += 1
    m = a / F(2) ** (e - 112)
    n, r = divmod(m.numerator, m.denominator)
    if 2 * r > m.denominator or (2 * r == m.denominator and n % 2 == 1):
        n += 1
    if n == 2 ** 113:
        n //= 2
        e += 1
    return (s, n, e)


def verified_triple(text):
    d = D(text)
    hi, mid, lo = split3(d)
    want = round_binary128(F(d))
    got = round_binary128(F(hi) + F(mid) + F(lo))
    if want != got:
        sys.exit(f"REFUSING to emit: triple for {text} rounds to a "
                 f"different binary128 than the source decimal")
    return hi, mid, lo


def parse(src):
    """name -> (bounds, [literal strings]), chunk refs expanded in order."""
    arrays = {}
    for name, bounds, body in DECL.findall(src):
        items = [t.strip() for t in body.replace("&", "").split(",")]
        items = [t for t in items if t]
        lits = []
        for it in items:
            if it.endswith("_qp"):
                lits.append(it[:-3])
            else:
                if it not in arrays:
                    sys.exit(f"chunk {it} referenced before definition")
                lits.extend(arrays[it][1])
        arrays[name] = (bounds, lits)
    return {n: v for n, v in arrays.items() if "*" not in v[0]}


def fmt_dp(x):
    return f"{x!r}_dp"


def emit_assign(name, idx, triple):
    args = ", ".join(fmt_dp(v) for v in triple)
    line = f"      {name}({idx}) = quad_from3({args})"
    if len(line) <= MAXCOL:
        return [line]
    a1, a2, a3 = (fmt_dp(v) for v in triple)
    return [f"      {name}({idx}) = quad_from3({a1}, &",
            f"         {a2}, {a3})"]


def main():
    src = open(sys.argv[1]).read()
    arrays = parse(src)
    if {n: len(v[1]) for n, v in arrays.items()} != EXPECTED:
        sys.exit(f"table shapes changed: "
                 f"{ {n: len(v[1]) for n, v in arrays.items()} }")

    fills = []          # (subroutine name, [statement lines])
    decls = []
    total = 0
    for name, (bounds, lits) in arrays.items():
        decls.append(f"   type(quad), protected :: {name}({bounds})")
        lo_bound = int(bounds.split(":")[0])
        for c0 in range(0, len(lits), CHUNK):
            sub = f"fill_{name.lower()}_{c0 // CHUNK + 1}"
            body = []
            for k, lit in enumerate(lits[c0:c0 + CHUNK]):
                body.extend(emit_assign(name, lo_bound + c0 + k,
                                        verified_triple(lit)))
                total += 1
            fills.append((sub, body))

    out = []
    w = out.append
    w("!")
    w("! Generated by scripts/c_tables_to_fortran.py -- do not edit.")
    w("! Re-emitted as three-double quad rebuilds by"
      " scripts/qp_tables_to_quad.py.")
    w("!")
    w("! Extended-precision Jacobi tables, from src/rys_wheeler.c. The C"
      " keeps two sets, long double and __float128; both become the quad"
      " ladder here, which PORT_TO_FORTRAN.md 3.1 argues for. The long"
      " double set is therefore MORE accurate here than in the C, where it"
      " is 80-bit x87, and so is deliberately not bitwise-comparable.")
    w("! Element order and values are those of the C tables.")
    w("!")
    w("! These are module VARIABLES, not parameters, and that is forced:"
      " the working type is `quad` -- binary128 reached through the C shim,"
      " for compilers with no real128 -- and rebuilding a quad is a call"
      " into that shim, which no initialiser may contain. So the tables"
      " are filled once at run time, element by element, from three-double"
      " splits whose sum the generator PROVED rounds to the same binary128"
      " as the 67-digit source value.")
    w("!")
    w("module cint_tab_jacobi_ext")
    w("   use cint_const, only: dp")
    w("   use cint_quad,  only: quad, quad_from3")
    w("   implicit none")
    w("   public")
    w("")
    out.extend(decls)
    w("")
    w("   logical, private :: tables_ready = .false.")
    priv = "   private :: " + ", ".join(s for s, _ in fills)
    while len(priv) > MAXCOL:
        cut = priv.rfind(",", 0, MAXCOL - 3)
        w(priv[:cut + 1] + " &")
        priv = "      " + priv[cut + 1:].lstrip()
    w(priv)
    w("")
    w("contains")
    w("")
    w("   subroutine cint_tab_jacobi_ext_init()")
    w("      !! Idempotent, and safe to call from inside a parallel region:")
    w("      !! the steady-state cost is one flag read, the fill happens")
    w("      !! once inside a named critical section, and the flag flips")
    w("      !! only after a flush, so a thread that sees .true. on the")
    w("      !! fast path also sees the filled tables.")
    w("      if (tables_ready) return")
    w("      !$omp critical (fint_tab_jacobi_ext_fill)")
    w("      if (.not. tables_ready) then")
    for sub, _ in fills:
        w(f"         call {sub}()")
    w("         !$omp flush")
    w("         tables_ready = .true.")
    w("      end if")
    w("      !$omp end critical (fint_tab_jacobi_ext_fill)")
    w("   end subroutine cint_tab_jacobi_ext_init")
    for sub, body in fills:
        w("")
        w(f"   subroutine {sub}()")
        out.extend(body)
        w(f"   end subroutine {sub}")
    w("")
    w("end module cint_tab_jacobi_ext")
    w("")
    sys.stdout.write("\n".join(out))
    sys.stderr.write(f"{total} constants emitted, every triple verified to "
                     f"round to the source's binary128\n")


main()
