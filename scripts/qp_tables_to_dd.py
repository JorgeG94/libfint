#!/usr/bin/env python3
"""Re-emit the extended-precision tables as double-double pairs.

Runs on the output of c_tables_to_fortran.py rather than on the C, because that
stage has already done the part that can lose bits -- it carries the C literal's
exact digits across and changes only the kind suffix. This stage reads those
digits and splits each into two doubles.

**Each `lo` is the remainder against the exact binary value of `hi`, never
against its printed form.** `repr(hi)` is the shortest decimal that round-trips,
which is not the same number; using it leaves every `lo` wrong by ~2e-17. Every
pair emitted here is verified against the source decimal before it is written.

A dd pair carries ~106 bits against binary128's 113, so this is a real if small
narrowing of the tables -- taken because nvfortran and LLVM Flang have no
real128 at all and cannot compile the qp module in any form.
"""
import re, sys
from decimal import Decimal as D, getcontext
getcontext().prec = 60

LIT = re.compile(r"([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)_qp")

def split(dec):
    hi = float(dec)
    lo = float(dec - D(hi))          # exact binary value of hi, not repr(hi)
    return hi, lo

def convert(text):
    worst = D(0)
    count = 0
    def repl(m):
        nonlocal worst, count
        d = D(m.group(1))
        hi, lo = split(d)
        if d != 0:
            rel = abs(d - (D(hi) + D(lo))) / abs(d)
            worst = max(worst, rel)
        count += 1
        return f"dd({hi!r}_dp, {lo!r}_dp)"
    out = LIT.sub(repl, text)
    return out, count, worst

def main():
    src = open(sys.argv[1]).read()
    out, count, worst = convert(src)
    out = out.replace("use cint_const, only: dp, qp", "use cint_const, only: dp\n   use cint_dd, only: dd")
    out = out.replace("real(qp), parameter", "type(dd), parameter")
    out = out.replace("-- do not edit.",
                      "-- do not edit.\n! Re-emitted as double-double by scripts/qp_tables_to_dd.py.")
    sys.stderr.write(f"{count} literals converted, worst pair error {float(worst):.2e}\n")
    if worst > D("1e-31"):
        sys.exit("a pair failed verification")
    sys.stdout.write(out)

main()
