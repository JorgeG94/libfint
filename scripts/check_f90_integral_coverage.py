#!/usr/bin/env python3
"""Fail if the C generator knows an integral the Fortran generator does not.

    python3 scripts/check_f90_integral_coverage.py

WHY THIS EXISTS
---------------
`scripts/auto_intor.cl` is the C's list of integrals.  `scripts/gen_f90_all.cl`
is the Fortran's, and it does not read the first one -- it carries its own copy
of every description.  So when upstream adds an integral, the Fortran generator
keeps producing exactly what it produced yesterday, and nothing anywhere says
so.  No error, no warning, no failing test: `catalogue_check` only checks what
was generated, so it passes happily while the port quietly stops being a
complete port.

That is what happened with #129, which added four second-derivative `int3c1e`
families.  They were noticed by hand, on the way to answering an unrelated
question.  Next time this should be a build failure instead.

The check is deliberately one-directional.  Names in `gen_f90_all.cl` that are
absent from `auto_intor.cl` are fine and expected: the "based on" families --
`cint1e_a.c` and `cint3c1e_a.c` -- are hand-written in the C with their
descriptions in comments, and the Fortran generates them from descriptions the
C never had in machine-readable form.

Matching is on the integral name only.  Comparing the descriptions themselves
would be better still, but the two files quote them differently and a
false alarm every release is worse than this.
"""

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# '("int3c1e_ip1"  ( nabla \, \,))  ->  int3c1e_ip1
NAME = re.compile(r"""^\s*'\("([A-Za-z_][A-Za-z0-9_]*)"\s""", re.M)


def names(path):
    with open(path) as fh:
        text = fh.read()
    # drop ;; comments so a commented-out description does not count as coverage
    text = re.sub(r";.*$", "", text, flags=re.M)
    return set(NAME.findall(text))


def main():
    # upstream's list, from the pinned submodule -- not a copy
    c_list = os.path.join(ROOT, "extern", "libcint", "scripts",
                          "auto_intor.cl")
    if not os.path.exists(c_list):
        sys.exit("extern/libcint is not checked out: "
                 "git submodule update --init")
    f_list = os.path.join(ROOT, "scripts", "gen_f90_all.cl")
    for p in (c_list, f_list):
        if not os.path.exists(p):
            sys.exit("missing %s" % p)

    in_c = names(c_list)
    in_f = names(f_list)
    if not in_c or not in_f:
        sys.exit("parsed no integral names -- the description syntax moved, "
                 "and this check is now silently vacuous")

    missing = sorted(in_c - in_f)
    print("auto_intor.cl: %d integrals   gen_f90_all.cl: %d   "
          "Fortran-only (expected): %d" % (len(in_c), len(in_f), len(in_f - in_c)))

    if missing:
        print()
        print("The C generator knows %d integral(s) the Fortran one does not:"
              % len(missing))
        for n in missing:
            print("    %s" % n)
        print()
        print("Copy each description from scripts/auto_intor.cl into the matching")
        print("gen-f90-cint block in scripts/gen_f90_all.cl, regenerate, then")
        print("re-run fortran/test/gen_catalogue_check.py so the new entry points")
        print("are actually verified against the C rather than merely built.")
        return 1

    print("OK: every integral in auto_intor.cl has a Fortran description")
    return 0


if __name__ == "__main__":
    sys.exit(main())
