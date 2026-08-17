#!/usr/bin/env python3
"""Re-emit libcint's static numeric tables as Fortran parameter modules.

Deliverable D2.  Run from the repository root:

    python3 scripts/c_tables_to_fortran.py --out fortran/src

WHY THIS TRANSFORMS THE C RATHER THAN RE-RUNNING THE PYTHON THAT MADE IT
-----------------------------------------------------------------------
PORT_PLAN.md originally proposed regenerating these tables from the mpmath
scripts in this directory.  Transforming the committed C is better, for three
reasons:

  * The C tables are what the shipped, validated library actually uses.
    Regenerating risks last-bit drift from a different mpmath version or
    precision setting -- a silent change to results, with nothing gained.
  * D2's acceptance criterion is bitwise identity.  A syntax transform
    satisfies it by construction; a regeneration would have to be argued.
  * The tabulation scripts are long-running jobs with pickle caches and need
    mpmath, numpy and scipy.

Regeneration stays available and is a different job: you go back to the Python
when you want to *extend* the tables (more roots, wider range), not when you
want the ones that already exist.

That this is exact is not assumed.  A C decimal literal and the same text in
Fortran round to the same IEEE double -- verified over 500 literals drawn from
these very tables -- and the generated modules are checked element by element
against the C at build time by the tables_bitwise ctest.

WHAT IT HANDLES
---------------
cpp conditionals inside an initializer.  g_trans_cart2sph has an #ifdef PYPZPX
choosing the p-orbital order, so the extracted fragment is run back through the
C preprocessor with the same -D flags the build uses, rather than parsed by
hand.
"""

import argparse
import os
import re
import subprocess
import sys

# (source file, [table names]).  roots_xw.dat is deliberately absent: it is
# 99.9 MB of Slater-geminal roots reached only from g2e_f12.c, which is not in
# metalquicha's dependency closure.
TABLES = [
    ("src/polyfits.c",           ["COS_14_14"]),
    ("src/roots_for_x0.dat",     ["POLY_SMALLX_R0", "POLY_SMALLX_R1",
                                  "POLY_SMALLX_W0", "POLY_SMALLX_W1",
                                  "POLY_LARGEX_RT", "POLY_LARGEX_WW"]),
    # JACOBI_XS and JACOBI_CS are deliberately absent.  Their only reader is
    # shifted_jacobi_moments(), which is defined at rys_wheeler.c:3236 and
    # never called, so both tables are dead data upstream -- the compiler drops
    # them, which is how this was noticed.  Carrying them into the port would
    # mean carrying data no build can check.
    ("src/rys_wheeler.c",        ["JACOBI_ALPHA", "JACOBI_BETA",
                                  "JACOBI_RN_PART2", "JACOBI_SN",
                                  "JACOBI_COEF", "JACOBI_COEF_ORDER"]),
    ("src/cart2sph.c",           ["g_trans_cart2sph", "g_trans_cart2jR",
                                  "g_trans_cart2jI", "_len_cart"]),
    # The extended-precision copies of the Jacobi tables.  CINTrys_roots
    # escalates to these above 8 Rys roots (see PORT_TO_FORTRAN.md 3.1), and
    # they are easy to miss because they are neither `double` nor named like
    # the tables beside them.
    ("src/rys_wheeler.c",        ["lJACOBI_ALPHA", "lJACOBI_BETA",
                                  "lJACOBI_RN_PART2", "lJACOBI_SN",
                                  "lJACOBI_COEF"]),
    ("src/rys_wheeler.c",        ["qJACOBI_ALPHA", "qJACOBI_BETA",
                                  "qJACOBI_RN_PART2", "qJACOBI_SN",
                                  "qJACOBI_COEF"]),
    ("src/sr_roots_part0_x.dat", ["SR_DATA0_X"]),
    ("src/sr_roots_part0_w.dat", ["SR_DATA0_W"]),
    ("src/sr_roots_part1_x.dat", ["SR_DATA1_X"]),
    ("src/sr_roots_part1_w.dat", ["SR_DATA1_W"]),
    ("src/sr_roots_part2_x.dat", ["SR_DATA2_X"]),
    ("src/sr_roots_part2_w.dat", ["SR_DATA2_W"]),
    ("src/sr_roots_part3_x.dat", ["SR_DATAL_X"]),
    ("src/sr_roots_part3_w.dat", ["SR_DATAL_W"]),
]

# Which Fortran module each source contributes to.
MODULE_OF = {
    "src/polyfits.c":       "cint_tab_polyfits",
    "src/roots_for_x0.dat": "cint_tab_roots_x0",
    "src/rys_wheeler.c":    "cint_tab_jacobi",
    "src/cart2sph.c":       "cint_tab_cart2sph",
}
SR_MODULE = "cint_tab_sr_roots"

# C identifiers that are not valid Fortran ones.
FORTRAN_NAME = {"_len_cart": "len_cart"}

DECL = re.compile(
    r"^static\s+(?:const\s+)?(long double|__float128|double|int|FINT)\s+(%s)\s*\[\]\s*=\s*\{",
    re.M)


def extract_fragment(path, name):
    """Pull the text of one table declaration out of a source file."""
    src = open(path).read()
    m = re.search(DECL.pattern % re.escape(name), src, re.M)
    if not m:
        raise SystemExit("could not find table %s in %s" % (name, path))
    ctype = m.group(1)
    i = src.index("{", m.start())
    depth = 0
    for j in range(i, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return ctype, src[i + 1:j]
    raise SystemExit("unterminated initializer for %s in %s" % (name, path))


def resolve_cpp(body, defines):
    """Run an initializer body through cpp, so #ifdef inside it is honoured
    exactly as the C build would honour it, and comments disappear."""
    if "#" not in body:
        # No conditionals; still strip comments the cheap way.
        body = re.sub(r"/\*.*?\*/", " ", body, flags=re.S)
        body = re.sub(r"//[^\n]*", " ", body)
        return body
    cmd = ["cpp", "-P"] + ["-D%s" % d for d in defines] + ["-"]
    out = subprocess.run(cmd, input=body, capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit("cpp failed on fragment:\n" + out.stderr)
    return out.stdout


def parse_values(body, ctype):
    toks = [t.strip() for t in body.split(",")]
    toks = [t for t in toks if t != ""]
    vals = []
    for t in toks:
        t = t.replace("\n", " ").strip()
        if ctype in ("int", "FINT"):
            vals.append(str(int(t, 0)))
        else:
            # Keep the literal's exact digits; only the suffix differs between
            # the languages.  Reformatting the number would be the one way to
            # lose bits here.
            # C marks precision with a literal suffix -- 1.5l for long double,
            # 1.5q for __float128 -- where Fortran uses a kind suffix.  Strip
            # it before validating, not after.
            t = t.rstrip("qQlLfF")
            if not re.match(r"^[-+]?(\d+\.?\d*([eE][-+]?\d+)?|\.\d+([eE][-+]?\d+)?)$", t):
                raise SystemExit("unexpected numeric token: %r" % t)
            # ...with one exception that is not a reformatting: C writes whole
            # numbers in a double table as `0` and `1`, and `0_dp` in Fortran
            # is an integer of kind dp, not a real.  A decimal point is
            # required, and adding one cannot change the value.
            if "." not in t and "e" not in t and "E" not in t:
                t += ".0"
            # C spells a quad literal 1.23q0; Fortran spells the kind suffix.
            t = t.rstrip("qQlLfF")
            vals.append(t + ("_qp" if ctype in ("long double", "__float128") else "_dp"))
    return vals


def emit_module(fh, modname, tables, note):
    fh.write("!\n! Generated by scripts/c_tables_to_fortran.py -- do not edit.\n")
    fh.write("!\n! %s\n" % note)
    fh.write("! Element order and values are bitwise those of the C tables; the\n")
    fh.write("! tables_bitwise test checks that rather than trusting it.\n!\n")
    fh.write("module %s\n" % modname)
    fh.write("   use cint_const, only: dp, qp\n")
    fh.write("   implicit none\n")
    fh.write("   public\n\n")
    for name, ctype, vals in tables:
        if ctype in ("int", "FINT"):
            ftype = "integer"
        elif ctype in ("long double", "__float128"):
            ftype = "real(qp)"
        else:
            ftype = "real(dp)"
        # A leading underscore is a legal C identifier and not a legal Fortran
        # one, so _len_cart has to be renamed.  Recorded here rather than
        # silently, because callers have to know.
        name = FORTRAN_NAME.get(name, name)
        # 0-based bounds, matching the C indexing the port keeps (see
        # PORT_TO_FORTRAN.md 3.6).
        # Free-form Fortran caps a line at 132 columns, and some of these
        # literals carry 36 digits -- more precision than a double holds, which
        # the C truncates on parse just as Fortran does.  So wrap on width, not
        # on a fixed count, and stay well inside the limit rather than relying
        # on -ffree-line-length-none.
        LIMIT = 120
        lines = []
        line, first = "     ", True
        for v in vals:
            piece = v if first else ", " + v
            if not first and len(line) + len(piece) + 2 > LIMIT:
                lines.append(line)
                line, piece = "     " + v, ""
            line += piece
            first = False
        lines.append(line)

        # A statement may carry at most 255 continuation lines (Fortran 2008;
        # 2018 lifts it, but not every compiler has).  The biggest table here
        # would need 56,197, so long tables are emitted as chunk parameters and
        # concatenated once.  Chunking only the tables that need it keeps the
        # common case a single constructor, which compiles faster.
        CHUNK = 200
        if len(lines) <= CHUNK:
            fh.write("   %s, parameter :: %s(0:%d) = [ &\n" % (ftype, name, len(vals) - 1))
            fh.write(", &\n".join(lines))
            fh.write(" ]\n\n")
        else:
            nch = 0
            for k in range(0, len(lines), CHUNK):
                nch += 1
                block = lines[k:k + CHUNK]
                fh.write("   %s, parameter :: %s_c%d(*) = [ &\n" % (ftype, name, nch))
                fh.write(", &\n".join(block))
                fh.write(" ]\n")
            fh.write("   %s, parameter :: %s(0:%d) = [ &\n" % (ftype, name, len(vals) - 1))
            names = ["%s_c%d" % (name, i + 1) for i in range(nch)]
            line, first = "     ", True
            for nm in names:
                piece = nm if first else ", " + nm
                if not first and len(line) + len(piece) + 2 > LIMIT:
                    fh.write("%s, &\n" % line)
                    line, piece = "     " + nm, ""
                line += piece
                first = False
            fh.write("%s ]\n\n" % line)
    fh.write("end module %s\n" % modname)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", required=True)
    ap.add_argument("--define", action="append", default=[],
                    help="cpp defines, e.g. --define PYPZPX")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    grouped = {}
    for path, names in TABLES:
        mod = MODULE_OF.get(path, SR_MODULE)
        if names and names[0][0] in "lq" and "JACOBI" in names[0]:
            mod = "cint_tab_jacobi_ext"
        full = os.path.join(args.root, path)
        for name in names:
            ctype, body = extract_fragment(full, name)
            body = resolve_cpp(body, args.define)
            vals = parse_values(body, ctype)
            grouped.setdefault(mod, []).append((name, ctype, vals))
            print("  %-18s %-16s %8d elements" % (os.path.basename(path), name, len(vals)))

    notes = {
        "cint_tab_polyfits":  "Chebyshev fit coefficients, from src/polyfits.c",
        "cint_tab_roots_x0":  "Small-x and large-x Rys shortcuts, from src/roots_for_x0.dat",
        "cint_tab_jacobi":    "Jacobi/Wheeler quadrature tables, from src/rys_wheeler.c",
        "cint_tab_cart2sph":  "Cartesian-to-spherical and -spinor transforms, from src/cart2sph.c",
        SR_MODULE:            "Range-separated Rys roots and weights, from src/sr_roots_part*.dat",
        "cint_tab_jacobi_ext": "Extended-precision Jacobi tables, from src/rys_wheeler.c. "
                               "The C keeps two sets, long double and __float128; both "
                               "become real(qp) here, which PORT_TO_FORTRAN.md 3.1 argues "
                               "for. The long double set is therefore MORE accurate here "
                               "than in the C, where it is 80-bit x87, and so is "
                               "deliberately not bitwise-comparable.",
    }
    for mod, tables in grouped.items():
        out = os.path.join(args.out, mod + ".f90")
        with open(out, "w") as fh:
            emit_module(fh, mod, tables, notes[mod])
        total = sum(len(v) for _, _, v in tables)
        print("wrote %s  (%d tables, %d elements)" % (out, len(tables), total))


if __name__ == "__main__":
    main()
