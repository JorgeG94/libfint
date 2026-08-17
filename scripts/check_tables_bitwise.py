#!/usr/bin/env python3
"""Check the generated Fortran tables against the C, element by element.

The C side is read out of the compiled object files -- the bytes gcc actually
emitted -- rather than by re-parsing the source.  Parsing the source twice
would only prove the parser agrees with itself.

    python3 scripts/check_tables_bitwise.py --objdir <build>/CMakeFiles/cint.dir/src \\
                                            --fortran-dump <dumpfile>

The Fortran side comes from fortran/test/dump_tables.f90, which writes every
table as an unformatted stream in the same order as TABLES below.
"""

import argparse
import os
import struct
import subprocess
import sys

# table name -> object file that defines it
TABLES = [
    ("COS_14_14",        "polyfits.c.o",        "d"),
    ("POLY_SMALLX_R0",   "rys_roots.c.o",       "d"),
    ("POLY_SMALLX_R1",   "rys_roots.c.o",       "d"),
    ("POLY_SMALLX_W0",   "rys_roots.c.o",       "d"),
    ("POLY_SMALLX_W1",   "rys_roots.c.o",       "d"),
    ("POLY_LARGEX_RT",   "rys_roots.c.o",       "d"),
    ("POLY_LARGEX_WW",   "rys_roots.c.o",       "d"),
    ("JACOBI_ALPHA",     "rys_wheeler.c.o",     "d"),
    ("JACOBI_BETA",      "rys_wheeler.c.o",     "d"),
    ("JACOBI_RN_PART2",  "rys_wheeler.c.o",     "d"),
    ("JACOBI_SN",        "rys_wheeler.c.o",     "d"),
    ("JACOBI_COEF",      "rys_wheeler.c.o",     "d"),
    ("JACOBI_COEF_ORDER","rys_wheeler.c.o",     "i"),
    ("g_trans_cart2sph", "cart2sph.c.o",        "d"),
    ("g_trans_cart2jR",  "cart2sph.c.o",        "d"),
    ("g_trans_cart2jI",  "cart2sph.c.o",        "d"),
    ("_len_cart",        "cart2sph.c.o",        "i"),
    ("SR_DATA0_X",       "sr_rys_polyfits.c.o", "d"),
    ("SR_DATA0_W",       "sr_rys_polyfits.c.o", "d"),
    ("SR_DATA1_X",       "sr_rys_polyfits.c.o", "d"),
    ("SR_DATA1_W",       "sr_rys_polyfits.c.o", "d"),
    ("SR_DATA2_X",       "sr_rys_polyfits.c.o", "d"),
    ("SR_DATA2_W",       "sr_rys_polyfits.c.o", "d"),
    ("SR_DATAL_X",       "sr_rys_polyfits.c.o", "d"),
    ("SR_DATAL_W",       "sr_rys_polyfits.c.o", "d"),
    # Extended precision.  __float128 is IEEE binary128 on both sides, so it
    # compares bitwise like the rest.  `long double` does not: the C stores
    # 80-bit x87 padded to 16 bytes, while Fortran's real128 is true
    # binary128, so the port is deliberately MORE accurate there and the two
    # cannot be bit-identical.  Those are checked by value instead.
    ("lJACOBI_ALPHA",    "rys_wheeler.c.o",     "l"),
    ("lJACOBI_BETA",     "rys_wheeler.c.o",     "l"),
    ("lJACOBI_RN_PART2", "rys_wheeler.c.o",     "l"),
    ("lJACOBI_SN",       "rys_wheeler.c.o",     "l"),
    ("lJACOBI_COEF",     "rys_wheeler.c.o",     "l"),
    ("qJACOBI_ALPHA",    "rys_wheeler.c.o",     "q"),
    ("qJACOBI_BETA",     "rys_wheeler.c.o",     "q"),
    ("qJACOBI_RN_PART2", "rys_wheeler.c.o",     "q"),
    ("qJACOBI_SN",       "rys_wheeler.c.o",     "q"),
    ("qJACOBI_COEF",     "rys_wheeler.c.o",     "q"),
]

ELEM_BYTES = {"d": 8, "i": 4, "q": 16, "l": 16}


def x87_to_float(b):
    """Decode an 80-bit x87 long double (padded to 16 bytes) to a Python float."""
    mant = int.from_bytes(b[0:8], "little")
    se = int.from_bytes(b[8:10], "little")
    sign = -1.0 if (se >> 15) else 1.0
    exp = se & 0x7FFF
    if exp == 0 and mant == 0:
        return 0.0 * sign
    return sign * mant * 2.0 ** (exp - 16383 - 63)


def binary128_to_float(b):
    """Decode an IEEE binary128 to a Python float (lossy to binary64)."""
    v = int.from_bytes(b, "little")
    sign = -1.0 if (v >> 127) else 1.0
    exp = (v >> 112) & 0x7FFF
    frac = v & ((1 << 112) - 1)
    if exp == 0 and frac == 0:
        return 0.0 * sign
    try:
        return sign * (1.0 + frac / float(1 << 112)) * 2.0 ** (exp - 16383)
    except OverflowError:
        return float("inf") * sign

SECTION_OF = {"d": ".data", "r": ".rodata", "b": ".bss", "D": ".data", "R": ".rodata"}


def symbol_info(obj, name):
    out = subprocess.run(["nm", "-S", obj], capture_output=True, text=True).stdout
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[3] == name:
            addr, size, typ = int(parts[0], 16), int(parts[1], 16), parts[2]
            return addr, size, typ
    raise SystemExit("symbol %s not found in %s" % (name, obj))


def section_bytes(obj, section):
    tmp = "/tmp/_sect.bin"
    r = subprocess.run(["objcopy", "-O", "binary", "--only-section=" + section, obj, tmp],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("objcopy failed on %s %s: %s" % (obj, section, r.stderr))
    data = open(tmp, "rb").read()
    os.unlink(tmp)
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--objdir", required=True)
    ap.add_argument("--fortran-dump", required=True)
    ap.add_argument("--no-polyfit", action="store_true",
                    help="WITH_POLYNOMIAL_FIT is off (the default), so neither "
                         "polyfits.c nor sr_rys_polyfits.c was compiled and "
                         "their tables are absent on both sides")
    args = ap.parse_args()

    fdump = open(args.fortran_dump, "rb").read()
    pos = 0
    failures = 0
    total_elems = 0

    tables = TABLES
    if args.no_polyfit:
        tables = [t for t in TABLES
                  if not t[0].startswith("SR_DATA") and t[0] != "COS_14_14"]

    for name, objfile, kind in tables:
        obj = os.path.join(args.objdir, objfile)
        if not os.path.exists(obj):
            raise SystemExit("missing object %s -- build libcint first" % obj)

        addr, size, typ = symbol_info(obj, name)
        section = SECTION_OF.get(typ)
        if section is None:
            raise SystemExit("%s: unhandled symbol type %r" % (name, typ))
        sect = section_bytes(obj, section)
        cbytes = sect[addr:addr + size]
        if len(cbytes) != size:
            raise SystemExit("%s: short read from %s" % (name, section))

        w = ELEM_BYTES[kind]
        n = size // w
        # The Fortran side stores long double as binary128, so its record for
        # an `l` table is the same width here (16 bytes) but a different
        # encoding.
        fbytes = fdump[pos:pos + n * w]
        pos += n * w
        total_elems += n

        if kind == "l":
            worst = 0.0
            bad = 0
            for e in range(n):
                cv = x87_to_float(cbytes[e * w:(e + 1) * w])
                fv = binary128_to_float(fbytes[e * w:(e + 1) * w])
                d = abs(cv - fv)
                rel = d / abs(cv) if cv != 0 else d
                worst = max(worst, rel)
                if rel > 1e-15:
                    bad += 1
                    if bad <= 3:
                        print("  %-20s VALUE MISMATCH element %d: C=%.17g F=%.17g" % (name, e, cv, fv))
            if bad:
                failures += 1
            else:
                print("  %-20s %8d elements  OK (promoted to binary128, worst rel %.1e)"
                      % (name, n, worst))
            continue

        if cbytes[:n * w] == fbytes:
            print("  %-20s %8d elements  OK" % (name, n))
            continue

        failures += 1
        # Report the first differing element concretely rather than "differs".
        fmt = {"d": "<d", "i": "<i", "q": None}[kind]
        shown = 0
        for e in range(n):
            cv = cbytes[e * w:(e + 1) * w]
            fv = fbytes[e * w:(e + 1) * w]
            if cv != fv:
                if fmt is None:
                    print("  %-20s MISMATCH at element %d: C=%.20g F=%.20g" % (
                        name, e, binary128_to_float(cv), binary128_to_float(fv)))
                else:
                    print("  %-20s MISMATCH at element %d: C=%r F=%r" % (
                        name, e, struct.unpack(fmt, cv)[0],
                        struct.unpack(fmt, fv)[0] if len(fv) == w else None))
                shown += 1
                if shown >= 5:
                    print("       ...")
                    break

    if pos != len(fdump):
        print("SIZE MISMATCH: Fortran dump has %d bytes, consumed %d" % (len(fdump), pos))
        failures += 1

    print()
    print("  %d tables, %d elements" % (len(tables), total_elems))
    if failures:
        print("  RESULT: FAIL (%d tables differ)" % failures)
        return 1
    print("  RESULT: PASS -- bitwise identical to the C, except the long double")
    print("          tables, which are promoted to binary128 and so agree by")
    print("          value rather than by bit pattern.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
