#!/usr/bin/env python3
"""Generate src/rotaxis/*.f90 (one module per class), src/cint_rotaxis_kernels.f90
(the dispatcher) and src/rotaxis/kernels.cmake.

    python3 scripts/rotaxis_mmd/generate.py [--lmax 2] [--no-L] [--root <checkout>]

Classes are the canonical ones over shell kinds s < p < L < d: the kinds
ascend within each pair, and the pairs ascend.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rotaxis_mmd.emit import emit_files
from rotaxis_mmd.derive import KIND_RANK, KIND_LMAX


def canonical_classes(lmax, with_L=True):
    kinds = [k for k in sorted(KIND_RANK, key=KIND_RANK.get)
             if KIND_LMAX[k] <= lmax and (with_L or k != "L")]
    pairs = [(a, b) for i, a in enumerate(kinds) for b in kinds[i:]]
    out = []
    for i, p in enumerate(pairs):
        for q in pairs[i:]:
            out.append(p + q)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=2)
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), help="libfint checkout to write into")
    ap.add_argument("--summary", action="store_true", help="print the per-class sizes only")
    ap.add_argument("--no-L", action="store_true", help="leave out the L-shell kernels")
    ap.add_argument("--unroll-limit", type=int, default=3000,
                    help="classes with more ket-level terms than this are table driven")
    args = ap.parse_args()
    classes = canonical_classes(args.lmax, not args.no_L)
    files, facts = emit_files(classes, args.unroll_limit)
    for f in facts:
        print(f.summary(), file=sys.stderr)
    if args.summary:
        return
    import glob
    for stale in glob.glob(os.path.join(args.root, "src", "rotaxis", "rotaxis_*.f90")):
        os.remove(stale)
    for rel, text in files.items():
        path = os.path.join(args.root, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)
        print(f"wrote {rel}: {text.count(chr(10))} lines", file=sys.stderr)


if __name__ == "__main__":
    main()
