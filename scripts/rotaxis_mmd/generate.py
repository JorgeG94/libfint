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

from rotaxis_mmd.emit import emit_files, emit_grad_files
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
    ap.add_argument("--grad", action="store_true",
                    help="emit the gradient kernels (d/dA, libcint int2e_ip1 layout)")
    ap.add_argument("--max-terms", type=int, default=0,
                    help="skip classes above this many level-2 terms; 0 = keep all. "
                         "The supported() predicates ask the dispatcher, so a "
                         "partial set is declined cleanly rather than trusted.")
    ap.add_argument("--unroll-limit", type=int, default=3000,
                    help="classes with more ket-level terms than this are table driven")
    args = ap.parse_args()
    if args.grad:
        # NO canonicalisation for gradients.  The energy driver is free to
        # permute a quartet into a canonical class because the integral is
        # symmetric under those swaps; d/dA is not -- it names one centre,
        # so whichever shell the caller put in slot 0 has to stay there.
        # Every ordered combination of kinds therefore gets its own kernel.
        kinds = [k for k in sorted(KIND_RANK, key=KIND_RANK.get)
                 if KIND_LMAX[k] <= args.lmax and (not args.no_L or k != "L")]
        classes = [(a, b, c, d) for a in kinds for b in kinds
                   for c in kinds for d in kinds]
    else:
        classes = canonical_classes(args.lmax, not args.no_L)
    if args.grad:
        files, facts = emit_grad_files(classes, args.unroll_limit, args.max_terms)
    else:
        files, facts = emit_files(classes, args.unroll_limit)
    for f in facts:
        print(f.summary(), file=sys.stderr)
    if args.summary:
        return
    import glob
    sub = "rotaxis_grad" if args.grad else "rotaxis"
    pat = "rotaxis_grad_*.f90" if args.grad else "rotaxis_*.f90"
    for stale in glob.glob(os.path.join(args.root, "src", sub, pat)):
        os.remove(stale)
    for rel, text in files.items():
        path = os.path.join(args.root, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)
        print(f"wrote {rel}: {text.count(chr(10))} lines", file=sys.stderr)


if __name__ == "__main__":
    main()
