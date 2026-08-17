#!/usr/bin/env python3
"""Translate libcint's closed-form rys_root1..5 from C to Fortran.

    python3 scripts/c_rysroot_to_fortran.py --out fortran/src/cint_rys_root_n.f90

These five functions are 1,376 lines of dense polynomial fits -- Chebyshev-like
expansions of the Rys roots and weights for one to five roots, branching on
ranges of X.  They are the one part of the Rys subsystem with no algorithm in
it, only arithmetic, and hand-copying 1,376 lines of coefficients is exactly
the kind of transcription that goes wrong silently.

The C they are written in is a narrow subset -- if/else-if chains on X, local
double scalars, straight-line arithmetic, and writes to roots[]/weights[] --
so a mechanical translation is both possible and checkable.  The check is that
the result reproduces the C to the last bit across a sweep, which
rys_root_n_check does.

What is deliberately NOT attempted: anything outside that subset.  If the
parser meets a construct it does not recognise it stops rather than guessing.
"""

import argparse
import re

FUNCS = ["rys_root1", "rys_root2", "rys_root3", "rys_root4", "rys_root5"]


def extract(src, name):
    """Pull one function body out of the C source."""
    m = re.search(r"^static int %s\(double X, double \*roots, double \*weights\)\s*\n\{" % name,
                  src, re.M)
    if not m:
        raise SystemExit("cannot find %s" % name)
    i = src.index("{", m.start())
    depth = 0
    for j in range(i, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[i+1:j]
    raise SystemExit("unterminated body for %s" % name)


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    text = re.sub(r"//[^\n]*", " ", text)
    return text


NUM = re.compile(r"""
    (?<![A-Za-z0-9_.])            # not part of an identifier
    (
      \d+\.\d*(?:[EeDd][-+]?\d+)? # 1.23e4, 1.
    | \.\d+(?:[EeDd][-+]?\d+)?    # .5
    | \d+[EeDd][-+]?\d+           # 1e4
    )
    """, re.X)


def kindify(expr):
    """Give every real literal an explicit _dp kind.

    Integer literals are left alone: they appear only as array subscripts and
    small counts here, and promoting them would change subscripts into reals.
    """
    def repl(m):
        t = m.group(1).replace("D", "E").replace("d", "e")
        return t + "_dp"
    return NUM.sub(repl, expr)


def convert_expr(e):
    # The C wraps long expressions across lines; Fortran needs its own
    # continuations, so flatten first and let wrap() re-break them.
    e = re.sub(r"\s+", " ", e).strip()
    e = re.sub(r"\broots\[([^\]]+)\]", r"roots(\1)", e)
    e = re.sub(r"\bweights\[([^\]]+)\]", r"weights(\1)", e)
    # C's pow(a,b) is Fortran's a**b.  exp and sqrt carry over unchanged, and
    # they are the only other functions these routines call.
    e = re.sub(r"\bpow\(([^,]+),\s*([^)]+)\)", r"(\1)**(\2)", e)
    return kindify(e)


def translate(body, name):
    body = strip_comments(body)
    out = []
    decls = set()
    indent = 2
    # Split into statements and braces while keeping order.
    tokens = re.split(r"(\{|\}|;)", body)
    buf = ""
    pending_else_if = False

    def emit(line):
        out.append("   " + "   " * (indent - 1) + line)

    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in ("{", "}", ";"):
            stmt = buf.strip()
            buf = ""
            if tok == ";":
                if not stmt:
                    i += 1
                    continue
                m = re.match(r"^double\s+(.*)$", stmt, re.S)
                if m:
                    # declaration, possibly with initialisers
                    for part in split_top(m.group(1)):
                        part = part.strip()
                        if "=" in part:
                            lhs, rhs = part.split("=", 1)
                            decls.add(lhs.strip())
                            emit("%s = %s" % (lhs.strip(), convert_expr(rhs)))
                        else:
                            decls.add(part)
                    i += 1
                    continue
                if re.match(r"^return\s+0$", stmt):
                    emit("err = 0")
                    emit("return")
                    i += 1
                    continue
                # A braceless `if (cond) stmt` -- the whole thing arrives as
                # one statement because there is no brace to split on.
                m = re.match(r"^(else\s+)?if\s*\(", stmt)
                if m:
                    cond, rest = split_cond(stmt)
                    if m.group(1) or pending_else_if:
                        emit("else if (%s) then" % cond)
                    else:
                        emit("if (%s) then" % cond)
                    indent += 1
                    emit_assign(emit, rest, decls, name)
                    indent -= 1
                    # A braceless if may still be followed by a braceless
                    # else; closing it here would orphan that else, and the
                    # error only shows up much later at the enclosing block.
                    if re.match(r"^else\b", "".join(tokens[i+1:i+2]).strip()):
                        pending_else_if = True
                    else:
                        emit("end if")
                        pending_else_if = False
                    i += 1
                    continue
                # ...and a braceless `else stmt`.
                m = re.match(r"^else\s+(?!if\b)(.*)$", stmt, re.S)
                if m:
                    emit("else")
                    indent += 1
                    emit_assign(emit, m.group(1).strip(), decls, name)
                    indent -= 1
                    emit("end if")
                    pending_else_if = False
                    i += 1
                    continue
                if "=" in stmt and not re.match(r"^(if|else)\b", stmt):
                    emit_assign(emit, stmt, decls, name)
                    i += 1
                    continue
                raise SystemExit("%s: unhandled statement %r" % (name, stmt))
            elif tok == "{":
                m = re.match(r"^(else\s+)?if\s*\((.*)\)$", stmt, re.S)
                if m:
                    cond = convert_expr(m.group(2)).replace("&&", ".and.").replace("||", ".or.")
                    if m.group(1) or pending_else_if:
                        emit("else if (%s) then" % cond)
                    else:
                        emit("if (%s) then" % cond)
                    indent += 1
                    pending_else_if = False
                elif re.match(r"^else$", stmt):
                    emit("else")
                    indent += 1
                    pending_else_if = False
                elif stmt == "":
                    indent += 1
                else:
                    raise SystemExit("%s: unhandled block opener %r" % (name, stmt))
            else:  # "}"
                indent -= 1
                # look ahead: an `else`/`else if` continues the same construct
                rest = "".join(tokens[i+1:i+3]).strip()
                if re.match(r"^else\b", rest):
                    pending_else_if = True
                else:
                    emit("end if")
                    pending_else_if = False
        else:
            buf += tok
        i += 1

    if buf.strip():
        raise SystemExit("%s: trailing text %r" % (name, buf.strip()))
    return decls, out


def split_cond(stmt):
    """Split `if (cond) rest` into (cond, rest) by matching parentheses."""
    i = stmt.index("(")
    depth = 0
    for j in range(i, len(stmt)):
        if stmt[j] == "(":
            depth += 1
        elif stmt[j] == ")":
            depth -= 1
            if depth == 0:
                cond = stmt[i+1:j]
                return (convert_expr(cond).replace("&&", ".and.").replace("||", ".or."),
                        stmt[j+1:].strip())
    raise SystemExit("unbalanced condition in %r" % stmt)


def emit_assign(emit, stmt, decls, name):
    """One assignment, expanding C's compound operators.

    Fortran has no `*=`, so `F1 *= .5` becomes `F1 = F1 * (.5)`.  The
    parentheses matter: the right-hand side may be a sum.
    """
    m = re.match(r"^(.*?)\s*([-+*/])=\s*(.*)$", stmt, re.S)
    if m and "==" not in stmt and "<=" not in stmt and ">=" not in stmt and "!=" not in stmt:
        lhs, op, rhs = m.group(1).strip(), m.group(2), m.group(3)
        if not re.match(r"^(roots|weights)\[", lhs):
            decls.add(lhs)
        emit("%s = %s %s (%s)" % (convert_expr(lhs), convert_expr(lhs), op, convert_expr(rhs)))
        return
    if "=" not in stmt:
        raise SystemExit("%s: unhandled statement %r" % (name, stmt))
    lhs, rhs = stmt.split("=", 1)
    lhs = lhs.strip()
    if not re.match(r"^(roots|weights)\[", lhs):
        decls.add(lhs)
    emit("%s = %s" % (convert_expr(lhs), convert_expr(rhs)))


def split_top(s):
    """Split on commas that are not inside parentheses."""
    parts, depth, cur = [], 0, ""
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return parts


def wrap(line, limit=118):
    """Break a long Fortran line onto continuations.

    Fortran allows a continuation anywhere, including inside parentheses, so
    the only real constraint is not splitting a token.  The one trap is the
    `-` inside an exponent: breaking `1.96E-01` after the E would silently
    change the number, so operators preceded by e/E/d/D are never break
    points.
    """
    if len(line) <= limit:
        return [line]
    indent = re.match(r"\s*", line).group(0) + "     "
    out = []
    while len(line) > limit:
        cut = -1
        for i in range(len(indent) + 4, min(limit, len(line))):
            ch = line[i]
            if ch in "+-*/" and line[i-1] not in "eEdD(+-*/":
                cut = i
        if cut < 0:
            # no safe break; leave it long rather than corrupt it
            break
        out.append(line[:cut].rstrip() + " &")
        line = indent + line[cut:].lstrip()
    out.append(line)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="src/rys_roots.c")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    src = open(args.src).read()
    with open(args.out, "w") as fh:
        fh.write("""!
! Closed-form Rys roots and weights for one to five roots.
!
! Generated by scripts/c_rysroot_to_fortran.py from src/rys_roots.c -- do not
! edit.  These are polynomial fits, not an algorithm: 1,376 lines of
! coefficients branching on ranges of X.  Translating them mechanically rather
! than by hand is the only way to be sure none of the digits moved, and the
! rys_root_n arm of the D3 tests checks that against the C bit for bit.
!
module cint_rys_root_n
   use cint_const, only: dp
   implicit none
   private

   public :: rys_root1, rys_root2, rys_root3, rys_root4, rys_root5

   real(dp), parameter :: PIE4 = 0.78539816339744827900_dp

contains
""")
        for name in FUNCS:
            body = extract(src, name)
            decls, lines = translate(body, name)
            fh.write("\n   function %s(x, roots, weights) result(err)\n" % name)
            fh.write("      real(dp), intent(in)  :: x\n")
            fh.write("      real(dp), intent(out) :: roots(0:), weights(0:)\n")
            fh.write("      integer :: err\n")
            drop = {"X", "x"}
            names = sorted(d for d in decls if d not in drop)
            if names:
                for k in range(0, len(names), 10):
                    fh.write("      real(dp) :: %s\n" % ", ".join(names[k:k+10]))
            fh.write("\n      err = 0\n")
            for l in lines:
                for w in wrap(l):
                    fh.write(w + "\n")
            fh.write("   end function %s\n" % name)
        fh.write("\nend module cint_rys_root_n\n")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
