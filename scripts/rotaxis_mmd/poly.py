"""Sparse multivariate polynomials with exact rational coefficients.

Deliberately dependency-free: sympy is not a build requirement of libfint
and the algebra needed here -- products, sums and grouping of monomials --
is small enough that a dictionary of exponent tuples does it exactly.

A Ring fixes the symbol order once; a monomial is a tuple of exponents in
that order and a Poly maps monomials to Fraction coefficients.  Every
coefficient that arises in the McMurchie-Davidson expansion is an integer
(the recursions only ever multiply by integers), but Fraction keeps the
normalisation in the factoriser exact if that ever changes.
"""

from fractions import Fraction


class Ring:
    def __init__(self, names):
        self.names = list(names)
        self.index = {n: i for i, n in enumerate(self.names)}
        self.n = len(self.names)

    def sym(self, name):
        e = [0] * self.n
        e[self.index[name]] = 1
        return Poly(self, {tuple(e): Fraction(1)})

    def const(self, c):
        c = Fraction(c)
        if c == 0:
            return Poly(self, {})
        return Poly(self, {tuple([0] * self.n): c})

    def zero(self):
        return Poly(self, {})

    def one(self):
        return self.const(1)

    def mono_str(self, mono, sep="*"):
        parts = []
        for name, e in zip(self.names, mono):
            if e == 1:
                parts.append(name)
            elif e > 1:
                parts.append(f"{name}^{e}")
        return sep.join(parts) if parts else "1"


class Poly:
    __slots__ = ("ring", "terms")

    def __init__(self, ring, terms):
        self.ring = ring
        self.terms = terms

    def is_zero(self):
        return not self.terms

    def copy(self):
        return Poly(self.ring, dict(self.terms))

    def __add__(self, other):
        if not isinstance(other, Poly):
            other = self.ring.const(other)
        t = dict(self.terms)
        for m, c in other.terms.items():
            v = t.get(m, 0) + c
            if v == 0:
                t.pop(m, None)
            else:
                t[m] = v
        return Poly(self.ring, t)

    __radd__ = __add__

    def __neg__(self):
        return Poly(self.ring, {m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        if not isinstance(other, Poly):
            other = self.ring.const(other)
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if not isinstance(other, Poly):
            c = Fraction(other)
            if c == 0:
                return self.ring.zero()
            return Poly(self.ring, {m: v * c for m, v in self.terms.items()})
        t = {}
        for m1, c1 in self.terms.items():
            for m2, c2 in other.terms.items():
                m = tuple(a + b for a, b in zip(m1, m2))
                v = t.get(m, 0) + c1 * c2
                if v == 0:
                    t.pop(m, None)
                else:
                    t[m] = v
        return Poly(self.ring, t)

    __rmul__ = __mul__

    def __pow__(self, n):
        r = self.ring.one()
        for _ in range(n):
            r = r * self
        return r

    def __eq__(self, other):
        if not isinstance(other, Poly):
            other = self.ring.const(other)
        return self.terms == other.terms

    def __hash__(self):
        return hash(frozenset(self.terms.items()))

    def __repr__(self):
        if not self.terms:
            return "0"
        out = []
        for m in sorted(self.terms):
            c = self.terms[m]
            s = self.ring.mono_str(m)
            if s == "1":
                out.append(f"{c}")
            elif c == 1:
                out.append(s)
            elif c == -1:
                out.append(f"-{s}")
            else:
                out.append(f"{c}*{s}")
        return " + ".join(out)
