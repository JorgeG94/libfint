"""Obara-Saika vertical recurrence and the Head-Gordon-Pople horizontal one.

The two together are HGP's split: build [e0|f0] -- all angular momentum on
A and C -- inside the primitive loop by the OS vertical recurrence,
contract those, and only then move momentum onto B and D by the horizontal
recurrence, which depends on A-B and C-D alone and so costs nothing per
primitive.

Unlike scripts/rotaxis_mmd, nothing here is expanded into a polynomial.
Each intermediate is a node computed once from earlier nodes, so the
generated code is the recurrence itself, straight-line, with common
subexpressions shared by construction.  Expanding [dd|dd] instead is what
produced the 96k-term tables that lost to Rys.

Node keys:
    ('v', e, f, m)    OS intermediate [e0|f0]^(m), e and f Cartesian
                      exponent triples, m the auxiliary index
    ('h', a, b, f)    bra-transferred (a b| f0), contracted
    ('g', a, b, c, d) fully transferred (a b| c d), contracted

Symbols an expression may carry, all per primitive quartet except the
last two, which are per shell quartet:
    PA0..2, WP0..2   P - A and W - P
    QC0..2, WQ0..2   Q - C and W - Q
    oo2p, oo2q       1/(2p), 1/(2q)
    rp, rq           rho/p, rho/q
    oo2pq            1/(2(p+q))
    AB0..2, CD0..2   A - B and C - D
"""

from .cart import cart_components, cart_index


def _first_dir(t):
    for i in range(3):
        if t[i] > 0:
            return i
    return None


def _dec(t, i, n=1):
    return tuple(t[j] - n if j == i else t[j] for j in range(3))


def _inc(t, i):
    return tuple(t[j] + 1 if j == i else t[j] for j in range(3))


class Graph:
    """A DAG of intermediates.  `add` returns a node key; `order` gives the
    nodes in an order where every node's inputs come first."""

    def __init__(self):
        self.expr = {}      # key -> list of (int coefficient, (symbols...), key)
        self._stack = []

    # ---- vertical: [e0|f0]^(m) ------------------------------------------

    def vrr(self, e, f, m):
        key = ('v', e, f, m)
        if key in self.expr:
            return key
        if sum(e) == 0 and sum(f) == 0:
            self.expr[key] = None          # a base: b(m)
            return key
        terms = []
        if sum(e) > 0:
            i = _first_dir(e)
            e1 = _dec(e, i)
            terms.append((1, (f'PA{i}',), self.vrr(e1, f, m)))
            terms.append((1, (f'WP{i}',), self.vrr(e1, f, m + 1)))
            if e[i] > 1:
                e2 = _dec(e, i, 2)
                terms.append((e[i] - 1, ('oo2p',), self.vrr(e2, f, m)))
                terms.append((-(e[i] - 1), ('oo2p', 'rp'), self.vrr(e2, f, m + 1)))
            if f[i] > 0:
                terms.append((f[i], ('oo2pq',), self.vrr(e1, _dec(f, i), m + 1)))
        else:
            i = _first_dir(f)
            f1 = _dec(f, i)
            terms.append((1, (f'QC{i}',), self.vrr(e, f1, m)))
            terms.append((1, (f'WQ{i}',), self.vrr(e, f1, m + 1)))
            if f[i] > 1:
                f2 = _dec(f, i, 2)
                terms.append((f[i] - 1, ('oo2q',), self.vrr(e, f2, m)))
                terms.append((-(f[i] - 1), ('oo2q', 'rq'), self.vrr(e, f2, m + 1)))
        self.expr[key] = terms
        return key

    # ---- horizontal, on contracted quantities ---------------------------

    def hrr_bra(self, a, b, f):
        """(a b| f0) from ((a+1) b'| f0) and (a b'| f0)."""
        if sum(b) == 0:
            return ('v', a, f, 0)          # a contracted vertical result
        key = ('h', a, b, f)
        if key in self.expr:
            return key
        i = _first_dir(b)
        b1 = _dec(b, i)
        self.expr[key] = [
            (1, (), self.hrr_bra(_inc(a, i), b1, f)),
            (1, (f'AB{i}',), self.hrr_bra(a, b1, f)),
        ]
        return key

    def hrr_ket(self, a, b, c, d):
        if sum(d) == 0:
            return self.hrr_bra(a, b, c)   # nothing left to move onto D
        key = ('g', a, b, c, d)
        if key in self.expr:
            return key
        i = _first_dir(d)
        d1 = _dec(d, i)
        self.expr[key] = [
            (1, (), self.hrr_ket(a, b, _inc(c, i), d1)),
            (1, (f'CD{i}',), self.hrr_ket(a, b, c, d1)),
        ]
        return key

    # ---- ordering --------------------------------------------------------

    def order(self, roots):
        seen, out = set(), []

        def walk(k):
            if k in seen:
                return
            seen.add(k)
            for _, _, dep in (self.expr.get(k) or []):
                walk(dep)
            out.append(k)

        for r in roots:
            walk(r)
        return out


def vrr_targets(la, lb, lc, ld):
    """The [e0|f0] the horizontal recurrence will consume."""
    out = []
    for e in range(la, la + lb + 1):
        for f in range(lc, lc + ld + 1):
            for ce in cart_components(e):
                for cf in cart_components(f):
                    out.append((ce, cf))
    return out


def build(la, lb, lc, ld):
    """Return (graph, vrr_roots, targets) for one plain class.

    `targets` is the list of ((a,b,c,d) exponent tuples, node key) in
    libcint's component order, i fastest.
    """
    g = Graph()
    # the vertical part, over every [e0|f0] the transfers will need
    vroots = [g.vrr(ce, cf, 0) for ce, cf in vrr_targets(la, lb, lc, ld)]
    # the transfers
    targets = []
    for cd_ in cart_components(ld):
        for cc_ in cart_components(lc):
            for cb in cart_components(lb):
                for ca in cart_components(la):
                    targets.append(((ca, cb, cc_, cd_),
                                    g.hrr_ket(ca, cb, cc_, cd_)))
    return g, vroots, targets
