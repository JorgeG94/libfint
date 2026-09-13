"""Shell kinds, as scripts/rotaxis_mmd has them: s, p, d and L, where an L
shell is one slot carrying s and p on shared exponents with a coefficient
column for each."""

from .cart import cart_components

# kind -> list of (l, coefficient type) in libcint's component order
KINDS = {
    "s": [(0, 0)],
    "p": [(1, 0)],
    "d": [(2, 0)],
    "L": [(0, 0), (1, 1)],
}
KIND_RANK = {"s": 0, "p": 1, "L": 2, "d": 3}
KIND_LMAX = {"s": 0, "p": 1, "d": 2, "L": 1}
KIND_NTYPE = {"s": 1, "p": 1, "d": 1, "L": 2}


def kind_blocks(k):
    """(l, type, component offset, ncart) per distinct l of the kind."""
    out, off = [], 0
    for l, t in KINDS[k]:
        n = len(cart_components(l))
        out.append((l, t, off, n))
        off += n
    return out


def kind_ncomp(k):
    return sum(n for _, _, _, n in kind_blocks(k))
