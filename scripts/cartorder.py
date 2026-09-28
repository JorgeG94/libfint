"""The Cartesian component order of a shell, defined once.

Two things have to agree about the order of Cartesian components within a
shell: the LIST that fixes the output layout, and the INVERSE that maps an
exponent triple back to its offset in that layout.  They used to be written
independently -- the list as a comprehension in two separate copies, the
inverse as a closed-form piece of arithmetic -- and nothing checked that the
two agreed.  They agree only for libcint's order, and nothing said so:
reordering the list without rewriting the arithmetic misplaces 5 of the 6 d
components and 9 of the 10 f components, with no error raised and
plausible-looking integrals coming out.

So here `index` is DERIVED from `components` by lookup.  The two cannot
disagree, and adding an order means adding one function to `ORDERS`.

`order` is an explicit parameter with a default, never module state.  A
build-wide ordering switch spanning several generators would otherwise
cross-contaminate silently the first time one process emitted two orders --
the failure would look like a miscompiled kernel rather than a stale global.
"""

from functools import lru_cache


def _libcint(l):
    """x-power descending, then y-power.  libcint's CINTcart_comp order.

    l=2 is XX XY XZ YY YZ ZZ.
    """
    return [(lx, ly, l - lx - ly)
            for lx in range(l, -1, -1)
            for ly in range(l - lx, -1, -1)]


def _gamess(l):
    """GAMESS's order: the pure powers first, then the mixed ones in
    libcint's relative order.

    l=2 is XX YY ZZ XY XZ YZ.  That was read out of gamess-libERI's own
    kernels rather than assumed: in rhf/rot_axis/int0022.F90 the 6x6 block's
    diagonal sits at eri_value(1, 8, 15), and eri_value(1) carries qx**4,
    eri_value(8) carries no geometry at all (the y direction has none in
    that frame) and eri_value(15) carries qz**4.  s and p coincide with
    libcint, so only d actually differs.

    l >= 3 is deliberately absent.  GAMESS's f convention has not been
    confirmed against the kernels, and guessing it would put a silent
    permutation into every f integral -- the exact failure this module
    exists to prevent.  Confirm it the same way d was confirmed, then add
    it here.
    """
    if l > 2:
        raise NotImplementedError(
            f"the gamess order is defined here only for l <= 2, not l={l}. "
            "Confirm the f ordering against gamess-libERI/rhf/rot_axis "
            "before adding it, the way d was confirmed from int0022.F90.")
    cs = _libcint(l)
    pure = [c for c in cs if sum(1 for e in c if e) <= 1]
    return pure + [c for c in cs if c not in pure]


ORDERS = {"libcint": _libcint, "gamess": _gamess}


@lru_cache(maxsize=None)
def components(l, order="libcint"):
    """The Cartesian exponent triples of angular momentum `l`, in `order`.

    Returns a tuple: it is cached, so a caller must not be able to mutate it.
    """
    try:
        f = ORDERS[order]
    except KeyError:
        raise ValueError(f"unknown Cartesian order {order!r}; "
                         f"have {sorted(ORDERS)}") from None
    return tuple(f(l))


@lru_cache(maxsize=None)
def _index_map(l, order):
    return {c: i for i, c in enumerate(components(l, order))}


def index(lx, ly, lz, order="libcint"):
    """The offset of (lx,ly,lz) within its shell, in `order`.

    A lookup into `components`, not an independent formula.  That is the
    whole point of this module.
    """
    return _index_map(lx + ly + lz, order)[(lx, ly, lz)]
