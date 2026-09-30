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


# Transcribed from GAMESS's own table rather than inferred from a rule.  The
# order turns out to be a grouping by NORMALISATION FACTOR -- see _gamess --
# in blocks of constant factor, which is why g puts the (2,2,0) family between
# the (3,1,0) and the (2,1,1) ones.  That is a useful consistency check on the
# transcription, but it is not how the table was obtained, and a rule fitted to
# d and f alone would still have been a silent permutation at g.  Source: gamess-jorge/source/int2a.src, subroutine
# SHELLS, the DATA LX/LY/LZ tables -- flat indices 5-10 for d, 11-20 for f,
# 21-35 for g, read across all three tables.  LX alone does not determine the
# triple wherever it is 0.
#
# d and f are corroborated by a second, independent implementation:
# gamess-libERI/rhf/rys/int3000_rysgen.F90 lines 69-116, whose ix/iy/iz
# offsets minus one reproduce the f triples exactly, and d additionally by
# the eri_value(1, 8, 15) diagonal of rhf/rot_axis/int0022.F90.  g rests on
# int2a.src alone -- which is GAMESS itself, so it is the authority rather
# than a derived copy, but it has not been cross-checked against a second
# implementation.
_GAMESS_TABLE = {
    0: [(0, 0, 0)],
    1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    2: [(2, 0, 0), (0, 2, 0), (0, 0, 2),
        (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    3: [(3, 0, 0), (0, 3, 0), (0, 0, 3),
        (2, 1, 0), (2, 0, 1), (1, 2, 0), (0, 2, 1), (1, 0, 2), (0, 1, 2),
        (1, 1, 1)],
    4: [(4, 0, 0), (0, 4, 0), (0, 0, 4),
        (3, 1, 0), (3, 0, 1), (1, 3, 0), (0, 3, 1), (1, 0, 3), (0, 1, 3),
        (2, 2, 0), (2, 0, 2), (0, 2, 2),
        (2, 1, 1), (1, 2, 1), (1, 1, 2)],
}


def _gamess(l):
    """GAMESS's Cartesian component order, l = 0..4.

    s and p coincide with libcint; d is XX YY ZZ XY XZ YZ; f is
    XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ.  See _GAMESS_TABLE for sources.

    NOTE, because it has bitten people: ORDER IS NOT THE ONLY DIFFERENCE.
    GAMESS also scales components individually within a shell, where libcint
    (and so libfint) normalise per l only.  int2a.src's GENRAL carries one
    running factor DUM1 across the component loop: its computed GO TO (near
    line 1251, indexed by the flat component number) sends the first component
    of each l to a label that SETS DUM1 from the contraction coefficient, some
    later components to labels that MULTIPLY it, and everything else to 220,
    which leaves it alone.  So the factor is cumulative and applies to a run of
    components, giving

        d   1 for XX YY ZZ            sqrt3 for XY XZ YZ
        f   1 for XXX YYY ZZZ         sqrt5 for the six XXY-type
            sqrt15 for XYZ            (sqrt5 then sqrt3, compounded)
        g   1 for XXXX YYYY ZZZZ      sqrt7 for the six XXXY-type
            sqrt(35/3) for XXYY XXZZ YYZZ    sqrt35 for XXYZ XYYZ XYZZ

    which is exactly the block structure of the order above -- the two are the
    same grouping.  Confirmed in gamess-libERI too: int0030_ericgen.F90's
    `angl` array is 1, sqrt3 ... and 1, sqrt5 ... sqrt15 by component type, in
    that code's own slot order, and shell_pair.F90 carries the d sqrt3.

    This module still does ONE thing: it says which component sits in which slot.
    It deliberately does not touch normalisation, because conflating a
    permutation with a diagonal rescale is how you get a bug that looks like
    a permutation bug.  Anyone comparing element by element against GAMESS
    needs both, and the rescale belongs at the boundary, as a per-AO diagonal.
    """
    if l not in _GAMESS_TABLE:
        raise NotImplementedError(
            f"the gamess order is transcribed here only for l <= 4, not l={l}. "
            "Read it out of int2a.src's DATA LX/LY/LZ tables (subroutine "
            "SHELLS) across all three arrays before adding it.")
    return list(_GAMESS_TABLE[l])


for _l, _t in _GAMESS_TABLE.items():
    assert sorted(_t) == sorted(_libcint(_l)), f"gamess table for l={_l} is not a permutation"
    assert len(set(_t)) == len(_t), f"gamess table for l={_l} repeats a component"

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
