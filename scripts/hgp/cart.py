"""The Cartesian component order -- now a thin shim over scripts/cartorder.py.

`cart_index` used to be an independent closed-form inverse of
`cart_components`'s ordering.  Nothing checked that the two agreed, and they
agree only for libcint's order; see cartorder.py.  Both now come from one
definition, so an added order cannot desynchronise them.
"""

from cartorder import components as _components, index as _index


def cart_components(l, order="libcint"):
    return _components(l, order)


def cart_index(lx, ly, lz, order="libcint"):
    return _index(lx, ly, lz, order)
