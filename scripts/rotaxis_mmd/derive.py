"""McMurchie-Davidson in the rotated frame, derived symbolically.

The frame (see doc/ROT_AXIS_MMD.md for the derivation in prose):

    origin at A, z along A->B, y along AB x CD, x = y x z.

    A = (0, 0, 0)          B = (0, 0, R_AB)
    C = (cx, cy, cz)       D = C + (Rs, 0, Rc)      Rs = R_CD sin g, Rc = R_CD cos g

For a primitive bra pair (a, b), p = a + b and the product centre P lies on
the z axis:  P = (0, 0, za) with za = b R_AB / p, and zb = za - R_AB is its
distance from B.  For a primitive ket pair (c, d), q = c + d, yc = c/q,
yd = d/q and Q = C + yd (Rs, 0, Rc) = (ax, cy, az).

The Hermite expansion of a bra pair therefore has X_PA = X_PB = Y_PA =
Y_PB = 0, and that of a ket pair has Y_QC = Y_QD = 0.  That is the whole
point of the rotation: the x and y expansions carry only the 1/(2p),
1/(2q) ladder terms, and every geometry-dependent quantity of the bra pair
is a scalar (za, zb) rather than a vector.

Symbols, grouped by what they depend on:

    bra primitive:   ip = 1/(2p), za, zb, zq = Z_PQ = za - az, B0..BL
                     (B_n = (-2 alpha)^n F_n(T) times the bra-pair weight)
    ket primitive:   iq = 1/(2q), yc, yd, ax
    quartet only:    Rs, Rc, cy

X_PQ = P_x - Q_x = -ax,  Y_PQ = -cy,  Z_PQ = zq.

The expansion is exactly Helgaker, Jorgensen and Olsen eq. 9.9.33,

    (ab|cd) = pref * sum_tuv E^ab_tuv sum_TUV (-1)^(T+U+V) E^cd_TUV R_{t+T,u+U,v+V}

with R^n_{000} = (-2 alpha)^n F_n and the recursions 9.9.18-20.  The
prefactor 2 pi^(5/2) / (p q sqrt(p+q)) K_ab K_cd is split between the pair
weights and the driver, see emit.py.
"""

from functools import lru_cache

from .poly import Ring


def cart_components(l):
    """libcint's Cartesian order: x-power descending, then y-power."""
    out = []
    for lx in range(l, -1, -1):
        for ly in range(l - lx, -1, -1):
            out.append((lx, ly, l - lx - ly))
    return out


# A slot of the quartet is a KIND, not just an angular momentum: an L shell
# (libcint's KAPPA_SP_SHELL, s and p on shared exponents with a coefficient
# column for each) is one slot of four components whose s component takes
# the s coefficient and whose p components take the p coefficient.  Each
# kind lists (Cartesian exponents, coefficient type) in libcint's output
# order; the type is the column block of the contraction coefficient.
KINDS = {
    "s": [((0, 0, 0), 0)],
    "p": [(c, 0) for c in cart_components(1)],
    "d": [(c, 0) for c in cart_components(2)],
    "f": [(c, 0) for c in cart_components(3)],
    "L": [((0, 0, 0), 0)] + [(c, 1) for c in cart_components(1)],
}
KIND_LMAX = {"s": 0, "p": 1, "d": 2, "f": 3, "L": 1}
KIND_NTYPES = {"s": 1, "p": 1, "d": 1, "f": 1, "L": 2}
# The canonical order of kinds within a pair and of pairs within a quartet.
KIND_RANK = {"s": 0, "p": 1, "L": 2, "d": 3, "f": 4}


class Derivation:
    def __init__(self, ka, kb, kc, kd):
        self.kinds = (ka, kb, kc, kd)
        self.l = tuple(KIND_LMAX[k] for k in self.kinds)
        self.L = sum(self.l)
        L = self.L
        self.bra_names = ["ip", "za", "zb", "zq"] + [f"B{n}" for n in range(L + 1)]
        self.ket_names = ["iq", "yc", "yd", "ax"]
        self.const_names = ["Rs", "Rc", "cy"]
        self.ring = Ring(self.bra_names + self.ket_names + self.const_names)
        R = self.ring
        self.s = {n: R.sym(n) for n in R.names}
        self._e_cache = {}
        self._r_cache = {}

    # ---- Hermite expansion coefficients, one Cartesian direction ---------

    def E(self, axis, side, i, j, t):
        """E^{ij}_t along `axis` ('x','y','z') for side 'bra' or 'ket'."""
        key = (axis, side, i, j, t)
        if key in self._e_cache:
            return self._e_cache[key]
        R = self.ring
        s = self.s
        if side == "bra":
            half = s["ip"]
            if axis == "z":
                xpa, xpb = s["za"], s["zb"]
            else:
                xpa, xpb = R.zero(), R.zero()
        else:
            half = s["iq"]
            if axis == "x":
                xpa, xpb = s["yd"] * s["Rs"], -s["yc"] * s["Rs"]
            elif axis == "z":
                xpa, xpb = s["yd"] * s["Rc"], -s["yc"] * s["Rc"]
            else:
                xpa, xpb = R.zero(), R.zero()

        if t < 0 or t > i + j:
            val = R.zero()
        elif i == 0 and j == 0:
            val = R.one() if t == 0 else R.zero()
        elif i > 0:
            val = (half * self.E(axis, side, i - 1, j, t - 1)
                   + xpa * self.E(axis, side, i - 1, j, t)
                   + (t + 1) * self.E(axis, side, i - 1, j, t + 1))
        else:
            val = (half * self.E(axis, side, i, j - 1, t - 1)
                   + xpb * self.E(axis, side, i, j - 1, t)
                   + (t + 1) * self.E(axis, side, i, j - 1, t + 1))
        self._e_cache[key] = val
        return val

    # ---- Hermite Coulomb integrals R^n_{tuv} ------------------------------

    def R(self, n, t, u, v):
        key = (n, t, u, v)
        if key in self._r_cache:
            return self._r_cache[key]
        s = self.s
        if t == 0 and u == 0 and v == 0:
            val = s[f"B{n}"]
        elif t > 0:
            val = -s["ax"] * self.R(n + 1, t - 1, u, v)
            if t > 1:
                val = val + (t - 1) * self.R(n + 1, t - 2, u, v)
        elif u > 0:
            val = -s["cy"] * self.R(n + 1, t, u - 1, v)
            if u > 1:
                val = val + (u - 1) * self.R(n + 1, t, u - 2, v)
        else:
            val = s["zq"] * self.R(n + 1, t, u, v - 1)
            if v > 1:
                val = val + (v - 1) * self.R(n + 1, t, u, v - 2)
        self._r_cache[key] = val
        return val

    # ---- the integral over one Cartesian quartet -------------------------

    @lru_cache(maxsize=None)
    def _pair_product(self, axis, ia, ib, ic, id_):
        """sum over t,T of E^ab_t E^cd_T (-1)^T, keyed by t+T."""
        out = {}
        for t in range(ia + ib + 1):
            eab = self.E(axis, "bra", ia, ib, t)
            if eab.is_zero():
                continue
            for T in range(ic + id_ + 1):
                ecd = self.E(axis, "ket", ic, id_, T)
                if ecd.is_zero():
                    continue
                term = eab * ecd
                if T % 2:
                    term = -term
                k = t + T
                out[k] = out[k] + term if k in out else term
        return {k: v for k, v in out.items() if not v.is_zero()}

    def integral(self, ca, cb, cc, cd):
        """The (ab|cd) polynomial for Cartesian exponents ca..cd (3-tuples)."""
        qx = self._pair_product("x", ca[0], cb[0], cc[0], cd[0])
        qy = self._pair_product("y", ca[1], cb[1], cc[1], cd[1])
        qz = self._pair_product("z", ca[2], cb[2], cc[2], cd[2])
        total = self.ring.zero()
        for t, px in qx.items():
            for u, py in qy.items():
                pxy = px * py
                for v, pz in qz.items():
                    total = total + pxy * pz * self.R(0, t, u, v)
        return total

    def all_components(self):
        """Every component of the class, in libcint's (i,j,k,l) column-major
        order (i fastest), as ((exponents, coefficient types), polynomial)."""
        ka, kb, kc, kd = (KINDS[k] for k in self.kinds)
        comps = []
        for cd, td in kd:
            for cc, tc in kc:
                for cb, tb in kb:
                    for ca, ta in ka:
                        comps.append((((ca, cb, cc, cd), (ta, tb, tc, td)),
                                      self.integral(ca, cb, cc, cd)))
        return comps
