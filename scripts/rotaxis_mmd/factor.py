"""Factorise the rotated-frame expansion into the three loop levels.

Every term of every Cartesian component is a product of

    (bra-primitive monomial) x (ket-primitive monomial) x (quartet constant)

and the kernel is a nest of two contraction loops, ket outside, bra inside.
So:

  level 1, inside the bra loop:  S[g,n,k] += kab * g(ip,za,zb) * B_n * zq^k
      one accumulator per distinct bra monomial.  Nothing here depends on
      which Cartesian component is being built, which is the reason the
      inner loop of a (pp|pp) class is a few dozen multiply-adds and not a
      few hundred.

  level 2, once per ket primitive after the bra loop:
      r[j] += kcd * P_j(iq, yc, yd, ax; S)
      one accumulator per distinct ket polynomial.  Two components that
      need the same combination -- and many do, because a ket p_x and a
      bra p_x pull the same R_{100} -- share it.

  level 3, once per quartet:  comp = sum  scale * (Rs^a Rc^b cy^c) * r[j]

The factoriser is exact: it groups monomials, it never approximates.  The
scale is a rational number, in practice always an integer.
"""

from fractions import Fraction


class Factorised:
    def __init__(self, derivation):
        self.d = derivation
        R = derivation.ring
        nb = len(derivation.bra_names)
        nk = len(derivation.ket_names)
        self.nb, self.nk = nb, nk
        self.L = derivation.L

        comps = derivation.all_components()
        self.components = [c for c, _ in comps]
        from .derive import KIND_NTYPES
        self.ntypes = tuple(KIND_NTYPES[k] for k in derivation.kinds)
        # bra type-pair index tt = ta + nta*tb, and the ket likewise
        self.ntt_bra = self.ntypes[0] * self.ntypes[1]
        self.ntt_ket = self.ntypes[2] * self.ntypes[3]

        # bra monomials -> S index.  A bra monomial is (g exponents over
        # ip,za,zb; n; k) and appears with exactly one B_n to the first power.
        self.s_index = {}
        self.s_list = []           # entries (g_exps, n, k, tt_bra)
        # ket monomials -> temporaries
        self.km_index = {}
        self.km_list = []
        # constant monomials -> temporaries
        self.cm_index = {}
        self.cm_list = []
        # ket-level accumulators
        self.r_index = {}
        self.r_list = []           # entries (terms, tt_ket); terms = tuple of ((km_id, s_id), coef)
        # assembly: per component, list of (scale, cm_id, r_id)
        self.assembly = []

        for (comp, types), poly in comps:
            tt_bra = types[0] + self.ntypes[0] * types[1]
            tt_ket = types[2] + self.ntypes[2] * types[3]
            groups = {}
            for mono, coef in poly.terms.items():
                bra = mono[:nb]
                ket = mono[nb:nb + nk]
                cst = mono[nb + nk:]
                s_id = self._s_id(bra, tt_bra)
                km_id = self._km_id(ket)
                cm_id = self._cm_id(cst)
                groups.setdefault(cm_id, {})
                key = (km_id, s_id)
                groups[cm_id][key] = groups[cm_id].get(key, 0) + coef
            terms = []
            for cm_id, sub in sorted(groups.items()):
                items = sorted((k, c) for k, c in sub.items() if c != 0)
                if not items:
                    continue
                lead = items[0][1]
                norm = (tuple((k, c / lead) for k, c in items), tt_ket)
                if norm not in self.r_index:
                    self.r_index[norm] = len(self.r_list)
                    self.r_list.append(norm)
                terms.append((lead, cm_id, self.r_index[norm]))
            self.assembly.append(terms)

        # the distinct g monomials (ip, za, zb) actually used at level 1
        self.g_index = {}
        self.g_list = []
        self.s_g = []
        for g, n, k, _ in self.s_list:
            if g not in self.g_index:
                self.g_index[g] = len(self.g_list)
                self.g_list.append(g)
            self.s_g.append(self.g_index[g])
        self.max_zq_power = max((k for _, _, k, _ in self.s_list), default=0)

    def _s_id(self, bra, tt_bra):
        ip, za, zb, zq = bra[:4]
        bs = bra[4:]
        nz = [i for i, e in enumerate(bs) if e]
        assert len(nz) == 1 and bs[nz[0]] == 1, "R must be linear in the Boys terms"
        key = ((ip, za, zb), nz[0], zq, tt_bra)
        if key not in self.s_index:
            self.s_index[key] = len(self.s_list)
            self.s_list.append(key)
        return self.s_index[key]

    def _km_id(self, ket):
        if ket not in self.km_index:
            self.km_index[ket] = len(self.km_list)
            self.km_list.append(ket)
        return self.km_index[ket]

    def _cm_id(self, cst):
        if cst not in self.cm_index:
            self.cm_index[cst] = len(self.cm_list)
            self.cm_list.append(cst)
        return self.cm_index[cst]

    # ---- reporting -------------------------------------------------------

    def summary(self):
        d = self.d
        name = "".join(d.kinds)
        nterms = sum(len(r) for r, _ in self.r_list)
        nasm = sum(len(a) for a in self.assembly)
        return (f"{name}: L={self.L} comps={len(self.components)} "
                f"S={len(self.s_list)} g={len(self.g_list)} "
                f"r={len(self.r_list)} (terms {nterms}) km={len(self.km_list)} "
                f"cm={len(self.cm_list)} assembly terms={nasm}")
