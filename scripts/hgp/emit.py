"""Emit the Obara-Saika/HGP recurrences as Fortran, one module per class.

The shape of a kernel:

    subroutine hgp_<kinds>(nbra, ncb, bp, kab, nket, nck, kp, kcd, ab, cd,
                           cutoff, res, any)

    bp(10, nbra)   p, 1/(2p), P-A (3), P (3), exp(-mu Rab^2)/p, mu Rab^2
    kab(ntb*ncb, nbra)   c_i c_j per coefficient-type pair and contraction
    kp(10, nket)   q, 1/(2q), Q-C (3), Q (3), exp(-nu Rcd^2)/q, nu Rcd^2
    kcd(ntk*nck, nket)   as kab
    ab(3), cd(3)   A - B and C - D, for the transfers
    res(ncomp, ncb*nck)  Cartesian integrals, i fastest, columns cA fastest

The vertical recurrence runs once per primitive quartet and its results are
contracted; the transfers run once per contracted quartet.  That split is
the whole point of HGP, and it is why the coefficients multiply in at the
contraction rather than riding in the primitive weight as they do in the
rotated-axis kernels.
"""

from .cart import cart_components, cart_index
from .kinds import KINDS, KIND_RANK, KIND_LMAX, KIND_NTYPE, kind_blocks, kind_ncomp
from .recur import Graph

MAX_LINE = 110


def cname(kinds):
    return "".join(kinds)


def class_code(kinds):
    return int("".join(str(KIND_RANK[k]) for k in kinds))


def wrap(lhs, terms, indent):
    out, cur = [], f"{indent}{lhs} = "
    for i, t in enumerate(terms):
        piece = t if i == 0 else (" - " + t[1:] if t.startswith("-") else " + " + t)
        if len(cur) + len(piece) > MAX_LINE and not cur.endswith("= "):
            out.append(cur + " &")
            cur = indent + "   " + piece.lstrip()
        else:
            cur += piece
    out.append(cur)
    return out


def term_text(coef, syms, ref):
    parts = []
    if abs(coef) != 1:
        parts.append(f"{abs(coef)}.0_dp")
    parts.extend(syms)
    parts.append(ref)
    body = "*".join(parts)
    return ("-" + body) if coef < 0 else body


class ClassPlan:
    """What one class needs: the shared vertical targets and one transfer
    block per pair of coefficient types."""

    def __init__(self, kinds):
        self.kinds = kinds
        self.blocks = [kind_blocks(k) for k in kinds]
        self.ncomp = [kind_ncomp(k) for k in kinds]
        self.ntb = KIND_NTYPE[kinds[0]] * KIND_NTYPE[kinds[1]]
        self.ntk = KIND_NTYPE[kinds[2]] * KIND_NTYPE[kinds[3]]

        self.bra_pairs = [(a, b) for a in self.blocks[0] for b in self.blocks[1]]
        self.ket_pairs = [(c, d) for c in self.blocks[2] for d in self.blocks[3]]

        # every (e, f) some transfer block consumes
        ef = set()
        for (la, _, _, _), (lb, _, _, _) in self.bra_pairs:
            for (lc, _, _, _), (ld, _, _, _) in self.ket_pairs:
                for e in range(la, la + lb + 1):
                    for f in range(lc, lc + ld + 1):
                        ef.add((e, f))
        self.targets = []
        for e, f in sorted(ef):
            for ce in cart_components(e):
                for cf in cart_components(f):
                    self.targets.append((ce, cf))
        self.tindex = {t: i for i, t in enumerate(self.targets)}

        self.g = Graph()
        vroots = [self.g.vrr(ce, cf, 0) for ce, cf in self.targets]
        self.vorder = self.g.order(vroots)
        self.vname = {k: i + 1 for i, k in enumerate(self.vorder)}

        # one transfer block per (bra type pair, ket type pair)
        self.hblocks = []
        for ib, ((la, ta, oa, _), (lb, tb, ob, _)) in enumerate(self.bra_pairs):
            for ik, ((lc, tc, oc, _), (ld, td, od, _)) in enumerate(self.ket_pairs):
                hg = Graph()
                roots = []
                for cd_ in cart_components(ld):
                    for cc_ in cart_components(lc):
                        for cb in cart_components(lb):
                            for ca in cart_components(la):
                                roots.append(((ca, cb, cc_, cd_),
                                              hg.hrr_ket(ca, cb, cc_, cd_)))
                order = hg.order([k for _, k in roots])
                self.hblocks.append(dict(
                    bcol=ta + KIND_NTYPE[kinds[0]] * tb + 1,
                    kcol=tc + KIND_NTYPE[kinds[2]] * td + 1,
                    offs=(oa, ob, oc, od), ls=(la, lb, lc, ld),
                    g=hg, order=order, roots=roots,
                    name={k: i + 1 for i, k in enumerate(order)}))

    def summary(self):
        nv = len(self.vorder)
        nh = sum(len(b["order"]) for b in self.hblocks)
        return (f"{cname(self.kinds)}: vrr={nv} carried={len(self.targets)} "
                f"hrr={nh} in {len(self.hblocks)} block(s) "
                f"comps={self.ncomp[0]*self.ncomp[1]*self.ncomp[2]*self.ncomp[3]}")


def emit_kernel(pl):
    k = pl.kinds
    L = sum(KIND_LMAX[x] for x in k)
    name = f"hgp_{cname(k)}"
    NV = len(pl.vorder)
    NT = len(pl.targets)
    NCOMP = pl.ncomp[0] * pl.ncomp[1] * pl.ncomp[2] * pl.ncomp[3]
    o = []
    o.append(f"   ! ({k[0]}{k[1]}|{k[2]}{k[3]})  {pl.summary()}")
    o.append(f"   subroutine {name}(nbra, ncb, bp, kab, nket, nck, kp, kcd, ab, cd, cutoff, res, any)")
    o.append("      integer,  intent(in)  :: nbra, ncb, nket, nck")
    o.append(f"      real(dp), intent(in)  :: bp(10, nbra), kab({pl.ntb}*ncb, nbra)")
    o.append(f"      real(dp), intent(in)  :: kp(10, nket), kcd({pl.ntk}*nck, nket)")
    o.append("      real(dp), intent(in)  :: ab(3), cd(3), cutoff")
    o.append(f"      real(dp), intent(out) :: res({NCOMP}, ncb*nck)")
    o.append("      logical,  intent(out) :: any")
    o.append(f"      real(dp) :: c({NT}, {pl.ntb}*ncb, {pl.ntk}*nck)")
    o.append(f"      real(dp) :: t({NV}), tv({NT}), f(0:{L}), b(0:{L})")
    o.append("      real(dp) :: p, q, pq, rho, oo2p, oo2q, oo2pq, rp, rq, tt, w, eab, ecd")
    o.append("      real(dp) :: PA0, PA1, PA2, QC0, QC1, QC2, WP0, WP1, WP2, WQ0, WQ1, WQ2")
    o.append("      real(dp) :: PQ0, PQ1, PQ2, AB0, AB1, AB2, CD0, CD1, CD2, wq_, wp_")
    o.append("      integer  :: kq, bq, bc, kc, cc, ck, n, col")
    o.append("")
    o.append("      AB0 = ab(1); AB1 = ab(2); AB2 = ab(3)")
    o.append("      CD0 = cd(1); CD1 = cd(2); CD2 = cd(3)")
    o.append("      any = .false.")
    o.append("      c = 0.0_dp")
    o.append("      do kq = 1, nket")
    o.append("         q = kp(1,kq); oo2q = kp(2,kq)")
    o.append("         QC0 = kp(3,kq); QC1 = kp(4,kq); QC2 = kp(5,kq)")
    o.append("         ecd = kp(10,kq)")
    o.append("         do bq = 1, nbra")
    o.append("            eab = bp(10,bq)")
    o.append("            if (eab + ecd > cutoff) cycle")
    o.append("            any = .true.")
    o.append("            p = bp(1,bq); oo2p = bp(2,bq)")
    o.append("            PA0 = bp(3,bq); PA1 = bp(4,bq); PA2 = bp(5,bq)")
    o.append("            PQ0 = bp(6,bq) - kp(6,kq)")
    o.append("            PQ1 = bp(7,bq) - kp(7,kq)")
    o.append("            PQ2 = bp(8,bq) - kp(8,kq)")
    o.append("            pq = p + q")
    o.append("            rho = p*q/pq")
    o.append("            oo2pq = 0.5_dp/pq")
    o.append("            rp = rho/p; rq = rho/q")
    o.append("            ! W - P = (q/(p+q)) (Q - P), W - Q = (p/(p+q)) (P - Q)")
    o.append("            wp_ = -q/pq; wq_ = p/pq")
    o.append("            WP0 = wp_*PQ0; WP1 = wp_*PQ1; WP2 = wp_*PQ2")
    o.append("            WQ0 = wq_*PQ0; WQ1 = wq_*PQ1; WQ2 = wq_*PQ2")
    o.append("            tt = rho*(PQ0*PQ0 + PQ1*PQ1 + PQ2*PQ2)")
    o.append(f"            call boys(f, tt, {L})")
    o.append("            w = TWO_PI_52*bp(9,bq)*kp(9,kq)/sqrt(pq)")
    o.append(f"            do n = 0, {L}")
    o.append("               b(n) = w*f(n)")
    o.append("            end do")
    # the vertical recurrence
    for key in pl.vorder:
        i = pl.vname[key]
        e = pl.g.expr[key]
        if e is None:
            o.append(f"            t({i}) = b({key[3]})")
        else:
            terms = [term_text(cf, sy, f"t({pl.vname[dep]})") for cf, sy, dep in e]
            o.extend(wrap(f"t({i})", terms, "            "))
    # gather the targets, then contract
    for j, tgt in enumerate(pl.targets):
        o.append(f"            tv({j+1}) = t({pl.vname[('v', tgt[0], tgt[1], 0)]})")
    o.append(f"            do kc = 1, {pl.ntk}*nck")
    o.append(f"               do bc = 1, {pl.ntb}*ncb")
    o.append(f"                  c(1:{NT},bc,kc) = c(1:{NT},bc,kc) + kab(bc,bq)*kcd(kc,kq)*tv(1:{NT})")
    o.append("               end do")
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    o.append("      res = 0.0_dp")
    o.append("      if (.not. any) return")
    # the transfers, once per contracted quartet
    o.append("      do ck = 1, nck")
    o.append("         do cc = 1, ncb")
    o.append("            col = cc + ncb*(ck-1)")
    for blk in pl.hblocks:
        bc = f"{blk['bcol']} + {pl.ntb}*(cc-1)" if pl.ntb > 1 else "cc"
        kc = f"{blk['kcol']} + {pl.ntk}*(ck-1)" if pl.ntk > 1 else "ck"
        o.append(f"            bc = {bc}; kc = {kc}")
        hg, nm = blk["g"], blk["name"]
        NH = len(blk["order"])
        for key in blk["order"]:
            i = nm[key]
            e = hg.expr.get(key)
            if e is None:  # a contracted vertical result
                j = pl.tindex[(key[1], key[2])] + 1
                o.append(f"            h({i}) = c({j},bc,kc)")
            else:
                terms = [term_text(cf, sy, f"h({nm[dep]})") for cf, sy, dep in e]
                o.extend(wrap(f"h({i})", terms, "            "))
        oa, ob, oc, od = blk["offs"]
        la, lb, lc, ld = blk["ls"]
        na, nb, nc = pl.ncomp[0], pl.ncomp[1], pl.ncomp[2]
        for (ca, cb, cc_, cd_), key in blk["roots"]:
            ia = oa + cart_index(*ca); ib = ob + cart_index(*cb)
            ic = oc + cart_index(*cc_); id_ = od + cart_index(*cd_)
            idx = ia + na*(ib + nb*(ic + nc*id_)) + 1
            o.append(f"            res({idx},col) = h({nm[key]})")
    o.append("         end do")
    o.append("      end do")
    o.append(f"   end subroutine {name}")
    return "\n".join(o)


def hmax(pl):
    return max((len(b["order"]) for b in pl.hblocks), default=1)


def emit_class_file(pl):
    name = f"hgp_{cname(pl.kinds)}"
    body = emit_kernel(pl)
    body = body.replace("      integer  :: kq, bq, bc, kc, cc, ck, n, col",
                        f"      real(dp) :: h({hmax(pl)})\n"
                        "      integer  :: kq, bq, bc, kc, cc, ck, n, col")
    o = ["! GENERATED by scripts/hgp -- do not edit.",
         f"! {pl.summary()}",
         f"module {name}_m",
         "   use cint_const, only: dp",
         "   ! the same guarded Boys the rotated-axis kernels use: libcint's",
         "   ! closed form cancels for F_1 at a tiny argument.",
         "   use cint_rotaxis_boys, only: boys",
         "   implicit none",
         "   private",
         f"   public :: {name}",
         "   real(dp), parameter :: PI = 3.14159265358979323846_dp",
         "   real(dp), parameter :: TWO_PI_52 = 2.0_dp*PI**2*sqrt(PI)",
         "contains",
         body,
         f"end module {name}_m"]
    return "\n".join(o) + "\n"
