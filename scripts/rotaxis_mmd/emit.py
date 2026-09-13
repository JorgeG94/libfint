"""Emit the factorised classes as one Fortran module, cint_rotaxis_kernels.

Every kernel has the same shape:

    subroutine rotaxis_LLLL(nbra, ncb, bp, kab, nket, nck, kp, kcd, gc, cutoff, res, any)

    bp(5, nbra)     bra pair table: p, 1/(2p), za, zb, mu R_AB^2
    kab(NTT*ncb, nbra)  bra pair weight c_i c_j exp(-mu R_AB^2) / p, one
                    per coefficient-type pair tt (an L shell has two types,
                    s and p; NTT is their product over the two slots) and
                    contraction pair cc, at tt + NTT*(cc-1) + 1
    kp(7, nket)     ket pair table: q, 1/(2q), c/q, d/q, ax, az, nu R_CD^2
    kcd(NTT*nck, nket)  ket pair weight c_k c_l exp(-nu R_CD^2) / q, laid
                    out as kab
    gc(3)           Rs, Rc, cy
    res(ncomp, ncb*nck)   local-frame Cartesian integrals, column
                    (ccb-1) + ncb*(cck-1) + 1, WITHOUT the 2 pi^(5/2) and
                    the per-l normalisation -- the driver applies those.

The 1/sqrt(p+q) of the prefactor is applied per primitive quartet, so that
B_n = (-2 alpha)^n F_n(T) / sqrt(p+q); with the 1/p in kab and the 1/q in
kcd that is the full 1/(p q sqrt(p+q)).
"""

from fractions import Fraction

from .derive import Derivation, KIND_RANK
from .factor import Factorised

MAX_LINE = 100


def lname(kinds):
    return "".join(kinds)


def class_code(kinds):
    """The integer the driver dispatches on: the kind ranks as digits."""
    return int("".join(str(KIND_RANK[k]) for k in kinds))


def fnum(c):
    c = Fraction(c)
    if c.denominator == 1:
        v = c.numerator
        return f"{v}.0_dp" if v >= 0 else f"({v}.0_dp)"
    return f"({c.numerator}.0_dp/{c.denominator}.0_dp)"


def wrap(lhs, terms, indent="      "):
    """`lhs = t1 + t2 + ...`, broken at blanks into continuation lines, and
    split into several statements when one would exceed the continuation
    limit (a term may itself be a long parenthesised sum)."""
    out = []
    chunks = [terms[i:i + 120] for i in range(0, len(terms), 120)]
    for si, ts in enumerate(chunks):
        expr = "" if si == 0 else lhs
        for i, t in enumerate(ts):
            if i == 0 and si == 0:
                expr = t
            elif t.startswith("-"):
                expr += " - " + t[1:]
            else:
                expr += " + " + t
        cur = f"{indent}{lhs} = "
        for tok in expr.split(" "):
            if len(cur) + len(tok) + 1 > MAX_LINE and not cur.endswith("= "):
                out.append(cur.rstrip() + " &")
                cur = indent + "   " + tok + " "
            else:
                cur += tok + " "
        out.append(cur.rstrip())
    return out


def mono_expr(names, exps, powers=None):
    """A monomial as a Fortran product, using precomputed power names."""
    parts = []
    for name, e in zip(names, exps):
        if e == 0:
            continue
        if powers and (name, e) in powers:
            parts.append(powers[(name, e)])
        elif e == 1:
            parts.append(name)
        else:
            parts.extend([name] * e)
    return "*".join(parts) if parts else None


def coef_term(c, expr):
    """c * expr with 1 and -1 folded."""
    c = Fraction(c)
    if expr is None:
        return fnum(c)
    if c == 1:
        return expr
    if c == -1:
        return f"-{expr}"
    return f"{fnum(c)}*{expr}"


def emit_kernel(f: Factorised):
    d = f.d
    k = d.kinds
    L = d.L
    name = f"rotaxis_{lname(k)}"
    NS, NR, NG = len(f.s_list), len(f.r_list), len(f.g_list)
    NKM, NCM = len(f.km_list), len(f.cm_list)
    NCOMP = len(f.components)
    NTTB, NTTK = f.ntt_bra, f.ntt_ket
    o = []
    o.append(f"   ! ({k[0]}{k[1]}|{k[2]}{k[3]})  "
             f"L={L}: {NS} bra accumulators, {NR} ket accumulators, {NCOMP} components")
    o.append(f"   subroutine {name}(nbra, ncb, bp, kab, nket, nck, kp, kcd, gc, cutoff, res, any)")
    o.append("      integer,  intent(in)  :: nbra, ncb, nket, nck")
    o.append(f"      real(dp), intent(in)  :: bp(5, nbra), kab({NTTB}*ncb, nbra)")
    o.append(f"      real(dp), intent(in)  :: kp(7, nket), kcd({NTTK}*nck, nket), gc(3), cutoff")
    o.append(f"      real(dp), intent(out) :: res({NCOMP}, ncb*nck)")
    o.append("      logical,  intent(out) :: any")
    o.append(f"      real(dp) :: s({max(NS,1)}, ncb), r({max(NR,1)}, ncb, nck), t({max(NR,1)})")
    o.append(f"      real(dp) :: f(0:{L}), b(0:{L})")
    o.append("      real(dp) :: p, ip, za, zb, eab, q, iq, yc, yd, ax, az, ecd")
    o.append("      real(dp) :: zq, alpha, tt, w, m2a, Rs, Rc, cy")
    decl = []
    if NG:
        decl.append(f"g({NG})")
    if NKM:
        decl.append(f"km({NKM})")
    if NCM:
        decl.append(f"cm({NCM})")
    if f.max_zq_power >= 2:
        decl.extend(f"zq{k}" for k in range(2, f.max_zq_power + 1))
    if decl:
        o.append("      real(dp) :: " + ", ".join(decl))
    o.append("      integer :: kq, bq, cc, ck, n")
    o.append("")
    o.append("      Rs = gc(1); Rc = gc(2); cy = gc(3)")
    o.append("      any = .false.")
    o.append("      r = 0.0_dp")
    o.append("      do kq = 1, nket")
    o.append("         q = kp(1,kq); iq = kp(2,kq); yc = kp(3,kq); yd = kp(4,kq)")
    o.append("         ax = kp(5,kq); az = kp(6,kq); ecd = kp(7,kq)")
    o.append("         s = 0.0_dp")
    o.append("         do bq = 1, nbra")
    o.append("            eab = bp(5,bq)")
    o.append("            if (eab + ecd > cutoff) cycle")
    o.append("            any = .true.")
    o.append("            p = bp(1,bq); ip = bp(2,bq); za = bp(3,bq); zb = bp(4,bq)")
    o.append("            zq = za - az")
    o.append("            alpha = p*q/(p + q)")
    o.append("            tt = alpha*(ax*ax + cy*cy + zq*zq)")
    o.append(f"            call boys(f, tt, {L})")
    o.append("            w = 1.0_dp/sqrt(p + q)")
    o.append("            m2a = -2.0_dp*alpha")
    o.append(f"            do n = 0, {L}")
    o.append("               b(n) = w*f(n)")
    o.append("               w = w*m2a")
    o.append("            end do")
    for k in range(2, f.max_zq_power + 1):
        prev = "zq" if k == 2 else f"zq{k-1}"
        o.append(f"            zq{k} = {prev}*zq")
    # the g monomials of the bra pair
    for gi, g in enumerate(f.g_list):
        e = mono_expr(["ip", "za", "zb"], g)
        o.append(f"            g({gi+1}) = {e if e else '1.0_dp'}")
    # the S values for this primitive quartet, then scattered over columns
    zqpow = {("zq", 1): "zq"}
    for k in range(2, f.max_zq_power + 1):
        zqpow[("zq", k)] = f"zq{k}"
    o.append("            do cc = 1, ncb")
    for si, (g, n, kz, tt) in enumerate(f.s_list):
        parts = []
        gi = f.g_index[g]
        if any(g):
            parts.append(f"g({gi+1})")
        parts.append(f"b({n})")
        if kz:
            parts.append(zqpow[("zq", kz)])
        kabi = "cc" if NTTB == 1 else f"{tt+1}+{NTTB}*(cc-1)"
        o.append(f"               s({si+1},cc) = s({si+1},cc) + kab({kabi},bq)*{'*'.join(parts)}")
    o.append("            end do")
    o.append("         end do")
    # ket level
    o.append("         if (.not. any) cycle")
    for ki, km in enumerate(f.km_list):
        e = mono_expr(["iq", "yc", "yd", "ax"], km)
        o.append(f"         km({ki+1}) = {e if e else '1.0_dp'}")
    o.append("         do cc = 1, ncb")
    for ri, (rterms, _) in enumerate(f.r_list):
        # group by km
        bykm = {}
        for (km_id, s_id), c in rterms:
            bykm.setdefault(km_id, []).append((c, s_id))
        terms = []
        for km_id in sorted(bykm):
            inner = [coef_term(c, f"s({s_id+1},cc)") for c, s_id in bykm[km_id]]
            kme = None if not any(f.km_list[km_id]) else f"km({km_id+1})"
            if len(inner) == 1:
                c, s_id = bykm[km_id][0]
                base = f"s({s_id+1},cc)" if kme is None else f"{kme}*s({s_id+1},cc)"
                terms.append(coef_term(c, base))
            else:
                body = " + ".join(inner).replace("+ -", "- ")
                terms.append(f"({body})" if kme is None else f"{kme}*({body})")
        o.extend(wrap(f"t({ri+1})", terms, indent="            "))
    o.append("            do ck = 1, nck")
    if NTTK == 1:
        o.append(f"               r(1:{NR},cc,ck) = r(1:{NR},cc,ck) + kcd(ck,kq)*t(1:{NR})")
    else:
        # group the accumulators by ket type so each group is one slice
        bytt = {}
        for ri, (_, tt) in enumerate(f.r_list):
            bytt.setdefault(tt, []).append(ri)
        for tt in sorted(bytt):
            ids = bytt[tt]
            # contiguous runs keep the slices short
            runs = []
            for ri in ids:
                if runs and runs[-1][1] == ri - 1:
                    runs[-1][1] = ri
                else:
                    runs.append([ri, ri])
            for a, b in runs:
                o.append(f"               r({a+1}:{b+1},cc,ck) = r({a+1}:{b+1},cc,ck) "
                         f"+ kcd({tt+1}+{NTTK}*(ck-1),kq)*t({a+1}:{b+1})")
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    # assembly
    o.append("      res = 0.0_dp")
    o.append("      if (.not. any) return")
    for ci, cm in enumerate(f.cm_list):
        e = mono_expr(["Rs", "Rc", "cy"], cm)
        o.append(f"      cm({ci+1}) = {e if e else '1.0_dp'}")
    o.append("      do ck = 1, nck")
    o.append("         do cc = 1, ncb")
    for comp_i, terms in enumerate(f.assembly):
        bycm = {}
        for scale, cm_id, r_id in terms:
            bycm.setdefault(cm_id, []).append((scale, r_id))
        parts = []
        for cm_id in sorted(bycm):
            inner = [coef_term(c, f"r({r_id+1},cc,ck)") for c, r_id in bycm[cm_id]]
            cme = None if not any(f.cm_list[cm_id]) else f"cm({cm_id+1})"
            if len(inner) == 1:
                c, r_id = bycm[cm_id][0]
                base = f"r({r_id+1},cc,ck)" if cme is None else f"{cme}*r({r_id+1},cc,ck)"
                parts.append(coef_term(c, base))
            else:
                body = " + ".join(inner).replace("+ -", "- ")
                parts.append(f"({body})" if cme is None else f"{cme}*({body})")
        if not parts:
            parts = ["0.0_dp"]
        o.extend(wrap(f"res({comp_i+1},cc+ncb*(ck-1))", parts, indent="            "))
    o.append("         end do")
    o.append("      end do")
    o.append(f"   end subroutine {name}")
    o.append("")
    return "\n".join(o)


UNROLL_LIMIT = 3000   # ket-level terms; above this a class is emitted table-driven


def data_table(name, values, per_line=16):
    """DATA statements for a module array, chunked so no statement is long.
    DATA rather than an array constructor because gfortran caps constructors
    at 65535 elements and the big d classes pass that."""
    o = []
    n = len(values)
    for i in range(0, n, per_line * 8):
        chunk = values[i:i + per_line * 8]
        lines = []
        for j in range(0, len(chunk), per_line):
            lines.append(", ".join(str(v) for v in chunk[j:j + per_line]))
        o.append(f"   data {name}({i+1}:{i+len(chunk)}) / &\n      " + ", &\n      ".join(lines) + " /")
    return o


def emit_kernel_tabled(f: Factorised):
    """The same three levels as emit_kernel, driven by index tables instead
    of unrolled statements.  Slower per term, but a (dd|dd) class has 96k
    ket-level terms and its unrolled form takes gfortran minutes."""
    d = f.d
    k = d.kinds
    L = d.L
    name = f"rotaxis_{lname(k)}"
    NS, NR, NG = len(f.s_list), len(f.r_list), len(f.g_list)
    NKM, NCM = len(f.km_list), len(f.cm_list)
    NCOMP = len(f.components)
    NTTB, NTTK = f.ntt_bra, f.ntt_ket
    KMAX = f.max_zq_power

    # tables
    s_g = [f.g_index[g] + 1 for g, _, _, _ in f.s_list]
    s_n = [n for _, n, _, _ in f.s_list]
    s_k = [kz for _, _, kz, _ in f.s_list]
    s_tt = [tt + 1 for _, _, _, tt in f.s_list]
    # coefficients are rationals (the ket polynomials are normalised by
    # their leading term), carried as numerator and denominator
    rt_r, rt_km, rt_s, rt_c, rt_d = [], [], [], [], []
    for ri, (terms, _) in enumerate(f.r_list):
        for (km_id, s_id), c in terms:
            c = Fraction(c)
            rt_r.append(ri + 1); rt_km.append(km_id + 1); rt_s.append(s_id + 1)
            rt_c.append(c.numerator); rt_d.append(c.denominator)
    r_tt = [tt + 1 for _, tt in f.r_list]
    at_comp, at_cm, at_r, at_c, at_d = [], [], [], [], []
    for ci, terms in enumerate(f.assembly):
        for scale, cm_id, r_id in terms:
            scale = Fraction(scale)
            at_comp.append(ci + 1); at_cm.append(cm_id + 1); at_r.append(r_id + 1)
            at_c.append(scale.numerator); at_d.append(scale.denominator)
    NT, NA = len(rt_r), len(at_r)
    g_e = [e for g in f.g_list for e in g]
    km_e = [e for m in f.km_list for e in m]
    cm_e = [e for m in f.cm_list for e in m]

    tab = []
    tab.append(f"   integer, save :: s_g({NS}), s_n({NS}), s_k({NS}), s_tt({NS})")
    tab.append(f"   integer, save :: rt_r({NT}), rt_km({NT}), rt_s({NT}), rt_c({NT}), rt_d({NT}), r_tt({NR})")
    tab.append(f"   integer, save :: at_comp({NA}), at_cm({NA}), at_r({NA}), at_c({NA}), at_d({NA})")
    tab.append(f"   integer, save :: g_e({3*NG}), km_e({4*NKM}), cm_e({3*NCM})   ! (3,NG), (4,NKM), (3,NCM) flat")
    for nm, vals in [("s_g", s_g), ("s_n", s_n), ("s_k", s_k), ("s_tt", s_tt),
                     ("rt_r", rt_r), ("rt_km", rt_km), ("rt_s", rt_s), ("rt_c", rt_c), ("rt_d", rt_d),
                     ("r_tt", r_tt),
                     ("at_comp", at_comp), ("at_cm", at_cm), ("at_r", at_r), ("at_c", at_c), ("at_d", at_d),
                     ("g_e", g_e), ("km_e", km_e), ("cm_e", cm_e)]:
        tab.extend(data_table(nm, vals))

    o = []
    o.append(f"   ! ({k[0]}{k[1]}|{k[2]}{k[3]})  L={L}: {NS} bra accumulators, {NR} ket accumulators "
             f"({NT} terms), {NCOMP} components -- table driven")
    o.append(f"   subroutine {name}(nbra, ncb, bp, kab, nket, nck, kp, kcd, gc, cutoff, res, any)")
    o.append("      integer,  intent(in)  :: nbra, ncb, nket, nck")
    o.append(f"      real(dp), intent(in)  :: bp(5, nbra), kab({NTTB}*ncb, nbra)")
    o.append(f"      real(dp), intent(in)  :: kp(7, nket), kcd({NTTK}*nck, nket), gc(3), cutoff")
    o.append(f"      real(dp), intent(out) :: res({NCOMP}, ncb*nck)")
    o.append("      logical,  intent(out) :: any")
    o.append(f"      real(dp) :: s({NS}, ncb), r({NR}, ncb, nck), t({NR})")
    o.append(f"      real(dp) :: f(0:{L}), b(0:{L}), g({NG}), km({NKM}), cm({NCM}), zqp(0:{KMAX}), v")
    o.append("      real(dp) :: p, ip, za, zb, eab, q, iq, yc, yd, ax, az, ecd")
    o.append("      real(dp) :: zq, alpha, tt, w, m2a, Rs, Rc, cy")
    o.append("      integer :: kq, bq, cc, ck, n, i, j, col")
    o.append("")
    o.append("      Rs = gc(1); Rc = gc(2); cy = gc(3)")
    o.append("      any = .false.")
    o.append("      r = 0.0_dp")
    o.append("      do kq = 1, nket")
    o.append("         q = kp(1,kq); iq = kp(2,kq); yc = kp(3,kq); yd = kp(4,kq)")
    o.append("         ax = kp(5,kq); az = kp(6,kq); ecd = kp(7,kq)")
    o.append("         s = 0.0_dp")
    o.append("         do bq = 1, nbra")
    o.append("            eab = bp(5,bq)")
    o.append("            if (eab + ecd > cutoff) cycle")
    o.append("            any = .true.")
    o.append("            p = bp(1,bq); ip = bp(2,bq); za = bp(3,bq); zb = bp(4,bq)")
    o.append("            zq = za - az")
    o.append("            alpha = p*q/(p + q)")
    o.append("            tt = alpha*(ax*ax + cy*cy + zq*zq)")
    o.append(f"            call boys(f, tt, {L})")
    o.append("            w = 1.0_dp/sqrt(p + q)")
    o.append("            m2a = -2.0_dp*alpha")
    o.append(f"            do n = 0, {L}")
    o.append("               b(n) = w*f(n)")
    o.append("               w = w*m2a")
    o.append("            end do")
    o.append("            zqp(0) = 1.0_dp")
    o.append(f"            do i = 1, {KMAX}")
    o.append("               zqp(i) = zqp(i-1)*zq")
    o.append("            end do")
    o.append(f"            do i = 1, {NG}")
    o.append("               g(i) = 1.0_dp")
    o.append("               do j = 1, g_e(3*i-2); g(i) = g(i)*ip; end do")
    o.append("               do j = 1, g_e(3*i-1); g(i) = g(i)*za; end do")
    o.append("               do j = 1, g_e(3*i); g(i) = g(i)*zb; end do")
    o.append("            end do")
    o.append(f"            do i = 1, {NS}")
    o.append("               v = g(s_g(i))*b(s_n(i))*zqp(s_k(i))")
    o.append("               do cc = 1, ncb")
    kabi = "cc" if NTTB == 1 else f"s_tt(i)+{NTTB}*(cc-1)"
    o.append(f"                  s(i,cc) = s(i,cc) + kab({kabi},bq)*v")
    o.append("               end do")
    o.append("            end do")
    o.append("         end do")
    o.append("         if (.not. any) cycle")
    o.append(f"         do i = 1, {NKM}")
    o.append("            km(i) = 1.0_dp")
    o.append("            do j = 1, km_e(4*i-3); km(i) = km(i)*iq; end do")
    o.append("            do j = 1, km_e(4*i-2); km(i) = km(i)*yc; end do")
    o.append("            do j = 1, km_e(4*i-1); km(i) = km(i)*yd; end do")
    o.append("            do j = 1, km_e(4*i); km(i) = km(i)*ax; end do")
    o.append("         end do")
    o.append("         do cc = 1, ncb")
    o.append("            t = 0.0_dp")
    o.append(f"            do i = 1, {NT}")
    o.append("               t(rt_r(i)) = t(rt_r(i)) + real(rt_c(i), dp)/real(rt_d(i), dp)*km(rt_km(i))*s(rt_s(i),cc)")
    o.append("            end do")
    o.append("            do ck = 1, nck")
    if NTTK == 1:
        o.append(f"               r(1:{NR},cc,ck) = r(1:{NR},cc,ck) + kcd(ck,kq)*t(1:{NR})")
    else:
        o.append(f"               do i = 1, {NR}")
        o.append(f"                  r(i,cc,ck) = r(i,cc,ck) + kcd(r_tt(i)+{NTTK}*(ck-1),kq)*t(i)")
        o.append("               end do")
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    o.append("      res = 0.0_dp")
    o.append("      if (.not. any) return")
    o.append(f"      do i = 1, {NCM}")
    o.append("         cm(i) = 1.0_dp")
    o.append("         do j = 1, cm_e(3*i-2); cm(i) = cm(i)*Rs; end do")
    o.append("         do j = 1, cm_e(3*i-1); cm(i) = cm(i)*Rc; end do")
    o.append("         do j = 1, cm_e(3*i); cm(i) = cm(i)*cy; end do")
    o.append("      end do")
    o.append("      do ck = 1, nck")
    o.append("         do cc = 1, ncb")
    o.append("            col = cc + ncb*(ck-1)")
    o.append(f"            do i = 1, {NA}")
    o.append("               res(at_comp(i),col) = res(at_comp(i),col) &")
    o.append("                  + real(at_c(i), dp)/real(at_d(i), dp)*cm(at_cm(i))*r(at_r(i),cc,ck)")
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    o.append(f"   end subroutine {name}")
    o.append("")
    return "\n".join(tab), "\n".join(o)


BOYS_MODULE = """! GENERATED by scripts/rotaxis_mmd -- do not edit.
!
! The Boys function F_0..F_m for the rotated-axis kernels.  gamma_inc_like
! takes its closed form above turnover_point(m), and that table is 0 for
! m = 0 and 1 -- the C never asks it for F_1 at a tiny argument, the Rys
! roots come from elsewhere -- so at t ~ 1e-33, which a primitive quartet
! with two equal pairs produces by rounding, the downward recursion
! (F_0 - exp(-t))/(2t) cancels to garbage.  The ascending series is stable
! there and cheap, so it is used below 1 as well.
module cint_rotaxis_boys
   use cint_const, only: dp
   use cint_fmt_tab, only: turnover_point
   use cint_fmt_dp, only: gamma_inc_like, fmt1_gamma_inc_like
   implicit none
   private
   public :: boys
   ! sqrt(pi)/2, the large-argument limit of sqrt(pi/(4t)) erf(sqrt(t)) t^(1/2)
   real(dp), parameter :: SQRT_PI_2 = 0.88622692545275801365_dp
contains
   subroutine boys(f, t, m)
      real(dp), intent(out) :: f(0:)
      real(dp), intent(in)  :: t
      integer,  intent(in)  :: m
      real(dp) :: b
      integer  :: i
      if (t > 50.0_dp) then
         ! THE FAR FIELD, and worth its own branch.  Above t = 50 the closed
         ! form's exp(-t) is under 2e-22 and its erf(sqrt(t)) is 1 to better
         ! than that, so both libm calls drop out and F_m is
         ! sqrt(pi)/(2 sqrt(t)) run up by (2i-1)/(2t).  Measured inside a
         ! Fock build, libm's erf and exp were 12 s of 101 s on this path
         ! against a Rys path that calls neither -- and a screened Fock
         ! build spends much of its time exactly here, on distant pairs.
         f(0) = SQRT_PI_2/sqrt(t)
         b = 0.5_dp/t
         do i = 1, m
            f(i) = b*real(2*i - 1, dp)*f(i-1)
         end do
      else if (t < max(turnover_point(m), 1.0_dp)) then
         call fmt1_gamma_inc_like(f, t, m)
      else
         call gamma_inc_like(f, t, m)
      end if
   end subroutine boys
end module cint_rotaxis_boys
"""


def emit_class_file(f, unroll_limit=UNROLL_LIMIT):
    """One module per class, so the classes compile in parallel and the big
    d kernels do not hold up the rest.  Small classes are unrolled; above
    `unroll_limit` ket-level terms the kernel is table driven."""
    name = f"rotaxis_{lname(f.d.kinds)}"
    nterms = sum(len(r) for r, _ in f.r_list)
    if nterms > unroll_limit:
        tables, body = emit_kernel_tabled(f)
    else:
        tables, body = "", emit_kernel(f)
    o = ["! GENERATED by scripts/rotaxis_mmd -- do not edit.",
         f"! {f.summary()}",
         f"module {name}_m",
         "   use cint_const, only: dp",
         "   use cint_rotaxis_boys, only: boys",
         "   implicit none",
         "   private",
         f"   public :: {name}"]
    if tables:
        o.append(tables)
    o += ["contains", body, f"end module {name}_m"]
    return "\n".join(o) + "\n"


def emit_files(classes, unroll_limit=UNROLL_LIMIT):
    """Return {relative path: text} for everything the generator writes."""
    from .derive import KIND_LMAX
    facts = [Factorised(Derivation(*k)) for k in classes]
    files = {"src/rotaxis/cint_rotaxis_boys.f90": BOYS_MODULE}
    for f in facts:
        files[f"src/rotaxis/rotaxis_{lname(f.d.kinds)}.f90"] = emit_class_file(f, unroll_limit)

    o = []
    o.append("! GENERATED by scripts/rotaxis_mmd -- do not edit.")
    o.append("!")
    o.append("! Rotated-axis McMurchie-Davidson kernels for the two-electron Coulomb")
    o.append("! integral: the dispatcher.  One module per canonical class of shell")
    o.append("! KINDS -- s, p, L (s and p on shared exponents, one four-component")
    o.append("! slot), d -- ordered by rank s < p < L < d within each pair and between")
    o.append("! the pairs, in src/rotaxis/.  The driver in cint_rotaxis_2e.f90 permutes")
    o.append("! any quartet into that order, builds the frame and the pair tables, and")
    o.append("! rotates the result back.  doc/ROT_AXIS_MMD.md derives the formulas;")
    o.append("! scripts/rotaxis_mmd/derive.py is the same derivation, executable.")
    o.append("!")
    for f in facts:
        o.append(f"!   {f.summary()}")
    o.append("module cint_rotaxis_kernels")
    o.append("   use cint_const, only: dp")
    for f in facts:
        n = f"rotaxis_{lname(f.d.kinds)}"
        o.append(f"   use {n}_m, only: {n}")
    o.append("   implicit none")
    o.append("   private")
    o.append("   public :: rotaxis_kernel, rotaxis_has_class")
    lmax = max(KIND_LMAX[x] for k in classes for x in k)
    o.append(f"   integer, parameter, public :: ROTAXIS_LMAX = {lmax}")
    o.append("   ! Kind ranks, the digits of the class code: s, p, L, d.")
    o.append("   integer, parameter, public :: ROTAXIS_KIND_S = 0, ROTAXIS_KIND_P = 1, &")
    o.append("                                 ROTAXIS_KIND_L = 2, ROTAXIS_KIND_D = 3")
    o.append("")
    o.append("contains")
    o.append("")
    o.append("   ! `code` is the four kind ranks as decimal digits, canonical order.")
    o.append("   pure logical function rotaxis_has_class(code) result(yes)")
    o.append("      integer, intent(in) :: code")
    o.append("      select case (code)")
    codes = [str(class_code(k)) for k in classes]
    o.append("      case (" + ", &\n            ".join(", ".join(codes[i:i+12]) for i in range(0, len(codes), 12)) + ")")
    o.append("         yes = .true.")
    o.append("      case default")
    o.append("         yes = .false.")
    o.append("      end select")
    o.append("   end function rotaxis_has_class")
    o.append("")
    o.append("   ! Dispatch on the canonical class.  `res` is (ncomp, ncb*nck); kab and")
    o.append("   ! kcd are (ntt*ncb, nbra) and (ntt*nck, nket), see the kernel head.")
    o.append("   subroutine rotaxis_kernel(code, nbra, ncb, bp, kab, nket, nck, kp, kcd, &")
    o.append("                             gc, cutoff, res, any)")
    o.append("      integer,  intent(in)  :: code, nbra, ncb, nket, nck")
    o.append("      real(dp), intent(in)  :: bp(5, nbra), kab(*)")
    o.append("      real(dp), intent(in)  :: kp(7, nket), kcd(*), gc(3), cutoff")
    o.append("      real(dp), intent(out) :: res(*)")
    o.append("      logical,  intent(out) :: any")
    o.append("      select case (code)")
    for k in classes:
        o.append(f"      case ({class_code(k)})")
        o.append(f"         call rotaxis_{lname(k)}(nbra, ncb, bp, kab, nket, nck, kp, kcd, gc, cutoff, res, any)")
    o.append("      case default")
    o.append("         error stop 'cint_rotaxis_kernels: no kernel for this class'")
    o.append("      end select")
    o.append("   end subroutine rotaxis_kernel")
    o.append("end module cint_rotaxis_kernels")
    files["src/cint_rotaxis_kernels.f90"] = "\n".join(o) + "\n"

    cm = ["# GENERATED by scripts/rotaxis_mmd -- do not edit.",
          "# The rotated-axis kernel sources, one per class, plus the shared Boys",
          "# wrapper; included from src/CMakeLists.txt.",
          "set(cintRotaxisSrc"]
    cm += [f"  rotaxis/{path.split('/')[-1]}" for path in files if path.startswith("src/rotaxis/")]
    cm.append(")")
    files["src/rotaxis/kernels.cmake"] = "\n".join(cm) + "\n"
    return files, facts
