"""Evaluate the recurrence graph numerically and compare against an
independent McMurchie-Davidson reference.  Catches an algebra error in the
recurrences before any of it reaches Fortran."""

import math
import random

from .cart import cart_components
from .recur import build


def boys(m, t):
    if t < 1.0e-14:
        return 1.0 / (2 * m + 1)
    if t < 25.0:
        term = 1.0 / (2 * m + 1)
        s = term
        k = 1
        while True:
            term *= (2 * t) / (2 * m + 2 * k + 1)
            s += term
            if term < 1e-18 * s:
                break
            k += 1
        return math.exp(-t) * s
    # downward from the asymptotic F_0
    f = 0.5 * math.sqrt(math.pi / t) * math.erf(math.sqrt(t))
    e = math.exp(-t)
    for i in range(1, m + 1):
        f = ((2 * i - 1) * f - e) / (2 * t)
    return f


def evaluate(g, order, geom):
    a, A, b, B, c, C, d, D, L = geom
    p, q = a + b, c + d
    P = [(a * A[i] + b * B[i]) / p for i in range(3)]
    Q = [(c * C[i] + d * D[i]) / q for i in range(3)]
    pq = p + q
    rho = p * q / pq
    W = [(p * P[i] + q * Q[i]) / pq for i in range(3)]
    sym = {}
    for i in range(3):
        sym[f'PA{i}'] = P[i] - A[i]
        sym[f'WP{i}'] = W[i] - P[i]
        sym[f'QC{i}'] = Q[i] - C[i]
        sym[f'WQ{i}'] = W[i] - Q[i]
        sym[f'AB{i}'] = A[i] - B[i]
        sym[f'CD{i}'] = C[i] - D[i]
    sym['oo2p'] = 0.5 / p
    sym['oo2q'] = 0.5 / q
    sym['oo2pq'] = 0.5 / pq
    sym['rp'] = rho / p
    sym['rq'] = rho / q
    rab2 = sum((A[i] - B[i]) ** 2 for i in range(3))
    rcd2 = sum((C[i] - D[i]) ** 2 for i in range(3))
    t = rho * sum((P[i] - Q[i]) ** 2 for i in range(3))
    pref = (2.0 * math.pi ** 2.5 / (p * q * math.sqrt(pq))
            * math.exp(-a * b / p * rab2) * math.exp(-c * d / q * rcd2))
    val = {}
    for k in order:
        e = g.expr.get(k)
        if e is None:
            val[k] = pref * boys(k[3], t)          # a base ('v', 0, 0, m)
        else:
            acc = 0.0
            for coef, syms, dep in e:
                x = float(coef)
                for s in syms:
                    x *= sym[s]
                acc += x * val[dep]
            val[k] = acc
    return val


# ---- the reference: McMurchie-Davidson in the lab frame -----------------

def e_coef(i, j, tt, xpa, xpb, half, kab, cache):
    key = (i, j, tt)
    if key in cache:
        return cache[key]
    if tt < 0 or tt > i + j:
        v = 0.0
    elif i == 0 and j == 0:
        v = kab
    elif i > 0:
        v = (half * e_coef(i - 1, j, tt - 1, xpa, xpb, half, kab, cache)
             + xpa * e_coef(i - 1, j, tt, xpa, xpb, half, kab, cache)
             + (tt + 1) * e_coef(i - 1, j, tt + 1, xpa, xpb, half, kab, cache))
    else:
        v = (half * e_coef(i, j - 1, tt - 1, xpa, xpb, half, kab, cache)
             + xpb * e_coef(i, j - 1, tt, xpa, xpb, half, kab, cache)
             + (tt + 1) * e_coef(i, j - 1, tt + 1, xpa, xpb, half, kab, cache))
    cache[key] = v
    return v


def r_int(n, t, u, v, alpha, R, cache):
    key = (n, t, u, v)
    if key in cache:
        return cache[key]
    if t == u == v == 0:
        val = (-2.0 * alpha) ** n * boys(n, alpha * sum(x * x for x in R))
    elif t > 0:
        val = R[0] * r_int(n + 1, t - 1, u, v, alpha, R, cache)
        if t > 1:
            val += (t - 1) * r_int(n + 1, t - 2, u, v, alpha, R, cache)
    elif u > 0:
        val = R[1] * r_int(n + 1, t, u - 1, v, alpha, R, cache)
        if u > 1:
            val += (u - 1) * r_int(n + 1, t, u - 2, v, alpha, R, cache)
    else:
        val = R[2] * r_int(n + 1, t, u, v - 1, alpha, R, cache)
        if v > 1:
            val += (v - 1) * r_int(n + 1, t, u, v - 2, alpha, R, cache)
    cache[key] = val
    return val


def reference(ca, cb, cc, cd, geom):
    a, A, b, B, c, C, d, D, _ = geom
    p, q = a + b, c + d
    P = [(a * A[i] + b * B[i]) / p for i in range(3)]
    Q = [(c * C[i] + d * D[i]) / q for i in range(3)]
    alpha = p * q / (p + q)
    R = [P[i] - Q[i] for i in range(3)]
    eb, ek, rc = [], [], {}
    for i in range(3):
        cb_ = {}
        kab = math.exp(-a * b / p * (A[i] - B[i]) ** 2)
        eb.append([e_coef(ca[i], cb[i], t, P[i] - A[i], P[i] - B[i],
                          0.5 / p, kab, cb_) for t in range(ca[i] + cb[i] + 1)])
        ck_ = {}
        kcd = math.exp(-c * d / q * (C[i] - D[i]) ** 2)
        ek.append([e_coef(cc[i], cd[i], t, Q[i] - C[i], Q[i] - D[i],
                          0.5 / q, kcd, ck_) for t in range(cc[i] + cd[i] + 1)])
    total = 0.0
    for t1, x1 in enumerate(eb[0]):
        for u1, y1 in enumerate(eb[1]):
            for v1, z1 in enumerate(eb[2]):
                for t2, x2 in enumerate(ek[0]):
                    for u2, y2 in enumerate(ek[1]):
                        for v2, z2 in enumerate(ek[2]):
                            sgn = (-1.0) ** (t2 + u2 + v2)
                            total += (x1 * y1 * z1 * x2 * y2 * z2 * sgn
                                      * r_int(0, t1 + t2, u1 + u2, v1 + v2,
                                              alpha, R, rc))
    return total * 2.0 * math.pi ** 2.5 / (p * q * math.sqrt(p + q))


def check(classes, seed=7):
    rng = random.Random(seed)
    worst, ncmp = 0.0, 0
    for cls in classes:
        g, vroots, targets = build(*cls)
        order = g.order([k for _, k in targets])
        for _ in range(2):
            geom = (rng.uniform(0.3, 3.0), [rng.uniform(-1.5, 1.5) for _ in range(3)],
                    rng.uniform(0.3, 3.0), [rng.uniform(-1.5, 1.5) for _ in range(3)],
                    rng.uniform(0.3, 3.0), [rng.uniform(-1.5, 1.5) for _ in range(3)],
                    rng.uniform(0.3, 3.0), [rng.uniform(-1.5, 1.5) for _ in range(3)],
                    sum(cls))
            val = evaluate(g, order, geom)
            scale = max(abs(val[k]) for _, k in targets)
            for comps, k in targets:
                ref = reference(*comps, geom)
                rel = abs(val[k] - ref) / max(scale, 1e-30)
                worst = max(worst, rel)
                ncmp += 1
        print(f"  {cls}  worst so far {worst:.2e}  ({ncmp} values)")
    return worst
