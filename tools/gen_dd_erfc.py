"""Emit double-double erfc coefficients for cint_dd.

Method and band structure from bitwise_adventures/tools/gen_erfc.py; the
changes for dd are more terms, the C/D boundary moved from 7 to 9 (an
asymptotic series has an error floor of its smallest term, 7.4e-22 at x=7 and
9.6e-36 at x=9), and a range-reduced cosine in the DCT.

Coefficients are split against the EXACT binary value of `hi`, never against
its printed form -- the difference is ~2e-17 and it is invisible until the
constant enters the arithmetic.
"""
from decimal import Decimal as D, getcontext
getcontext().prec = 120
PI = D('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651')
TWOPI = 2*PI
SQRTPI = PI.sqrt(); TWO_OVER_SQRTPI = 2/SQRTPI

def dcos(x):
    x = x - int(x/TWOPI)*TWOPI          # range reduction: see module docstring
    s, term, n = D(1), D(1), 0
    while abs(term) > D('1e-115'):
        n += 2; term *= -x*x/(n*(n-1)); s += term
    return s

def derfc(x):
    s, term, n = D(0), x, 0
    while True:
        s += term/(2*n+1); n += 1; term *= -x*x/n
        if abs(term) < D('1e-110') and n > 5: break
        if n > 2000: raise RuntimeError('no convergence')
    return 1 - TWO_OVER_SQRTPI*s

def derfcx(x): return derfc(x)*(x*x).exp()

def cheb(f, a, b, N, keep):
    a, b = D(a), D(b); nodes, fv = [], []
    for k in range(N):
        th = PI*(2*D(k)+1)/(2*N)
        nodes.append(th); fv.append(f((a+b)/2 + (b-a)/2*dcos(th)))
    out = []
    for j in range(keep):
        s = sum((fv[k]*dcos(j*nodes[k]) for k in range(N)), D(0))
        c = 2*s/N
        out.append(c/2 if j == 0 else c)
    return out

def split(d):
    hi = float(d)
    lo = float(d - D(hi))               # EXACT binary value of hi, not repr(hi)
    return hi, lo

def emit(name, vals, lo_index=0):
    ok = True
    for v in vals:
        hi, lo = split(v)
        if v != 0 and abs(v - (D(hi)+D(lo)))/abs(v) > D('1e-31'): ok = False
    print(f"   ! {name}: {len(vals)} terms, every pair verified to <1e-31")
    if not ok: print("   ! *** A PAIR FAILED VERIFICATION ***")
    print(f"   type(dd), parameter :: {name}({lo_index}:{lo_index+len(vals)-1}) = [ &")
    for i, v in enumerate(vals):
        hi, lo = split(v)
        tail = ", &" if i < len(vals)-1 else " ]"
        print(f"      dd({hi!r}_dp, {lo!r}_dp){tail}")

# band A: erf Maclaurin on |x| <= 0.46875, erf = x*P(x^2)
mac, fact = [], D(1)
for n in range(0, 40):
    if n: fact *= n
    t = TWO_OVER_SQRTPI * (-1)**n / (fact*(2*n+1))
    mac.append(t)
    if abs(t) * D('0.46875')**(2*n) < D('1e-36') and n > 4: break
emit('ERF_MAC', mac)
print()
emit('CXB', cheb(derfcx,'0.46875','2',220,35))
print()
emit('CXC', cheb(derfcx,'2','9',220,55))
print()
# band D: asymptotic in s = 1/x^2, x > 9
asy, df = [D(1)], D(1)
for n in range(1, 200):
    df *= (2*n-1)
    t = (-1)**n * df / D(2)**n
    asy.append(t)
    if abs(t) / D(81)**n < D('1e-36'): break
emit('ASY_D', asy)
