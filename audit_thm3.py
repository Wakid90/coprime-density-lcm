"""
Audit of Theorem 3: Sharp form at fixed M.

Claims:
  (C1) q**(M) < q*(M) strictly for M >= 3.
  (C2) R_k(M) <= (q**)^k * exp(Lambda(M)).
  (C3) Sigma_k(M) >= (q*)^k / (M-1).
  (C4) For M >= 11, q** = (M-3)/(M-1), q* = (M-2)/(M-1), ratio = (M-3)/(M-2).
  (C5) Therefore R_k/Sigma_k <= (M-1) exp(Lambda) * (q**/q*)^k, decays exponentially.
"""
import math
from itertools import combinations
from sympy import Rational, primerange

def primes_up_to(n): return list(primerange(2, n+1))
def q_p_rat(p, M): return 1 - Rational(M // p, M - 1)
def beta_T_rat(T, M):
    numerator = 0
    for r in range(len(T)+1):
        for S in combinations(T, r):
            prod = 1
            for p in S: prod *= p
            if prod > M:
                if r == 0: numerator += M
                continue
            numerator += ((-1)**r) * (M // prod)
    return Rational(numerator - 1, M - 1)
def Sigma_k_rat(k, M):
    return sum(q_p_rat(p, M)**k * Rational(1, p-1) for p in primes_up_to(M))
def R_k_rat(k, M):
    primes = primes_up_to(M)
    R = Rational(0)
    for r in range(2, len(primes)+1):
        for T in combinations(primes, r):
            b = beta_T_rat(T, M)
            w = Rational(1)
            for p in T: w *= Rational(1, p-1)
            R += b**k * w
    return R
def Lambda_M(M):
    return sum(Rational(1, p-1) for p in primes_up_to(M))

print("="*90)
print("AUDIT: Theorem 3 — Sharp form at fixed M (the new theorem)")
print("="*90)

# =====================================================================
# (C1): q** < q* strictly for M >= 3. Exhaustive up to M = 50.
# =====================================================================
print("(C1) q**(M) < q*(M) strictly for all M in [3, 50]. EXHAUSTIVE CHECK.")
violations_c1 = []
for M in range(3, 51):
    primes = primes_up_to(M)
    q_star = max(q_p_rat(p, M) for p in primes) if primes else Rational(0)
    q_star_star = Rational(0)
    if len(primes) >= 2:
        for r in range(2, len(primes)+1):
            for T in combinations(primes, r):
                b = beta_T_rat(T, M)
                if b > q_star_star:
                    q_star_star = b
    if q_star_star >= q_star and len(primes) >= 2:
        violations_c1.append((M, q_star, q_star_star))

if violations_c1:
    print(f"  !!!! VIOLATIONS !!!! ")
    for M, q_s, q_ss in violations_c1:
        print(f"    M={M}: q*={q_s}, q**={q_ss}")
else:
    print(f"  PASS: q** < q* strictly at every M in [3, 50].")

print()

# =====================================================================
# (C4): q** = (M-3)/(M-1) for M >= 11.
# =====================================================================
print("(C4) For M >= 11, q**(M) = (M-3)/(M-1). EXHAUSTIVE up to M = 60.")
violations_c4 = []
for M in range(11, 61):
    primes = primes_up_to(M)
    q_star_star = Rational(0)
    for r in range(2, len(primes)+1):
        for T in combinations(primes, r):
            b = beta_T_rat(T, M)
            if b > q_star_star:
                q_star_star = b
    predicted = Rational(M-3, M-1)
    if q_star_star != predicted:
        violations_c4.append((M, q_star_star, predicted))

if violations_c4:
    print(f"  !!!! VIOLATIONS !!!! ")
    for M, qss, pred in violations_c4[:5]:
        print(f"    M={M}: q**={qss}, predicted={pred}")
else:
    print(f"  PASS: q**(M) = (M-3)/(M-1) exactly for all M in [11, 60].")

print()

# =====================================================================
# (C2) and (C3): the individual bounds in the theorem.
# =====================================================================
print("(C2) R_k(M) <= (q**)^k * exp(Lambda) EXACT check.")
print(f"  {'M':>4} {'k':>4} {'R_k':>14} {'(q**)^k exp(Lam)':>20} {'ratio':>8} {'OK?':>5}")
for M in [12, 20, 30, 40, 50]:
    primes = primes_up_to(M)
    q_star_star = Rational(0)
    for r in range(2, len(primes)+1):
        for T in combinations(primes, r):
            b = beta_T_rat(T, M)
            if b > q_star_star: q_star_star = b
    Lam = float(Lambda_M(M))
    for k in [10, 50, 100]:
        Rk = float(R_k_rat(k, M))
        bnd = (float(q_star_star))**k * math.exp(Lam)
        ratio = Rk/bnd if bnd > 0 else 0
        ok = "Y" if Rk <= bnd + 1e-300 else "NO"
        print(f"  {M:>4} {k:>4} {Rk:>14.4e} {bnd:>20.4e} {ratio:>8.4f} {ok:>5}")

print()

print("(C3) Sigma_k(M) >= (q*)^k / (M-1) EXACT check.")
print(f"  {'M':>4} {'k':>4} {'Sigma_k':>14} {'(q*)^k/(M-1)':>18} {'ratio':>8} {'OK?':>5}")
for M in [12, 20, 30, 40, 50]:
    primes = primes_up_to(M)
    q_star = max(q_p_rat(p, M) for p in primes)
    for k in [10, 50, 100]:
        Sk = float(Sigma_k_rat(k, M))
        lo = float(q_star**k) / (M-1)
        ok = "Y" if Sk >= lo - 1e-300 else "NO"
        ratio = Sk/lo if lo > 0 else 0
        print(f"  {M:>4} {k:>4} {Sk:>14.4e} {lo:>18.4e} {ratio:>8.4f} {ok:>5}")

print()

# =====================================================================
# (C5): The composite ratio R_k/Sigma_k and the theoretical rate.
# =====================================================================
print("(C5) R_k/Sigma_k vs. theoretical exponential rate. Push k much higher.")
for M in [12, 20, 30, 40, 50]:
    primes = primes_up_to(M)
    q_star = max(q_p_rat(p, M) for p in primes)
    q_star_star = Rational(0)
    for r in range(2, len(primes)+1):
        for T in combinations(primes, r):
            b = beta_T_rat(T, M)
            if b > q_star_star: q_star_star = b
    rate = float(q_star_star / q_star)
    Lam = float(Lambda_M(M))
    prefactor = (M-1) * math.exp(Lam)
    print(f"\n  M={M}, q**/q*={rate:.5f}, prefactor (M-1)exp(Lam) = {prefactor:.2f}")
    print(f"  {'k':>4} {'R_k/Sigma_k':>14} {'bound':>14} {'ratio':>8}")
    for k in [20, 50, 100, 200, 500, 1000]:
        Rk = float(R_k_rat(k, M))
        Sk = float(Sigma_k_rat(k, M))
        obs = Rk/Sk if Sk > 0 else 0
        bnd = prefactor * rate**k
        r = obs/bnd if bnd > 0 else 0
        # The bound should hold: obs <= bnd
        ok = "Y" if obs <= bnd else "NO"
        print(f"  {k:>4} {obs:>14.4e} {bnd:>14.4e} {r:>8.4f} {ok}")

print()
print("VERDICT: Theorem 3 — VALID if all checks pass above.")
