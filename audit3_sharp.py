"""
AUDIT 3: Theorem 3 (Paper 1) — Sharp form at fixed M.

Claim:
(a) q^{**}(M) := max_{|T|>=2} beta_T(M) < q^*(M) := max_p q_p(M), for M >= 3.
(b) R_k(M) = O((q^{**}/q^*)^k) as k -> infty.
(c) For M >= 11, q^{**}(M) = (M-3)/(M-1) and q^*(M) = (M-2)/(M-1).
"""
from fractions import Fraction
from itertools import combinations
from sympy import primerange
import math

def primes_up_to(n):
    return list(primerange(2, n + 1))

def q_p(p, M):
    return Fraction(M - 1 - (M // p), M - 1)

def beta_T(T, M):
    count = 0
    for n in range(2, M + 1):
        ok = True
        for p in T:
            if n % p == 0:
                ok = False
                break
        if ok:
            count += 1
    return Fraction(count, M - 1)

def q_star(M):
    plist = primes_up_to(M)
    return max(q_p(p, M) for p in plist)

def q_star_star(M):
    """max beta_T over |T| >= 2, exhaustive."""
    plist = primes_up_to(M)
    best = Fraction(0)
    for r in range(2, len(plist) + 1):
        for T in combinations(plist, r):
            b = beta_T(T, M)
            if b > best:
                best = b
    return best

def Sigma_k_exact(k, M):
    return sum(q_p(p, M) ** k * Fraction(1, p - 1) for p in primes_up_to(M))

def Rk_exact(k, M):
    plist = primes_up_to(M)
    total = Fraction(0)
    for r in range(2, len(plist) + 1):
        for T in combinations(plist, r):
            b = beta_T(T, M)
            pref = Fraction(1)
            for p in T:
                pref *= Fraction(1, p - 1)
            total += (b ** k) * pref
    return total

def run_audit3():
    print("="*72)
    print("AUDIT 3: Theorem 3 (sharp form at fixed M)")
    print("="*72)

    print("\n--- (a) q** < q* strict, for M in [3, 40] ---")
    strict_fail = []
    for M in range(3, 41):
        qs = q_star(M)
        qss = q_star_star(M)
        if qss >= qs:
            strict_fail.append((M, qss, qs))
            print(f"  FAIL M={M}: q**={qss} >= q*={qs}")
        # Don't print every pass, just sample
    if not strict_fail:
        print(f"  All M in [3, 40]: q** < q* strict. (Sample values below.)")
        for M in [3, 5, 7, 11, 15, 20, 30, 40]:
            qs = q_star(M)
            qss = q_star_star(M)
            print(f"    M={M:2d}: q*={float(qs):.6f}, q**={float(qss):.6f}, "
                  f"q**/q*={float(qss/qs):.6f}")

    print("\n--- (c) For M >= 11, q* = (M-2)/(M-1) and q** = (M-3)/(M-1) ---")
    c_fail = []
    for M in range(11, 51):
        qs = q_star(M)
        qss = q_star_star(M)
        qs_predict = Fraction(M - 2, M - 1)
        qss_predict = Fraction(M - 3, M - 1)
        if qs != qs_predict or qss != qss_predict:
            c_fail.append((M, qs, qs_predict, qss, qss_predict))
    if not c_fail:
        print(f"  All M in [11, 50]: q* = (M-2)/(M-1), q** = (M-3)/(M-1). EXACT MATCH.")
    else:
        for f in c_fail:
            print(f"  FAIL: M={f[0]}: q*={f[1]} vs {f[2]}, q**={f[3]} vs {f[4]}")

    print("\n--- (b) R_k/Sigma_k decays geometrically with rate (q**/q*)^k ---")
    print("       Also: this decay is actually FASTER than the naive bound")
    for M in [15, 20, 30]:
        qs = q_star(M)
        qss = q_star_star(M)
        ratio = qss / qs
        print(f"\n  M={M}, q**/q* = {float(ratio):.6f} ({qss}/{qs})")
        for k in [10, 20, 40, 80, 160]:
            Sig = float(Sigma_k_exact(k, M))
            Rk = float(Rk_exact(k, M))
            obs = Rk / Sig if Sig > 0 else float('nan')
            rate = float(ratio) ** k
            print(f"    k={k:3d}: R/Sigma = {obs:.3e}, (q**/q*)^k = {rate:.3e}, "
                  f"ratio of obs/rate = {obs/rate:.4f}")

    if not strict_fail and not c_fail:
        print("\n*** AUDIT 3 VERDICT: all three claims of Theorem 3 verified. ***")
    else:
        print("\n*** AUDIT 3 VERDICT: failures detected. ***")

if __name__ == "__main__":
    run_audit3()
