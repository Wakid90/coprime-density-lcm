"""
Test the sharp form at fixed M: does R_k(M) / Sigma_k(M) -> 0 as k -> infty?

Compute exactly for M = 10, 15, 20 for a range of k.
"""
import math
from itertools import combinations
from sympy import primerange

def primes_up_to(n): return list(primerange(2, n+1))
def q_p(p, M): return 1.0 - (M // p) / (M - 1)

def Sigma_k(k, M):
    return sum(q_p(p, M)**k / (p-1) for p in primes_up_to(M))

def beta_T(T, M):
    total = 0
    for r in range(len(T) + 1):
        for S in combinations(T, r):
            prod = 1; stop = False
            for p in S:
                prod *= p
                if prod > M: stop = True; break
            if stop and r >= 1: continue
            total += ((-1)**r) * (M // prod)
    return (total - 1) / (M - 1)

def rho_over_PM(k, M):
    primes = primes_up_to(M)
    total = 0.0
    for m in range(len(primes) + 1):
        for T in combinations(primes, m):
            b = beta_T(T, M)
            w = 1.0
            for p in T: w /= (p-1)
            total += b**k * w
    return total

def R_k_exact(k, M):
    return rho_over_PM(k, M) - 1 - Sigma_k(k, M)


# Also find q_max and the primes attaining it
def analyze_M(M):
    primes = primes_up_to(M)
    qs = [(q_p(p, M), p) for p in primes]
    qs.sort(reverse=True)
    return qs

# Test
for M in [10, 15, 20, 30, 40]:
    print(f"\n=== M = {M} ===")
    qs = analyze_M(M)
    print(f"Top 5 q_p values: {[(f'{q:.4f}', p) for q, p in qs[:5]]}")
    # Identify "top stratum" (all primes with q close to max)
    q_max = qs[0][0]
    top = [p for q, p in qs if q == q_max]
    print(f"Primes achieving q_max = {q_max:.4f}: {top}")
    
    print(f"{'k':>6} {'Sigma_k':>12} {'R_k':>12} {'R_k/Sigma_k':>14} {'Sigma_(k/2)^2/Sigma_k':>23}")
    for k in [20, 40, 60, 80, 120, 160, 200, 300, 500]:
        Sk = Sigma_k(k, M)
        Rk = R_k_exact(k, M)
        Sk2 = Sigma_k(k//2, M)
        ratio1 = Rk / Sk if Sk > 0 else 0
        ratio2 = (Sk2**2)/Sk if Sk > 0 else 0
        print(f"{k:>6} {Sk:>12.4e} {Rk:>12.4e} {ratio1:>14.6f} {ratio2:>23.6f}")
