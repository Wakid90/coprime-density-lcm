"""
Verify that q** := max_{|T|>=2} beta_T(M) is strictly less than q* := max_p q_p(M),
for a range of M.

Also compute the ratio (q**/q*)^k to see the exponential decay rate.
"""
import math
from itertools import combinations
from sympy import primerange

def primes_up_to(n): return list(primerange(2, n+1))
def q_p(p, M): return 1.0 - (M // p) / (M - 1)
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

# Compute q*, q**, and their ratio for M in range
print(f"{'M':>4} {'q*':>10} {'q**':>10} {'q**/q*':>10} {'T achieving q**':>30}")
for M in range(5, 51):
    primes = primes_up_to(M)
    q_star = max(q_p(p, M) for p in primes)
    
    # Compute beta_T for all T with |T|>=2, find max
    q_star_star = 0.0
    best_T = None
    for r in range(2, len(primes) + 1):
        for T in combinations(primes, r):
            b = beta_T(T, M)
            if b > q_star_star:
                q_star_star = b
                best_T = T
    
    ratio = q_star_star / q_star if q_star > 0 else 0
    T_str = str(best_T) if best_T and len(best_T) <= 5 else f"|T|={len(best_T) if best_T else 0}"
    print(f"{M:>4} {q_star:>10.6f} {q_star_star:>10.6f} {ratio:>10.6f} {T_str:>30}")
