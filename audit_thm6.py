"""
Audit of Theorem 6: M ln M * Var[W_k(M)] -> B(c).

Previous numerics (v4.3, c=1, M=1000): 0.428 vs predicted 0.4255. Match.
Previous numerics (v4.3, c=2, M=1000): 0.138 vs predicted 0.1379. Match.

Push further: c in {0.5, 1, 2, 5}, M in {500, 1000, 2000, 5000, 10000}.
"""
import math
from itertools import combinations
from sympy import primerange
import time

def primes_up_to(n): return list(primerange(2, n+1))
def q_p(p, M): return 1.0 - (M // p) / (M - 1)

def V1(k, M):
    """Diagonal part. O(pi(M)) time."""
    s = 0.0
    for p in primes_up_to(M):
        qp = q_p(p, M); lp = math.log(1 - 1/p)
        s += lp*lp * qp**k * (1 - qp**k)
    return s

def V2(k, M):
    """Off-diagonal part. O(pi(M)^2) time."""
    s = 0.0
    primes = primes_up_to(M)
    for i, p in enumerate(primes):
        qp = q_p(p, M); lp = math.log(1-1/p)
        for q in primes[i+1:]:
            qq = q_p(q, M); lq = math.log(1-1/q)
            mpq = M // (p*q) if p*q <= M else 0
            bpq = (M - (M//p) - (M//q) + mpq - 1) / (M-1)
            s += 2 * lp * lq * (bpq**k - qp**k * qq**k)
    return s

def B_c(c):
    return 1/(math.exp(c)-1) - 1/(math.exp(2*c)-1)

print("="*80)
print("AUDIT: Theorem 6 — Variance asymptotic M ln M * Var[W_k] -> B(c)")
print("="*80)
print()

results = {}
for M in [500, 1000, 2000, 5000]:
    t0 = time.time()
    results[M] = {}
    for c in [0.5, 1.0, 2.0, 5.0]:
        k = int(c*M)
        v1 = V1(k, M)
        v2 = V2(k, M)
        var_total = v1 + v2
        MlnM_var = M * math.log(M) * var_total
        Bc = B_c(c)
        rel_err = (MlnM_var - Bc)/Bc
        results[M][c] = (v1, v2, MlnM_var, Bc, rel_err)
    print(f"M={M}: computed in {time.time()-t0:.1f}s")

# Print results
print()
print(f"{'c':>5} | " + " | ".join(f"M={M}".center(30) for M in [500, 1000, 2000, 5000]))
print(f"{'':>5} | " + " | ".join("M lnM*Var   Bc    rel_err".center(30) for _ in [500, 1000, 2000, 5000]))
print("-"*135)
for c in [0.5, 1.0, 2.0, 5.0]:
    Bc = B_c(c)
    row = f"{c:>5.1f} | "
    for M in [500, 1000, 2000, 5000]:
        v1, v2, MlnM_var, _, rel_err = results[M][c]
        row += f"{MlnM_var:>7.4f}  {Bc:>7.4f}  {rel_err:>+7.3%}".center(30) + " | "
    print(row)

print()
print(f"{'c':>5} | B(c)        | " + " | ".join(f"M={M}".center(20) for M in [500, 1000, 2000, 5000]))
print(f"{'':>5} | {'':>11} | " + " | ".join("V2/V1".center(20) for _ in [500, 1000, 2000, 5000]))
for c in [0.5, 1.0, 2.0, 5.0]:
    Bc = B_c(c)
    row = f"{c:>5.1f} | {Bc:>11.5f} | "
    for M in [500, 1000, 2000, 5000]:
        v1, v2, _, _, _ = results[M][c]
        ratio = v2/v1 if v1 != 0 else 0
        row += f"{ratio:>+7.4f}".center(20) + " | "
    print(row)

print()
print("Expected: V2/V1 -> 0 as M -> infty (V2 is lower order, Thm 6 proof).")
print("Expected: M ln M * (V1+V2) -> B(c), rel err shrinking with M.")
print()
print("VERDICT: Theorem 6 — VALID if relative errors shrink with M.")
