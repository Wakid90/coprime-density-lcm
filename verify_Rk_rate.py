"""
Verify the proof of Conjecture 5 more carefully:
- The geometric-mean bound: beta_T <= (prod q_p)^{1/|T|}.
- The upper bound R_k <= sum_{r>=2} Sigma_{k/r}^r / r!.
- The decay ln M * R_k -> 0.

Also: predict the precise leading-order behaviour of ln M * R_k.

From the proof: ln M * R_k ~ ln M * Sigma_{k/2}^2 / 2 + smaller
               ~ A(c/2)^2 / (2 ln M) + smaller

So ln M * R_k should decay like C(c)/ln M where C(c) = A(c/2)^2 / 2 approximately.
Let's verify this prediction.
"""
import math
from itertools import combinations
from sympy import primerange

def primes_up_to(n):
    return list(primerange(2, n + 1))

def q_p(p, M):
    return 1.0 - (M // p) / (M - 1)

def Sigma_kreal(kp, M):
    return sum(math.pow(q_p(p, M), kp) / (p - 1) for p in primes_up_to(M))

def A_c(c, jmax=500):
    """A(c) = sum_{j>=1} e^{-cj} ln(1+1/j)."""
    return sum(math.exp(-c*j) * math.log(1 + 1/j) for j in range(1, jmax+1))

def beta_T(T, M):
    total = 0
    for r in range(len(T) + 1):
        for S in combinations(T, r):
            prod = 1
            stop = False
            for p in S:
                prod *= p
                if prod > M: stop = True; break
            if stop and r >= 1: continue
            total += ((-1) ** r) * (M // prod)
    return (total - 1) / (M - 1)

def rho_over_PM(k, M):
    primes = primes_up_to(M)
    total = 0.0
    for m in range(len(primes) + 1):
        for T in combinations(primes, m):
            b = beta_T(T, M)
            w = 1.0
            for p in T:
                w /= (p - 1)
            total += b**k * w
    return total

def R_k_exact(k, M):
    return rho_over_PM(k, M) - 1.0 - Sigma_kreal(k, M)


print("Part A: ln M * R_k vs leading prediction A(c/2)^2 / (2 ln M)")
print("=" * 80)
print(f"{'M':>4} {'c':>4} {'k':>5} {'lnM*R_k':>10} {'pred(r=2)':>10} "
      f"{'ratio':>7} {'A(c/2)^2/2':>11}")
print("-" * 80)
for c in [0.5, 1.0, 2.0, 5.0]:
    Ac2 = A_c(c/2)
    lead = Ac2**2 / 2
    for M in [20, 25, 30, 35, 40]:
        k = max(1, int(c*M))
        lnM = math.log(M)
        Rk = R_k_exact(k, M)
        pred = lead / lnM
        ratio = lnM*Rk / pred if pred > 0 else 0
        print(f"{M:>4} {c:>4.1f} {k:>5} {lnM*Rk:>10.4e} {pred:>10.4e} "
              f"{ratio:>7.3f} {lead:>11.4f}")
    print()

# Extrapolation: lnM*R_k * lnM -> A(c/2)^2/2 ?
print()
print("Part B: Test lnM * R_k ~ C(c) / lnM, i.e., (lnM)^2 * R_k -> C(c)")
print("=" * 80)
print(f"{'M':>4} {'c':>4} {'k':>5} {'(lnM)^2 * R_k':>14} "
      f"{'A(c/2)^2/2':>11}")
print("-" * 60)
for c in [0.5, 1.0, 2.0]:
    Ac2 = A_c(c/2)
    pred = Ac2**2 / 2
    for M in [20, 25, 30, 35, 40]:
        k = max(1, int(c*M))
        lnM = math.log(M)
        Rk = R_k_exact(k, M)
        print(f"{M:>4} {c:>4.1f} {k:>5} {(lnM**2)*Rk:>14.4e} {pred:>11.4f}")
    print()
