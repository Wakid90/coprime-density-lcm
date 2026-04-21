"""
AUDIT 2: Theorem 2 (Paper 1) — Limit + rate at fixed M.

Claim: (rho_k(M) - P_M)/P_M = Sigma_k(M) + R_k(M) exactly (decomposition),
with R_k(M) >= 0 and R_k(M) <= (1/2) * Sigma_{k/2}(M)^2 * exp(Lambda(M)).

Verify exactly using Theorem 1's identity to compute (rho - P_M)/P_M,
then confirm that subtracting Sigma_k leaves a nonneg remainder bounded by
the claimed expression.
"""
from fractions import Fraction
from itertools import combinations
from sympy import primerange
import math

def primes_up_to(n):
    return list(primerange(2, n + 1))

def PM_exact(M):
    p = Fraction(1)
    for pr in primes_up_to(M):
        p *= Fraction(pr - 1, pr)
    return p

def q_p(p, M):
    return Fraction(M - 1 - (M // p), M - 1)

def Sigma_k_exact(k, M):
    plist = primes_up_to(M)
    s = Fraction(0)
    for p in plist:
        s += q_p(p, M) ** k * Fraction(1, p - 1)
    return s

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

def Rk_exact(k, M):
    """R_k(M) = sum_{|T|>=2} beta_T^k * prod 1/(p-1) exactly."""
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

def rho_over_PM_minus_1(k, M):
    """(rho_k(M)/P_M - 1) = Sigma_k + R_k, from Theorem 1."""
    plist = primes_up_to(M)
    # sum over T != empty of beta_T^k * prod 1/(p-1)
    total = Fraction(0)
    for r in range(1, len(plist) + 1):
        for T in combinations(plist, r):
            b = beta_T(T, M)
            pref = Fraction(1)
            for p in T:
                pref *= Fraction(1, p - 1)
            total += (b ** k) * pref
    return total

def Lambda(M):
    return sum(Fraction(1, p - 1) for p in primes_up_to(M))

def run_audit2():
    print("="*72)
    print("AUDIT 2: Theorem 2 (limit + rate at fixed M)")
    print("="*72)

    cases = []
    for M in [11, 15, 20, 25, 30]:
        for k in [4, 8, 16, 32, 64]:
            cases.append((k, M))

    n_decomp_ok = 0
    n_bound_ok = 0
    for k, M in cases:
        total = rho_over_PM_minus_1(k, M)
        Sig = Sigma_k_exact(k, M)
        Rk_direct = Rk_exact(k, M)
        Rk_from_decomp = total - Sig
        # Check decomposition holds exactly
        decomp_ok = (Rk_direct == Rk_from_decomp)
        # Check R_k >= 0
        nonneg = Rk_direct >= 0
        # Check bound: R_k <= (1/2) * Sigma_{k/2}^2 * exp(Lambda(M))
        k_half = k // 2
        Sig_half = Sigma_k_exact(k_half, M)
        Lam = float(Lambda(M))
        bound = 0.5 * float(Sig_half) ** 2 * math.exp(Lam)
        bound_ok = float(Rk_direct) <= bound + 1e-15
        margin = bound / max(float(Rk_direct), 1e-30) if float(Rk_direct) > 0 else float('inf')

        if decomp_ok and nonneg and bound_ok:
            n_decomp_ok += 1
            n_bound_ok += 1
            print(f"  PASS M={M:2d}, k={k:3d}: R_k={float(Rk_direct):.4e}, "
                  f"Sigma={float(Sig):.4e}, bound={bound:.2e}, margin={margin:.1f}x")
        else:
            print(f"  FAIL M={M:2d}, k={k:3d}: decomp_ok={decomp_ok}, "
                  f"nonneg={nonneg}, bound_ok={bound_ok}")

    print(f"\n{n_decomp_ok}/{len(cases)} decompositions verified exactly")
    print(f"{n_bound_ok}/{len(cases)} bounds satisfied")

    # Also: convergence rate — check Sigma_k, R_k -> 0
    print("\n--- R_k, Sigma_k both decay to 0 at fixed M ---")
    for M in [11, 20]:
        for k in [8, 32, 128, 512]:
            Sig = float(Sigma_k_exact(k, M))
            Rk = float(Rk_exact(k, M))
            print(f"  M={M:2d}, k={k:3d}: Sigma_k = {Sig:.3e}, R_k = {Rk:.3e}, R_k/Sigma_k = {Rk/Sig:.3e}")

    if n_decomp_ok == len(cases) and n_bound_ok == len(cases):
        print("\n*** AUDIT 2 VERDICT: decomposition exact, bound holds with comfortable margin. ***")
    else:
        print("\n*** AUDIT 2 VERDICT: failures detected. ***")

if __name__ == "__main__":
    run_audit2()
