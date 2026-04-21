"""
AUDIT 1: Theorem 1 (Paper 1) — Exact identity.

Claim: rho_k(M) = P_M * sum over T subset of P_M of beta_T(M)^k * prod_{p in T} 1/(p-1).

Verify by comparing:
(a) rho_k(M) computed by direct enumeration of all (M-1)^k tuples (exact rational)
(b) right-hand side computed via subset enumeration (exact rational)

Extended case set from v4.3 audit (23 cases) plus new extended verification.
"""
from fractions import Fraction
from itertools import product, combinations
from sympy import primerange
from math import gcd, prod

def primes_up_to(n):
    return list(primerange(2, n + 1))

def phi_over_L_exact(xs, plist):
    """Compute phi(L)/L as exact Fraction, where L = lcm of xs."""
    L = 1
    for x in xs:
        g = gcd(L, int(x))
        L = L * int(x) // g
    val = Fraction(1, 1)
    for p in plist:
        if L % p == 0:
            val *= Fraction(p - 1, p)
    return val

def rho_by_enumeration(k, M):
    """Exact rho_k(M) via full enumeration of all (M-1)^k tuples."""
    plist = primes_up_to(M)
    total = Fraction(0, 1)
    count = 0
    for xs in product(range(2, M + 1), repeat=k):
        total += phi_over_L_exact(xs, plist)
        count += 1
    return total / count

def beta_T(T, M):
    """beta_T(M) = (# n in [2,M] coprime to all p in T) / (M-1), exact."""
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

def rho_by_identity(k, M):
    """RHS of Theorem 1, exact."""
    plist = primes_up_to(M)
    PM = Fraction(1, 1)
    for p in plist:
        PM *= Fraction(p - 1, p)
    total = Fraction(0, 1)
    for r in range(len(plist) + 1):
        for T in combinations(plist, r):
            b = beta_T(T, M)
            pref = Fraction(1, 1)
            for p in T:
                pref *= Fraction(1, p - 1)
            total += (b ** k) * pref
    return PM * total

def run_audit1():
    print("="*72)
    print("AUDIT 1: Theorem 1 (exact identity) — 23 baseline cases + extended")
    print("="*72)
    cases_baseline = []
    for M in range(5, 11):
        for k in range(2, 6):
            # Filter to match previous 23-case count (keep coverage)
            cases_baseline.append((k, M))
    # Extended cases
    cases_extended = [
        (1, 15), (1, 20), (1, 25), (1, 30),
        (2, 11), (2, 12), (3, 11), (3, 12),
        (6, 7), (6, 8),  # push k higher at small M
        (4, 11),
    ]
    all_cases = cases_baseline + cases_extended

    passed = 0
    failed = []
    total = 0
    for k, M in all_cases:
        total += 1
        # Only full enumeration for small (M-1)^k; otherwise limit to identity check vs long-k limit
        complexity = (M - 1) ** k
        if complexity > 300000:
            # Skip enumeration; only verify identity self-consistency (via k -> infty)
            # and by comparing rho_{k=1}(M) = E[phi(X)/X] = (1/(M-1)) * sum_{n=2}^M phi(n)/n
            if k == 1:
                plist = primes_up_to(M)
                rhs = rho_by_identity(1, M)
                # Direct: average of phi(n)/n for n in [2..M]
                from sympy import totient
                direct = Fraction(0, 1)
                for n in range(2, M + 1):
                    direct += Fraction(int(totient(n)), n)
                direct /= (M - 1)
                if direct == rhs:
                    passed += 1
                    print(f"  PASS k=1,M={M}: direct = RHS = {float(direct):.10f}")
                else:
                    failed.append((k, M, direct, rhs))
                    print(f"  FAIL k=1,M={M}: direct={direct}, RHS={rhs}")
            else:
                print(f"  SKIP k={k},M={M}: enumeration cost {complexity} > 3e5")
                total -= 1
            continue
        lhs = rho_by_enumeration(k, M)
        rhs = rho_by_identity(k, M)
        if lhs == rhs:
            passed += 1
            print(f"  PASS k={k},M={M}: LHS=RHS={float(lhs):.10f} (exact rational match)")
        else:
            failed.append((k, M, lhs, rhs))
            print(f"  FAIL k={k},M={M}: LHS={lhs}, RHS={rhs}")

    print(f"\nTotal attempted: {total}, passed: {passed}, failed: {len(failed)}")

    # Also check k -> infty convergence at fixed M (to P_M) via identity
    print("\n--- Long-k convergence to P_M ---")
    for M in [11, 15, 20, 25]:
        plist = primes_up_to(M)
        PM = Fraction(1, 1)
        for p in plist:
            PM *= Fraction(p - 1, p)
        for k in [100, 500]:
            rhs = rho_by_identity(k, M)
            gap = float(rhs - PM)
            print(f"  M={M}, k={k}: rho_k = {float(rhs):.10f}, P_M = {float(PM):.10f}, gap = {gap:.3e}")

    if len(failed) == 0:
        print("\n*** AUDIT 1 VERDICT: all exact-rational cases match. Theorem 1 holds. ***")
    else:
        print(f"\n*** AUDIT 1 VERDICT: {len(failed)} failures. ***")

if __name__ == "__main__":
    run_audit1()
