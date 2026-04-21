"""
Test the proposed geometric-mean bound strategy for Conjecture 5.

Claim (Generalised Lemma 10):
    For any T subseteq P_M with |T| = r >= 1, beta_T^k <= prod_{p in T} q_p^{k/r}.

This follows from Lemma 9 (beta_T <= q_p for each p in T) applied |T| times:
    beta_T^r <= prod_{p in T} beta_T <= prod q_p.
    beta_T <= (prod q_p)^{1/r}.
    beta_T^k <= prod q_p^{k/r}.   [for k real, q_p in (0,1)]

Consequence: R_k(M) <= sum_{r >= 2} Sigma_{k/r}(M)^r / r!
where Sigma_{k'}(M) = sum_p q_p(M)^{k'} / (p-1).

In the joint scaling regime k = floor(cM), Theorem 4 applied at rate c/r
gives ln M * Sigma_{k/r}(M) -> A(c/r), so:
    ln M * R_k(M) <= sum_{r>=2} A(c/r)^r / (r! * (ln M)^{r-1}) + o(1)
    ...which --> 0 term by term, and uniformly by DCT.

This test:
1. Verifies the generalised Lemma 10 bound holds numerically.
2. Compares the bound sum_r Sigma_{k/r}^r/r! to actual R_k(M).
3. Compares ln M * (bound) to ln M * R_k.
"""
import math
from itertools import combinations
from sympy import primerange

def primes_up_to(n):
    return list(primerange(2, n + 1))

def q_p(p, M):
    return 1.0 - (M // p) / (M - 1)

def Sigma_kreal(kp, M):
    """Sigma_{k'} with real-valued k'. q_p^{k'} = exp(k' * ln q_p)."""
    total = 0.0
    for p in primes_up_to(M):
        qp = q_p(p, M)
        total += math.pow(qp, kp) / (p - 1)
    return total

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
    """R_k = sum_{|T|>=2} beta_T^k * prod 1/(p-1)."""
    primes = primes_up_to(M)
    S = Sigma_kreal(k, M)
    return rho_over_PM(k, M) - 1.0 - S

def R_k_upper_bound(k, M, Rmax=None):
    """Compute the upper bound sum_{r=2..Rmax} Sigma_{k/r}^r / r!."""
    primes = primes_up_to(M)
    if Rmax is None:
        Rmax = len(primes)
    total = 0.0
    log_fact = 0.0
    for r in range(2, Rmax + 1):
        log_fact += math.log(r)
        Srk = Sigma_kreal(k / r, M)
        total += Srk**r / math.exp(log_fact)
    return total

# ----------------------------------------------------------------------
# Test 1: Verify generalised Lemma 10 bound beta_T <= (prod q_p)^{1/|T|}
# ----------------------------------------------------------------------
def test_gen_lemma10():
    print("TEST 1: beta_T <= (prod q_p)^{1/|T|}")
    print("=" * 70)
    violations = 0
    total = 0
    max_ratio = 0.0
    worst_case = None
    for M in range(5, 40):
        primes = primes_up_to(M)
        for r in range(1, min(len(primes) + 1, 5)):
            for T in combinations(primes, r):
                b = beta_T(T, M)
                prod_q = 1.0
                for p in T:
                    prod_q *= q_p(p, M)
                rhs = math.pow(prod_q, 1.0 / r) if prod_q > 0 else 0
                total += 1
                if b > rhs + 1e-12:
                    violations += 1
                if b > 0 and rhs > 0:
                    ratio = b / rhs
                    if ratio > max_ratio:
                        max_ratio = ratio
                        worst_case = (M, T, b, rhs)
    print(f"Total (M, T) pairs tested: {total}")
    print(f"Violations: {violations}")
    print(f"Max ratio beta_T / (prod q_p)^{{1/|T|}}: {max_ratio:.6f}")
    if worst_case:
        M, T, b, rhs = worst_case
        print(f"Worst case: M={M}, T={T}, beta={b:.6f}, bound={rhs:.6f}")
    print()

# ----------------------------------------------------------------------
# Test 2: Upper bound vs actual R_k
# ----------------------------------------------------------------------
def test_upper_bound():
    print("TEST 2: Upper bound sum_{r>=2} Sigma_{k/r}^r/r!  vs  actual R_k")
    print("=" * 70)
    print(f"{'M':>4} {'c':>5} {'k':>5} {'R_k':>10} {'UB':>10} {'UB/R_k':>8} "
          f"{'lnM*R_k':>9} {'lnM*UB':>9}")
    print("-" * 70)
    for c in [0.5, 1.0, 2.0, 5.0]:
        for M in [20, 25, 30, 35]:
            k = max(1, int(c * M))
            Rk = R_k_exact(k, M)
            UB = R_k_upper_bound(k, M)
            lnM = math.log(M)
            ratio = UB / Rk if Rk > 0 else 0
            print(f"{M:>4} {c:>5.1f} {k:>5} {Rk:>10.3e} {UB:>10.3e} "
                  f"{ratio:>8.2f} {lnM*Rk:>9.4f} {lnM*UB:>9.4f}")
        print()

# ----------------------------------------------------------------------
# Test 3: Decomposition by |T| = r
# ----------------------------------------------------------------------
def R_k_by_size(k, M):
    """Return dict: |T| -> contribution to R_k."""
    primes = primes_up_to(M)
    result = {}
    for r in range(2, len(primes) + 1):
        s = 0.0
        for T in combinations(primes, r):
            b = beta_T(T, M)
            w = 1.0
            for p in T:
                w /= (p - 1)
            s += b**k * w
        result[r] = s
    return result

def test_decomposition():
    print("TEST 3: Contributions to R_k by |T|")
    print("=" * 70)
    M = 30; c = 1.0; k = int(c * M)
    contribs = R_k_by_size(k, M)
    total = sum(contribs.values())
    print(f"M={M}, c={c}, k={k}, R_k={total:.4e}")
    print(f"{'r=|T|':>6} {'contribution':>14} {'fraction':>10} {'bound_r':>14}")
    print("-" * 55)
    log_fact = 0.0
    for r in sorted(contribs.keys()):
        log_fact += math.log(r) if r > 0 else 0
    log_fact = 0.0
    for r in sorted(contribs.keys()):
        log_fact = math.lgamma(r + 1)
        bnd = Sigma_kreal(k/r, M)**r / math.exp(log_fact)
        frac = contribs[r]/total if total > 0 else 0
        print(f"{r:>6} {contribs[r]:>14.4e} {frac:>10.4%} {bnd:>14.4e}")


if __name__ == "__main__":
    test_gen_lemma10()
    test_upper_bound()
    test_decomposition()
