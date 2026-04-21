"""
Verify: for M large enough, q** = beta_{p1,p2}(M) where p1, p2 are the two
largest primes <= M, and both satisfy m_p = 1 (i.e., p > M/2).

Claim: for M >= 5 with Bertrand guaranteeing a prime in (M/2, M],
  q** = 1 - 2/(M-1) = (M-3)/(M-1).

Compare this to q* = 1 - 1/(M-1) = (M-2)/(M-1).

Then q**/q* = (M-3)/(M-2).

Verify at several M.
"""
import math
from itertools import combinations
from sympy import primerange

def primes_up_to(n): return list(primerange(2, n+1))

# Check that the maximizer is always two primes > M/2
# and that q** = (M-3)/(M-1)
print(f"{'M':>4} {'q** (brute)':>12} {'(M-3)/(M-1)':>14} {'match?':>8}")
all_match = True
for M in range(5, 60):
    primes = primes_up_to(M)
    q_ss_brute = 0.0
    for r in range(2, len(primes) + 1):
        for T in combinations(primes, r):
            # beta_T via inclusion-exclusion
            total = 0
            for rr in range(len(T) + 1):
                for S in combinations(T, rr):
                    prod = 1
                    stop = False
                    for p in S:
                        prod *= p
                        if prod > M: stop = True; break
                    if stop and rr >= 1: continue
                    total += ((-1)**rr) * (M // prod)
            b = (total - 1) / (M-1)
            if b > q_ss_brute:
                q_ss_brute = b
    predicted = (M-3)/(M-1)
    match = abs(q_ss_brute - predicted) < 1e-10
    if not match: all_match = False
    print(f"{M:>4} {q_ss_brute:>12.8f} {predicted:>14.8f} {str(match):>8}")

print(f"\nAll match: {all_match}")
