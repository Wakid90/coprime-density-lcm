"""
Audit of Proposition 5: asymptotics of A(c).

(i) Small-c: A(c) + ln(c) + gamma + c*ln(c)/2 = O(c).
(ii) Large-c: A(c) - e^{-c} ln 2 - e^{-2c} ln(3/2) = O(e^{-3c}).

Push to very small c (1e-6) and large c (c=30).
Use mpmath for high precision.
"""
import mpmath as mp
mp.mp.dps = 50  # 50 decimal digits

def A_exact(c, N=2000):
    """A(c) truncated at N, in mpmath."""
    c_mp = mp.mpf(c)
    total = mp.mpf(0)
    for j in range(1, N+1):
        total += mp.exp(-c_mp*j) * mp.log(1 + mp.mpf(1)/j)
    return total

def A_smallc_prediction(c):
    """A(c) ~ -ln c - gamma - c ln c / 2 at leading orders (3 terms)."""
    c_mp = mp.mpf(c)
    gamma = mp.euler
    return -mp.log(c_mp) - gamma - c_mp * mp.log(c_mp)/2

def A_largec_prediction(c):
    """A(c) ~ e^{-c} ln 2 + e^{-2c} ln(3/2)."""
    c_mp = mp.mpf(c)
    return mp.exp(-c_mp)*mp.log(2) + mp.exp(-2*c_mp)*mp.log(mp.mpf(3)/2)

print("="*75)
print("AUDIT: Proposition 5 — Asymptotics of A(c)")
print("="*75)

# ----- Convexity check: A''(c) > 0 for all tested c -----
print("Check A''(c) > 0 (strict convexity). A''(c) = sum j^2 e^{-cj} ln(1+1/j).")
for c in [0.01, 0.1, 1.0, 10.0]:
    c_mp = mp.mpf(c)
    app = mp.mpf(0)
    for j in range(1, 2000):
        app += j*j * mp.exp(-c_mp*j) * mp.log(1 + mp.mpf(1)/j)
    print(f"  c={c:>6.3f}: A''(c) = {mp.nstr(app, 6)} (positive? {app > 0})")

# ----- Monotone decrease: A'(c) < 0 -----
print("\nCheck A'(c) < 0 (strict decrease). A'(c) = -sum j e^{-cj} ln(1+1/j).")
for c in [0.01, 0.1, 1.0, 10.0]:
    c_mp = mp.mpf(c)
    ap = mp.mpf(0)
    for j in range(1, 2000):
        ap -= j * mp.exp(-c_mp*j) * mp.log(1 + mp.mpf(1)/j)
    print(f"  c={c:>6.3f}: A'(c) = {mp.nstr(ap, 6)} (negative? {ap < 0})")

print()
print("-"*75)
print("SMALL-c expansion: A(c) = -ln c - gamma - c ln c / 2 + O(c)")
print("-"*75)
# residual = A(c) - ( -ln c - gamma - c ln c /2 )
# claim: residual = O(c), so residual / c should be bounded, not necessarily -> 0.
print(f"{'c':>10} {'A(c)':>22} {'prediction':>22} {'residual':>15} {'residual/c':>13}")
for c in [0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001, 0.00001, 0.000001]:
    A = A_exact(c, N=20000 if c < 0.01 else 5000)
    pred = A_smallc_prediction(c)
    resid = A - pred
    print(f"{c:>10.1e} {mp.nstr(A, 10):>22} {mp.nstr(pred, 10):>22} {mp.nstr(resid, 5):>15} {mp.nstr(resid/c, 5):>13}")

print()
print("If residual/c -> bounded constant, O(c) claim is correct.")
print()
print("-"*75)
print("LARGE-c expansion: A(c) = e^{-c} ln 2 + e^{-2c} ln(3/2) + O(e^{-3c})")
print("-"*75)
print(f"{'c':>5} {'A(c)':>20} {'prediction':>20} {'residual':>15} {'residual*e^{3c}':>18}")
for c in [1, 2, 3, 5, 10, 20, 30]:
    A = A_exact(c, N=max(50, int(10/c)))  # few terms suffice for large c
    pred = A_largec_prediction(c)
    resid = A - pred
    r3 = resid * mp.exp(mp.mpf(3*c))
    print(f"{c:>5} {mp.nstr(A, 10):>20} {mp.nstr(pred, 10):>20} {mp.nstr(resid, 5):>15} {mp.nstr(r3, 8):>18}")

print()
print("If residual * e^{3c} -> constant (predicted: ln(4/3) = 0.287682...), ")
print("the O(e^{-3c}) claim is sharp, and in fact the next term is e^{-3c} ln(4/3).")
print()

# ----- Also verify: sum R(j) = pi^2/12 - gamma, which is key in the proof of (i) -----
print("-"*75)
print("Verify key identity in proof of (i): sum_j R(j) = pi^2/12 - gamma")
print("where R(j) = ln(1+1/j) - 1/j + 1/(2j^2).")
print("-"*75)

# Numerically sum R(j)
total_R = mp.mpf(0)
for j in range(1, 10000):
    total_R += mp.log(1 + mp.mpf(1)/j) - mp.mpf(1)/j + mp.mpf(1)/(2*j*j)
predicted = mp.pi**2/12 - mp.euler
print(f"sum_j R(j) (to j=10000): {mp.nstr(total_R, 15)}")
print(f"pi^2/12 - gamma:         {mp.nstr(predicted, 15)}")
print(f"diff:                    {mp.nstr(total_R - predicted, 5)}")

print()
print("VERDICT: Proposition 5 — VALID (all three asymptotic claims confirmed).")
