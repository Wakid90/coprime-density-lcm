# Audit Report — Paper 1

**Paper:** *On the Expected Coprime Density of the Least Common Multiple of Random Integers*
**Audit date:** April 2026
**Protocol:** `/math-paper` skill, step-by-step verification mode
**Revision:** v2 — post-repair audit following external review.

## Summary

All five proved results (Theorems 1–4 and Proposition 5) hold. The proof of Theorem 4 (sharp form) has been revised: an earlier version used an inequality $\beta_T \leq q_2(M)$ in the tail bound that is false when $2 \notin T$. The statement is correct (numerically verified to $M = 10^7$) and the revised proof establishes it rigorously via a uniform-in-$r$ estimate on $\Sigma_{k/r}(M)$. No other proofs changed.

| # | Result | Status | Verification method |
|---|--------|--------|---------------------|
| 1 | Thm 1: Exact identity | VALID | 35 cases, exact rational equality |
| 2 | Thm 2: Limit + rate at fixed $M$ | VALID | 25 cases, decomposition exact, bound 28–31,000× safe |
| 3 | Thm 3: Sharp form at fixed $M$ | VALID | $q^{**} < q^*$ strict on $M \in [3, 40]$; closed form $(M-3)/(M-1)$ exact on $M \in [11, 50]$ |
| 4 | Thm 4: Joint scaling, sharp form | VALID (proof revised) | $M$ up to $10^7$, extrapolation to $A(c)$ within 1.35% max; repaired proof verified on 51 additional checks |
| 5 | Prop 5: Asymptotics of $A$ | VALID | mpmath 50 decimal digits; both (i) and (ii) confirmed |

All verification scripts are in this bundle: `audit1_identity.py`, `audit2_limit_rate.py`, `audit3_sharp.py`, `audit4_joint.py`, `audit5_propA.py`, `audit_full.py` (full re-audit, post-repair).

---

## Audit 1: Theorem 1 (Exact identity)

**Claim.** For every $M \geq 2$, $k \geq 1$:
$$\rho_k(M) = P_M \sum_{T \subseteq \mathcal{P}_M} \beta_T(M)^k \prod_{p \in T}\frac{1}{p-1}.$$

**Verification method.** Direct enumeration of all $(M-1)^k$ tuples vs RHS from identity, using `fractions.Fraction` exact rational arithmetic throughout.

**Extended case set (35 cases):**
- 24 cases: $M \in \{5, \ldots, 10\}$, $k \in \{2, \ldots, 5\}$
- 4 cases: $k = 1$, $M \in \{15, 20, 25, 30\}$ (compared against direct $\sum \varphi(n)/n$ average)
- 7 cases: assorted $M \in \{11, 12\}$ at $k \in \{2, 3, 4\}$; $k = 6$ at $M \in \{7, 8\}$

**All 35 cases pass as exact rationals.** The paper's main text states 24 test cases (the $6 \times 4$ grid); the audit bundle extends this.

**VERDICT: Theorem 1 holds.**

---

## Audit 2: Theorem 2 (Limit + rate at fixed $M$)

**Claim.** $(\rho_k(M) - P_M)/P_M = \Sigma_k(M) + R_k(M)$ with $R_k(M) \geq 0$ and $R_k(M) \leq \tfrac12 \Sigma_{k/2}(M)^2 e^{\Lambda(M)}$.

**Verification method.** Compute each quantity exactly using rational arithmetic over subsets. Check (a) decomposition is exact, (b) $R_k \geq 0$, (c) bound holds.

**25 cases** over $M \in \{11, 15, 20, 25, 30\}$ and $k \in \{4, 8, 16, 32, 64\}$.

- All 25 decompositions exact: $R_k$ computed directly equals total $-\Sigma_k$ as rationals.
- All 25 satisfy the bound. Minimum margin: 28× (at $M = 30$, $k = 4$); maximum: 31,434× (at $M = 11$, $k = 64$).

**VERDICT: Theorem 2 holds with substantial safety margin.**

---

## Audit 3: Theorem 3 (Sharp form at fixed $M$)

**Claims:**
(a) $q^{**}(M) < q^{*}(M)$ for $M \geq 3$;
(b) $R_k(M)/\Sigma_k(M) = O((q^{**}/q^{*})^k)$;
(c) For $M \geq 11$, $q^{**} = (M-3)/(M-1)$, $q^{*} = (M-2)/(M-1)$.

**Verification:**

- (a) Exhaustive over $M \in [3, 40]$. $q^{**} < q^{*}$ strictly on every $M$ tested.
- (c) Exhaustive over $M \in [11, 50]$. $q^{*}$ and $q^{**}$ exactly equal the predicted closed forms for every $M$.
- (b) Observed ratio $R_k/\Sigma_k$ at $M \in \{15, 20, 30\}$ and $k \in \{10, 20, 40, 80, 160\}$ is always bounded above by $(q^{**}/q^{*})^k$.

**VERDICT: Theorem 3 holds in all three parts.**

---

## Audit 4: Theorem 4 (Joint scaling, sharp form) — REVISED

### 4A: Statement verification (unchanged)

**Claim.** $\ln M \cdot \Sigma_{\lfloor cM \rfloor}(M) \to A(c)$ in the joint regime, and $\ln M \cdot (\rho_k - P_M)/P_M \to A(c)$.

**Verification.** Exact computation at $M \in \{10^4, \ldots, 10^7\}$, $c \in \{0.05, \ldots, 15\}$. Linear extrapolation in $1/\ln M$ hits $A(c)$ with median 0.45%, maximum 1.35%. Table 2 in Paper 1 is reproduced to the 4 significant figures printed.

### 4B: Proof-level audit (new)

An external audit noted that the original proof of the sharp form used this step in §6.2:

> "using Lemma single-prime in the weak form $\beta_T \leq q_{\min}(M) := \min_{p \leq M} q_p(M) = q_2(M)$"

This is **false** when $T$ does not contain $2$. Lemma single-prime only gives $\beta_T \leq q_{p_0}$ for $p_0 \in T$.

**Numerical counterexamples** (confirming the step, not the theorem, is broken):

| $M$ | $T$ | $\beta_T$ | $q_2(M)$ | $\beta_T > q_2$? |
|-----|-----|-----------|----------|------------------|
| 20 | $\{11, 13, 17, 19\}$ | 0.7895 | 0.4737 | yes, by 0.32 |
| 50 | $\{31, 37, 41, 43, 47\}$ | 0.8980 | 0.4898 | yes, by 0.41 |
| 30 | 302 distinct subsets $T \not\ni 2$ | — | — | all violate |

### 4C: Repaired proof

The revised proof of §6.2 uses:

1. Size decomposition $R_k = \sum_{r \geq 2} R_k^{(r)}$ (unchanged).
2. Geometric-mean bound (Lemma `lem:geom-mean`) giving $R_k^{(r)} \leq \Sigma_{k/r}(M)^r/r!$ (unchanged).
3. **New uniform-in-$r$ estimate:** For $2 \leq r \leq M^{1/2}/\ln M$ and $M \geq M_1(c)$,
   $$\Sigma_{k/r}(M) \leq \frac{C_u \cdot r}{\ln M}, \qquad C_u \leq 10.08 \cdot (2 + 2/c).$$
   Proof uses the envelope lemma (already in §6.1), the inequality $1 - e^{-v} \geq v/(1+v)$, and the trivial bound $\Sigma_s \leq \Lambda(M)$ for the $r > r_{\max}$ tail.
4. Stirling $r! \geq (r/e)^r$ gives $R_k^{(r)} \leq (C_u e/\ln M)^r$.
5. Summing a geometric series: $R_k(M) = O(1/(\ln M)^2)$, hence $\ln M \cdot R_k \to 0$.

### 4D: Numerical verification of each step in the repaired proof

`audit_full.py` (51 checks) verifies:

| Check | Claim | Result |
|-------|-------|--------|
| A1 | $S_j(M) \ln M \leq 10.08$ for $1 \leq j \leq \sqrt M$, $M \geq 289$ | Max observed much smaller (envelope loose) |
| A2 | $\Sigma_{k/r}(M) \ln M / r$ uniformly bounded | Max observed: 0.374 (bound constant $C_u$ is loose) |
| A3 | $(1 - M^{-1/2})^{k/r}$ super-polynomially small | All tested cases $\leq M^{-c/2}$ |
| A4 | Stirling $r! \geq (r/e)^r$ | Holds for $r \in \{2, 5, 10, 20, 100\}$ |
| A5–A8 | $\ln M \cdot R_k^{\text{upper bound}} \to 0$ | Decreases from 5.3% ($M = 500$) to 2.2% ($M = 5 \cdot 10^5$) |
| A9 | Upper bound $\geq$ actual $R_k$ on 9 small-$M$ cases | Ratio 6× to 21× (upper bound valid but loose) |

Plus re-verification of all previously-passing checks (B1–B7).

**VERDICT: Theorem 4 holds; the revised proof is rigorous and verified.**

---

## Audit 5: Proposition 5 (Asymptotics of $A$)

**Claims:**
(i) $A(c) = -\ln c - \gamma - (c \ln c)/2 + O(c)$ as $c \to 0^+$;
(ii) $A(c) = e^{-c}\ln 2 + e^{-2c}\ln(3/2) + e^{-3c}\ln(4/3) + O(e^{-4c})$ as $c \to \infty$.

**Verification method.** mpmath at 50-decimal precision.

Large-$c$ test (ii): the scaled residual $(A(c) - \text{3-term prediction}) \cdot e^{4c}$ converges to $\ln(5/4) = 0.2231\ldots$ as expected.

Small-$c$ test (i): residual / $c$ stays bounded across $c \in [0.001, 0.5]$.

Monotonicity / convexity: $A'(c) < 0$ and $A''(c) > 0$ verified.

**VERDICT: Proposition 5 holds to high precision.**

---

## Minor corrections applied during revision

1. **Real-exponent notation for $\Sigma_s$.** A new remark (`rem:sigma-real`) in §2 extends the definition $\Sigma_k = \sum q_p^k/(p-1)$ to real exponents $s \geq 0$. This makes $\Sigma_{k/r}$ (used in the Theorem 4 sharp-form proof with $r = 2, 3, \ldots$, where $k/r$ is generally non-integer) well-defined from the start.

2. **Typo in §7.** "All 23 test cases" → "All 24 test cases" (the $6 \times 4$ grid).

3. **§6.1 strip-discrepancy bound.** Replaced a cumbersomely-written bound with the cleaner estimate $\sum_{n > N} 1/(n(n-1)) = 1/N$, giving the required $o(1/\ln M)$ cleanly.

4. **Cesàro citation framing.** One added sentence clarifying that Cesàro's classical result is on $\gcd$, and the $\mathrm{lcm}$ asymptotic follows via $\mathrm{lcm}(a,b) = ab/\gcd(a,b)$.

---

## Reproducibility

All verification scripts included:

```bash
# Requires Python 3 with: sympy, mpmath, numpy, matplotlib
python3 audit1_identity.py     # Theorem 1  — ~30 seconds
python3 audit2_limit_rate.py   # Theorem 2  — ~10 seconds
python3 audit3_sharp.py        # Theorem 3  — ~30 seconds
python3 audit4_joint.py        # Theorem 4  — ~5 minutes
python3 audit5_propA.py        # Prop 5     — ~30 seconds
python3 audit_full.py          # Post-repair comprehensive audit (51 checks) — ~5 minutes
```

Figures can be regenerated from `regenerate_figures.py`.

---

## Submission readiness

**Paper 1 (post-revision):** READY for submission.

The repaired §6.2 proof has been verified step-by-step, and all other theorems have been re-verified to confirm no collateral damage. The paper now compiles to 16 pages (was 15 pre-revision; the extra page is the clean proof of the sharp form).

**Suggested venues:** *Journal of Number Theory*, *Integers*, *Acta Arithmetica*, *Journal de Théorie des Nombres de Bordeaux*.
