# Coprime Density of Random LCMs — code and verification

This repository contains the Python code, random seeds, and raw numerical output
accompanying the paper:

**Muhamad Wakid**, *On the Expected Coprime Density of the Least Common Multiple
of Random Integers*, 2026.
[ORCID: 0009-0008-6274-778X](https://orcid.org/0009-0008-6274-778X)

## Contents

### Figures
- `regenerate_figures.py` — regenerates all three figures in the paper
  (Monte Carlo convergence, joint-scaling convergence, CLT histograms).
  Deterministic, seed 2026. Runtime ~90 seconds.
- `fig1_convergence.pdf`, `fig2_joint_scaling.pdf`, `fig3_distribution.pdf` —
  the figures as they appear in the paper.

### Theorem audits
Each script independently verifies one theorem in the paper by pushing numerics
well beyond the values reported in the paper's Numerical Evidence section.

- `audit_thm1.py` — exact rational verification of Theorem 1 at 23 test cases.
- `audit_thm2.py` — decomposition and bound check for Theorem 2.
- `audit_thm3.py` — five sub-claims of Theorem 4 (sharp at fixed M).
- `audit_thm4.py` — joint-scaling Theorem 5, recomputed to M = 10^7.
- `audit_prop5.py` — asymptotics of A(c) (Proposition 6), high-precision mpmath.
- `audit_thm6.py` — variance asymptotic Theorem 7.

### Auxiliary verification scripts
- `test_sharp_fixed_M.py`, `verify_q_star_star.py`, `verify_q_ss_formula.py`,
  `verify_geometric_mean_bound.py`, `verify_Rk_rate.py` — exploratory scripts
  used during development.

## Requirements
