# Coprime Density of Random LCMs

Source code, figures, and audit scripts for the paper:

**M. Wakid**, *On the Expected Coprime Density of the Least Common Multiple of Random Integers* (April 2026).

## Files

- `Coprime_Density_LCM.pdf` — the paper (15 pages)
- `Coprime_Density_LCM.tex` — LaTeX source
- `fig1_convergence.pdf`, `fig2_joint_scaling.pdf` — figures
- `regenerate_figures.py` — regenerates the figures (seed 2026)
- `audit1_identity.py` through `audit5_propA.py` — numerical verification of Theorems 1–4 and Proposition 5
- `AUDIT_REPORT.md` — consolidated verification report
- `paper1_page1_current.png` — preview of the paper's first page

## Reproducing the figures and audits

Requires Python 3 with `sympy`, `mpmath`, `numpy`, `matplotlib`:
