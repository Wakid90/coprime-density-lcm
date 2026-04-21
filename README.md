# Coprime Density of Random LCMs

Source code, figures, and audit scripts for the paper:

**M. Wakid**, *On the Expected Coprime Density of the Least Common Multiple of Random Integers* (April 2026).

## Files

- `Coprime_Density_LCM.pdf` — the paper (16 pages, v2 post-revision)
- `Coprime_Density_LCM.tex` — LaTeX source
- `fig1_convergence.pdf`, `fig2_joint_scaling.pdf` — figures
- `regenerate_figures.py` — regenerates the figures (seed 2026)
- `audit1_identity.py` through `audit5_propA.py` — numerical verification of Theorems 1–4 and Proposition 5
- `audit_full.py` — comprehensive post-revision audit (51 checks: repair verification + re-verification of all theorems)
- `AUDIT_REPORT.md` — consolidated verification report
- `paper1_page1_current.png` — preview of the paper's first page

## Reproducing the figures and audits

Requires Python 3 with `sympy`, `mpmath`, `numpy`, `matplotlib`:

```
python3 regenerate_figures.py
python3 audit1_identity.py
python3 audit2_limit_rate.py
python3 audit3_sharp.py
python3 audit4_joint.py
python3 audit5_propA.py
python3 audit_full.py
```

## Compiling the paper

Requires a LaTeX installation (TeX Live or MiKTeX):

```
pdflatex Coprime_Density_LCM.tex
pdflatex Coprime_Density_LCM.tex
pdflatex Coprime_Density_LCM.tex
```

## Author

Muhamad Wakid — Independent researcher
ORCID: [0009-0008-6274-778X](https://orcid.org/0009-0008-6274-778X)
