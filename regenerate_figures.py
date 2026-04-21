"""
Regenerate figures for v4.3. Identical to the regenerator from v4.2.
Deterministic seed: 2026.

Fig 1: Monte Carlo convergence of rho_k(M) to P_M, M in {100, 300, 1000}
Fig 2: ln M * Sigma_{cM}(M) vs A(c) for M in {1e3, 1e4, 1e5, 1e6}, c in grid
Fig 3: Standardised distribution histograms at k=M for M in {100, 300, 1000}
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange
from math import gcd

np.random.seed(2026)
SEED = 2026

# ---------------- Utilities ----------------
def primes_up_to(n):
    return list(primerange(2, n + 1))

def P_M(M):
    p = 1.0
    for pr in primes_up_to(M):
        p *= (1.0 - 1.0/pr)
    return p

def q_p(p, M):
    return 1.0 - (M // p) / (M - 1)

def Sigma_k_exact(k, M):
    return sum(q_p(pr, M)**k / (pr - 1) for pr in primes_up_to(M))

def A_c(c, jmax=200):
    return sum(math.exp(-c*j) * math.log(1 + 1/j) for j in range(1, jmax + 1))

def phi_over_L_single(X_samples, primes_list):
    """phi(L)/L for L = lcm of X_samples."""
    L = 1
    for x in X_samples:
        g = gcd(L, int(x))
        L = L * int(x) // g
    prod = 1.0
    for pr in primes_list:
        if L % pr == 0:
            prod *= (1.0 - 1.0/pr)
    return prod

def mc_rho_k(k, M, N, seed):
    rng = np.random.default_rng(seed)
    primes_list = primes_up_to(M)
    vals = np.zeros(N)
    for i in range(N):
        xs = rng.integers(2, M+1, size=k)
        vals[i] = phi_over_L_single(xs, primes_list)
    return vals

# ---------------- Figure 1 ----------------
def figure1():
    print("Figure 1: convergence plot")
    Ms = [100, 300, 1000]
    N = 2000
    k_list = [5, 8, 12, 20, 35, 50, 100, 200, 500, 1000, 2000]
    # For M=1000 we need up to k=2000 say; for M=100 same k range

    fig, ax = plt.subplots(figsize=(7.5, 5))
    colors = {100: 'tab:blue', 300: 'tab:red', 1000: 'tab:green'}
    markers = {100: 'o', 300: 's', 1000: '^'}
    for M in Ms:
        P = P_M(M)
        means, stderr, ks = [], [], []
        for k in k_list:
            vals = mc_rho_k(k, M, N, seed=SEED + M*1000 + k)
            m = vals.mean()
            s = vals.std(ddof=1) / np.sqrt(len(vals))
            means.append(m); stderr.append(s); ks.append(k)
        ax.errorbar(ks, means, yerr=stderr, marker=markers[M], color=colors[M],
                    markersize=6, capsize=2, linewidth=1.2,
                    label=f'$M={M}$')
        ax.axhline(P, linestyle='--', color=colors[M], linewidth=0.8, alpha=0.6)
        ax.text(max(ks)*1.05, P, f'$\\mathcal{{M}}({M})$', va='center', fontsize=9,
                color=colors[M])
    ax.set_xscale('log')
    ax.set_xlabel('$k$ (number of samples)')
    ax.set_ylabel(r'$\rho_k(M) = \mathbb{E}[\varphi(L_k)/L_k]$')
    ax.set_title('Convergence of coprime density to Mertens product')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('fig1_convergence.pdf', bbox_inches='tight')
    plt.close()

# ---------------- Figure 2 ----------------
def figure2():
    print("Figure 2: joint scaling plot")
    Ms = [1000, 10000, 100000, 1000000]
    cs = np.concatenate([
        np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75]),
        np.arange(1.0, 5.01, 0.25)
    ])

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    markers = {1000: 'o', 10000: 's', 100000: '^', 1000000: 'd'}
    colors = {1000: '#b0b0ff', 10000: '#d0c060', 100000: '#90c090', 1000000: '#e8b070'}

    for M in Ms:
        xs, ys = [], []
        for c in cs:
            k = int(c * M)
            if k < 1: continue
            S = Sigma_k_exact(k, M)
            xs.append(c)
            ys.append(math.log(M) * S)
        ax.plot(xs, ys, marker=markers[M], color=colors[M],
                markersize=7, linestyle='', alpha=0.9,
                label=f'$M={M}$')

    # Smooth A(c) curve
    c_fine = np.linspace(0.05, 5.0, 500)
    A_vals = [A_c(c) for c in c_fine]
    ax.plot(c_fine, A_vals, '-', color='black', linewidth=1.4,
            label=r'$A(c) = \sum_{j \geq 1} e^{-cj} \ln(1 + 1/j)$')

    ax.set_xlabel(r'$c = k/M$')
    ax.set_ylabel(r'$\ln M \cdot (\rho_{\lfloor cM \rfloor}(M) - \mathcal{M}(M))/\mathcal{M}(M)$')
    ax.set_title(r'Joint scaling: empirical convergence to $A(c)$')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 5)
    ax.set_ylim(bottom=-0.05)
    plt.tight_layout()
    plt.savefig('fig2_joint_scaling.pdf', bbox_inches='tight')
    plt.close()

# ---------------- Figure 3 ----------------
def figure3():
    print("Figure 3: CLT distribution panels")
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))

    specs = [(100, 100, 3000), (300, 300, 3000), (1000, 1000, 1000)]
    for idx, (M, k, N) in enumerate(specs):
        print(f"  Sampling M={M}, k={k}, N={N}...")
        rng = np.random.default_rng(SEED + M)
        primes_list = primes_up_to(M)
        vals = np.zeros(N)
        for i in range(N):
            xs = rng.integers(2, M+1, size=k)
            vals[i] = phi_over_L_single(xs, primes_list)
        mean = vals.mean()
        std = vals.std(ddof=1)
        z = (vals - mean) / std
        # Skewness and kurtosis
        n = len(z)
        m3 = np.mean((z)**3)
        m4 = np.mean((z)**4)
        skew = m3
        excess_kurt = m4 - 3

        ax = axs[idx]
        ax.hist(z, bins=35, range=(-4, 4), density=True, color='lightblue',
                edgecolor='steelblue', alpha=0.75)
        xs_plot = np.linspace(-4, 4, 400)
        ax.plot(xs_plot, np.exp(-xs_plot**2/2)/np.sqrt(2*np.pi), 'r-', linewidth=1.5,
                label=r'$\mathcal{N}(0, 1)$')
        ax.set_title(f'$M={M}$, $k={k}$\n'
                     f'skew = {skew:.2f}, excess kurt = {excess_kurt:.2f}',
                     fontsize=10)
        ax.set_xlabel(r'$(\varphi(L_k)/L_k - \hat\rho_k)/\hat\sigma_k$')
        if idx == 0:
            ax.set_ylabel('density')
        ax.set_xlim(-4, 4)
        ax.set_ylim(0, 0.55)
        ax.legend(fontsize=9, loc='upper left')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('fig3_distribution.pdf', bbox_inches='tight')
    plt.close()

# ---------------- Run ----------------
if __name__ == "__main__":
    import sys, os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    figure1()
    figure2()
    figure3()
    print("All figures generated.")
