"""
Post-process LEAP-UCD-2017 sloping-ground results and write figures to report/.

Reads the recorder output written by ``input/model_2D_slope.py`` and produces:
  * base vs surface acceleration time history
  * excess pore-pressure ratio (ru) time histories at selected depths
  * peak ru profile with depth
  * lateral-displacement time history (surface) and final profile

ru is computed two ways and both are reported:
  ru_sigma = 1 - sigma'_v(t) / sigma'_v0     (from element effective stress)
  ru_pwp   = (u(t) - u0) / sigma'_v0         (from nodal pore pressure)

Run:  python postprocess/plot_results.py
"""

from __future__ import annotations

import os
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "output")
REPORT = os.path.join(ROOT, "report")
INPUT = os.path.join(ROOT, "input")
os.makedirs(REPORT, exist_ok=True)

# geometry must match the model
SOIL_DEPTH = 4.0


def _load(name):
    path = os.path.join(OUT, name)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found. Run input/model_2D_slope.py first.")
    return np.loadtxt(path)


def main():
    sig0 = np.loadtxt(os.path.join(OUT, "sigma_v0.out"))
    sigma_v0 = sig0[:, 1]
    n_ele = sigma_v0.size
    n_node = n_ele + 1
    dy = SOIL_DEPTH / n_ele

    stress = _load("stress.out")
    disp = _load("disp.out")
    pwp = _load("porepressure.out")

    t = stress[:, 0]
    ncomp = (stress.shape[1] - 1) // n_ele       # stress components per element

    # ---- effective vertical stress and ru_sigma per element ----
    syy = np.empty((t.size, n_ele))
    for j in range(n_ele):
        syy[:, j] = stress[:, 1 + j * ncomp + 1]   # +1 -> sigma_yy
    ru_sigma = 1.0 - np.abs(syy) / sigma_v0[None, :]
    ru_sigma = np.clip(ru_sigma, 0.0, 1.2)

    # ---- ru from pore pressure (excess relative to start of shaking) ----
    # pwp columns: time + one column (dof3) per node, base..surface
    u = pwp[:, 1:1 + n_node]
    u_exc = u - u[0:1, :]                          # excess pp
    # map node excess pp to element (use lower node of each element)
    ru_pwp = u_exc[:, :n_ele] / sigma_v0[None, :]
    ru_pwp = np.clip(ru_pwp, -0.2, 1.2)

    depth_ele = SOIL_DEPTH - (np.arange(n_ele) + 0.5) * dy

    # ---- surface acceleration (relative) + base motion -> absolute ----
    accel = _load("accel.out")
    surf_ax_rel = accel[:, 1 + (n_node - 1) * 2]   # surface node ux accel
    base = np.loadtxt(os.path.join(INPUT, "base_motion.txt"))
    bt, bg = base[:, 0], base[:, 1]                # base time [s], accel [g]
    base_interp = np.interp(t, bt, bg) * 9.81
    surf_ax_abs = surf_ax_rel + base_interp

    # ---- lateral displacement ----
    ux = disp[:, 1:1 + n_node * 2:2]               # ux for every node
    surf_ux = ux[:, -1]
    final_ux = ux[-1, :]
    y_nodes = np.arange(n_node) * dy

    # ========================================================================
    # FIGURES
    # ========================================================================
    # 1. acceleration
    fig, ax = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
    ax[0].plot(t, base_interp / 9.81, lw=0.7, color="tab:gray")
    ax[0].set_ylabel("base accel [g]")
    ax[0].grid(alpha=0.3)
    ax[1].plot(t, surf_ax_abs / 9.81, lw=0.7, color="tab:blue")
    ax[1].set_ylabel("surface accel [g]")
    ax[1].set_xlabel("time [s]")
    ax[1].grid(alpha=0.3)
    ax[0].set_title("Input (base) vs surface acceleration")
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT, "acceleration.png"), dpi=130)
    plt.close(fig)

    # 2. ru time histories at a few depths
    fig, ax = plt.subplots(figsize=(9, 4))
    picks = [0, n_ele // 4, n_ele // 2, 3 * n_ele // 4]
    for j in picks:
        ax.plot(t, ru_sigma[:, j], lw=0.9,
                label=f"depth {depth_ele[j]:.1f} m")
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xlabel("time [s]")
    ax.set_ylabel(r"$r_u = 1 - \sigma'_v/\sigma'_{v0}$")
    ax.set_title("Excess pore-pressure ratio (from effective stress)")
    ax.set_ylim(0, 1.15)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT, "excess_pwp_ratio.png"), dpi=130)
    plt.close(fig)

    # 3. peak ru profile
    fig, ax = plt.subplots(figsize=(4.5, 5))
    ax.plot(ru_sigma.max(axis=0), depth_ele, "-o", ms=3, label=r"$r_u$ (stress)")
    ax.plot(np.clip(ru_pwp.max(axis=0), 0, 1.2), depth_ele, "-s", ms=3,
            label=r"$r_u$ (pwp)")
    ax.invert_yaxis()
    ax.set_xlabel("peak $r_u$")
    ax.set_ylabel("depth [m]")
    ax.set_title("Peak liquefaction profile")
    ax.axvline(1.0, color="k", ls="--", lw=0.8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT, "ru_profile.png"), dpi=130)
    plt.close(fig)

    # 4. lateral displacement
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
    ax[0].plot(t, surf_ux, lw=0.9, color="tab:red")
    ax[0].set_xlabel("time [s]")
    ax[0].set_ylabel("surface lateral disp [m]")
    ax[0].set_title("Down-slope surface displacement")
    ax[0].grid(alpha=0.3)
    ax[1].plot(final_ux, y_nodes, "-o", ms=3, color="tab:red")
    ax[1].set_xlabel("lateral disp [m]")
    ax[1].set_ylabel("elevation [m]")
    ax[1].set_title("Final lateral-disp profile")
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT, "lateral_displacement.png"), dpi=130)
    plt.close(fig)

    # ---- summary ----
    summary = os.path.join(REPORT, "summary.txt")
    with open(summary, "w") as f:
        f.write("LEAP-UCD-2017 sloping-ground - results summary\n")
        f.write("=" * 48 + "\n")
        f.write(f"duration analysed        : {t[-1]:.2f} s\n")
        f.write(f"peak base accel          : {np.max(np.abs(base_interp))/9.81:.3f} g\n")
        f.write(f"peak surface accel       : {np.max(np.abs(surf_ax_abs))/9.81:.3f} g\n")
        f.write(f"max ru (any depth)       : {ru_sigma.max():.3f}\n")
        f.write(f"residual surface disp    : {surf_ux[-1]:.4f} m\n")
        f.write(f"peak surface disp        : {np.max(np.abs(surf_ux)):.4f} m\n")
    print(open(summary).read())
    print(f"Figures written to {REPORT}")


if __name__ == "__main__":
    main()
