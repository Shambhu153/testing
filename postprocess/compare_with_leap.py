"""
Overlay the simulated response on the recorded LEAP-UCD-2017 data.

Requires:
  * ``data/leap_recorded.csv``  -> produced by ``input/leap_csv_to_motion.py``
  * simulated recorder output in ``output/`` (run ``input/model_2D_slope.py``)

Produces figures in ``report/``:
  * compare_surface_accel.png  -> recorded surface accel vs simulated surface accel
  * compare_porepressure.png   -> recorded PPT excess pore pressure vs simulated

These comparisons are qualitative: matching a recorded sensor to a model depth
requires the LEAP sensor-layout sheet (PRJ-1843). Set the channel names and the
matching model depths below to line them up quantitatively.

Run:  python postprocess/compare_with_leap.py
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
DATA = os.path.join(ROOT, "data")
REPORT = os.path.join(ROOT, "report")
os.makedirs(REPORT, exist_ok=True)

SOIL_DEPTH = 4.0

# --- which recorded channels to compare (edit to match the sensor layout) ---
SURFACE_ACCEL_CHANNEL = "AH4 (g)"            # near-surface central-array sensor
PPT_CHANNELS = ["P1 (kpa)", "P2 (kpa)", "P9 (kpa)"]


def _load_recorded():
    path = os.path.join(DATA, "leap_recorded.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found. Run input/leap_csv_to_motion.py first.")
    header = open(path).readline().strip().split(",")
    arr = np.genfromtxt(path, delimiter=",", skip_header=1)
    return {name: arr[:, i] for i, name in enumerate(header)}, arr[:, 0]


def _load_sim():
    sig0 = np.loadtxt(os.path.join(OUT, "sigma_v0.out"))
    sigma_v0 = sig0[:, 1]
    n_ele = sigma_v0.size
    n_node = n_ele + 1
    accel = np.loadtxt(os.path.join(OUT, "accel.out"))
    pwp = np.loadtxt(os.path.join(OUT, "porepressure.out"))
    t = accel[:, 0]
    surf_rel = accel[:, 1 + (n_node - 1) * 2]
    base = np.loadtxt(os.path.join(ROOT, "input", "base_motion.txt"))
    base_a = np.interp(t, base[:, 0], base[:, 1]) * 9.81
    surf_abs = (surf_rel + base_a) / 9.81
    u = pwp[:, 1:1 + n_node]
    u_exc = u - u[0:1, :]
    return dict(t=t, surf_abs=surf_abs, sigma_v0=sigma_v0,
                u_exc=u_exc, n_ele=n_ele)


def main():
    rec, rt = _load_recorded()
    sim = _load_sim()

    # ---- surface acceleration ----
    fig, ax = plt.subplots(figsize=(9, 4))
    if SURFACE_ACCEL_CHANNEL in rec:
        ax.plot(rt, rec[SURFACE_ACCEL_CHANNEL], lw=0.7, color="k",
                label=f"recorded {SURFACE_ACCEL_CHANNEL}")
    ax.plot(sim["t"], sim["surf_abs"], lw=0.7, color="tab:blue", alpha=0.8,
            label="simulated surface")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("accel [g]")
    ax.set_title("Surface acceleration: recorded vs simulated")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT, "compare_surface_accel.png"), dpi=130)
    plt.close(fig)

    # ---- pore pressure (excess) ----
    fig, ax = plt.subplots(figsize=(9, 4))
    for ch in PPT_CHANNELS:
        if ch in rec:
            ax.plot(rt, rec[ch] - rec[ch][0], lw=0.8, label=f"recorded {ch}")
    # simulated excess pp at a few node depths (qualitative)
    for j in [1, sim["n_ele"] // 2]:
        ax.plot(sim["t"], sim["u_exc"][:, j], lw=0.9, ls="--",
                label=f"sim excess pp node {j}")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("excess pore pressure [kPa]")
    ax.set_title("Excess pore pressure: recorded PPTs vs simulated")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT, "compare_porepressure.png"), dpi=130)
    plt.close(fig)

    print("Wrote report/compare_surface_accel.png and "
          "report/compare_porepressure.png")
    print("Note: comparisons are qualitative until sensor depths are mapped "
          "to model elements (see the LEAP sensor-layout sheet).")


if __name__ == "__main__":
    main()
