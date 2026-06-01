"""
Validate the simulation against BOTH kinds of recorded LEAP-UCD-2017 data:
the central-array **accelerations** (AH1..AH4) and the **pore pressures**
(P1..P4), mapped to model depth.

Requires:
  * ``data/leap_recorded.csv``  -> produced by ``input/leap_csv_to_motion.py``
  * simulated recorder output in ``output/`` (run ``input/model_2D_slope.py``)

Produces figures in ``report/``:
  * compare_accel.png        -> recorded vs simulated acceleration at matched depths
  * compare_ru_time.png      -> recorded vs simulated ru(t) at matched depths
  * compare_ru_profile.png   -> peak ru with depth, recorded (P1..P4) vs simulated

Sensor depths (prototype, below ground surface) follow the LEAP-UCD-2017
specification (Kutter et al., 2017/2018): the central array P1..P4 / AH1..AH4
sit at roughly 4, 3, 2, 1 m, with initial sigma'_v0 ~ 40, 30, 20, 10 kPa. The
recorded pore-pressure channels are excess pore pressure, so
``ru_recorded = P / sigma'_v0(depth)``. Adjust the depths below to the as-built
values from the sensor-layout sheet for your specific test.
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
G = 9.81

# central-array sensor depths below surface [m]  (LEAP-UCD-2017 standard)
ACC_DEPTH = {"AH1 (g)": 4.0, "AH2 (g)": 3.0, "AH3 (g)": 2.0, "AH4 (g)": 1.0}
PPT_DEPTH = {"P1 (kpa)": 4.0, "P2 (kpa)": 3.0, "P3 (kpa)": 2.0, "P4 (kpa)": 1.0}


def load_recorded():
    path = os.path.join(DATA, "leap_recorded.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found. Run input/leap_csv_to_motion.py first.")
    header = open(path).readline().strip().split(",")
    arr = np.genfromtxt(path, delimiter=",", skip_header=1)
    rec = {name: arr[:, i] for i, name in enumerate(header)}
    return rec, arr[:, 0]


def load_sim():
    sig0 = np.loadtxt(os.path.join(OUT, "sigma_v0.out"))
    sigma_v0 = sig0[:, 1]
    n_ele = sigma_v0.size
    n_node = n_ele + 1
    dy = SOIL_DEPTH / n_ele

    stress = np.loadtxt(os.path.join(OUT, "stress.out"))
    accel = np.loadtxt(os.path.join(OUT, "accel.out"))
    t = stress[:, 0]
    ncomp = (stress.shape[1] - 1) // n_ele
    syy = np.column_stack([stress[:, 1 + j * ncomp + 1] for j in range(n_ele)])
    ru = np.clip(1.0 - np.abs(syy) / sigma_v0[None, :], 0.0, 1.2)

    base = np.loadtxt(os.path.join(ROOT, "input", "base_motion.txt"))
    base_a = np.interp(t, base[:, 0], base[:, 1]) * G   # m/s^2

    return dict(t=t, ru=ru, sigma_v0=sigma_v0, accel=accel, base_a=base_a,
                n_ele=n_ele, n_node=n_node, dy=dy)


def ele_at_depth(sim, depth):
    """Element index whose centre is closest to `depth` below surface."""
    centres = SOIL_DEPTH - (np.arange(sim["n_ele"]) + 0.5) * sim["dy"]
    return int(np.argmin(np.abs(centres - depth)))


def node_at_depth(sim, depth):
    j = int(round((SOIL_DEPTH - depth) / sim["dy"]))
    return max(0, min(sim["n_node"] - 1, j))


def sim_abs_accel(sim, node_j):
    rel = sim["accel"][:, 1 + node_j * 2]
    return (rel + sim["base_a"]) / G            # [g], absolute


def main():
    rec, rt = load_recorded()
    sim = load_sim()

    # ---------------------------------------------------------------- accel
    chans = [c for c in ACC_DEPTH if c in rec]
    fig, axs = plt.subplots(len(chans), 1, figsize=(9, 2.1 * len(chans)),
                            sharex=True, squeeze=False)
    for ax, ch in zip(axs[:, 0], chans):
        d = ACC_DEPTH[ch]
        nj = node_at_depth(sim, d)
        ax.plot(rt, rec[ch], lw=0.7, color="k", label=f"recorded {ch}")
        ax.plot(sim["t"], sim_abs_accel(sim, nj), lw=0.7, color="tab:blue",
                alpha=0.8, label=f"sim @ {d:.0f} m")
        ax.set_ylabel("a [g]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2, loc="upper right")
    axs[-1, 0].set_xlabel("time [s]")
    axs[0, 0].set_title("Central-array acceleration: recorded vs simulated")
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT, "compare_accel.png"), dpi=130)
    plt.close(fig)

    # ----------------------------------------------------------- ru time hist
    ppts = [c for c in PPT_DEPTH if c in rec]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(ppts)))
    for ch, c in zip(ppts, colors):
        d = PPT_DEPTH[ch]
        ej = ele_at_depth(sim, d)
        sv0 = sim["sigma_v0"][ej]
        ru_rec = rec[ch] / sv0                  # recorded P is excess pp
        ax.plot(rt, np.clip(ru_rec, -0.2, 1.4), lw=0.9, color=c,
                label=f"rec {ch} @ {d:.0f} m")
        ax.plot(sim["t"], sim["ru"][:, ej], lw=0.9, ls="--", color=c,
                label=f"sim @ {d:.0f} m")
    ax.axhline(1.0, color="k", ls=":", lw=0.8)
    ax.set_xlabel("time [s]")
    ax.set_ylabel(r"$r_u$")
    ax.set_title("Excess pore-pressure ratio: recorded (solid) vs simulated (dashed)")
    ax.set_ylim(-0.1, 1.4)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT, "compare_ru_time.png"), dpi=130)
    plt.close(fig)

    # ----------------------------------------------------------- ru profile
    rec_d, rec_ru = [], []
    for ch in ppts:
        d = PPT_DEPTH[ch]
        ej = ele_at_depth(sim, d)
        rec_d.append(d)
        rec_ru.append(np.clip(rec[ch] / sim["sigma_v0"][ej], 0, 1.4).max())
    sim_depth = SOIL_DEPTH - (np.arange(sim["n_ele"]) + 0.5) * sim["dy"]
    fig, ax = plt.subplots(figsize=(4.8, 5))
    ax.plot(sim["ru"].max(axis=0), sim_depth, "-", color="tab:blue",
            label="simulated")
    ax.plot(rec_ru, rec_d, "o", color="k", ms=7, label="recorded P1..P4")
    ax.invert_yaxis()
    ax.set_xlabel("peak $r_u$")
    ax.set_ylabel("depth below surface [m]")
    ax.set_title("Peak liquefaction profile")
    ax.axvline(1.0, color="k", ls=":", lw=0.8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT, "compare_ru_profile.png"), dpi=130)
    plt.close(fig)

    # ---- console summary ----
    print("Recorded vs simulated peak ru (central array):")
    print(f"{'depth[m]':>9}{'recorded':>10}{'simulated':>11}")
    for ch in ppts:
        d = PPT_DEPTH[ch]
        ej = ele_at_depth(sim, d)
        rr = np.clip(rec[ch] / sim["sigma_v0"][ej], 0, 1.4).max()
        ss = sim["ru"][:, ej].max()
        print(f"{d:>9.0f}{rr:>10.2f}{ss:>11.2f}")
    print("\nFigures: report/compare_accel.png, compare_ru_time.png, "
          "compare_ru_profile.png")


if __name__ == "__main__":
    main()
