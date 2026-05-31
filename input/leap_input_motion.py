"""
LEAP-UCD-2017 base input motion generator.

The LEAP-UCD-2017 destructive motion is a tapered (ramped) ~1 Hz sine wave with a
superimposed high-frequency component. The real recorded achieved base motions are
available from DesignSafe (PRJ-1843). Until you drop a recorded record into
``input/base_motion.txt`` (two columns: time[s]  accel[g], or a single column of
accel with a known dt), this module synthesises a representative ramped-sine motion
so the model is runnable end-to-end.

Run directly to (re)generate ``input/base_motion.txt`` and a quick-look plot:

    python input/leap_input_motion.py
"""

from __future__ import annotations

import os
import math

# ----------------------------------------------------------------------------
# Motion parameters (representative LEAP-UCD-2017 "destructive" motion, prototype)
# ----------------------------------------------------------------------------
DT = 0.01            # time step [s]  (100 points / cycle at 1 Hz)
DURATION = 20.0      # total duration [s]
FREQ = 1.0           # dominant frequency [Hz]
PGA_G = 0.15         # target peak ground acceleration [g] (tune per test)
RAMP_UP = 1.0        # cosine taper up [s]
RAMP_DOWN = 1.0      # cosine taper down [s]
HF_FREQ = 5.0        # superimposed high-frequency component [Hz]
HF_FRACTION = 0.08   # amplitude of HF component as fraction of PGA

G = 9.81             # gravity [m/s^2]

HERE = os.path.dirname(os.path.abspath(__file__))
MOTION_FILE = os.path.join(HERE, "base_motion.txt")


def _taper(t: float, total: float) -> float:
    """Cosine taper window: 0 -> 1 over RAMP_UP, 1 -> 0 over RAMP_DOWN."""
    if t < RAMP_UP:
        return 0.5 * (1.0 - math.cos(math.pi * t / RAMP_UP))
    if t > total - RAMP_DOWN:
        return 0.5 * (1.0 - math.cos(math.pi * (total - t) / RAMP_DOWN))
    return 1.0


def generate(write: bool = True):
    """Return (times, accel_g). Optionally write to MOTION_FILE."""
    n = int(round(DURATION / DT)) + 1
    times, accel_g = [], []
    for i in range(n):
        t = i * DT
        env = _taper(t, DURATION)
        main = math.sin(2.0 * math.pi * FREQ * t)
        hf = HF_FRACTION * math.sin(2.0 * math.pi * HF_FREQ * t)
        a = PGA_G * env * (main + hf)
        times.append(t)
        accel_g.append(a)

    if write:
        with open(MOTION_FILE, "w") as f:
            f.write("# LEAP-UCD-2017 synthetic base motion\n")
            f.write("# col1: time [s]   col2: acceleration [g]\n")
            for t, a in zip(times, accel_g):
                f.write(f"{t:.5f}\t{a:.6f}\n")
        print(f"Wrote {len(times)} samples to {MOTION_FILE}")
    return times, accel_g


def load(path: str = MOTION_FILE):
    """Load a base motion file. Returns (dt, accel_ms2_list).

    Accepts either two columns (time, accel[g]) or a single accel column.
    Acceleration is converted from g to m/s^2.
    """
    if not os.path.exists(path):
        generate(write=True)

    times, accel_g = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) >= 2:
                times.append(float(parts[0]))
                accel_g.append(float(parts[1]))
            else:
                accel_g.append(float(parts[0]))

    dt = (times[1] - times[0]) if len(times) > 1 else DT
    accel_ms2 = [a * G for a in accel_g]
    return dt, accel_ms2


if __name__ == "__main__":
    t, a = generate(write=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        report_dir = os.path.join(os.path.dirname(HERE), "report")
        os.makedirs(report_dir, exist_ok=True)
        plt.figure(figsize=(9, 3))
        plt.plot(t, a, lw=0.8)
        plt.xlabel("time [s]")
        plt.ylabel("base accel [g]")
        plt.title(f"LEAP base motion (PGA={PGA_G} g, {FREQ} Hz)")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        out = os.path.join(report_dir, "base_motion.png")
        plt.savefig(out, dpi=130)
        print(f"Saved plot to {out}")
    except Exception as exc:  # pragma: no cover - plotting is optional
        print(f"(plot skipped: {exc})")
