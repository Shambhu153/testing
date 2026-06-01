"""
Convert a LEAP-UCD-2017 *Processed* CSV (DesignSafe PRJ-1843) into the input
files the model uses.

The processed LEAP files (e.g. ``CU2_Motion1_Processed_4222.csv``) contain, at
*prototype* scale:
  Time (sec), AH1..AH12 (g)  -> horizontal accelerometers,
              AV1..AV2 (g)   -> vertical accelerometers,
              P1..P10  (kpa) -> pore-pressure transducers (PPTs).

They are sampled at a high rate during shaking (dt ~ 0.0067 s) followed by a
coarse post-shaking consolidation tail. This script:

  1. reads the CSV,
  2. keeps the high-rate dynamic segment and trims to [T_START, T_END],
  3. resamples the chosen base-input channel onto a uniform dt,
  4. writes ``input/base_motion.txt``  (time[s]  accel[g])  -> used by the model,
  5. writes ``data/leap_recorded.csv`` (trimmed, uniform) -> used by
     ``postprocess/plot_results.py`` to overlay the measured response.

IMPORTANT - choose the base channel
-----------------------------------
``BASE_CHANNEL`` must be the accelerometer that recorded the *achieved base /
container motion*. In the standard LEAP-UCD-2017 layout the base/container
accelerometers are AH11/AH12 (they carry the raw, unfiltered table content);
the soil-embedded central-array sensors (AH2..AH9) show the filtered ~1 Hz
motion. ALWAYS confirm against the sensor-layout sheet that ships with the test
in PRJ-1843. Run with ``--list`` to print a channel summary (PGA + dominant
frequency) to help identify it.

Usage
-----
    python input/leap_csv_to_motion.py --list
    python input/leap_csv_to_motion.py --base "AH11 (g)" --t0 0 --t1 25
"""

from __future__ import annotations

import os
import argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")

# ----------------------------------------------------------------------------
# Defaults (override on the command line)
# ----------------------------------------------------------------------------
CSV_DEFAULT = os.path.join(DATA, "CU2_Motion1_Processed_4222.csv")
BASE_CHANNEL = "AH11 (g)"  # container/base accel (VERIFY vs sensor-layout sheet!)
T_START = 0.0              # trim window start [s]
T_END = 25.0              # trim window end   [s] (covers ramp + strong motion)
MAX_DYN_DT = 0.05         # rows with dt above this belong to the slow tail

MOTION_FILE = os.path.join(HERE, "base_motion.txt")
RECORDED_FILE = os.path.join(DATA, "leap_recorded.csv")


def read_csv(path):
    with open(path) as f:
        header = [h.strip() for h in f.readline().strip().strip(",").split(",")]
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    data = data[:, : len(header)]
    return header, data


def dominant_freq(x, dt):
    x = np.nan_to_num(x - np.nanmean(x))
    spec = np.abs(np.fft.rfft(x))
    fr = np.fft.rfftfreq(len(x), dt)
    return fr[1:][np.argmax(spec[1:])] if len(fr) > 1 else 0.0


def summarize(header, data):
    t = data[:, 0]
    dt = float(np.median(np.diff(t)))
    print(f"file: dt~{dt:.4f}s  fs~{1/dt:.0f}Hz  rows={len(t)}  end={t[-1]:.1f}s")
    print(f"{'channel':<12}{'PGA[g]/max':>12}{'fdom[Hz]':>10}")
    for i, name in enumerate(header[1:], start=1):
        col = data[:, i]
        if name.upper().startswith("A"):
            print(f"{name:<12}{np.nanmax(np.abs(col)):>12.4f}"
                  f"{dominant_freq(col, dt):>10.2f}")
        else:
            print(f"{name:<12}{np.nanmax(col):>12.2f}{'-':>10}")


def convert(path, base_channel, t0, t1):
    header, data = read_csv(path)
    if base_channel not in header:
        raise SystemExit(
            f"channel '{base_channel}' not found. Available:\n  "
            + "\n  ".join(header))

    t = data[:, 0]
    dt = np.concatenate(([np.median(np.diff(t))], np.diff(t)))
    keep = (dt < MAX_DYN_DT) & (t >= t0) & (t <= t1)   # high-rate, in-window
    if keep.sum() < 10:
        raise SystemExit("window too small / no high-rate samples in [t0,t1]")

    tt = t[keep]
    target_dt = float(np.median(np.diff(tt)))
    uni_t = np.arange(tt[0], tt[-1], target_dt)

    bi = header.index(base_channel)
    base = np.interp(uni_t, tt, data[keep, bi])
    base = base - base[0]                              # remove DC offset

    with open(MOTION_FILE, "w") as f:
        f.write(f"# LEAP base motion from {os.path.basename(path)} "
                f"channel '{base_channel}', prototype scale\n")
        f.write("# col1: time [s]   col2: acceleration [g]\n")
        for ti, ai in zip(uni_t - uni_t[0], base):
            f.write(f"{ti:.5f}\t{ai:.6f}\n")
    print(f"wrote {MOTION_FILE}  ({base.size} samples, dt={target_dt:.4f}s, "
          f"PGA={np.max(np.abs(base)):.3f} g)")

    cols = [uni_t - uni_t[0]]
    for i in range(1, len(header)):
        cols.append(np.interp(uni_t, tt, data[keep, i]))
    rec = np.column_stack(cols)
    np.savetxt(RECORDED_FILE, rec, delimiter=",",
               header=",".join(["Time (sec)"] + header[1:]), comments="")
    print(f"wrote {RECORDED_FILE}  ({rec.shape[0]} rows, {rec.shape[1]} cols) "
          f"for validation overlay")


def main():
    ap = argparse.ArgumentParser(description="LEAP CSV -> model input")
    ap.add_argument("csv", nargs="?", default=CSV_DEFAULT)
    ap.add_argument("--list", action="store_true",
                    help="print channel summary and exit")
    ap.add_argument("--base", default=BASE_CHANNEL)
    ap.add_argument("--t0", type=float, default=T_START)
    ap.add_argument("--t1", type=float, default=T_END)
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        raise SystemExit(f"CSV not found: {args.csv}\n"
                         f"Put your LEAP file in data/ or pass its path.")

    header, data = read_csv(args.csv)
    if args.list:
        summarize(header, data)
        print("\nPick the base channel with --base \"AHx (g)\" "
              "(verify against the sensor layout).")
        return
    convert(args.csv, args.base, args.t0, args.t1)


if __name__ == "__main__":
    main()
