# opensees_leap — LEAP-UCD-2017 sloping-ground validation (OpenSeesPy)

A small, runnable OpenSeesPy mini-project that simulates the **LEAP-UCD-2017**
mildly **sloping ground** liquefaction / lateral-spreading experiment and
validates it against the DesignSafe **PRJ-1843** centrifuge data, using **both**
the recorded **accelerations** and **pore pressures**.

The deposit (saturated, medium-dense Ottawa F-65 sand, ~5° slope) is modelled
as an **effective-stress, fully-coupled (u-p)** infinite slope using:

- **Element:** `quadUP` — the *FourNodeQuadUP* u-p element (3 DOF/node:
  `ux`, `uy`, pore pressure `p`). This replaces the plain `quad` element, whose
  2-DOF nodes are **DOF-incompatible** with the 3-DOF pore-pressure formulation
  (the error you hit earlier).
- **Material:** `PressureDependMultiYield02` (PDMY02), calibrated for
  medium-dense Ottawa F-65 sand (Dr ≈ 65 %).
- **Slope:** represented by **tilting gravity** by the slope angle so a static
  driving shear stress develops — on liquefaction this drives the down-slope
  lateral spreading LEAP measures.

## Why `FourNodeQuadUP` (`quadUP`) and not `quad`

The standard `quad` element has 2 DOF per node (`ux`, `uy`). A two-phase
effective-stress (u-p) analysis needs a **third DOF for pore pressure** at each
node. Mixing a 2-DOF element into a 3-DOF (`-ndf 3`) model triggers the DOF
compatibility error. `quadUP` (FourNodeQuadUP) is the u-p version with 3 DOF per
node, so it is the correct element for saturated, liquefiable soil.

> In OpenSees**Py** the element keyword is `quadUP` (the C++ class is
> `FourNodeQuadUP`). Using the string `"FourNodeQuadUP"` raises
> *"element type FourNodeQuadUP is unknown"*.

## Project layout

```
opensees_leap/
├── input/
│   ├── model_2D_slope.py      # main effective-stress model (gravity + dynamic)
│   ├── leap_input_motion.py   # synthetic base-motion generator / loader
│   ├── leap_csv_to_motion.py  # convert a LEAP PRJ-1843 CSV -> base_motion.txt
│   └── base_motion.txt        # (generated) base acceleration record
├── output/                    # (generated) recorder time histories *.out
├── postprocess/
│   ├── plot_results.py        # figures: accel, ru, lateral displacement
│   └── compare_with_leap.py   # recorded vs simulated: accel + pore pressure
├── report/                    # (generated) PNG figures + summary.txt
├── data/                      # LEAP-UCD-2017 data (DesignSafe PRJ-1843)
├── requirements.txt
└── README.md
```

## Quick start (Windows / PowerShell, VS Code terminal)

From `C:\opensees_leap` with your virtual environment activated:

```powershell
# 1. (once) install dependencies
pip install -r requirements.txt

# 2a. EITHER generate a synthetic ramped-sine base motion ...
python input\leap_input_motion.py

# 2b. ... OR build the base motion from your downloaded LEAP CSV (see below)
python input\leap_csv_to_motion.py --list      # inspect channels
python input\leap_csv_to_motion.py             # base = mean(AH11, AH12)

# 3. run the effective-stress model (gravity + dynamic stages)
python input\model_2D_slope.py

# 4. post-process -> figures and summary in report\
python postprocess\plot_results.py

# 5. (if you used real data) validate against recorded accel + pore pressure
python postprocess\compare_with_leap.py
```

On macOS/Linux use forward slashes (`python input/model_2D_slope.py`).

## Using your downloaded LEAP-UCD-2017 data

The processed LEAP CSV (e.g. `CU2_Motion1_Processed_4222.csv`) has columns:

```
Time (sec), AH1..AH12 (g), AV1..AV2 (g), P1..P10 (kpa)
```

- **AHx (g):** horizontal accelerometers — central array `AH1..AH4`, container
  base `AH11`/`AH12`.
- **Px (kpa):** pore-pressure transducers (excess pore pressure) — central
  array `P1..P4`, base corners `P9`/`P10`.

Key facts (verified on `CU2_Motion1`):

- **Scale:** the data are at **prototype scale** (dominant frequency ≈ 1 Hz,
  PGA ≈ 0.2 g). **No centrifuge (1/N) scaling is required** — feed it directly.
  (If you ever use a *model-scale* file, convert first: prototype time = N·model
  time, prototype accel = model accel / N, stresses are 1:1.)
- **Sampling:** a high-rate dynamic segment (~149 Hz, dt ≈ 0.0067 s) followed by
  a coarse post-shaking **consolidation tail** out to ~200 s.
- **Base input motion:** per the LEAP-UCD-2017 spec, it is the **average of
  `AH11` and `AH12`** (the two container-base accelerometers). `leap_csv_to_motion.py`
  does this by default.
- **Central array (for validation):** `P1..P4` and `AH1..AH4` sit at prototype
  depths ≈ 4, 3, 2, 1 m, with initial σ'v0 ≈ 40, 30, 20, 10 kPa — matching the
  model's stress profile. Recorded `Px` is *excess* pore pressure, so
  `ru_recorded = Px / σ'v0(depth)`.

Steps:

1. Put the CSV in `data/`.
2. `python input\leap_csv_to_motion.py --list` to inspect channels (PGA +
   dominant frequency). Confirm channel names against the PRJ-1843
   **sensor-layout sheet** for your test.
3. `python input\leap_csv_to_motion.py` writes `input/base_motion.txt`
   (base = mean of `AH11`,`AH12`) and `data/leap_recorded.csv`.
4. Run the model, then `compare_with_leap.py` to overlay **both** the recorded
   accelerations (`AH1..AH4`) and pore-pressure ratios (`P1..P4`) at matched
   depths.

> As-built sensor depths vary slightly by facility. The standard LEAP depths are
> pre-set in `postprocess/compare_with_leap.py` (`ACC_DEPTH`, `PPT_DEPTH`) — edit
> them to the exact values from the layout sheet for a precise match.

## Analysis stages

1. **Gravity, elastic** (`updateMaterialStage … -stage 0`) — establishes the
   correct initial effective-stress profile (linear, buoyant unit weight).
2. **Gravity, plastic** (`-stage 1`) — engages PDMY02 plasticity at the
   in-situ stresses.
3. **Dynamic** — base motion applied as a `UniformExcitation`, integrated with
   the L-stable **TRBDF2** integrator (good for the near-singular liquefied
   tangent; suppresses spurious acceleration spikes), `BandGeneral` solver and
   an adaptive sub-stepping recovery loop.

## Outputs

- `output/disp.out`, `accel.out`, `porepressure.out` — nodal time histories.
- `output/stress.out`, `strain.out` — element (gauss-point) time histories.
- `output/sigma_v0.out` — initial effective vertical stress per element.
- `report/acceleration.png`, `excess_pwp_ratio.png`, `ru_profile.png`,
  `lateral_displacement.png` — simulated response.
- `report/compare_accel.png`, `compare_ru_time.png`, `compare_ru_profile.png` —
  recorded-vs-simulated validation (accelerations and pore pressures).
- `report/summary.txt` — key scalar results.

Excess pore-pressure ratio is reported as `rᵤ = 1 − σ'v(t)/σ'v0` (from element
effective stress) and cross-checked against the recorded PPT excess.

## Example validation (CU2, Motion 1)

Driven by the recorded base motion (mean of `AH11`,`AH12`, PGA ≈ 0.21 g), the
model reproduces near-full liquefaction through the central array:

| depth (m) | recorded peak rᵤ | simulated peak rᵤ |
|----------:|-----------------:|------------------:|
| 4 | 1.07 | 0.97 |
| 3 | 1.01 | 1.00 |
| 2 | 1.03 | 1.00 |
| 1 | 1.19 | 0.98 |

(Recorded rᵤ slightly exceeds 1 due to transient PPT spikes and the
approximate sensor depths.)

## Tuning notes

- `SLOPE_DEG`, `SOIL_DEPTH`, `N_ELE_Y` and the `PDMY` dict in
  `input/model_2D_slope.py` are the main knobs.
- The PDMY02 contraction/dilation and permeability parameters control the rate
  of pore-pressure build-up and the magnitude of lateral spreading; calibrate
  them to your specific Dr and the measured LEAP response.
