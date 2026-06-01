# opensees_leap — LEAP-UCD-2017 sloping-ground validation (OpenSeesPy)

A small, runnable OpenSeesPy mini-project that simulates the **LEAP-UCD-2017**
mildly **sloping ground** liquefaction / lateral-spreading experiment and
post-processes the results for validation against the DesignSafe **PRJ-1843**
data.

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
│   └── compare_with_leap.py   # overlay recorded LEAP data vs simulation
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
python input\leap_csv_to_motion.py --list
python input\leap_csv_to_motion.py --base "AH11 (g)" --t0 0 --t1 25

# 3. run the effective-stress model (gravity + dynamic stages)
python input\model_2D_slope.py

# 4. post-process -> figures and summary in report\
python postprocess\plot_results.py

# 5. (if you used real data) overlay recorded vs simulated response
python postprocess\compare_with_leap.py
```

On macOS/Linux use forward slashes (`python input/model_2D_slope.py`).

## Using your downloaded LEAP-UCD-2017 data

The processed LEAP CSV (e.g. `CU2_Motion1_Processed_4222.csv`) has columns:

```
Time (sec), AH1..AH12 (g), AV1..AV2 (g), P1..P10 (kpa)
```

- **AHx (g):** horizontal accelerometers (soil array + base/container).
- **AVx (g):** vertical accelerometers.
- **Px (kpa):** pore-pressure transducers (PPTs).

Key facts about these files (verified on `CU2_Motion1`):

- **Scale:** the data are at **prototype scale** (dominant frequency ≈ 1 Hz,
  PGA ≈ 0.15–0.23 g). **No centrifuge (1/N) scaling is required** — feed it to
  the model directly. (If you ever use a *model-scale* file, convert first:
  prototype time = N·model time, prototype accel = model accel / N, stresses
  are 1:1.)
- **Sampling:** a high-rate dynamic segment (~149 Hz, dt ≈ 0.0067 s) up to
  ~96 s, followed by a coarse post-shaking **consolidation tail** out to ~200 s.
- **Strong motion:** ~6.8–14.9 s; PPTs build up during shaking then dissipate —
  i.e. real liquefaction + reconsolidation.

Steps to use it:

1. Put the CSV in `data/`.
2. `python input\leap_csv_to_motion.py --list` and identify the **base /
   container** channel from the PRJ-1843 **sensor-layout sheet**. In the
   standard layout the base accelerometers are `AH11`/`AH12` (raw table content,
   ~3 Hz), while the soil-embedded sensors (`AH2..AH9`) show the filtered ~1 Hz
   motion. **Confirm before trusting the default `AH11`.**
3. `python input\leap_csv_to_motion.py --base "AH11 (g)" --t0 0 --t1 25` writes
   `input/base_motion.txt` (used by the model) and `data/leap_recorded.csv`
   (used by the comparison plot).
4. Run the model and `compare_with_leap.py`.

> Quantitative PPT/accel comparison needs each sensor's depth/location mapped to
> a model element — take those positions from the sensor-layout sheet and set
> them in `postprocess/compare_with_leap.py`.

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
- `report/*.png` — acceleration, excess pore-pressure ratio (rᵤ), rᵤ-with-depth
  profile, lateral-displacement, and recorded-vs-simulated comparison figures.
- `report/summary.txt` — key scalar results.

Excess pore-pressure ratio is reported as `rᵤ = 1 − σ'v(t)/σ'v0` (from element
effective stress) and cross-checked against the nodal pore-pressure excess.

## Tuning notes

- `SLOPE_DEG`, `SOIL_DEPTH`, `N_ELE_Y` and the `PDMY` dict in
  `input/model_2D_slope.py` are the main knobs.
- The PDMY02 contraction/dilation and permeability parameters control the rate
  of pore-pressure build-up and the magnitude of lateral spreading; calibrate
  them to your specific Dr and the measured LEAP response.
