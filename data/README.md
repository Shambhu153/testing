# data/

Place the **LEAP-UCD-2017** experimental data here for validation.

- Source: DesignSafe project **PRJ-1843** (LEAP-UCD-2017).
- Dataset used: **LEAP-2017 sloping ground** (mildly sloping, ~5°, saturated
  Ottawa F-65 sand, rigid container, ramped ~1 Hz base motion).

What to download for validation:

- The **achieved base motion** for the specific centrifuge test you are
  validating against. Save it as `input/base_motion.txt` (two columns:
  `time[s]  accel[g]`) so `input/model_2D_slope.py` uses the real record
  instead of the synthetic ramped sine.
- The recorded **acceleration**, **pore-pressure transducer (PPT)** and
  **surface lateral-displacement** time histories, to overlay against the
  simulated results from `postprocess/plot_results.py`.

Large raw files (`*.csv`, `*.zip`, `*.h5`) are git-ignored; keep them local.

> Note: LEAP data are at model (centrifuge) scale. Convert to prototype scale
> (consistent with the model's metres / seconds / kPa units) before comparing.
