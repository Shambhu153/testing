"""
LEAP-UCD-2017 sloping-ground effective-stress model (OpenSeesPy).

2D plane-strain, fully-coupled (u-p) effective-stress analysis of a mildly
sloping, fully saturated Ottawa F-65 sand deposit, idealised as an infinite
slope. This is the standard simplified model used to validate LEAP-UCD-2017
lateral-spreading / liquefaction response.

Element : quadUP  (the FourNodeQuadUP u-p element; 3 DOF/node = ux, uy, pore p)
Material: PressureDependMultiYield02 (PDMY02), calibrated for medium-dense
          Ottawa F-65 sand (Dr ~ 65%).

The slope is represented by tilting gravity by ``SLOPE_DEG`` so a static
driving shear stress develops; on liquefaction this produces the down-slope
lateral spreading that LEAP measures. Lateral boundaries are periodic
(equalDOF at each elevation) to reproduce 1D-style shear of an infinite slope.

Stages
------
1. Gravity, elastic   (updateMaterialStage stage 0)   -> initial stresses
2. Gravity, plastic   (updateMaterialStage stage 1)   -> let PDMY02 engage
3. Dynamic            (UniformExcitation base motion)  -> shaking response

The dynamic stage uses the TRBDF2 integrator (L-stable, with high-frequency
numerical dissipation) which is well suited to the near-singular tangent of a
liquefied soil and suppresses spurious acceleration spikes.

Run:  python input/model_2D_slope.py
"""

from __future__ import annotations

import os
import sys
import math

import openseespy.opensees as ops

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "output")
os.makedirs(OUT, exist_ok=True)

# Make the input-motion helper importable regardless of CWD.
sys.path.insert(0, HERE)
import leap_input_motion as motion  # noqa: E402

# ============================================================================
# CONFIGURATION  (units: m, s, kN, kPa, ton(=Mg))
# ============================================================================
G = 9.81             # gravity [m/s^2]

# --- geometry (prototype scale) ---
SOIL_DEPTH = 4.0     # vertical soil depth [m]
N_ELE_Y = 20         # number of elements over the depth
ELE_WIDTH = 1.0      # element width [m] (infinite-slope: 1 column)
THICK = 1.0          # out-of-plane thickness [m]
SLOPE_DEG = 5.0      # ground slope [deg]

# --- water / two-phase ---
RHO_SOLID = 2.0      # saturated soil mass density [ton/m^3]
RHO_WATER = 1.0      # fluid mass density [ton/m^3]
FLUID_BULK = 2.2e6   # combined undrained fluid bulk modulus [kPa]
PERM = 1.0e-5        # permeability coefficient [m/s] (Ottawa F-65 ~ 1e-5)

# --- PDMY02 (medium-dense Ottawa F-65 sand, Dr ~ 65%) ---
PDMY = dict(
    nd=2,
    rho=RHO_SOLID,
    refShearModul=9.0e4,    # Gr at refPress [kPa]
    refBulkModul=2.2e5,     # Br at refPress [kPa]
    frictionAng=33.0,       # [deg]
    peakShearStra=0.1,
    refPress=80.0,          # [kPa]
    pressDependCoe=0.5,
    PTAng=26.0,             # phase-transformation angle [deg]
    contrac1=0.067,
    contrac3=0.23,
    dilat1=0.06,
    dilat3=0.27,
)

# --- analysis control ---
GRAV_STEPS_ELASTIC = 100
GRAV_STEPS_PLASTIC = 100
GRAV_DT = 5.0e2          # large pseudo-time for static-like consolidation
DAMP_RATIO = 0.03        # target Rayleigh damping ratio
MAT_TAG = 1


# ============================================================================
# MODEL BUILDERS
# ============================================================================
def _node_left(j):
    return 100 + j


def _node_right(j):
    return 200 + j


def _ele(j):
    return 1000 + j


def build_model(slope_deg=SLOPE_DEG, perm=PERM):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    dy = SOIL_DEPTH / N_ELE_Y
    theta = math.radians(slope_deg)

    # nodes: two vertical columns (left x=0, right x=ELE_WIDTH)
    for j in range(N_ELE_Y + 1):
        y = j * dy
        ops.node(_node_left(j), 0.0, y)
        ops.node(_node_right(j), ELE_WIDTH, y)

    # base: fully fixed in translation (pore pressure free -> impermeable base)
    ops.fix(_node_left(0), 1, 1, 0)
    ops.fix(_node_right(0), 1, 1, 0)

    # periodic lateral boundary: tie left & right translation at every level
    for j in range(1, N_ELE_Y + 1):
        ops.equalDOF(_node_left(j), _node_right(j), 1, 2)

    # drained free surface: pore pressure = 0 at the top
    ops.fix(_node_left(N_ELE_Y), 0, 0, 1)
    ops.fix(_node_right(N_ELE_Y), 0, 0, 1)

    # material
    ops.nDMaterial("PressureDependMultiYield02", MAT_TAG, *PDMY.values())

    # body force: gravity tilted by the slope angle -> static driving shear
    b1 = G * math.sin(theta)
    b2 = -G * math.cos(theta)

    for j in range(N_ELE_Y):
        ops.element("quadUP", _ele(j),
                    _node_left(j), _node_right(j),
                    _node_right(j + 1), _node_left(j + 1),
                    THICK, MAT_TAG, FLUID_BULK, RHO_WATER, perm, perm, b1, b2)

    return dy, b1, b2


def run_gravity():
    """Two-stage gravity: elastic then plastic, with high numerical damping."""
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.test("NormDispIncr", 1.0e-5, 30, 0)
    ops.algorithm("Newton")
    gamma, beta = 1.5, (1.5 + 0.5) ** 2 / 4.0   # large numerical damping
    ops.integrator("Newmark", gamma, beta)
    ops.analysis("Transient")

    ops.updateMaterialStage("-material", MAT_TAG, "-stage", 0)
    if ops.analyze(GRAV_STEPS_ELASTIC, GRAV_DT) != 0:
        raise RuntimeError("gravity (elastic) failed")

    ops.updateMaterialStage("-material", MAT_TAG, "-stage", 1)
    if ops.analyze(GRAV_STEPS_PLASTIC, GRAV_DT) != 0:
        raise RuntimeError("gravity (plastic) failed")

    ops.setTime(0.0)
    ops.wipeAnalysis()


def initial_effective_sigma_v():
    """Effective vertical stress (kPa, +compression) at each element gauss pt."""
    return {_ele(j): abs(ops.eleResponse(_ele(j), "material", 1, "stress")[1])
            for j in range(N_ELE_Y)}


def add_recorders():
    """Recorders write plain-text time histories to output/."""
    node_ids = [_node_left(j) for j in range(N_ELE_Y + 1)]
    ele_ids = [_ele(j) for j in range(N_ELE_Y)]

    ops.recorder("Node", "-file", os.path.join(OUT, "disp.out"), "-time",
                 "-node", *node_ids, "-dof", 1, 2, "disp")
    ops.recorder("Node", "-file", os.path.join(OUT, "accel.out"), "-time",
                 "-node", *node_ids, "-dof", 1, 2, "accel")
    # pore pressure is the value of DOF 3 -> record as 'disp'
    ops.recorder("Node", "-file", os.path.join(OUT, "porepressure.out"), "-time",
                 "-node", *node_ids, "-dof", 3, "disp")
    ops.recorder("Element", "-file", os.path.join(OUT, "stress.out"), "-time",
                 "-ele", *ele_ids, "material", 1, "stress")
    ops.recorder("Element", "-file", os.path.join(OUT, "strain.out"), "-time",
                 "-ele", *ele_ids, "material", 1, "strain")


def _set_integrator(kind="TRBDF2"):
    if kind == "TRBDF2":
        ops.integrator("TRBDF2")
    else:
        ops.integrator("Newmark", 0.6, 0.3025)


def _setup_dynamic_solver():
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")        # robust for the near-singular liquefied tangent
    ops.test("NormDispIncr", 1.0e-4, 30, 0)
    ops.algorithm("Newton")
    _set_integrator("TRBDF2")        # L-stable, damps spurious high-freq accel
    ops.analysis("Transient")


def _adaptive_step(dt):
    """Advance one base step of size dt, recovering with sub-stepping/algorithms."""
    if ops.analyze(1, dt) == 0:
        return True

    for algo in (("Newton",), ("KrylovNewton",),
                 ("ModifiedNewton",), ("NewtonLineSearch",)):
        for nsub in (2, 4, 8, 16):
            ops.algorithm(*algo)
            ops.test("NormDispIncr", 5.0e-4, 100, 0)
            ok = ops.analyze(nsub, dt / nsub)
            if ok == 0:
                ops.algorithm("Newton")
                ops.test("NormDispIncr", 1.0e-4, 30, 0)
                return True
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1.0e-4, 30, 0)
    return False


def run_dynamic(dt, accel_ms2):
    """Apply base motion as a UniformExcitation and integrate (adaptive)."""
    nstep = len(accel_ms2)

    ops.timeSeries("Path", 2, "-dt", dt, "-values", *accel_ms2)
    ops.pattern("UniformExcitation", 2, 1, "-accel", 2)

    # Rayleigh damping ~DAMP_RATIO at two representative frequencies
    f1, f2 = 0.5, 8.0
    w1, w2 = 2 * math.pi * f1, 2 * math.pi * f2
    a0 = DAMP_RATIO * 2.0 * w1 * w2 / (w1 + w2)
    a1 = DAMP_RATIO * 2.0 / (w1 + w2)
    ops.rayleigh(a0, a1, 0.0, 0.0)

    _setup_dynamic_solver()

    fails = 0
    for step in range(nstep):
        if not _adaptive_step(dt):
            fails += 1
        if (step + 1) % 400 == 0:
            print(f"  dynamic step {step + 1}/{nstep}  "
                  f"t={ops.getTime():.2f}s  (give-up steps: {fails})")
    if fails:
        print(f"  note: {fails} step(s) marched without full convergence")


def main():
    print("Building LEAP sloping-ground model ...")
    build_model(slope_deg=SLOPE_DEG, perm=PERM)

    print("Running gravity (elastic + plastic) ...")
    run_gravity()
    sig_v0 = initial_effective_sigma_v()
    mid = _ele(N_ELE_Y // 2)
    print(f"  initial sigma'v: bottom = {sig_v0[_ele(0)]:.2f} kPa, "
          f"mid = {sig_v0[mid]:.2f} kPa")

    with open(os.path.join(OUT, "sigma_v0.out"), "w") as f:
        f.write("# eleTag\tsigma_v0_eff[kPa]\n")
        for j in range(N_ELE_Y):
            f.write(f"{_ele(j)}\t{sig_v0[_ele(j)]:.6f}\n")

    print("Loading base motion ...")
    dt, accel = motion.load()
    print(f"  {len(accel)} samples, dt={dt}s, "
          f"PGA={max(abs(a) for a in accel) / G:.3f} g")

    add_recorders()
    print("Running dynamic stage ...")
    run_dynamic(dt, accel)

    top = _node_left(N_ELE_Y)
    print("Done.")
    print(f"  residual surface lateral disp = {ops.nodeDisp(top, 1):.4f} m")
    ops.wipe()


if __name__ == "__main__":
    main()
