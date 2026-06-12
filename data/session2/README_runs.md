# Session 2 Data — Run Log
**Date:** 20 November 2025  
**Sample:** NU1783 — GaAs/AlGaAs two-dimensional electron gas, Hall bar geometry  
**Cryostat:** Helium insert cryostat, DEWAR 01  
**Magnet:** Superconducting solenoid, calibration B = 1.3445 × V_B (T/V)  
**Excitation:** SR830 lock-in amplifiers, f = 67 Hz, I ≈ 0.84 µA (plateau-derived)  

All quantitative results in the final report derive exclusively from this session.

---

## CSV File Schema

Every file in this directory was produced by `src/acquisition/qhe_data_acquisition.py`
and shares the following column structure:

| Column | Header | Description |
|--------|--------|-------------|
| 0 | `t_s` | Elapsed time since acquisition start (s) |
| 1 | `Vxx_X_V` | Longitudinal voltage in-phase component, V_xx X (V) |
| 2 | `Vxx_Y_V` | Longitudinal voltage quadrature component (V) |
| 3 | `Vxx_R_V` | Longitudinal voltage magnitude V_xx R (V) ← used in analysis |
| 4 | `Vxx_theta_deg` | Longitudinal lock-in phase angle (°) |
| 5 | `Vxy_X_V` | Hall voltage in-phase component V_xy X (V) ← used in analysis |
| 6 | `Vxy_Y_V` | Hall voltage quadrature component (V) |
| 7 | `Vxy_R_V` | Hall voltage magnitude (V) |
| 8 | `Vxy_theta_deg` | Hall lock-in phase angle (°) |
| 9 | `VB_V` | Magnet power supply monitor voltage (V) |

**To convert to physical quantities:**
```python
B     = 1.3445 * VB_V              # magnetic field (Tesla)
Vxy_X = Vxy_X_V * (1.0 / 0.95)    # corrected Hall voltage (5% gain offset, see below)
Rxx   = Vxx_R_V / I_used           # longitudinal resistance (Ω)
Rxy   = Vxy_X   / I_used           # Hall resistance (Ω)
```

**Important — 5% Vxy gain correction:**  
A gain offset of approximately 5% was identified on the upper SR830 lock-in
(V_xy channel) at the end of Session 1, caused by an instrument configuration
setting. All V_xy values must be multiplied by `1/0.95 ≈ 1.0526` before
analysis. This correction is applied automatically inside `load_qhe_run()`
in both analysis scripts. See Appendix A.2 of the report for details.

---

## Instrument Test Files (pre-measurement, excluded from analysis)

These six files were recorded during instrument checks before measurements
began at approximately 15:00. All have V_B near zero and are not used
in any analysis.

| Filename | Rows | Notes |
|----------|------|-------|
| `QHE_mergedDATA_20251120_103609.csv` | 26 | Morning instrument check |
| `QHE_mergedDATA_20251120_103825.csv` | 20 | Morning instrument check |
| `QHE_mergedDATA_20251120_104103.csv` | 120 | Morning instrument check |
| `QHE_mergedDATA_20251120_145719.csv` | 25 | Afternoon pre-run check |
| `QHE_mergedDATA_20251120_145758.csv` | 36 | Afternoon pre-run check |
| `QHE_mergedDATA_20251120_145905.csv` | 130 | Afternoon pre-run check |

---

## Primary Measurement Runs

### Run 1 — Zero-field stability hold
**File:** `QHE_mergedDATA_20251120_150023.csv`  
**Rows:** 1488  
**Temperature:** 3 K  
**Field:** B ≈ 0 T (stability hold)  
**Purpose:** Extract the zero-field longitudinal resistance for carrier mobility
determination. The stable analysis window is t = 60–300 s, discarding an
initial transient. Mean R_xx in the stable window: 711.749 Ω (σ = 0.93 Ω),
confirming statistical stability of the measurement.  
**Used in:** Section 2 of `qhe_analysis_pipeline.py`; Stage 3G of
`qhe_final_analysis.py`; Fig12.

---

### Run 2 — Short sweep 0 to −2 T
**File:** `QHE_mergedDATA_20251120_151803.csv`  
**Rows:** 390  
**Temperature:** 3 K  
**Field:** 0 → −2 T  
**Purpose:** Diagnostic sweep. Hall voltage phase approximately −162°
(phase-flipped relative to positive-field runs). Weak localisation visible
near B = 0. **Not used in the main QHE analysis.**

---

### Run 3 — Primary bipolar sweep ⭐
**File:** `QHE_mergedDATA_20251120_152653.csv`  
**Rows:** 1500  
**Temperature:** 3 K  
**Field:** −2 → +2 T (ramp rate 0.02 A/s)  
**Purpose:** Primary run for carrier density extraction and tensor analysis.
Provides the SdH oscillation data (positive field 0.2–1.2 T) for Δ(1/B)
and Landau fan analysis, the low-field Hall slope (|B| < 0.25 T), and the
full bipolar sweep for the resistivity/conductivity tensor figures.  
**Used in:** Sections 1, 3 of `qhe_analysis_pipeline.py`; Stages 3D, 3F, 3H
of `qhe_final_analysis.py`; Figs 10, 11, 15, 16 and temperature overlay Fig13.

---

### Run 4 — High-field plateau sweep ⭐
**File:** `QHE_mergedDATA_20251120_155700.csv`  
**Rows:** 1113  
**Temperature:** 3 K  
**Field:** +2 → +5 → +2 T (up-sweep then down-sweep)  
**Purpose:** Primary run for ν = 2 Hall plateau characterisation. The
turn-around point near 5 T separates the up-sweep (2→5 T) from the
down-sweep (5→2 T), enabling reproducibility and hysteresis checks.
The plateau-derived excitation current I_used = 8.351×10⁻⁷ A is extracted
from this run and used for all downstream resistance conversions.  
**Used in:** Sections 1, 3 of `qhe_analysis_pipeline.py`; Stages 3A–3C, 3H
of `qhe_final_analysis.py`; Figs 05, 06, 07, 08, 09, 10, 11.

---

### Run 5 — Plateau hold at ν = 2
**File:** `QHE_mergedDATA_20251120_162504.csv`  
**Rows:** 450  
**Temperature:** 3 K  
**Field:** Fixed at B ≈ 3.224 T (V_B ≈ 2.3979 V)  
**Purpose:** Fixed-field hold at the ν = 2 plateau centre for high-precision
plateau statistics. A single sharp spike early in the run is removed in
analysis. After spike removal, V_xy clusters within a relative spread of
0.006%, demonstrating excellent plateau flatness.  
**Used in:** Section 3 of `qhe_analysis_pipeline.py` (single-point
conductivity on the plateau).

---

### ⚠️ Aborted Run 6 attempt (do not use)
**File:** `QHE_mergedDATA_20251120_163959.csv`  
**Rows:** 19  
**Temperature:** 3 K  
**Field:** V_B ≈ 1.51–1.58 V  
**Purpose:** Aborted first attempt at the 4 K sweep. **This file should not
be used for temperature dependence analysis.** It is present in the directory
for completeness only.  

> **Known bug:** Early versions of the analysis scripts (including some
> intermediate consolidation scripts) incorrectly assigned this file as
> `filename_run6`. This bug is fixed in both `qhe_analysis_pipeline.py` and
> `qhe_final_analysis.py`, which correctly use `164654.csv` for Run 6.

---

### Run 6 — 4 K sweep ⭐
**File:** `QHE_mergedDATA_20251120_164654.csv`  
**Rows:** ~400  
**Temperature:** 4 K  
**Field:** 0 → +2 T  
**Purpose:** Temperature dependence comparison. SdH oscillations are weaker
than the 3 K run, consistent with thermal broadening of Landau levels.  
**Used in:** Section 4 of `qhe_analysis_pipeline.py`; Stage 3I of
`qhe_final_analysis.py`; Figs 13, 14.

---

### Run 7 — 5 K sweep ⭐
**File:** `QHE_mergedDATA_20251120_165538.csv`  
**Rows:** 336  
**Temperature:** 5 K  
**Field:** +2 → 0 T  
**Purpose:** Temperature dependence comparison at the highest measured
temperature. SdH oscillations are almost entirely suppressed, confirming
strong thermal broadening. The sweep direction is reversed relative to
Run 6 (sorted by B in analysis).  
**Used in:** Section 4 of `qhe_analysis_pipeline.py`; Stage 3I of
`qhe_final_analysis.py`; Figs 13, 14.

---

## Quick Reference — Primary Runs

| Run | File (timestamp) | T (K) | Field range | Primary purpose |
|-----|-----------------|--------|-------------|-----------------|
| 1 | 150023 | 3 | 0 T hold | Mobility |
| 3 | 152653 | 3 | −2 → +2 T | SdH, Hall slope, tensors |
| 4 | 155700 | 3 | 2 → 5 → 2 T | ν = 2 plateau, current cal. |
| 5 | 162504 | 3 | Hold ~3.22 T | Plateau statistics |
| 6 | 164654 | 4 | 0 → +2 T | Temperature dependence |
| 7 | 165538 | 5 | 2 → 0 T | Temperature dependence |

---

## Experimental Notes

**Cryostat leak:** A leak in the cryostat vacuum space prevented stable
operation at the nominal base temperature of ~2 K. Active thermal
stabilisation was used throughout Session 2, limiting the lowest
reproducible temperature to approximately 3 K. This restricted the
temperature dependence study to T = 3, 4, 5 K and reduced visibility
of higher filling-factor features, but did not prevent observation of
the ν = 2 quantum Hall signatures.

**Contact configuration:** Current injection through device pads 2 and 9;
V_xy through pads 5 and 10; V_xx through pads 6 and 4. Pads 8 and 11
had bonding failures and were not used. See Fig. 4 of the report and
Appendix D.1 for the Hall-bar geometry.

**Sample details:** GaAs/AlGaAs heterostructure, MBE growth sheet NU1783
(provided by Dr C. Mellor, University of Nottingham). Layer structure:
17 nm GaAs cap / 40 nm n-Al₀.₃₃Ga₀.₆₇As (n = 1.3×10¹⁸ cm⁻³) /
40 nm undoped AlGaAs spacer / 2DEG interface / 500 nm GaAs /
250 nm GaAs-AlGaAs superlattice / 1000 nm GaAs buffer /
semi-insulating GaAs (100) substrate.
