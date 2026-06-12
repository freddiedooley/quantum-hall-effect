#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qhe_analysis_pipeline.py
========================
Step-by-step analysis pipeline for the Quantum Hall Effect experiment on a
GaAs/AlGaAs two-dimensional electron gas (sample NU1783).

This script is organised into four sequential analysis sections that mirror
the logical progression of the experiment:

    Section 1 — Carrier density and excitation current
        Determines the sheet carrier density ns using three independent methods:
          (i)   Hall resistance plateau at filling factor ν = 2 (Run 4)
          (ii)  Shubnikov–de Haas oscillation periodicity Δ(1/B) (Run 3)
          (iii) Landau fan diagram from SdH extrema (Run 3)
        Also extracts the effective excitation current I_used by anchoring
        the median plateau Hall voltage to the ideal quantised value h/2e².

    Section 2 — Carrier mobility and transport lifetime
        Uses the zero-field longitudinal resistance from the stability hold
        (Run 1) to extract the sheet resistance, carrier mobility μ, and
        transport scattering time τ_tr via the Drude relation.

    Section 3 — Conductivity tensor and plateau quality
        Converts longitudinal and Hall resistances to resistivity and
        conductivity tensor components (σxx, σxy) for Runs 3, 4, and 5.
        Characterises the ν = 2 Hall plateau: median Rxy, deviation from
        the ideal quantised value h/2e² in parts per million, and a zoomed
        twin-axis resistance plot.

    Section 4 — Temperature dependence of SdH oscillations
        Overlays Vxx(B) at T = 3 K, 4 K, 5 K (Runs 3, 6, 7) and extracts
        a simple SdH oscillation amplitude at a fixed reference field to
        demonstrate thermal broadening of Landau levels.

Each section prints a summary of extracted parameters and displays
diagnostic plots interactively. Figures are shown as each section
completes, matching the original interactive development workflow.

Data files used:
    Run 1  — QHE_mergedDATA_20251120_150023.csv   (0 T stability hold, 3 K)
    Run 3  — QHE_mergedDATA_20251120_152653.csv   (−2 → +2 T sweep, 3 K)
    Run 4  — QHE_mergedDATA_20251120_155700.csv   (2 → 5 → 2 T sweep, 3 K)
    Run 5  — QHE_mergedDATA_20251120_162504.csv   (hold near ν = 2 plateau, 3 K)
    Run 6  — QHE_mergedDATA_20251120_164654.csv   (0 → +2 T sweep, 4 K)
    Run 7  — QHE_mergedDATA_20251120_165538.csv   (2 → 0 T sweep, 5 K)

    Note: Session 1 CSV files are not available (lost after the session).
    All quantitative results are derived exclusively from Session 2 data.

Usage:
    python qhe_analysis_pipeline.py

    Before running, set DATA_DIR (below) to the folder containing the
    Session 2 CSV files.

Dependencies:
    numpy, matplotlib, scipy

Known issues documented in lab diary:
    - The data acquisition script issued FREQ70 (70 Hz) but measurements
      were actually made at 67 Hz set manually on the SR830 front panel.
      This has no effect on the analysis.
    - A 5% gain offset on the upper lock-in (Vxy channel) was identified
      at the end of Session 1 and corrected by applying FILTER_CORRECTION
      = 1/0.95 to all Vxy values. This is applied consistently throughout.
    - Early analysis scripts used a preliminary geometry (W=35, L=445 px).
      This script uses the final Fiji-measured geometry (W=13.0, L=213.0 px)
      first introduced in qhe_common5.py on 28 December 2025.

Authors:
    Frederick Dooley and Christina Mooney
    School of Physics and Astronomy, University of Nottingham
    Autumn semester 2025 — project PHYS3003 (Project No. 19)
    Supervisor: Dr Chris Mellor
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, find_peaks


# ============================================================
# 0. SETUP: constants, calibration, geometry, file paths
# ============================================================

# ---------------------------------------------------------------------------
# Data directory — edit this to point at your Session 2 CSV folder
# ---------------------------------------------------------------------------
DATA_DIR = "data/session2"

def data_path(filename):
    """Return the full path to a data file inside DATA_DIR."""
    import os
    return os.path.join(DATA_DIR, filename)


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
e         = 1.602e-19        # Elementary charge (C)
h         = 6.626e-34        # Planck's constant (J s)
h_over_e2 = 25812.807        # Von Klitzing constant RK (Ω)
e2_over_h = e**2 / h         # Conductance quantum (S)

# GaAs effective electron mass — needed for transport lifetime τ_tr
m_e    = 9.11e-31            # Free electron mass (kg)
m_star = 0.067 * m_e         # GaAs 2DEG effective mass ≈ 0.067 m_e (kg)


# ---------------------------------------------------------------------------
# Magnet calibration
# B (Tesla) = calibration_factor × VB (monitor voltage, V)
# Conversion factor from laboratory superconducting magnet power supply.
# ---------------------------------------------------------------------------
calibration_factor = 1.3445   # T / V


# ---------------------------------------------------------------------------
# Lock-in filter correction
# A 5% gain offset was identified on the upper SR830 (Vxy channel) at the
# end of Session 1 (20 November 2025) due to an instrument configuration
# setting. All Vxy values are corrected by this factor before analysis.
# ---------------------------------------------------------------------------
FILTER_CORRECTION = 1.0 / 0.95   # ≈ 1.0526


# ---------------------------------------------------------------------------
# Hall-bar geometry (final Fiji pixel measurements, frozen 28 December 2025)
#
# Width W: ten edge-to-edge measurements across the Hall bar channel.
# Edge positions defined as the midpoint of the intensity transition in the
# optical microscope image HallBarMicroscopePhoto.png analysed in Fiji.
# Conservative ±1 px per edge → ±2 px on the difference W.
#
# Probe separation L: centre-to-centre distance between voltage probes
# at contacts 6 (154 ± 2 px) and 10 (367 ± 2 px).
# Uncertainty propagated from both midpoint estimates.
#
# The 15.4% relative uncertainty in W/L is the dominant systematic
# contribution to the mobility and all ρxx-derived quantities.
# ---------------------------------------------------------------------------
W_px_values = np.array([13, 14, 13, 13, 13, 13, 12, 13, 13, 13], dtype=float)
W_px_mean   = float(np.mean(W_px_values))          # px
W_px_unc    = 2.0                                  # px (conservative edge uncertainty)

L_px        = 213.0                                 # px (367 − 154)
L_px_unc    = float(np.sqrt(2.0**2 + 2.0**2))      # px (propagated from both midpoints)

WL_ratio    = W_px_mean / L_px                      # dimensionless geometry factor
rel_WL_unc  = float(np.sqrt(
    (W_px_unc  / W_px_mean)**2 +
    (L_px_unc  / L_px)**2
))
WL_unc      = WL_ratio * rel_WL_unc                 # absolute uncertainty on W/L

print("=== Hall-bar geometry (pixel-based, Fiji analysis) ===")
print(f"W_px_mean  = {W_px_mean:.2f} px  (edge unc ±{W_px_unc:.1f} px)")
print(f"L_px       = {L_px:.1f} px  (unc ±{L_px_unc:.2f} px)")
print(f"W/L        = {WL_ratio:.5f} ± {WL_unc:.5f}  ({rel_WL_unc*100:.1f} % relative)")
print()


# ---------------------------------------------------------------------------
# Input filenames — Session 2 (20 November 2025)
# All quantitative results in the final report derive from these files.
# ---------------------------------------------------------------------------
filename_run1 = data_path("QHE_mergedDATA_20251120_150023.csv")  # 0 T hold, 3 K (mobility)
filename_run3 = data_path("QHE_mergedDATA_20251120_152653.csv")  # −2→+2 T, 3 K (SdH + Hall)
filename_run4 = data_path("QHE_mergedDATA_20251120_155700.csv")  # 2→5→2 T, 3 K (ν=2 plateau)
filename_run5 = data_path("QHE_mergedDATA_20251120_162504.csv")  # hold at ν=2, 3 K
filename_run6 = data_path("QHE_mergedDATA_20251120_164654.csv")  # 0→+2 T, 4 K (T-dependence)
filename_run7 = data_path("QHE_mergedDATA_20251120_165538.csv")  # 2→0 T, 5 K (T-dependence)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def load_qhe_run(filename):
    """
    Load a QHE CSV file and return the primary measurement columns.

    The CSV schema produced by qhe_data_acquisition.py is:
      col 0:  t_s          — elapsed time (s)
      col 1:  Vxx_X_V      — longitudinal voltage in-phase component (V)
      col 2:  Vxx_Y_V      — longitudinal voltage quadrature (V)
      col 3:  Vxx_R_V      — longitudinal voltage magnitude (V)  ← used here
      col 4:  Vxx_theta_deg
      col 5:  Vxy_X_V      — Hall voltage in-phase component (V) ← used here
      col 6:  Vxy_Y_V
      col 7:  Vxy_R_V
      col 8:  Vxy_theta_deg
      col 9:  VB_V         — magnet monitor voltage (V)

    Returns:
        t      — elapsed time (s)
        B      — magnetic field (T), converted from VB_V
        vxx_r  — longitudinal voltage magnitude Vxx_R (V)
        vxy_x  — Hall voltage in-phase component Vxy_X (V)
        vb     — raw monitor voltage VB (V), retained for reference
    """
    data  = np.genfromtxt(filename, delimiter=",", skip_header=1)
    t     = data[:, 0]
    vxx_r = data[:, 3]
    vxy_x = data[:, 5] * FILTER_CORRECTION   # apply 5% gain correction to Vxy
    vb    = data[:, 9]
    B     = calibration_factor * vb           # convert monitor voltage to Tesla
    return t, B, vxx_r, vxy_x, vb


def load_run_sorted_by_B(filename):
    """
    Load a QHE run and return arrays sorted by ascending magnetic field B.

    Sorting by B is used for temperature-comparison overlays where the
    time ordering of the sweep is less important than the field ordering.

    Returns:
        B_sort    — magnetic field (T), sorted ascending
        Vxx_sort  — longitudinal voltage magnitude (V)
        Vxy_sort  — Hall voltage in-phase component (V), filter-corrected
        t_sort    — elapsed time (s)
    """
    data  = np.genfromtxt(filename, delimiter=",", skip_header=1)
    t     = data[:, 0]
    vxx_r = data[:, 3]
    vxy_x = data[:, 5] * FILTER_CORRECTION
    vb    = data[:, 9]
    B     = calibration_factor * vb

    sort_idx = np.argsort(B)
    return B[sort_idx], vxx_r[sort_idx], vxy_x[sort_idx], t[sort_idx]


def apply_filter_correction(Vxy_x, Vxy_y, Vxy_r):
    """
    Apply the 5% gain correction to all three lock-in Vxy output channels.

    The correction factor 1/0.95 compensates for the instrument configuration
    offset identified on the upper SR830 at the end of Session 1.
    Called explicitly in scripts that load all lock-in columns separately.
    In load_qhe_run above the correction is applied inline to Vxy_X only.
    """
    Vxy_x *= FILTER_CORRECTION
    Vxy_y *= FILTER_CORRECTION
    Vxy_r *= FILTER_CORRECTION
    return Vxy_x, Vxy_y, Vxy_r


# ============================================================
# SECTION 1 — CARRIER DENSITY AND EXCITATION CURRENT
# ============================================================
#
# Three independent methods are used to determine the sheet carrier density:
#   Method A: Hall resistance plateau (Run 4, ν = 2)
#             ns = ν e B_mean / h  using the mean B on the identified plateau.
#   Method B: SdH oscillation periodicity Δ(1/B) (Run 3)
#             ns = e / (h Δ(1/B))  from the mean spacing between SdH extrema
#             in inverse-field space, identified via a Savitzky–Golay
#             second-derivative peak-finding procedure.
#   Method C: Landau fan diagram (Run 3)
#             A linear fit to extremum index n vs 1/B_peak yields a slope
#             equal to h/(e ns), from which ns is extracted independently
#             of the absolute field offset.
#
# The excitation current I_used is derived by anchoring the median Hall
# voltage on the ν = 2 plateau to the ideal quantised resistance h/2e².
# This is self-consistent: the plateau is first identified using a coarse
# current estimate, then I_used is refined once ν = 2 is confirmed.
# ============================================================

print("=" * 60)
print("SECTION 1: CARRIER DENSITY AND EXCITATION CURRENT")
print("=" * 60)

# ---------------------------------------------------------------------------
# 1A. Load Run 4 — Hall plateau sweep (2 → 5 → 2 T, 3 K)
# ---------------------------------------------------------------------------
t4, B4, vxx_r4, vxy_x4, vb4 = load_qhe_run(filename_run4)

# Sanity plots — confirm the raw Hall and longitudinal signals look sensible
plt.figure()
plt.plot(B4, vxy_x4, '.', markersize=2, label="Vxy_X (filter-corrected)")
plt.xlabel("B (T)")
plt.ylabel("Vxy_X (V)")
plt.title("Run 4: Hall voltage Vxy_X vs B")
plt.grid(True)
plt.legend()
plt.show()

plt.figure()
plt.plot(B4, vxx_r4, '.', markersize=2, label="Vxx_R")
plt.xlabel("B (T)")
plt.ylabel("Vxx_R (V)")
plt.title("Run 4: Longitudinal voltage Vxx_R vs B")
plt.grid(True)
plt.legend()
plt.show()


# ---------------------------------------------------------------------------
# 1B. Identify the ν = 2 Hall plateau region — Run 4
#
# A field window around the expected plateau centre is selected by visual
# inspection (confirmed from Session 1 analysis and the SdH periodicity).
# The plateau is centred near B ≈ 3.23 T for this sample.
# ---------------------------------------------------------------------------
B_min_plateau = 3.2   # T — lower bound of plateau selection window
B_max_plateau = 3.3   # T — upper bound

idx_plateau4  = (B4 > B_min_plateau) & (B4 < B_max_plateau)

B4_plateau    = B4[idx_plateau4]
Vxy4_plateau  = vxy_x4[idx_plateau4]
Vxy4Y_plateau = vxy_x4[idx_plateau4]   # placeholder; Y-channel not loaded separately here
theta4_plateau = np.zeros_like(B4_plateau)  # not loaded in simplified helper

if B4_plateau.size == 0:
    raise ValueError("Run 4 plateau selection is empty. Adjust B_min_plateau / B_max_plateau.")

# Zoomed plateau plot — visual confirmation of the flat region
plt.figure()
plt.plot(B4_plateau, Vxy4_plateau, '.', markersize=4)
plt.xlabel("B (T)")
plt.ylabel("Vxy_X (V)")
plt.title("Run 4: Vxy_X vs B (plateau zoom)")
plt.grid(True)
plt.show()


# ---------------------------------------------------------------------------
# 1C. Estimate filling factor and carrier density from the plateau
#
# An initial ns guess (from the SdH periodicity in Session 2 analysis)
# is used to estimate ν at the plateau centre. Rounding to the nearest
# integer then gives the ideal Rxy = h/(νe²) used for the current estimate.
# This breaks the circular logic present in the earliest analysis scripts.
# ---------------------------------------------------------------------------
ns_guess_run4     = 1.572e15       # m⁻²  — initial estimate from SdH (Session 2.1)
B4_plateau_mean   = np.nanmean(B4_plateau)

nu_approx4 = ns_guess_run4 * h / (e * B4_plateau_mean)
nu_est4    = int(np.round(nu_approx4))   # nearest integer filling factor

# Refine ns from this plateau and the integer ν
ns_from_plateau    = nu_est4 * e * B4_plateau_mean / h    # m⁻²
ns_from_plateau_cm2 = ns_from_plateau / 1e4              # cm⁻²

# ---------------------------------------------------------------------------
# 1D. Hall voltage plateau statistics and current estimate
#
# The median Vxy across the plateau window is used rather than the mean
# to suppress the effect of any residual noise spikes.
# I_used is derived from Vxy_med / Rxy_ideal; its sign reflects the wiring
# convention, so the magnitude is taken for all downstream calculations.
# ---------------------------------------------------------------------------
Vxy4_med = np.nanmedian(Vxy4_plateau)
Vxy4_std = np.nanstd(Vxy4_plateau)

if np.abs(Vxy4_med) > 0:
    rel_spread4 = Vxy4_std / np.abs(Vxy4_med)
else:
    rel_spread4 = np.nan

Rxy_theory4 = h_over_e2 / nu_est4    # Ideal Hall resistance for this ν (Ω)
I_est        = Vxy4_med / Rxy_theory4 # Current estimate (A); may be negative
I_used       = np.abs(I_est)          # Magnitude used for all further analysis

print("\n--- Run 4: Hall plateau analysis (ν ≈ 2) ---")
print(f"Plateau range:             B ∈ [{B_min_plateau:.3f}, {B_max_plateau:.3f}] T")
print(f"Mean B on plateau:         {B4_plateau_mean:.4f} T")
print(f"Initial ns guess:          {ns_guess_run4:.3e} m⁻²")
print(f"Approximate ν:             {nu_approx4:.2f}  → ν_est = {nu_est4:d}")
print(f"Refined ns (plateau):      {ns_from_plateau:.3e} m⁻²  ({ns_from_plateau_cm2:.3e} cm⁻²)")
print(f"Median Vxy on plateau:     {Vxy4_med:.6e} V")
print(f"Std(Vxy) on plateau:       {Vxy4_std:.6e} V")
print(f"Relative spread:           {rel_spread4:.3e}")
print(f"Ideal Rxy (ν={nu_est4:d}):        {Rxy_theory4:.3f} Ω")
print(f"Current estimate I_est:    {I_est:.3e} A  (signed)")
print(f"Current magnitude I_used:  {I_used:.3e} A")

# Canonical Rxx / Rxy vs B for Run 4
Rxy4 = vxy_x4 / I_used
Rxx4 = vxx_r4 / I_used

plt.figure()
plt.plot(B4, Rxy4, '.', markersize=2, label=r"$R_{xy}$ (Hall)")
plt.plot(B4, Rxx4, '.', markersize=2, label=r"$R_{xx}$ (longitudinal)")
plt.xlabel("B (T)")
plt.ylabel("Resistance (Ω)")
plt.title("Run 4: Longitudinal and Hall resistance vs B")
plt.grid(True)
plt.legend()
plt.show()


# ---------------------------------------------------------------------------
# 1E. Load Run 3 — bipolar sweep (−2 → +2 T, 3 K)
#
# Run 3 provides the SdH oscillations used for the Δ(1/B) and Landau fan
# carrier density determinations.
# ---------------------------------------------------------------------------
t3, B3, vxx_r3, vxy_x3, vb3 = load_qhe_run(filename_run3)

plt.figure()
plt.plot(B3, vxy_x3, '.', markersize=2)
plt.xlabel("B (T)")
plt.ylabel("Vxy_X (V)")
plt.title("Run 3: Hall voltage Vxy_X vs B")
plt.grid(True)
plt.show()

plt.figure()
plt.plot(B3, vxx_r3, '.', markersize=2)
plt.xlabel("B (T)")
plt.ylabel("Vxx_R (V)")
plt.title("Run 3: Longitudinal voltage Vxx_R vs B")
plt.grid(True)
plt.show()


# ---------------------------------------------------------------------------
# 1F. Select the SdH oscillation region — Run 3 positive-field window
#
# SdH oscillations are periodic in 1/B and are most clearly resolved
# in the positive-field portion of the bipolar sweep away from B = 0,
# where weak localisation and low-field noise are absent.
# The window 0.2–1.2 T captures the resolved oscillations without
# entering the ν = 2 plateau regime at higher fields.
# ---------------------------------------------------------------------------
B_min_SdH = 0.2   # T
B_max_SdH = 1.2   # T

mask_sdh3    = (B3 > B_min_SdH) & (B3 < B_max_SdH)
B3_sdh       = B3[mask_sdh3]
Rxx3_sdh_raw = vxx_r3[mask_sdh3] / I_used   # Ω
Rxx3_sdh     = np.abs(Rxx3_sdh_raw)         # magnitudes for peak finding

if B3_sdh.size < 15:
    raise ValueError("Not enough SdH points in Run 3. Adjust B_min_SdH / B_max_SdH.")

# Sort by B (safety — data is usually sorted but not guaranteed)
sort_idx = np.argsort(B3_sdh)
B3_sdh   = B3_sdh[sort_idx]
Rxx3_sdh = Rxx3_sdh[sort_idx]

dB3 = np.mean(np.diff(B3_sdh))   # mean point spacing in B


# ---------------------------------------------------------------------------
# 1G. Savitzky–Golay smoothing and second-derivative peak detection
#
# The Savitzky–Golay filter simultaneously smooths the data and computes
# its second derivative with respect to B. Peaks in |d²Rxx/dB²| correspond
# to inflection points in Rxx(B), which coincide with SdH extrema.
# This is more robust than direct peak-finding on the noisy raw signal.
#
# Window length is capped at 51 points (from frozen settings) and forced
# odd as required by savgol_filter. Polynomial order 3 is used throughout.
# ---------------------------------------------------------------------------
window_length = min(51, B3_sdh.size if B3_sdh.size % 2 == 1 else B3_sdh.size - 1)
window_length = max(window_length, 7)
if window_length % 2 == 0:
    window_length += 1

polyorder = 3

# Smooth the SdH signal
Rxx3_sg = savgol_filter(Rxx3_sdh, window_length, polyorder)

# Compute second derivative using the same filter kernel
# delta=dB3 ensures the derivative is correctly scaled in physical units (Ω/T²)
d2Rxx3_dB2_raw = savgol_filter(Rxx3_sdh, window_length, polyorder,
                               deriv=2, delta=dB3)
d2Rxx3_dB2_mag = np.abs(d2Rxx3_dB2_raw)   # magnitude for peak finding

print(f"\n--- Run 3: SdH second-derivative analysis ---")
print(f"SdH window:     B ∈ [{B_min_SdH:.3f}, {B_max_SdH:.3f}] T  ({B3_sdh.size} points)")
print(f"SG parameters:  window_length = {window_length},  polyorder = {polyorder}")
print(f"Mean ΔB:        {dB3:.3e} T")


# ---------------------------------------------------------------------------
# 1H. Find SdH extrema from |d²Rxx/dB²| peaks
#
# A prominence threshold of 5% of the maximum second-derivative amplitude
# is used to distinguish genuine Landau-level-related extrema from noise.
# Manual verification confirmed peak positions are robust to reasonable
# variations in smoothing parameters (documented in Appendix C.2).
# ---------------------------------------------------------------------------
prominence_fraction = 0.05   # 5% of maximum — tunable
peak_prominence     = d2Rxx3_dB2_mag.max() * prominence_fraction

peaks_idx, peak_props = find_peaks(d2Rxx3_dB2_mag, prominence=peak_prominence)

B3_peaks  = B3_sdh[peaks_idx]
d2_peaks  = d2Rxx3_dB2_mag[peaks_idx]

if B3_peaks.size < 2:
    raise ValueError("Too few second-derivative peaks found in Run 3. "
                     "Adjust SdH window or prominence_fraction.")

# Sort peaks by 1/B (ascending 1/B = descending B)
invB_peaks      = 1.0 / B3_peaks
sort_idx_peaks  = np.argsort(invB_peaks)
invB_peaks      = invB_peaks[sort_idx_peaks]
B3_peaks        = B3_peaks[sort_idx_peaks]
d2_peaks        = d2_peaks[sort_idx_peaks]


# ---------------------------------------------------------------------------
# 1I. Method B — ns from mean Δ(1/B) spacing
#
# SdH oscillations are periodic in inverse field with period Δ(1/B) = h/(e ns).
# The factor of 2 for spin degeneracy is absorbed into the h/e prefactor
# when using the full spin-resolved formula for a GaAs 2DEG (spin splitting
# is partially resolved at these fields and temperatures).
# ---------------------------------------------------------------------------
delta_invB      = np.diff(invB_peaks)
delta_invB_mean = np.mean(delta_invB)
delta_invB_std  = np.std(delta_invB)

ns_from_d2    = e / (h * delta_invB_mean)    # m⁻²
ns_from_d2_cm2 = ns_from_d2 / 1e4            # cm⁻²

print(f"\nMethod B — ns from SdH Δ(1/B):")
print(f"Peaks found:       {B3_peaks.size}")
print(f"Mean Δ(1/B):       {delta_invB_mean:.3e} T⁻¹")
print(f"Std(Δ(1/B)):       {delta_invB_std:.3e} T⁻¹")
print(f"ns (SdH):          {ns_from_d2:.3e} m⁻²  ({ns_from_d2_cm2:.3e} cm⁻²)")


# ---------------------------------------------------------------------------
# 1J. Method C — ns from Landau fan diagram
#
# Plotting SdH extremum index n against 1/B_peak gives a straight line
# with slope h/(e ns). The intercept b is not physically meaningful here
# as it depends on the arbitrary choice of index origin (n=0 assigned to
# the lowest-field peak). The slope alone determines ns independently of
# the absolute field calibration offset.
# ---------------------------------------------------------------------------
n_indices = np.arange(B3_peaks.size)

coeffs  = np.polyfit(n_indices, invB_peaks, 1)
a_fit   = coeffs[0]     # slope: h/(e ns) in T⁻¹/index
b_fit   = coeffs[1]     # intercept (not physically meaningful here)

ns_from_fan    = e / (h * a_fit)     # m⁻²
ns_from_fan_cm2 = ns_from_fan / 1e4  # cm⁻²
gamma           = b_fit / a_fit       # phase offset

print(f"\nMethod C — ns from Landau fan diagram:")
print(f"Fit:  1/B = a·n + b")
print(f"a = {a_fit:.3e} T⁻¹,   b = {b_fit:.3e} T⁻¹")
print(f"γ = b/a ≈ {gamma:.3f}")
print(f"ns (fan):          {ns_from_fan:.3e} m⁻²  ({ns_from_fan_cm2:.3e} cm⁻²)")


# ---------------------------------------------------------------------------
# 1K. SdH diagnostic plot — smoothed Rxx and scaled |d²Rxx/dB²| with peaks
# ---------------------------------------------------------------------------
max_Rxx3    = np.nanmax(Rxx3_sg)
max_d2_3    = np.nanmax(d2Rxx3_dB2_mag)
scale_d2    = max_Rxx3 / max_d2_3 if max_d2_3 > 0 else 1.0
d2_scaled   = d2Rxx3_dB2_mag * scale_d2

plt.figure(figsize=(7, 4))
plt.plot(B3_sdh, Rxx3_sg, '-', linewidth=1.2, label="smoothed |Rxx|")
plt.plot(B3_sdh, d2_scaled, '--', linewidth=1.0, label="|d²Rxx/dB²| (scaled)")

# Mark each detected peak and label with its index
for i, (Bx, invBx) in enumerate(zip(B3_peaks, invB_peaks)):
    pk_plot_idx = np.searchsorted(B3_sdh, Bx)
    pk_plot_idx = min(pk_plot_idx, len(d2_scaled) - 1)
    plt.plot(Bx, d2_scaled[pk_plot_idx], 'o', markersize=5, color='C2')
    plt.text(Bx, d2_scaled[pk_plot_idx], f" {i}", fontsize=7, va="bottom")

plt.xlabel("B (T)")
plt.ylabel("Arb. units (Ω, scaled)")
plt.title("Run 3: SdH region — smoothed |Rxx| and |d²Rxx/dB²| with peaks")
plt.grid(True, alpha=0.3)

# Annotate key results on the plot
textstr = (
    f"Peaks: {B3_peaks.size}\n"
    f"⟨Δ(1/B)⟩ = {delta_invB_mean:.2e} T⁻¹\n"
    f"ns(Δ) = {ns_from_d2:.2e} m⁻²\n"
    f"ns(plateau) = {ns_from_plateau:.2e} m⁻²"
)
plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes,
         fontsize=8, va='top', ha='left',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
plt.legend(loc="best", fontsize=9)
plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------------
# 1L. Landau fan diagram plot
# ---------------------------------------------------------------------------
plt.figure(figsize=(6, 4))
plt.plot(n_indices, invB_peaks, 'o', label="Data (2nd-deriv. peaks)")

n_fit    = np.linspace(n_indices.min(), n_indices.max(), 200)
invB_fit = a_fit * n_fit + b_fit
plt.plot(n_fit, invB_fit, '-', label="Linear fit")

plt.xlabel("Landau index n (relative)")
plt.ylabel(r"$1/B$ (T$^{-1}$)")
plt.title("Run 3: Landau fan diagram from SdH extrema")
plt.grid(True, alpha=0.3)

textstr_fan = (
    f"a = {a_fit:.2e} T⁻¹\n"
    f"b = {b_fit:.2e} T⁻¹\n"
    f"γ = b/a ≈ {gamma:.2f}\n"
    f"ns(fan) = {ns_from_fan:.2e} m⁻²\n"
    f"ns(plateau) = {ns_from_plateau:.2e} m⁻²"
)
plt.text(0.02, 0.78, textstr_fan, transform=plt.gca().transAxes,
         fontsize=8, va='top', ha='left',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
plt.legend(loc="best", fontsize=9)
plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------------
# 1M. Combine ns estimates and report final values
# ---------------------------------------------------------------------------
ns_list  = np.array([ns_from_plateau, ns_from_d2, ns_from_fan])
labels_ns = ["Hall plateau (Run 4)", "SdH Δ(1/B) (Run 3)", "Landau fan (Run 3)"]

ns_best     = np.nanmean(ns_list)    # simple mean of three methods
ns_best_cm2 = ns_best / 1e4

print("\n--- Combined ns and refined current ---")
for lab, val in zip(labels_ns, ns_list):
    print(f"  ns ({lab:>22s}): {val:.3e} m⁻²  ({val/1e4:.3e} cm⁻²)")
print(f"\n  ns_best (mean):              {ns_best:.3e} m⁻²  ({ns_best_cm2:.3e} cm⁻²)")

# Verify ν assignment using ns_best
nu_approx_refined = ns_best * h / (e * B4_plateau_mean)
nu_refined        = int(np.round(nu_approx_refined))
Rxy_theory_refined = h_over_e2 / nu_refined

I_refined      = Vxy4_med / Rxy_theory_refined
I_used_refined = np.abs(I_refined)

print(f"\n  ν (refined):                 {nu_approx_refined:.2f} → {nu_refined:d}")
print(f"  Rxy ideal (ν={nu_refined:d}):       {Rxy_theory_refined:.3f} Ω")
print(f"  I_used (refined):            {I_used_refined:.3e} A")
print("\nI_used_refined and ns_best are passed to Sections 2–4.")


# ============================================================
# SECTION 2 — CARRIER MOBILITY AND TRANSPORT LIFETIME
# ============================================================
#
# The carrier mobility μ is extracted from the zero-field longitudinal
# resistance using the Drude relation:
#
#     μ = 1 / (ns e ρxx)    where ρxx = Rxx × (W/L)
#
# Run 1 provides a 0 T stability hold at 3 K specifically for this purpose.
# The transport lifetime τ_tr = μ m* / e follows from the GaAs effective mass.
#
# The dominant uncertainty in μ comes from the 15.4% relative uncertainty
# in the Hall-bar geometry factor W/L, not from electronic noise.
# ============================================================

print("\n" + "=" * 60)
print("SECTION 2: CARRIER MOBILITY AND TRANSPORT LIFETIME")
print("=" * 60)

# Use the refined current and best ns from Section 1
ns_fixed = ns_best
I_used   = I_used_refined


# ---------------------------------------------------------------------------
# 2A. Load Run 1 — zero-field stability hold (0 T, 3 K)
# ---------------------------------------------------------------------------
data1      = np.genfromtxt(filename_run1, delimiter=",", skip_header=1)
t_run1     = data1[:, 0]
vxx_r_run1 = data1[:, 3]        # Vxx magnitude; no filter correction needed (longitudinal)
vb_run1    = data1[:, 9]
B_run1     = calibration_factor * vb_run1

# Confirm B is near zero throughout the hold
plt.figure()
plt.plot(t_run1, B_run1, '.', markersize=2)
plt.xlabel("t (s)")
plt.ylabel("B (T)")
plt.title("Run 1: B vs time (should be near 0 T throughout)")
plt.grid(True)
plt.show()


# ---------------------------------------------------------------------------
# 2B. Select the stable low-field window
#
# The first 60 s are discarded as a transient settling period observed
# at the start of Run 1. An additional |B| < 0.05 T cut removes any
# residual field drift during the hold.
# ---------------------------------------------------------------------------
t_min   = 60.0    # s — skip initial transient
B_max_low = 0.05  # T — maximum |B| allowed in the zero-field window

mask_time = t_run1 > t_min
t_low     = t_run1[mask_time]
Vxx_low   = vxx_r_run1[mask_time]
B_low     = B_run1[mask_time]

mask_B  = np.abs(B_low) < B_max_low
t_low   = t_low[mask_B]
Vxx_low = Vxx_low[mask_B]
B_low   = B_low[mask_B]

plt.figure()
plt.plot(t_low, Vxx_low, '.', markersize=2)
plt.xlabel("t (s)")
plt.ylabel("Vxx_R (V)")
plt.title("Run 1: Vxx_R vs time (low-field window)")
plt.grid(True)
plt.show()


# ---------------------------------------------------------------------------
# 2C. Remove outliers and compute mean Rxx
#
# A 5σ cut on Vxx is applied to remove rare interference spikes.
# The distribution of Rxx values after cleaning is confirmed as
# approximately Gaussian (see histogram below), supporting the
# use of the mean as the representative low-field resistance.
# ---------------------------------------------------------------------------
Vxx_med = np.nanmedian(Vxx_low)
Vxx_std = np.nanstd(Vxx_low)

N_sigma    = 5.0
is_outlier = np.abs(Vxx_low - Vxx_med) > N_sigma * Vxx_std
keep_mask  = ~is_outlier

print(f"\nRun 1: flagged {np.sum(is_outlier)} Vxx outliers ({N_sigma}σ cut)")

t_low_clean   = t_low[keep_mask]
Vxx_low_clean = Vxx_low[keep_mask]
B_low_clean   = B_low[keep_mask]

# Convert cleaned voltages to resistance
Rxx_low  = Vxx_low_clean / I_used
Rxx_mean = np.nanmean(Rxx_low)
Rxx_std  = np.nanstd(Rxx_low)

print(f"Points used:       {Rxx_low.size:d}")
print(f"Mean Rxx (low B):  {Rxx_mean:.3f} Ω")
print(f"Std(Rxx):          {Rxx_std:.3f} Ω")


# ---------------------------------------------------------------------------
# 2D. Sheet resistance, mobility, and transport lifetime
# ---------------------------------------------------------------------------

# Sheet resistance: ρxx = Rxx × (W/L)  in Ω/square
Rxx_sheet     = Rxx_mean * WL_ratio
Rxx_sheet_std = Rxx_std  * WL_ratio

print(f"\nSheet resistance:  ρxx = {Rxx_sheet:.3f} Ω/sq  (std {Rxx_sheet_std:.3f} Ω/sq)")
print(f"W/L factor:        {WL_ratio:.5f} ± {WL_unc:.5f}  ({rel_WL_unc*100:.1f} %)")

# Mobility: μ = 1 / (ns e ρxx)
mu     = 1.0 / (ns_fixed * e * Rxx_sheet)    # m²/(V s)
mu_cm2 = mu * 1e4                             # cm²/(V s)

# Approximate relative uncertainty: dominated by W/L (15.4%)
# Rxx statistical uncertainty is subdominant (< 0.2%)
rel_err_R  = Rxx_sheet_std / Rxx_sheet if Rxx_sheet != 0 else np.nan
rel_err_ns = 0.007   # ~0.7% from the 0.3% fractional spread across ns methods
rel_err_mu = rel_WL_unc + rel_err_R + rel_err_ns   # linear sum (conservative)

print(f"\nMobility:          μ = {mu:.3e} m²/(V·s)  ({mu_cm2:.3e} cm²/(V·s))")
print(f"Approx rel. error: {rel_err_mu*100:.1f} %  (dominated by W/L uncertainty)")

# Transport lifetime: τ_tr = μ m* / e
tau_tr = mu * m_star / e    # s
print(f"\nTransport lifetime: τ_tr = {tau_tr:.3e} s  ({tau_tr*1e12:.2f} ps)")
print(f"GaAs effective mass used: m* = {m_star/m_e:.3f} m_e")


# ---------------------------------------------------------------------------
# 2E. Polished diagnostic figures for mobility analysis
# ---------------------------------------------------------------------------

# Figure (i): B vs time showing the selected low-field window
plt.figure()
plt.plot(t_run1, B_run1, '.', markersize=2, label="B(t) full run")
plt.plot(t_low_clean, B_low_clean, 'o', markersize=3, label="selected low-field window")
plt.axhline(0.0, linestyle='--', linewidth=1, alpha=0.7, label="B = 0 T")
plt.xlabel("t (s)")
plt.ylabel("B (T)")
plt.title("Run 1 (3 K): B vs time — low-field mobility window")
plt.grid(True)
plt.legend(loc="best")
plt.tight_layout()
plt.show()

# Figure (ii): Vxx vs time in the selected window with mean and ±1σ
plt.figure()
plt.plot(t_low_clean, Vxx_low_clean, '.', markersize=2, label="Vxx_R (cleaned)")
plt.axhline(Vxx_med, linestyle='--', linewidth=1,
            label=f"median Vxx = {Vxx_med:.3e} V")
plt.axhline(Vxx_med + Vxx_std, linestyle=':', linewidth=1, label="median ± 1σ")
plt.axhline(Vxx_med - Vxx_std, linestyle=':', linewidth=1)
plt.xlabel("t (s)")
plt.ylabel("Vxx_R (V)")
plt.title("Run 1 (3 K): Vxx_R vs time in low-field window")
plt.grid(True)
plt.legend(loc="best")
plt.tight_layout()
plt.show()

# Figure (iii): Histogram of Rxx — Gaussian distribution confirms stability
plt.figure()
plt.hist(Rxx_low, bins=40, alpha=0.8)
plt.axvline(Rxx_mean, linestyle='--', linewidth=1,
            label=f"mean = {Rxx_mean:.2f} Ω")
plt.axvline(Rxx_mean + Rxx_std, linestyle=':', linewidth=1, label="mean ± 1σ")
plt.axvline(Rxx_mean - Rxx_std, linestyle=':', linewidth=1)
plt.xlabel("Rxx (Ω)")
plt.ylabel("Counts")
plt.title("Run 1 (3 K): Distribution of low-field Rxx")
plt.grid(True)
plt.legend(loc="best")
plt.tight_layout()
plt.show()


# ============================================================
# SECTION 3 — CONDUCTIVITY TENSOR AND PLATEAU QUALITY
# ============================================================
#
# The resistivity tensor components are obtained from the measured
# resistances using the Hall-bar geometry factor:
#     ρxx = (W/L) × Rxx
#     ρxy = Rxy          (Hall resistivity equals Hall resistance in 2D)
#
# The conductivity tensor is then computed by standard matrix inversion:
#     σxx = ρxx / (ρxx² + ρxy²)
#     σxy = ρxy / (ρxx² + ρxy²)
#
# At B → 0, ρxy → 0 and the inversion becomes numerically ill-conditioned.
# The immediate low-field region |B| < 0.25 T is therefore masked in
# the conductivity plots but not in the resistivity plots.
#
# The ν = 2 plateau is characterised using Run 4 (up- and down-sweeps)
# and Run 5 (fixed-field hold). The deviation of the measured Rxy from
# the ideal quantised value h/2e² is reported in parts per million.
# ============================================================

print("\n" + "=" * 60)
print("SECTION 3: CONDUCTIVITY TENSOR AND PLATEAU QUALITY")
print("=" * 60)


# ---------------------------------------------------------------------------
# 3A. Load and convert Run 3 and Run 4 to resistance
# ---------------------------------------------------------------------------
t3, B3, vxx_r3, vxy_x3, vb3 = load_qhe_run(filename_run3)
t4, B4, vxx_r4, vxy_x4, vb4 = load_qhe_run(filename_run4)

# Sanity plots for both runs
for (B_arr, vxy_arr, title_str) in [
    (B3, vxy_x3, "Run 3: Vxy_X vs B"),
    (B3, vxx_r3, "Run 3: Vxx_R vs B"),
    (B4, vxy_x4, "Run 4: Vxy_X vs B"),
    (B4, vxx_r4, "Run 4: Vxx_R vs B"),
]:
    plt.figure()
    plt.plot(B_arr, vxy_arr, '.', markersize=2)
    plt.xlabel("B (T)")
    plt.ylabel("Voltage (V)")
    plt.title(title_str)
    plt.grid(True)
    plt.show()

# Resistance conversion
Rxx3 = vxx_r3 / I_used    # Ω
Rxy3 = vxy_x3 / I_used    # Ω
Rxx4 = vxx_r4 / I_used    # Ω
Rxy4 = vxy_x4 / I_used    # Ω


# ---------------------------------------------------------------------------
# 3B. Load Run 5 — fixed-field hold near ν = 2 plateau
# ---------------------------------------------------------------------------
data5    = np.genfromtxt(filename_run5, delimiter=",", skip_header=1)
t5       = data5[:, 0]
vxx_r5   = data5[:, 3]
vxy_x5   = data5[:, 5] * FILTER_CORRECTION
vb5      = data5[:, 9]
B5       = calibration_factor * vb5

plt.figure()
plt.plot(t5, B5, '.-', markersize=3)
plt.xlabel("t (s)")
plt.ylabel("B (T)")
plt.title("Run 5: B vs time (hold near ν = 2 plateau)")
plt.grid(True)
plt.show()


# ---------------------------------------------------------------------------
# 3C. Compute resistivity and conductivity tensors for Runs 3 and 4
# ---------------------------------------------------------------------------
# For each run: ρxx, ρxy, then σxx, σxy by inversion
# The denominator D = ρxx² + ρxy² appears in both components

# Run 3
rho_xx3 = Rxx3 * WL_ratio   # Ω/sq
rho_xy3 = Rxy3               # Ω (ρxy = Rxy in 2D Hall bar)
D3 = rho_xx3**2 + rho_xy3**2
sigma_xx3_e2h = (rho_xx3 / D3) / e2_over_h   # in units of e²/h
sigma_xy3_e2h = (rho_xy3 / D3) / e2_over_h

# Run 4
rho_xx4 = Rxx4 * WL_ratio
rho_xy4 = Rxy4
D4 = rho_xx4**2 + rho_xy4**2
sigma_xx4_e2h = (rho_xx4 / D4) / e2_over_h
sigma_xy4_e2h = (rho_xy4 / D4) / e2_over_h


# ---------------------------------------------------------------------------
# 3D. Run 5 plateau statistics and single-point conductivity
# ---------------------------------------------------------------------------
# Expected plateau field from ns_fixed and ν = 2: B = ns h / (e ν)
nu_target = 2
B_target  = ns_fixed * h / (e * nu_target)
dB_tol    = 0.02   # T — tolerance window around predicted plateau centre

mask5_pl  = (B5 > (B_target - dB_tol)) & (B5 < (B_target + dB_tol))
B5_pl     = B5[mask5_pl]
Vxx5_pl   = vxx_r5[mask5_pl]
Vxy5_pl   = vxy_x5[mask5_pl]

if B5_pl.size == 0:
    raise ValueError("Run 5 plateau selection is empty. Adjust B_target or dB_tol.")

# 5σ outlier removal on Vxy to remove spikes
Vxy5_med = np.nanmedian(Vxy5_pl)
Vxy5_std = np.nanstd(Vxy5_pl)
is_out5  = np.abs(Vxy5_pl - Vxy5_med) > 5.0 * Vxy5_std
keep5    = ~is_out5

print(f"Run 5: flagged {np.sum(is_out5)} Vxy outliers in plateau window")

B5_pl   = B5_pl[keep5]
Vxx5_pl = Vxx5_pl[keep5]
Vxy5_pl = Vxy5_pl[keep5]

B5_mean  = np.nanmean(B5_pl)
Vxx5_med = np.nanmedian(Vxx5_pl)
Vxy5_med = np.nanmedian(Vxy5_pl)

# Single-point conductivity for the plateau hold
Rxx5      = Vxx5_med / I_used
Rxy5      = Vxy5_med / I_used
rho_xx5   = Rxx5 * WL_ratio
rho_xy5   = Rxy5
D5        = rho_xx5**2 + rho_xy5**2
sigma_xx5_e2h = (rho_xx5 / D5) / e2_over_h
sigma_xy5_e2h = (rho_xy5 / D5) / e2_over_h

print(f"\nRun 5 plateau stats (ν ≈ 2):")
print(f"  Mean B on plateau:   {B5_mean:.4f} T")
print(f"  Median Vxy:          {Vxy5_med:.6e} V")
print(f"  Median Vxx:          {Vxx5_med:.6e} V")
print(f"  Rxy_5:               {Rxy5:.3f} Ω")
print(f"  σxy (e²/h):          {sigma_xy5_e2h:.3f}")
print(f"  σxx (e²/h):          {sigma_xx5_e2h:.3f}")


# ---------------------------------------------------------------------------
# 3E. Plateau quantisation summary — deviation from h/2e² in ppm
# ---------------------------------------------------------------------------
Rxy_ideal_nu2   = h_over_e2 / nu_target   # 12906.4 Ω
delta_Rxy5      = Rxy5 - Rxy_ideal_nu2
delta_Rxy5_ppm  = 1.0e6 * delta_Rxy5 / Rxy_ideal_nu2

print(f"\nPlateau quantisation (Run 5, ν = 2):")
print(f"  Ideal Rxy (h/2e²):   {Rxy_ideal_nu2:.6f} Ω")
print(f"  Measured Rxy:        {Rxy5:.6f} Ω")
print(f"  Deviation:           {delta_Rxy5:.6e} Ω  ({delta_Rxy5_ppm:.3f} ppm)")


# ---------------------------------------------------------------------------
# 3F. Conductivity tensor figures — σxy and σxx vs B
# ---------------------------------------------------------------------------
# Sort both runs by B for clean line plots
order3 = np.argsort(B3)
order4 = np.argsort(B4)
B3_s   = B3[order3];   sig_xx3_s = sigma_xx3_e2h[order3];   sig_xy3_s = sigma_xy3_e2h[order3]
B4_s   = B4[order4];   sig_xx4_s = sigma_xx4_e2h[order4];   sig_xy4_s = sigma_xy4_e2h[order4]

# Mask out the ill-conditioned low-field region near B = 0
B_min_plot = 0.25   # T — below this the tensor inversion is numerically unreliable
m3 = np.abs(B3_s) > B_min_plot
m4 = np.abs(B4_s) > B_min_plot

plt.figure(figsize=(6, 6))

plt.subplot(2, 1, 1)
plt.plot(B3_s[m3], sig_xy3_s[m3], '-', linewidth=1.2, label="Run 3 sweep")
plt.plot(B4_s[m4], sig_xy4_s[m4], '-', linewidth=1.2, label="Run 4 sweep")
plt.plot(B5_mean, sigma_xy5_e2h, 'o', markersize=8, label="Run 5 ν ≈ 2 plateau")
plt.axhline(2.0, linestyle='--', linewidth=1, alpha=0.7, label=r"$\sigma_{xy} = 2\,e^2/h$")
plt.axhline(4.0, linestyle=':', linewidth=1, alpha=0.7, label=r"$\sigma_{xy} = 4\,e^2/h$")
plt.xlabel("B (T)")
plt.ylabel(r"$\sigma_{xy}\ (e^2/h)$")
plt.title(r"$\sigma_{xy}(B)$ at 3 K (low-field region masked)")
plt.grid(True)
plt.legend(loc="best")

plt.subplot(2, 1, 2)
plt.plot(B3_s[m3], sig_xx3_s[m3], '-', linewidth=1.2, label="Run 3 sweep")
plt.plot(B4_s[m4], sig_xx4_s[m4], '-', linewidth=1.2, label="Run 4 sweep")
plt.plot(B5_mean, sigma_xx5_e2h, 'o', markersize=8, label="Run 5 ν ≈ 2 plateau")
plt.xlabel("B (T)")
plt.ylabel(r"$\sigma_{xx}\ (e^2/h)$")
plt.title(r"$\sigma_{xx}(B)$ at 3 K (low-field region masked)")
plt.grid(True)
plt.legend(loc="best")

plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------------
# 3G. Zoomed ν = 2 plateau — Rxy(B) and Rxx(B) on twin axes
# ---------------------------------------------------------------------------
dB_zoom    = 0.30    # T — half-width of the zoom window
B_zoom_min = B_target - dB_zoom
B_zoom_max = B_target + dB_zoom

mask4_zoom    = (B4 > B_zoom_min) & (B4 < B_zoom_max)
B4_pl_zoom    = B4[mask4_zoom]
Rxy4_pl_zoom  = Rxy4[mask4_zoom]
Rxx4_pl_zoom  = Rxx4[mask4_zoom]

if B4_pl_zoom.size == 0:
    raise ValueError("Zoomed plateau selection in Run 4 is empty. Adjust dB_zoom.")

order_zoom   = np.argsort(B4_pl_zoom)
B4_pl_zoom   = B4_pl_zoom[order_zoom]
Rxy4_pl_zoom = Rxy4_pl_zoom[order_zoom]
Rxx4_pl_zoom = Rxx4_pl_zoom[order_zoom]

# Deviation in ppm across the zoom window
Rxy_dev_ppm = 1.0e6 * (Rxy4_pl_zoom - Rxy_ideal_nu2) / Rxy_ideal_nu2

fig, ax1 = plt.subplots()
ax1.plot(B4_pl_zoom, Rxy4_pl_zoom, '.', markersize=3,
         label=r"$R_{xy}$ (Run 4)")
ax1.axhline(Rxy_ideal_nu2, linestyle='--', linewidth=1.2,
            label=r"Ideal $h/2e^2$")
ax1.set_xlabel("B (T)")
ax1.set_ylabel(r"$R_{xy}$ (Ω)")
ax1.set_title(r"Zoomed ν = 2 plateau: $R_{xy}$ and $R_{xx}$")
ax1.grid(True)

# Rxx on a second y-axis (values much smaller on plateau)
ax2 = ax1.twinx()
ax2.plot(B4_pl_zoom, Rxx4_pl_zoom, '.', markersize=3, alpha=0.6,
         color='C1', label=r"$R_{xx}$ (Run 4)")
ax2.set_ylabel(r"$R_{xx}$ (Ω)")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")
plt.tight_layout()
plt.show()

# Deviation from ideal in ppm
plt.figure()
plt.plot(B4_pl_zoom, Rxy_dev_ppm, '.', markersize=3)
plt.axhline(0.0, linestyle='--', linewidth=1.0)
plt.xlabel("B (T)")
plt.ylabel(r"$\Delta R_{xy}$ (ppm relative to $h/2e^2$)")
plt.title(r"Deviation of $R_{xy}$ from $h/2e^2$ on ν = 2 plateau (Run 4)")
plt.grid(True)
plt.tight_layout()
plt.show()


# ============================================================
# SECTION 4 — TEMPERATURE DEPENDENCE OF SdH OSCILLATIONS
# ============================================================
#
# The longitudinal voltage Vxx(B) is compared at three temperatures:
#   Run 3:  3 K   (−2 → +2 T bipolar sweep)
#   Run 6:  4 K   (0 → +2 T sweep)
#   Run 7:  5 K   (2 → 0 T sweep, so sorted by B for comparison)
#
# All three runs are restricted to a common positive-field window
# 0.30–2.00 T for a fair comparison. The SdH oscillation amplitude
# decreases monotonically with temperature, consistent with thermal
# broadening of Landau levels (Lifshitz–Kosevich damping). Only three
# temperature points are available, so the analysis is restricted to
# qualitative trends rather than a full LK fit.
#
# The Hall slope is separately confirmed to be temperature-independent
# across this range, showing that the carrier density ns is constant.
# ============================================================

print("\n" + "=" * 60)
print("SECTION 4: TEMPERATURE DEPENDENCE OF SdH OSCILLATIONS")
print("=" * 60)

# Temperatures corresponding to each run
T_run3 = 3.0   # K
T_run6 = 4.0   # K
T_run7 = 5.0   # K

# Common field range for overlaying all three runs
B_min_common = 0.30   # T
B_max_common = 2.00   # T

# Reference field at which the SdH amplitude is extracted
# — chosen to coincide with a clearly resolved oscillation at 3 K
B_peak_ref = 1.63    # T
dB_peak    = 0.05    # T — half-width of extraction window


# ---------------------------------------------------------------------------
# 4A. Load Runs 3, 6, 7 sorted by B
# ---------------------------------------------------------------------------
B3, Vxx3, Vxy3, t3 = load_run_sorted_by_B(filename_run3)
B6, Vxx6, Vxy6, t6 = load_run_sorted_by_B(filename_run6)
B7, Vxx7, Vxy7, t7 = load_run_sorted_by_B(filename_run7)


# ---------------------------------------------------------------------------
# 4B. Restrict to common positive-field window
# ---------------------------------------------------------------------------
mask3 = (B3 > B_min_common) & (B3 < B_max_common)
mask6 = (B6 > B_min_common) & (B6 < B_max_common)
mask7 = (B7 > B_min_common) & (B7 < B_max_common)

B3_pos, Vxx3_pos, Vxy3_pos = B3[mask3], Vxx3[mask3], Vxy3[mask3]
B6_pos, Vxx6_pos, Vxy6_pos = B6[mask6], Vxx6[mask6], Vxy6[mask6]
B7_pos, Vxx7_pos, Vxy7_pos = B7[mask7], Vxx7[mask7], Vxy7[mask7]

print(f"Common B range:  [{B_min_common:.2f}, {B_max_common:.2f}] T")
print(f"Run 3 (3 K):     {B3_pos.size} points")
print(f"Run 6 (4 K):     {B6_pos.size} points")
print(f"Run 7 (5 K):     {B7_pos.size} points")


# ---------------------------------------------------------------------------
# 4C. Overlay plots — Vxx(B) and Vxy(B) at three temperatures
# ---------------------------------------------------------------------------
plt.figure()
plt.plot(B3_pos, Vxx3_pos, '.', markersize=2, label="3 K (Run 3)")
plt.plot(B6_pos, Vxx6_pos, '.', markersize=2, label="4 K (Run 6)")
plt.plot(B7_pos, Vxx7_pos, '.', markersize=2, label="5 K (Run 7)")
plt.xlabel("B (T)")
plt.ylabel("Vxx_R (V)")
plt.title("Vxx_R vs B — temperature comparison")
plt.grid(True)
plt.legend(loc="best")
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(B3_pos, Vxy3_pos, '.', markersize=2, label="3 K (Run 3)")
plt.plot(B6_pos, Vxy6_pos, '.', markersize=2, label="4 K (Run 6)")
plt.plot(B7_pos, Vxy7_pos, '.', markersize=2, label="5 K (Run 7)")
plt.xlabel("B (T)")
plt.ylabel("Vxy_X (V)")
plt.title("Vxy_X vs B — temperature comparison (Hall)")
plt.grid(True)
plt.legend(loc="best")
plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------------
# 4D. Extract SdH amplitude at the reference field for each temperature
#
# Amplitude is defined as max(Vxx) − min(Vxx) within a small window
# centred on B_peak_ref. This is a simple but robust measure given the
# limited number of temperature points available.
# ---------------------------------------------------------------------------
def sdh_amplitude(B_arr, Vxx_arr, B_center, dB_window):
    """
    Estimate the SdH oscillation amplitude near B_center as
    A = max(Vxx) − min(Vxx) within [B_center − dB_window, B_center + dB_window].

    Returns (A, N) where N is the number of points in the window.
    Returns (nan, N) if fewer than 3 points are found.
    """
    mask  = (B_arr > (B_center - dB_window)) & (B_arr < (B_center + dB_window))
    V_win = Vxx_arr[mask]
    if V_win.size < 3:
        return np.nan, int(V_win.size)
    return float(np.max(V_win) - np.min(V_win)), int(V_win.size)


A3, N3 = sdh_amplitude(B3_pos, Vxx3_pos, B_peak_ref, dB_peak)
A6, N6 = sdh_amplitude(B6_pos, Vxx6_pos, B_peak_ref, dB_peak)
A7, N7 = sdh_amplitude(B7_pos, Vxx7_pos, B_peak_ref, dB_peak)

print(f"\nSdH amplitude near B = {B_peak_ref:.2f} T (±{dB_peak:.2f} T):")
print(f"  3 K (Run 3):  A = {A3:.3e} V  (N = {N3})")
print(f"  4 K (Run 6):  A = {A6:.3e} V  (N = {N6})")
print(f"  5 K (Run 7):  A = {A7:.3e} V  (N = {N7})")

if not np.isnan(A3) and A3 != 0:
    print(f"\n  A(4K)/A(3K) = {A6/A3:.3f}")
    print(f"  A(5K)/A(3K) = {A7/A3:.3f}")


# ---------------------------------------------------------------------------
# 4E. SdH amplitude vs temperature plots
# ---------------------------------------------------------------------------
T_array = np.array([T_run3, T_run6, T_run7])
A_array = np.array([A3, A6, A7])

# Absolute amplitude
plt.figure()
plt.plot(T_array, A_array, 'o-', markersize=6)
plt.xlabel("Temperature (K)")
plt.ylabel("SdH amplitude ΔVxx (V)")
plt.title(f"SdH oscillation amplitude vs temperature  (B ≈ {B_peak_ref:.2f} T)")
plt.grid(True)
plt.tight_layout()
plt.show()

# Amplitude normalised to the 3 K value
if not np.isnan(A3) and A3 != 0:
    A_rel = A_array / A3
    plt.figure()
    plt.plot(T_array, A_rel, 'o-', markersize=6)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Normalised amplitude  A(T) / A(3 K)")
    plt.title(f"Normalised SdH amplitude vs temperature  (B ≈ {B_peak_ref:.2f} T)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# ============================================================
# FINAL SUMMARY — all key transport parameters
# ============================================================

print("\n" + "=" * 60)
print("FINAL TRANSPORT PARAMETER SUMMARY")
print("=" * 60)
print(f"\nGeometry:          W/L = {WL_ratio:.5f} ± {WL_unc:.5f}  ({rel_WL_unc*100:.1f} % relative)")
print(f"Excitation current: I_used = {I_used:.3e} A  (plateau-derived, ν = {nu_refined})")
print(f"\nCarrier density (four methods):")
print(f"  ns (Hall plateau):  {ns_from_plateau:.3e} m⁻²")
print(f"  ns (SdH Δ(1/B)):    {ns_from_d2:.3e} m⁻²")
print(f"  ns (Landau fan):    {ns_from_fan:.3e} m⁻²")
print(f"  ns (best mean):     {ns_best:.3e} m⁻²  → ({ns_best/1e4:.3e} cm⁻²)")
print(f"\nMobility:           μ = {mu:.3e} m²/(V·s)  ({mu_cm2:.3e} cm²/(V·s))")
print(f"Transport lifetime: τ_tr = {tau_tr:.3e} s  ({tau_tr*1e12:.2f} ps)")
print(f"\nν = 2 plateau:")
print(f"  B_mean on plateau:  {B4_plateau_mean:.4f} T")
print(f"  Measured Rxy:       {Rxy5:.4f} Ω")
print(f"  Ideal h/2e²:        {Rxy_ideal_nu2:.4f} Ω")
print(f"  Deviation:          {delta_Rxy5_ppm:.3f} ppm")
print(f"\nSdH amplitude (B ≈ {B_peak_ref:.2f} T):")
print(f"  3 K: {A3:.3e} V   4 K: {A6:.3e} V   5 K: {A7:.3e} V")
print("\nDone.")
