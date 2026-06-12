#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qhe_final_analysis.py
=====================
Definitive one-pass analysis and figure export script for the Quantum Hall
Effect experiment on a GaAs/AlGaAs two-dimensional electron gas (sample NU1783).

This script reproduces every figure in the submitted laboratory report
(QHE_Report_Submission.pdf, January 2026) in a single run. It is the
authoritative version of the analysis: all input parameters are frozen at
the top of the script, and all downstream calculations follow from them
without manual intervention.

Pipeline overview:
------------------
The script runs as a linear pipeline through eight numbered stages, populating
two shared dictionaries:
    results  — all extracted scalar quantities and statistics
    clean    — cleaned and converted array data for each run

Stage 3A  Load Run 4 (2→5→2 T, 3 K); coarse current estimate for plateau search.
Stage 3B  Despike Run 4; split into up- and down-sweeps; find the ν = 2 plateau
          window in each sweep direction using an objective moving-slope criterion.
Stage 3C  Derive the final excitation current I_used and ns from the plateau.
Stage 3D  Load Run 3 (−2→+2 T, 3 K); Savitzky–Golay second-derivative peak
          detection for SdH extrema; compute ns from Δ(1/B) and Landau fan.
Stage 3E  Combine ns estimates; verify ν assignment; store ns_best.
Stage 3F  Low-field Hall slope from Run 3 (|B| < 0.25 T, ρ-space); cross-check ns.
Stage 3G  Run 1 (0 T hold, 3 K); clean; extract sheet resistance, mobility, τ_tr.
Stage 3H  Full cleaned arrays and resistivity/conductivity tensors for Runs 3 and 4.
Stage 3I  Temperature dependence: Runs 3, 6, 7 (3 K, 4 K, 5 K); SdH amplitudes.

Figures produced:
-----------------
Main report figures (exported to report_figures/):
    Fig05_Overview_Rxy_Rxx_T3K
    Fig06_CarrierDensity_Comparison_BarChart
    Fig07_Nu2_Rxy_Zoom
    Fig08_Nu2_Rxx_Suppression_Zoom
    Fig09_Nu2_UpDown_Deviation_From_Quantised
    Fig10_ResistivityTensor_rho_xx_rho_xy
    Fig11_ConductivityTensor_sigma_xx_sigma_xy_Masked
    Fig12_Mobility_LowField_Diagnostics
    Fig13_Vxx_Overlay_T3_T4_T5K
    Fig14_SdH_Amplitude_vs_Temperature
    Fig15_LandauFan_SdH_Index_vs_InvB
    Fig16_Rxx_Raw_vs_SecondDerivative_PeakDetection

Appendix figures (exported to appendix_figures/):
    FigA1_SdH_window

All figures are saved as both .png (300 dpi) and .pdf.

Data files used:
----------------
    Run 1  — QHE_mergedDATA_20251120_150023.csv   (0 T stability hold, 3 K)
    Run 3  — QHE_mergedDATA_20251120_152653.csv   (−2 → +2 T sweep, 3 K)
    Run 4  — QHE_mergedDATA_20251120_155700.csv   (2 → 5 → 2 T sweep, 3 K)
    Run 6  — QHE_mergedDATA_20251120_164654.csv   (0 → +2 T sweep, 4 K)
    Run 7  — QHE_mergedDATA_20251120_165538.csv   (2 → 0 T sweep, 5 K)

    Note: Run 5 (plateau hold, 163959.csv) was an aborted attempt and is
    not used here. The correct Run 6 file is 164654.csv — see lab diary
    entry 20 November 2025 for the full run log.

Usage:
------
    python qhe_final_analysis.py

    Set DATA_DIR (below) to the folder containing the Session 2 CSV files.
    Output directories report_figures/ and appendix_figures/ are created
    automatically.

Dependencies:
    numpy, matplotlib, scipy

Known issues documented in lab diary:
--------------------------------------
    - The acquisition script issued FREQ70 (70 Hz) but measurements were
      made at 67 Hz set manually on the SR830 front panel. No effect on
      the analysis.
    - A 5% gain offset on the upper SR830 (Vxy channel) was identified at
      the end of Session 1. All Vxy values are corrected by FILTER_CORRECTION
      = 1/0.95 before any analysis. See Appendix A.2 of the report.
    - The low-field conductivity inversion is masked for |B| < 0.25 T where
      ρxy → 0 causes numerical instability. See Appendix A.3 of the report.
    - Early analysis scripts used a preliminary geometry (W=35, L=445 px).
      This script uses the final Fiji-measured geometry (W=13.0, L=213.0 px)
      introduced in qhe_common5.py on 28 December 2025.

Authors:
    Frederick Dooley and Christina Mooney
    School of Physics and Astronomy, University of Nottingham
    Autumn semester 2025 — project PHYS3003 (Project No. 19)
    Supervisor: Dr Chris Mellor
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import savgol_filter, find_peaks, medfilt


# ============================================================
# 0) SETUP: constants, calibration, file paths, plot style
# ============================================================

# ---------------------------------------------------------------------------
# Data directory — edit this to point at your Session 2 CSV folder
# ---------------------------------------------------------------------------
DATA_DIR = "data/session2"

def data_path(filename):
    """Return the full path to a data file inside DATA_DIR."""
    return os.path.join(DATA_DIR, filename)


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
e         = 1.602e-19        # Elementary charge (C)
h         = 6.626e-34        # Planck's constant (J s)
h_over_e2 = 25812.807        # Von Klitzing constant RK (Ω)
e2_over_h = e**2 / h         # Conductance quantum (S)
m_e       = 9.11e-31         # Free electron mass (kg)
m_star    = 0.067 * m_e      # GaAs 2DEG effective mass ≈ 0.067 m_e (kg)


# ---------------------------------------------------------------------------
# Magnet calibration
# B (Tesla) = calibration_factor × VB (monitor voltage, V)
# ---------------------------------------------------------------------------
calibration_factor = 1.3445   # T / V


# ---------------------------------------------------------------------------
# Lock-in filter correction
# 5% gain offset identified on the upper SR830 (Vxy channel) at end of
# Session 1 (20 November 2025). Corrected uniformly before all analysis.
# ---------------------------------------------------------------------------
FILTER_CORRECTION = 1.0 / 0.95   # ≈ 1.0526


# ---------------------------------------------------------------------------
# Input filenames — Session 2 (20 November 2025)
# ---------------------------------------------------------------------------
filename_run4 = data_path("QHE_mergedDATA_20251120_155700.csv")  # 2→5→2 T, 3 K (plateau)
filename_run3 = data_path("QHE_mergedDATA_20251120_152653.csv")  # −2→+2 T, 3 K (SdH + tensors)
filename_run1 = data_path("QHE_mergedDATA_20251120_150023.csv")  # 0 T hold, 3 K (mobility)
filename_run6 = data_path("QHE_mergedDATA_20251120_164654.csv")  # 0→+2 T, 4 K (T-dep)
filename_run7 = data_path("QHE_mergedDATA_20251120_165538.csv")  # 2→0 T, 5 K (T-dep)


# ---------------------------------------------------------------------------
# Plot style — applied globally for all report figures
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "font.size"        : 10,
    "axes.linewidth"   : 1.0,
    "axes.grid"        : True,
    "grid.linestyle"   : "--",
    "grid.linewidth"   : 0.6,
    "grid.alpha"       : 0.7,
    "xtick.direction"  : "in",
    "ytick.direction"  : "in",
    "xtick.top"        : True,
    "ytick.right"      : True,
    "figure.dpi"       : 150,
})


# ============================================================
# 1) FROZEN INPUTS
# ============================================================
# All analysis parameters are declared here. Edit ONLY this block
# to adjust the analysis; nothing else in the script should require
# manual changes for a rerun on the same dataset.

frozen = {}

# ---- Hall-bar geometry (pixel-based, Fiji analysis, frozen 28 Dec 2025) ----
# Width W: ten edge-to-edge measurements across the Hall bar channel.
# Edge defined as midpoint of intensity transition in the optical image.
# Conservative ±1 px per edge → ±2 px uncertainty on the difference W.
# Probe separation L: centre-to-centre of contacts 6 (154 px) and 10 (367 px).
W_px_values         = np.array([13, 14, 13, 13, 13, 13, 12, 13, 13, 13], dtype=float)
frozen["W_px_mean"] = float(np.mean(W_px_values))
frozen["W_px_std"]  = float(np.std(W_px_values, ddof=1))
frozen["W_px_unc"]  = 2.0                                   # px (edge-picking unc)
frozen["L_px"]      = 213.0                                  # px (367 − 154)
frozen["L_px_unc"]  = float(np.sqrt(2.0**2 + 2.0**2))       # px (propagated)

frozen["WL_ratio"]   = frozen["W_px_mean"] / frozen["L_px"]
frozen["WL_rel_unc"] = float(np.sqrt(
    (frozen["W_px_unc"] / frozen["W_px_mean"])**2 +
    (frozen["L_px_unc"] / frozen["L_px"])**2
))
frozen["WL_unc"] = frozen["WL_ratio"] * frozen["WL_rel_unc"]

# ---- Plateau detection settings ----
# The ν = 2 plateau window is found objectively using a moving-slope
# criterion: the longest contiguous field interval in which |dRxy/dB|
# remains below slope_threshold, allowing small gaps of gap_tol_pts points.
frozen["nu_plateau"]       = 2
frozen["B_center_guess"]   = 3.249    # T — approximate plateau centre
frozen["search_halfwidth"] = 0.25     # T — search window half-width
frozen["slope_win_pts"]    = 21       # number of points for local slope fit
frozen["slope_threshold"]  = 200.0   # Ω/T — maximum allowed local slope
frozen["gap_tol_pts"]      = 3        # allowed gap in flat region (points)
frozen["min_pts_plateau"]  = 8        # minimum plateau width (points)

# ---- SdH analysis settings ----
frozen["B_min_SdH"]      = 0.20   # T — lower bound of SdH window
frozen["B_max_SdH"]      = 1.20   # T — upper bound
frozen["sg_window_max"]  = 51     # maximum SG filter window length (points)
frozen["sg_polyorder"]   = 3      # SG polynomial order
frozen["peak_prom_frac"] = 0.05   # peak prominence as fraction of max |d²Rxx/dB²|

# ---- Low-field Hall window for cross-check ns extraction ----
frozen["B0_lowfield"] = 0.25   # T — |B| < B0 used for low-field Hall slope

# ---- Run 1 (mobility) settings ----
frozen["t_min_run1"]    = 60.0   # s — skip initial transient
frozen["B_max_low_run1"] = 0.05  # T — maximum |B| for zero-field selection
frozen["N_sigma_run1"]  = 5.0    # outlier rejection threshold (σ)

# ---- Temperature dependence settings ----
frozen["B_min_common_T"] = 0.30    # T — common positive-field window lower bound
frozen["B_max_common_T"] = 2.00    # T — upper bound
frozen["B_peak_ref"]     = 1.63    # T — reference field for SdH amplitude extraction
frozen["dB_peak"]        = 0.05    # T — half-width of amplitude extraction window
frozen["T_run3"]         = 3.0     # K
frozen["T_run6"]         = 4.0     # K
frozen["T_run7"]         = 5.0     # K

print("=== FROZEN INPUTS FOR FINAL RESULTS ===")
print(f"W/L used:   {frozen['WL_ratio']:.5f} ± {frozen['WL_unc']:.5f}  ({frozen['WL_rel_unc']*100:.1f} %)")
print(f"ν used:     {frozen['nu_plateau']}")
print(
    f"Slope rule: win={frozen['slope_win_pts']} pts, "
    f"|dRxy/dB| ≤ {frozen['slope_threshold']} Ω/T, "
    f"gap_tol = {frozen['gap_tol_pts']} pts"
)
print()


# ============================================================
# 2) HELPERS: I/O, cleaning, analysis utilities
# ============================================================

def ensure_dir(path: str) -> None:
    """Create directory (and parents) if it does not already exist."""
    if not os.path.exists(path):
        os.makedirs(path)


def savefig_both(fig, outdir: str, name_base: str, dpi: int = 300) -> None:
    """
    Save a matplotlib figure as both PNG and PDF in outdir.
    Calls tight_layout() before saving, then closes the figure.
    """
    ensure_dir(outdir)
    png_path = os.path.join(outdir, f"{name_base}.png")
    pdf_path = os.path.join(outdir, f"{name_base}.pdf")
    fig.tight_layout()
    fig.savefig(png_path, dpi=dpi)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


def require_keys(d: dict, keys: list, name: str = "dict") -> None:
    """Raise KeyError if any of keys are missing from dict d."""
    missing = [k for k in keys if k not in d]
    if missing:
        raise KeyError(f"Missing keys in {name}: {missing}")


def load_qhe_run(filename: str):
    """
    Load a QHE CSV file and return the primary measurement columns.

    CSV schema (produced by qhe_data_acquisition.py):
      col 0  t_s          — elapsed time (s)
      col 1  Vxx_X_V      — longitudinal in-phase (V)
      col 2  Vxx_Y_V      — longitudinal quadrature (V)
      col 3  Vxx_R_V      — longitudinal magnitude (V)   ← used
      col 4  Vxx_theta_deg
      col 5  Vxy_X_V      — Hall in-phase (V)            ← used (filter-corrected)
      col 6  Vxy_Y_V
      col 7  Vxy_R_V
      col 8  Vxy_theta_deg
      col 9  VB_V         — magnet monitor voltage (V)   ← used

    Returns: t, B, vxx_r, vxy_x, vb
    """
    data  = np.genfromtxt(filename, delimiter=",", skip_header=1)
    t     = data[:, 0]
    vxx_r = data[:, 3]
    vxy_x = data[:, 5] * FILTER_CORRECTION   # apply 5% gain correction inline
    vb    = data[:, 9]
    B     = calibration_factor * vb
    return t, B, vxx_r, vxy_x, vb


def load_run_sorted(filename: str):
    """
    Load a QHE run and return (B, Vxx_R, Vxy_X, t) sorted by ascending B.
    Used for temperature-comparison overlays where field ordering matters.
    """
    data  = np.genfromtxt(filename, delimiter=",", skip_header=1)
    t     = data[:, 0]
    vxx_r = data[:, 3]
    vxy_x = data[:, 5] * FILTER_CORRECTION
    vb    = data[:, 9]
    B     = calibration_factor * vb
    idx   = np.argsort(B)
    return B[idx], vxx_r[idx], vxy_x[idx], t[idx]


def despike(arr, kernel=21, nsig=5.0):
    """
    Median-filter despiking: removes isolated transient spikes from an array.

    Method: subtract a running median (window = kernel points) to obtain
    residuals; flag points where |residual| > nsig × std(residuals) as spikes
    and replace them with NaN.

    Returns: clean (NaN at spike locations), keep_mask (True where kept).
    Passes through unchanged if sigma is zero or non-finite.
    """
    arr = np.asarray(arr, dtype=float)
    if kernel % 2 == 0:
        kernel += 1   # medfilt requires odd kernel

    med   = medfilt(arr, kernel_size=kernel)
    resid = arr - med
    sigma = np.nanstd(resid)

    if (not np.isfinite(sigma)) or sigma <= 0:
        return arr.copy(), np.ones_like(arr, dtype=bool)

    keep  = np.abs(resid) < nsig * sigma
    clean = arr.copy()
    clean[~keep] = np.nan
    return clean, keep


def resistances_to_tensors(Rxx, Rxy, WL_ratio):
    """
    Convert measured Hall-bar resistances to resistivity and conductivity
    tensor components using the standard 2D relations:

        ρxx = (W/L) × Rxx
        ρxy = Rxy               (Hall resistivity equals Hall resistance in 2D)

    Conductivity tensor by inversion of the resistivity tensor:
        σxx = ρxx / (ρxx² + ρxy²)
        σxy = ρxy / (ρxx² + ρxy²)

    The denominator is guarded against division by zero using a small epsilon.
    Near B = 0, ρxy → 0 and the inversion becomes ill-conditioned; those
    points are handled by masking in the figure functions rather than here.

    Returns: rho_xx, rho_xy, sigma_xx, sigma_xy  (all in SI units)
    """
    rho_xx = np.asarray(Rxx, dtype=float) * float(WL_ratio)
    rho_xy = np.asarray(Rxy, dtype=float)

    denom      = rho_xx**2 + rho_xy**2
    eps        = 1e-30
    denom_safe = np.where(denom > eps, denom, np.nan)   # avoid exact zero division

    sigma_xx = rho_xx / denom_safe
    sigma_xy = rho_xy / denom_safe
    return rho_xx, rho_xy, sigma_xx, sigma_xy


def sdh_amplitude(B_arr, Vxx_arr, B_center, dB_window):
    """
    Estimate SdH amplitude near B_center as max(Vxx) − min(Vxx) within
    [B_center − dB_window, B_center + dB_window].

    Returns (A, N) where N is the number of points in the window.
    Returns (nan, N) if fewer than 3 points are found.
    """
    B_arr   = np.asarray(B_arr,   dtype=float)
    Vxx_arr = np.asarray(Vxx_arr, dtype=float)
    mask    = (B_arr > (B_center - dB_window)) & (B_arr < (B_center + dB_window))
    V_win   = Vxx_arr[mask]
    if V_win.size < 3:
        return np.nan, int(V_win.size)
    return float(np.nanmax(V_win) - np.nanmin(V_win)), int(V_win.size)


def plateau_deviation_ppm(B_arr, Rxy_arr, Bmin, Bmax, Rxy_ideal):
    """
    Compute the median Hall resistance and its deviation from the ideal
    quantised value within the plateau window [Bmin, Bmax].

    Returns: (dev_frac, dev_ppm, N, Rxy_med, Rxy_std)
    Returns NaN values if fewer than 5 finite points are found.
    """
    B_arr   = np.asarray(B_arr,   dtype=float)
    Rxy_arr = np.asarray(Rxy_arr, dtype=float)

    mask    = (B_arr >= Bmin) & (B_arr <= Bmax) & np.isfinite(Rxy_arr)
    N       = int(np.sum(mask))
    if N < 5:
        return np.nan, np.nan, N, np.nan, np.nan

    Rxy_med = float(np.nanmedian(Rxy_arr[mask]))
    Rxy_std = float(np.nanstd(Rxy_arr[mask]))
    dev     = (Rxy_med - Rxy_ideal) / Rxy_ideal
    return float(dev), float(dev * 1e6), N, Rxy_med, Rxy_std


def moving_slope_linear_fit(B, Y, win_pts=21):
    """
    Compute a local linear slope dY/dB at each point using a sliding
    window of win_pts points fitted by least squares.

    Points with fewer than max(6, win_pts//2) finite neighbours are left
    as NaN. Window is forced odd and at least 5 points wide.
    Used by plateau_window_from_slope to identify the flat plateau region.
    """
    B = np.asarray(B, dtype=float)
    Y = np.asarray(Y, dtype=float)

    win_pts = max(win_pts, 5)
    if win_pts % 2 == 0:
        win_pts += 1

    half  = win_pts // 2
    slope = np.full_like(Y, np.nan, dtype=float)

    for i in range(half, len(B) - half):
        Bi = B[i-half : i+half+1]
        Yi = Y[i-half : i+half+1]
        ok = np.isfinite(Bi) & np.isfinite(Yi)
        if np.sum(ok) < max(6, win_pts // 2):
            continue
        m, _c    = np.polyfit(Bi[ok], Yi[ok], 1)
        slope[i] = m

    return slope


def longest_segment_allow_gaps(mask, gap_tol_pts=3):
    """
    Find the longest contiguous True segment in a boolean mask, allowing
    short False runs of up to gap_tol_pts consecutive points.

    Returns (i0, i1) indices of the best segment, or (None, None) if none found.
    Used by plateau_window_from_slope to tolerate small noise excursions
    within the otherwise flat plateau region.
    """
    mask = np.asarray(mask, dtype=bool)
    if np.sum(mask) == 0:
        return None, None

    best_len         = 0
    best_i0, best_i1 = None, None
    i0               = 0
    false_run        = 0

    for i in range(mask.size):
        if mask[i]:
            false_run = 0
        else:
            false_run += 1
            if false_run > gap_tol_pts:
                i1      = i - false_run
                seg_len = i1 - i0 + 1
                if seg_len > best_len:
                    best_len         = seg_len
                    best_i0, best_i1 = i0, i1
                i0        = i
                false_run = 0

    # Check final segment
    i1      = mask.size - 1
    seg_len = i1 - i0 + 1
    if seg_len > best_len:
        best_i0, best_i1 = i0, i1

    return best_i0, best_i1


def plateau_window_from_slope(B, Rxy, B_center, search_halfwidth=0.25,
                              win_pts=21, slope_thresh=200.0,
                              gap_tol_pts=3, min_pts=8):
    """
    Objectively identify the ν = 2 Hall plateau window using a local-slope
    criterion applied to Rxy(B).

    Method:
      1. Restrict to B ∈ [B_center ± search_halfwidth].
      2. Compute the local slope dRxy/dB using a sliding linear fit
         of width win_pts.
      3. Flag points where |slope| ≤ slope_thresh as 'flat'.
      4. Find the longest flat segment, tolerating gaps of gap_tol_pts.
      5. Return the field bounds Bmin, Bmax of that segment.

    Returns a dict with keys 'ok' (bool), and if ok:
      Bmin, Bmax, N, slope_median, slope_rms.
    """
    B   = np.asarray(B,   dtype=float)
    Rxy = np.asarray(Rxy, dtype=float)

    mask_search = (
        (B   >= (B_center - search_halfwidth)) &
        (B   <= (B_center + search_halfwidth)) &
        np.isfinite(Rxy)
    )
    if np.sum(mask_search) < max(min_pts, win_pts):
        return {"ok": False, "reason": "Not enough points in search region."}

    Bs  = B[mask_search]
    Ys  = Rxy[mask_search]
    idx = np.argsort(Bs)
    Bs  = Bs[idx];  Ys = Ys[idx]

    slope = moving_slope_linear_fit(Bs, Ys, win_pts=win_pts)
    good  = np.isfinite(slope) & (np.abs(slope) <= slope_thresh)

    i0, i1 = longest_segment_allow_gaps(good, gap_tol_pts=gap_tol_pts)
    if i0 is None or i1 is None:
        return {"ok": False, "reason": "No flat segment found."}

    seg_B    = Bs[i0 : i1+1]
    seg_s    = slope[i0 : i1+1]
    seg_good = np.isfinite(seg_s) & (np.abs(seg_s) <= slope_thresh)

    if np.sum(seg_good) < min_pts:
        return {"ok": False, "reason": "Flat segment found but too short."}

    out = {"ok": True}
    out["Bmin"]        = float(np.nanmin(seg_B[seg_good]))
    out["Bmax"]        = float(np.nanmax(seg_B[seg_good]))
    out["N"]           = int(np.sum(seg_good))
    out["slope_median"] = float(np.nanmedian(seg_s[seg_good]))
    out["slope_rms"]   = float(np.sqrt(np.nanmean(seg_s[seg_good]**2)))
    return out


# ============================================================
# 3) MAIN PIPELINE: compute everything once into results + clean
# ============================================================

results = {"frozen": frozen.copy()}   # all scalar outputs go here
clean   = {}                          # all cleaned array outputs go here

# ---------------------------------------------------------------------------
# 3A) RUN 4: coarse current estimate for plateau search
#
# A median Vxy from a broad window around the expected plateau centre is
# divided by the ideal Rxy = h/2e² to give a coarse current I_coarse.
# This is used only for the initial despiking and plateau detection; the
# final I_used is derived from the objectively-defined plateau window in 3C.
# ---------------------------------------------------------------------------
t4, B4, Vxx4, Vxy4, vb4 = load_qhe_run(filename_run4)

nu_plateau      = frozen["nu_plateau"]
Rxy_ideal_nu2   = h_over_e2 / nu_plateau    # 12906.404 Ω for ν = 2

B_center_guess = frozen["B_center_guess"]
B_coarse_hw    = frozen["search_halfwidth"]

mask_coarse = (B4 >= (B_center_guess - B_coarse_hw)) & \
              (B4 <= (B_center_guess + B_coarse_hw))
if int(np.sum(mask_coarse)) < 10:
    raise RuntimeError("Run 4 coarse selection too small. "
                       "Check B_center_guess / search_halfwidth.")

Vxy4_coarse_med = float(np.nanmedian(Vxy4[mask_coarse]))
I_coarse        = np.abs(Vxy4_coarse_med) / Rxy_ideal_nu2

results["run4"] = {
    "Vxy_med_coarse": Vxy4_coarse_med,
    "I_coarse":       I_coarse,
}

# ---------------------------------------------------------------------------
# 3B) RUN 4: despike, split up/down, find slope-derived plateau windows
#
# The run is split at its field maximum (the turn-around point near 5 T)
# into an up-sweep (2→5 T) and a down-sweep (5→2 T). The plateau window
# is found independently in each sweep direction using the moving-slope
# criterion. The common intersection of the two windows is used for all
# subsequent plateau statistics, ensuring results are not biased by either
# sweep direction.
# ---------------------------------------------------------------------------
# Coarse resistance conversion for despiking (uses I_coarse, not I_used)
Rxy4_full_coarse = Vxy4 / I_coarse
Rxx4_full_coarse = np.abs(Vxx4) / I_coarse

# Median-filter despiking: kernel=61 for Rxy (broader, Hall is smoother)
# kernel=21 for Rxx (shorter, longitudinal oscillates more rapidly)
Rxx4_clean, keep_xx4 = despike(Rxx4_full_coarse, kernel=21, nsig=5.0)
Rxy4_clean, keep_xy4 = despike(Rxy4_full_coarse, kernel=61, nsig=5.0)
keep4  = keep_xx4 & keep_xy4   # keep only points clean in both channels

B4g    = B4[keep4]
Rxx4g  = Rxx4_clean[keep4]
Rxy4g  = Rxy4_clean[keep4]

# Split at the field turning point
idx_max4    = int(np.nanargmax(B4g))
B4_up,  Rxx_up,  Rxy_up  = B4g[:idx_max4+1], Rxx4g[:idx_max4+1], Rxy4g[:idx_max4+1]
B4_dn,  Rxx_dn,  Rxy_dn  = B4g[idx_max4:],   Rxx4g[idx_max4:],   Rxy4g[idx_max4:]

# Find plateau window in each sweep direction
res_up = plateau_window_from_slope(
    B4_up, Rxy_up,
    B_center       = B_center_guess,
    search_halfwidth = frozen["search_halfwidth"],
    win_pts        = frozen["slope_win_pts"],
    slope_thresh   = frozen["slope_threshold"],
    gap_tol_pts    = frozen["gap_tol_pts"],
    min_pts        = frozen["min_pts_plateau"],
)
res_dn = plateau_window_from_slope(
    B4_dn, Rxy_dn,
    B_center       = B_center_guess,
    search_halfwidth = frozen["search_halfwidth"],
    win_pts        = frozen["slope_win_pts"],
    slope_thresh   = frozen["slope_threshold"],
    gap_tol_pts    = frozen["gap_tol_pts"],
    min_pts        = frozen["min_pts_plateau"],
)

if (not res_up.get("ok", False)) or (not res_dn.get("ok", False)):
    raise RuntimeError("Slope-derived plateau failed for up or down sweep. "
                       "Check cleaning / slope settings.")

# Common plateau window = intersection of up and down windows
Bmin_common = max(res_up["Bmin"], res_dn["Bmin"])
Bmax_common = min(res_up["Bmax"], res_dn["Bmax"])

results["run4"]["plateau_up"]            = res_up
results["run4"]["plateau_dn"]            = res_dn
results["run4"]["Bmin_common"]           = float(Bmin_common)
results["run4"]["Bmax_common"]           = float(Bmax_common)
results["run4"]["plateau_width_common"]  = float(Bmax_common - Bmin_common)

# ---------------------------------------------------------------------------
# 3C) RUN 4: final current I_used and ns from the common plateau window
#
# The median Vxy within the objectively-defined common plateau window is
# divided by the ideal h/2e² to give the final excitation current I_used.
# This is anchored to fundamental constants and is independent of the
# lock-in readback (which was found to be ~17% lower; see Appendix B.3).
# I_used is then held fixed for all downstream resistance conversions.
# ---------------------------------------------------------------------------
mask_pl_common = (B4 >= Bmin_common) & (B4 <= Bmax_common) & np.isfinite(Vxy4)
if int(np.sum(mask_pl_common)) < frozen["min_pts_plateau"]:
    raise RuntimeError("Not enough raw Run 4 points in common plateau window. "
                       "Check Bmin_common / Bmax_common.")

B4_plateau_mean = float(np.nanmean(B4[mask_pl_common]))
Vxy4_med        = float(np.nanmedian(Vxy4[mask_pl_common]))

I_used              = np.abs(Vxy4_med) / Rxy_ideal_nu2
ns_from_plateau     = nu_plateau * e * B4_plateau_mean / h   # m⁻²

results["run4"]["B_plateau_mean"]  = B4_plateau_mean
results["run4"]["Vxy_med_plateau"] = Vxy4_med
results["run4"]["I_used"]          = I_used
results["run4"]["ns_from_plateau"] = ns_from_plateau

# ---------------------------------------------------------------------------
# 3D) RUN 3: SdH Landau fan → ns from Δ(1/B) and fan slope
#
# Second-derivative peak detection is used to locate SdH extrema robustly:
#   1. Apply Savitzky–Golay filter to compute d²Rxx/dB² simultaneously
#      with smoothing (avoids separate differentiation noise).
#   2. Find peaks in |d²Rxx/dB²| using a prominence threshold.
#   3. Convert peak B-positions to 1/B and compute mean spacing Δ(1/B).
#   4. Fit a line to (index n, 1/B_peak) to get the Landau fan slope.
#
# Note: peak detection uses |d²Rxx/dB²| but Fig16 plots the SIGNED
# second derivative for visual clarity of the oscillation structure.
# ---------------------------------------------------------------------------
t3, B3, Vxx3, Vxy3, vb3 = load_qhe_run(filename_run3)

B_min_SdH = frozen["B_min_SdH"]
B_max_SdH = frozen["B_max_SdH"]
mask_sdh3 = (B3 > B_min_SdH) & (B3 < B_max_SdH)

B3_sdh   = B3[mask_sdh3]
Rxx3_sdh = np.abs(Vxx3[mask_sdh3] / I_used)   # |Rxx| in SdH window

if B3_sdh.size < 15:
    raise RuntimeError("Not enough SdH points in Run 3; "
                       "adjust B_min_SdH / B_max_SdH.")

idx      = np.argsort(B3_sdh)
B3_sdh   = B3_sdh[idx]
Rxx3_sdh = Rxx3_sdh[idx]
dB3      = float(np.mean(np.diff(B3_sdh)))   # mean point spacing in B

# SG window: cap at sg_window_max, force odd, minimum 7 points
window_length = min(
    frozen["sg_window_max"],
    B3_sdh.size if B3_sdh.size % 2 == 1 else B3_sdh.size - 1
)
window_length = max(window_length, 7)
if window_length % 2 == 0:
    window_length += 1

polyorder    = frozen["sg_polyorder"]

# Signed second derivative (used for Fig16 visualisation)
d2Rxx3_dB2 = savgol_filter(Rxx3_sdh, window_length, polyorder,
                            deriv=2, delta=dB3)

# Peak picking on |d²Rxx/dB²| only — magnitude avoids sign ambiguity
d2_mag            = np.abs(d2Rxx3_dB2)
peak_prominence   = float(np.nanmax(d2_mag) * frozen["peak_prom_frac"])
peaks_idx, _props = find_peaks(d2_mag, prominence=peak_prominence)

B3_peaks = B3_sdh[peaks_idx]
if B3_peaks.size < 2:
    raise RuntimeError("Too few SdH peaks; tweak prominence or SdH window.")

# Store diagnostic arrays for Fig16
results["run3_diag"] = {
    "B_sdh"          : B3_sdh,
    "Rxx_sdh"        : Rxx3_sdh,
    "d2Rxx_dB2"      : d2Rxx3_dB2,   # signed — used in Fig16
    "peaks_idx"      : peaks_idx,
    "window_length"  : int(window_length),
    "polyorder"      : int(polyorder),
    "peak_prominence": float(peak_prominence),
}

# Landau fan: sort peaks by ascending 1/B
invB_peaks    = 1.0 / B3_peaks
order         = np.argsort(invB_peaks)
invB_peaks    = invB_peaks[order]
B3_peaks      = B3_peaks[order]

delta_invB_mean = float(np.mean(np.diff(invB_peaks)))
ns_from_d2      = e / (h * delta_invB_mean)   # m⁻²

n_indices    = np.arange(invB_peaks.size)
a_fit, b_fit = np.polyfit(n_indices, invB_peaks, 1)
ns_from_fan  = e / (h * a_fit)   # m⁻²

results["run3"] = {
    "SdH_window"      : (B_min_SdH, B_max_SdH),
    "n_peaks"         : int(B3_peaks.size),
    "delta_invB_mean" : delta_invB_mean,
    "ns_from_d2"      : ns_from_d2,
    "a_fit"           : float(a_fit),
    "b_fit"           : float(b_fit),
    "ns_from_fan"     : ns_from_fan,
}

# ---------------------------------------------------------------------------
# 3E) Combine ns estimates → ns_best; verify ν assignment
# ---------------------------------------------------------------------------
ns_list          = np.array([ns_from_plateau, ns_from_d2, ns_from_fan], dtype=float)
ns_best          = float(np.nanmean(ns_list))
nu_approx_refined = ns_best * h / (e * B4_plateau_mean)
nu_refined        = int(np.round(nu_approx_refined))

results["combined"] = {
    "ns_best"    : ns_best,
    "nu_approx"  : float(nu_approx_refined),
    "nu_refined" : int(nu_refined),
}

ns_fixed = ns_best   # alias used throughout remaining stages

# ---------------------------------------------------------------------------
# 3F) Low-field Hall slope from Run 3 (ρ-space, |B| < B0_lowfield)
#
# A linear fit to ρxy(B) at low fields (classical Hall regime) provides a
# fourth independent ns estimate and a cross-check on the Hall slope value.
# The slope a = 1/(ns e) in Ω/T; ns is extracted as 1/(e|a|).
# The low-field ρxx(0) also provides an alternative mobility estimate.
# ---------------------------------------------------------------------------
B0           = frozen["B0_lowfield"]
Rxx3_full    = np.abs(Vxx3) / I_used
Rxy3_full    = Vxy3 / I_used

rho_xx3_full = Rxx3_full * frozen["WL_ratio"]
rho_xy3_full = Rxy3_full

mask_low    = (np.abs(B3) <= B0) & np.isfinite(rho_xx3_full) & np.isfinite(rho_xy3_full)
B_low       = B3[mask_low]
rho_xy_low  = rho_xy3_full[mask_low]
rho_xx_low  = rho_xx3_full[mask_low]

a_hall = b_hall = ns_hall = rho_xx0 = mu_hall = np.nan

if B_low.size >= 10:
    a_hall, b_hall = np.polyfit(B_low, rho_xy_low, 1)
    ns_hall        = 1.0 / (e * a_hall)              # signed (negative = electrons)
    rho_xx0        = float(np.nanmean(rho_xx_low))
    mu_hall        = 1.0 / (ns_hall * e * rho_xx0)   # signed

results["lowfield"] = {
    "B0"           : float(B0),
    "a_hall"       : float(a_hall)       if np.isfinite(a_hall)  else np.nan,
    "b_hall"       : float(b_hall)       if np.isfinite(b_hall)  else np.nan,
    "ns_hall_abs"  : float(np.abs(ns_hall)) if np.isfinite(ns_hall) else np.nan,
    "rho_xx0"      : float(rho_xx0)      if np.isfinite(rho_xx0) else np.nan,
    "mu_hall_abs"  : float(np.abs(mu_hall)) if np.isfinite(mu_hall) else np.nan,
}

# ---------------------------------------------------------------------------
# 3G) RUN 1: mobility and transport lifetime
#
# The zero-field longitudinal resistance from the stability hold is used to
# extract the sheet resistance ρxx(0) = Rxx × (W/L), from which the mobility
# follows: μ = 1/(ns e ρxx). The dominant uncertainty is the 15.4% relative
# error in W/L.
# ---------------------------------------------------------------------------
data1   = np.genfromtxt(filename_run1, delimiter=",", skip_header=1)
t1      = data1[:, 0]
Vxx1    = data1[:, 3]
B1      = calibration_factor * data1[:, 9]

# Apply time and field cuts to isolate the stable zero-field window
mask    = (t1 > frozen["t_min_run1"]) & (np.abs(B1) < frozen["B_max_low_run1"])
t1c     = t1[mask]
Vxx1c   = Vxx1[mask]

# 5σ outlier removal on Vxx
Vxx_med = np.nanmedian(Vxx1c)
Vxx_std = np.nanstd(Vxx1c)
keep    = np.abs(Vxx1c - Vxx_med) <= frozen["N_sigma_run1"] * Vxx_std

t1c    = t1c[keep]
Vxx1c  = Vxx1c[keep]

Rxx1     = Vxx1c / I_used       # Ω (longitudinal resistance at B = 0)
Rxx_mean = float(np.nanmean(Rxx1))

R_sheet  = Rxx_mean * frozen["WL_ratio"]   # Ω/sq (sheet resistance)
mu_run1  = 1.0 / (ns_fixed * e * R_sheet)  # m²/(V s)
tau_tr   = mu_run1 * m_star / e            # s

results["mobility"] = {
    "Rxx_mean" : float(Rxx_mean),
    "R_sheet"  : float(R_sheet),
    "mu"       : float(mu_run1),    # m²/(V s)
    "tau_tr"   : float(tau_tr),     # s
}

# ---------------------------------------------------------------------------
# 3H) Cleaned arrays and full tensor conversion for figures
#
# Run 3: full bipolar sweep cleaned and converted.
# Run 4: up-sweep only used for ρ/σ tensor figures (the down-sweep retraces
#        the same transport without adding new information for these plots).
#
# Hall sign consistency: if the median Rxy on the plateau is negative
# (depending on wiring convention), both Rxy3 and Rxy4 are flipped to
# ensure positive-valued Hall resistance for conventional plotting.
# ---------------------------------------------------------------------------
Rxx4_full = np.abs(Vxx4) / I_used
Rxy4_full = Vxy4 / I_used

# Determine sign of Hall resistance from the plateau region
mask_pl_r  = (B4 >= Bmin_common) & (B4 <= Bmax_common) & np.isfinite(Rxy4_full)
hall_sign  = np.sign(np.nanmedian(Rxy4_full[mask_pl_r]))
if hall_sign < 0:
    Rxy3_full = -Rxy3_full
    Rxy4_full = -Rxy4_full

# Despike both channels of Runs 3 and 4
Rxx3_clean, keep_xx3 = despike(Rxx3_full, kernel=21, nsig=5.0)
Rxy3_clean, keep_xy3 = despike(Rxy3_full, kernel=61, nsig=5.0)
keep3  = keep_xx3 & keep_xy3
B3g    = B3[keep3];   Rxx3g = Rxx3_clean[keep3];   Rxy3g = Rxy3_clean[keep3]

Rxx4_clean, keep_xx4 = despike(Rxx4_full, kernel=21, nsig=5.0)
Rxy4_clean, keep_xy4 = despike(Rxy4_full, kernel=61, nsig=5.0)
keep4  = keep_xx4 & keep_xy4
B4g2   = B4[keep4];   Rxx4g2 = Rxx4_clean[keep4];   Rxy4g2 = Rxy4_clean[keep4]

# Split Run 4 into up/down for plateau deviation figure
idx_max4b      = int(np.nanargmax(B4g2))
B4_up2, Rxx_up2, Rxy_up2 = B4g2[:idx_max4b+1], Rxx4g2[:idx_max4b+1], Rxy4g2[:idx_max4b+1]
B4_dn2, Rxx_dn2, Rxy_dn2 = B4g2[idx_max4b:],   Rxx4g2[idx_max4b:],   Rxy4g2[idx_max4b:]

# Tensor conversion: Run 3 full sweep, Run 4 up-sweep
rho_xx3, rho_xy3, sigma_xx3, sigma_xy3       = resistances_to_tensors(Rxx3g,   Rxy3g,   frozen["WL_ratio"])
rho_xx4_up, rho_xy4_up, sigma_xx4_up, sigma_xy4_up = resistances_to_tensors(Rxx_up2, Rxy_up2, frozen["WL_ratio"])

# Plateau statistics on the common window (up and down separately)
_dev_up, ppm_up, N_up, Rxy_med_up, Rxy_std_up = plateau_deviation_ppm(
    B4_up2, Rxy_up2, Bmin_common, Bmax_common, Rxy_ideal_nu2)
_dev_dn, ppm_dn, N_dn, Rxy_med_dn, Rxy_std_dn = plateau_deviation_ppm(
    B4_dn2, Rxy_dn2, Bmin_common, Bmax_common, Rxy_ideal_nu2)

results["plateau_stats"] = {
    "Rxy_ideal_nu2" : float(Rxy_ideal_nu2),
    "Bmin_common"   : float(Bmin_common),
    "Bmax_common"   : float(Bmax_common),
    "up" : {"N": int(N_up), "Rxy_med": float(Rxy_med_up),
             "Rxy_std": float(Rxy_std_up), "ppm": float(ppm_up)},
    "dn" : {"N": int(N_dn), "Rxy_med": float(Rxy_med_dn),
             "Rxy_std": float(Rxy_std_dn), "ppm": float(ppm_dn)},
}

# Stash cleaned arrays for figure functions
clean["run3"] = {
    "B": B3g, "Rxx": Rxx3g, "Rxy": Rxy3g,
    "rho_xx": rho_xx3, "rho_xy": rho_xy3,
    "sigma_xx": sigma_xx3, "sigma_xy": sigma_xy3,
}
clean["run4"] = {
    "B": B4g2, "Rxx": Rxx4g2, "Rxy": Rxy4g2,
    "up": {
        "B": B4_up2, "Rxx": Rxx_up2, "Rxy": Rxy_up2,
        "rho_xx": rho_xx4_up, "rho_xy": rho_xy4_up,
        "sigma_xx": sigma_xx4_up, "sigma_xy": sigma_xy4_up,
    },
    "dn": {"B": B4_dn2, "Rxx": Rxx_dn2, "Rxy": Rxy_dn2},
}

# ---------------------------------------------------------------------------
# 3I) Temperature dependence (qualitative, Vxx-based)
#
# Vxx(B) overlays at 3 K, 4 K, 5 K and a simple peak-to-peak amplitude
# extraction at a fixed reference field. Only three temperature points are
# available, so a full Lifshitz–Kosevich fit is not attempted.
# ---------------------------------------------------------------------------
B3_T, Vxx3_T, _Vxy3_T, _t3_T = load_run_sorted(filename_run3)
B6_T, Vxx6_T, _Vxy6_T, _t6_T = load_run_sorted(filename_run6)
B7_T, Vxx7_T, _Vxy7_T, _t7_T = load_run_sorted(filename_run7)

mask3T = (B3_T > frozen["B_min_common_T"]) & (B3_T < frozen["B_max_common_T"])
mask6T = (B6_T > frozen["B_min_common_T"]) & (B6_T < frozen["B_max_common_T"])
mask7T = (B7_T > frozen["B_min_common_T"]) & (B7_T < frozen["B_max_common_T"])

B3_pos, Vxx3_pos = B3_T[mask3T], Vxx3_T[mask3T]
B6_pos, Vxx6_pos = B6_T[mask6T], Vxx6_T[mask6T]
B7_pos, Vxx7_pos = B7_T[mask7T], Vxx7_T[mask7T]

A3, N3 = sdh_amplitude(B3_pos, Vxx3_pos, frozen["B_peak_ref"], frozen["dB_peak"])
A6, N6 = sdh_amplitude(B6_pos, Vxx6_pos, frozen["B_peak_ref"], frozen["dB_peak"])
A7, N7 = sdh_amplitude(B7_pos, Vxx7_pos, frozen["B_peak_ref"], frozen["dB_peak"])

results["temperature"] = {
    "B_peak_ref": float(frozen["B_peak_ref"]),
    "dB_peak"   : float(frozen["dB_peak"]),
    "A3": float(A3), "N3": int(N3),
    "A4": float(A6), "N4": int(N6),
    "A5": float(A7), "N5": int(N7),
}


# ============================================================
# 4) RESULTS SUMMARY PRINT
# ============================================================

print("\n" + "=" * 60)
print("FINAL TRANSPORT PARAMETERS SUMMARY (RESULTS-READY)")
print("=" * 60)

print("\n[Geometry]")
print(f"W/L used:              {frozen['WL_ratio']:.5f} ± {frozen['WL_unc']:.5f}  ({frozen['WL_rel_unc']*100:.1f} %)")

print("\n[Plateau]")
print(f"ν used:                {nu_plateau:d}")
print(f"Plateau window:        B ∈ [{Bmin_common:.4f}, {Bmax_common:.4f}] T")
print(f"Plateau width:         {(Bmax_common - Bmin_common):.4f} T")
print(f"B_plateau_mean:        {B4_plateau_mean:.4f} T")

print("\n[Current + density]")
print(f"I_used:                {I_used:.3e} A  (plateau-derived)")
print(f"n_s (plateau):         {ns_from_plateau:.3e} m⁻²")
print(f"n_s (SdH Δ(1/B)):      {ns_from_d2:.3e} m⁻²")
print(f"n_s (Landau fan):      {ns_from_fan:.3e} m⁻²")
print(f"n_s (best mean):       {ns_best:.3e} m⁻²  ({ns_best/1e4:.3e} cm⁻²)")
print(f"ν check:               {nu_approx_refined:.2f} → ν_refined = {nu_refined:d}")

print("\n[Low-field Hall cross-check (Run 3)]")
if np.isfinite(results["lowfield"]["ns_hall_abs"]):
    print(f"Hall slope a:          {results['lowfield']['a_hall']:.6g} Ω/T")
    print(f"|ns| (Hall slope):     {results['lowfield']['ns_hall_abs']:.3e} m⁻²")
    print(f"ρxx(0):                {results['lowfield']['rho_xx0']:.3f} Ω/sq")
    print(f"|μ| (Hall+ρxx):        {results['lowfield']['mu_hall_abs']:.3e} m²/(V·s)  "
          f"({results['lowfield']['mu_hall_abs']*1e4:.3e} cm²/(V·s))")
else:
    print("Low-field Hall fit:    FAILED / insufficient points")

print("\n[Mobility (Run 1 + ns_best)]")
print(f"R_sheet:               {results['mobility']['R_sheet']:.3f} Ω/sq")
print(f"μ:                     {results['mobility']['mu']:.3e} m²/(V·s)  "
      f"({results['mobility']['mu']*1e4:.3e} cm²/(V·s))")
print(f"τ_tr:                  {results['mobility']['tau_tr']:.3e} s  "
      f"({results['mobility']['tau_tr']*1e12:.3e} ps)")

print("\n[Plateau quantisation (common window)]")
print(f"Rxy ideal (ν=2):       {results['plateau_stats']['Rxy_ideal_nu2']:.3f} Ω")
print(f"Up:  median ± std:     {results['plateau_stats']['up']['Rxy_med']:.3f} ± "
      f"{results['plateau_stats']['up']['Rxy_std']:.3f} Ω   "
      f"({results['plateau_stats']['up']['ppm']:.1f} ppm)")
print(f"Down: median ± std:    {results['plateau_stats']['dn']['Rxy_med']:.3f} ± "
      f"{results['plateau_stats']['dn']['Rxy_std']:.3f} Ω   "
      f"({results['plateau_stats']['dn']['ppm']:.1f} ppm)")

print("\n[Temperature dependence (SdH amplitude)]")
print(f"At B ≈ {frozen['B_peak_ref']:.2f} T (±{frozen['dB_peak']:.2f} T):")
print(f"  3 K: A = {A3:.3e} V  (N = {N3})")
print(f"  4 K: A = {A6:.3e} V  (N = {N6})")
print(f"  5 K: A = {A7:.3e} V  (N = {N7})")

print("\n" + "=" * 60)
print("DONE: analysis pipeline complete. Running figure export...")
print("=" * 60)


# ============================================================
# 5) REPORT FIGURE FUNCTIONS
# ============================================================
# Each function takes the pre-computed clean and results dicts and
# produces exactly one output file (PNG + PDF). Figure names match
# the LaTeX \includegraphics filenames in the submitted report.
# require_keys() guards each function against missing upstream data.

def fig05_overview_Rxy_Rxx_T3K(clean: dict, results: dict, outdir: str) -> None:
    """
    Fig05: overview magnetotransport at T = 3 K.
    Top panel: Rxy(B) — Hall resistance staircase.
    Bottom panel: Rxx(B) — longitudinal resistance with SdH oscillations.
    Shaded region marks the common ν = 2 plateau window.
    Data source: Run 4 up-sweep (2→5 T, 3 K).
    """
    require_keys(results, ["plateau_stats"], "results")
    require_keys(results["plateau_stats"], ["Bmin_common", "Bmax_common"], "results['plateau_stats']")
    require_keys(clean, ["run4"], "clean")
    require_keys(clean["run4"], ["up"], "clean['run4']")
    require_keys(clean["run4"]["up"], ["B", "Rxx", "Rxy"], "clean['run4']['up']")

    Bmin = results["plateau_stats"]["Bmin_common"]
    Bmax = results["plateau_stats"]["Bmax_common"]
    B    = clean["run4"]["up"]["B"]
    Rxx  = clean["run4"]["up"]["Rxx"]
    Rxy  = clean["run4"]["up"]["Rxy"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True)

    ax1.plot(B, Rxy, linewidth=1.0)
    ax1.axvspan(Bmin, Bmax, alpha=0.15)
    ax1.set_ylabel(r"$R_{xy}$ (Ω)")
    ax1.set_title(r"Overview magnetotransport at $T=3\,\mathrm{K}$ (Run 4 up-sweep)")

    ax2.plot(B, Rxx, linewidth=1.0)
    ax2.axvspan(Bmin, Bmax, alpha=0.15)
    ax2.set_xlabel(r"$B$ (T)")
    ax2.set_ylabel(r"$R_{xx}$ (Ω)")

    savefig_both(fig, outdir, "Fig05_Overview_Rxy_Rxx_T3K")


def fig06_carrier_density_comparison(results: dict, outdir: str) -> None:
    """
    Fig06: bar chart comparing ns extracted by all four independent methods.
    Methods: Hall plateau, SdH Δ(1/B), Landau fan, low-field Hall slope.
    The mean (ns_best) is shown as a fifth bar for reference.
    """
    require_keys(results, ["run3", "run4", "combined", "lowfield"], "results")
    require_keys(results["run4"],     ["ns_from_plateau"],          "results['run4']")
    require_keys(results["run3"],     ["ns_from_d2", "ns_from_fan"],"results['run3']")
    require_keys(results["lowfield"], ["ns_hall_abs"],               "results['lowfield']")
    require_keys(results["combined"], ["ns_best"],                   "results['combined']")

    ns_plateau = float(results["run4"]["ns_from_plateau"])
    ns_sdh     = float(results["run3"]["ns_from_d2"])
    ns_fan     = float(results["run3"]["ns_from_fan"])
    ns_hall    = float(results["lowfield"]["ns_hall_abs"])
    ns_best_l  = float(results["combined"]["ns_best"])

    labels = ["Hall plateau", "SdH Δ(1/B)", "Landau fan", "Low-field Hall", "Mean (best)"]
    values = np.array([ns_plateau, ns_sdh, ns_fan, ns_hall, ns_best_l], dtype=float) / 1e15
    x      = np.arange(values.size)

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.bar(x, values)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel(r"$n_s$ ($10^{15}\,\mathrm{m^{-2}}$)")
    ax.set_title("Carrier density from independent methods")

    savefig_both(fig, outdir, "Fig06_CarrierDensity_Comparison_BarChart")


def fig07_nu2_rxy_zoom(clean: dict, results: dict, outdir: str) -> None:
    """
    Fig07: zoomed view of Rxy(B) around the ν = 2 Hall plateau.
    Shaded region = common plateau window. Dashed line = ideal h/2e².
    """
    require_keys(results["plateau_stats"],
                 ["Bmin_common", "Bmax_common", "Rxy_ideal_nu2"], "results['plateau_stats']")
    require_keys(clean["run4"]["up"], ["B", "Rxy"], "clean['run4']['up']")

    B       = clean["run4"]["up"]["B"]
    Rxy     = clean["run4"]["up"]["Rxy"]
    Bmin    = results["plateau_stats"]["Bmin_common"]
    Bmax    = results["plateau_stats"]["Bmax_common"]
    Rxy_ideal = results["plateau_stats"]["Rxy_ideal_nu2"]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(B, Rxy, linewidth=1.0)
    ax.axvspan(Bmin, Bmax, alpha=0.20, label="Common plateau window")
    ax.axhline(Rxy_ideal, linestyle="--", linewidth=1.0, label=r"$h/2e^2$")
    ax.set_xlim(Bmin - 0.15, Bmax + 0.15)
    ax.set_xlabel(r"$B$ (T)")
    ax.set_ylabel(r"$R_{xy}$ (Ω)")
    ax.set_title(r"Zoom: $\nu=2$ Hall plateau (Run 4 up-sweep)")
    ax.legend(loc="best")

    savefig_both(fig, outdir, "Fig07_Nu2_Rxy_Zoom")


def fig08_nu2_rxx_suppression_zoom(clean: dict, results: dict, outdir: str) -> None:
    """
    Fig08: zoomed Rxx(B) showing strong suppression over the ν = 2 plateau.
    Shaded region = common plateau window, matching Fig07 for direct comparison.
    """
    require_keys(results["plateau_stats"], ["Bmin_common", "Bmax_common"],
                 "results['plateau_stats']")
    require_keys(clean["run4"]["up"], ["B", "Rxx"], "clean['run4']['up']")

    B    = clean["run4"]["up"]["B"]
    Rxx  = clean["run4"]["up"]["Rxx"]
    Bmin = results["plateau_stats"]["Bmin_common"]
    Bmax = results["plateau_stats"]["Bmax_common"]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(B, Rxx, linewidth=1.0)
    ax.axvspan(Bmin, Bmax, alpha=0.20, label="Common plateau window")
    ax.set_xlim(Bmin - 0.15, Bmax + 0.15)
    ax.set_xlabel(r"$B$ (T)")
    ax.set_ylabel(r"$R_{xx}$ (Ω)")
    ax.set_title(r"Suppression of $R_{xx}$ on the $\nu=2$ plateau (Run 4 up-sweep)")
    ax.legend(loc="best")

    savefig_both(fig, outdir, "Fig08_Nu2_Rxx_Suppression_Zoom")


def fig09_nu2_updown_deviation(clean: dict, results: dict, outdir: str) -> None:
    """
    Fig09: deviation of Rxy from ideal h/2e² in ppm for up- and down-sweeps.
    Demonstrates reproducibility and absence of hysteresis across the plateau.
    Only points within the common plateau window are plotted.
    """
    require_keys(results["plateau_stats"],
                 ["Bmin_common", "Bmax_common", "Rxy_ideal_nu2"], "results['plateau_stats']")
    require_keys(clean["run4"]["up"], ["B", "Rxy"], "clean['run4']['up']")
    require_keys(clean["run4"]["dn"], ["B", "Rxy"], "clean['run4']['dn']")

    Bmin      = results["plateau_stats"]["Bmin_common"]
    Bmax      = results["plateau_stats"]["Bmax_common"]
    Rxy_ideal = results["plateau_stats"]["Rxy_ideal_nu2"]

    Bu, Rxyu = clean["run4"]["up"]["B"], clean["run4"]["up"]["Rxy"]
    Bd, Rxyd = clean["run4"]["dn"]["B"], clean["run4"]["dn"]["Rxy"]

    def ppm_curve(B, Rxy):
        """Compute ppm deviation from Rxy_ideal; NaN outside plateau window."""
        ppm  = np.full_like(Rxy, np.nan, dtype=float)
        mask = (B >= Bmin) & (B <= Bmax) & np.isfinite(Rxy)
        ppm[mask] = (Rxy[mask] - Rxy_ideal) / Rxy_ideal * 1e6
        return ppm

    ppm_u = ppm_curve(Bu, Rxyu)
    ppm_d = ppm_curve(Bd, Rxyd)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(Bu, ppm_u, linewidth=1.0, label="Up-sweep")
    ax.plot(Bd, ppm_d, linewidth=1.0, label="Down-sweep")
    ax.axvspan(Bmin, Bmax, alpha=0.20)
    ax.set_xlim(Bmin - 0.10, Bmax + 0.10)
    ax.set_xlabel(r"$B$ (T)")
    ax.set_ylabel(r"Deviation from $h/2e^2$ (ppm)")
    ax.set_title(r"Plateau deviation from $h/2e^2$: up- vs down-sweep")
    ax.legend(loc="best")

    savefig_both(fig, outdir, "Fig09_Nu2_UpDown_Deviation_From_Quantised")


def fig10_resistivity_tensor(clean: dict, results: dict, outdir: str) -> None:
    """
    Fig10: resistivity tensor components ρxy(B) and ρxx(B).
    Combined two-panel figure showing both components together.
    A linear classical Hall trend ρxy = B/(ns e) is overlaid on ρxy
    to show the crossover from classical to quantised behaviour.
    """
    require_keys(results["combined"],    ["ns_best"],                     "results['combined']")
    require_keys(results["plateau_stats"],["Bmin_common","Bmax_common"],   "results['plateau_stats']")
    require_keys(clean["run3"],          ["B","rho_xx","rho_xy"],          "clean['run3']")
    require_keys(clean["run4"]["up"],    ["B","rho_xx","rho_xy"],          "clean['run4']['up']")

    Bmin = results["plateau_stats"]["Bmin_common"]
    Bmax = results["plateau_stats"]["Bmax_common"]
    B3   = clean["run3"]["B"]
    B4   = clean["run4"]["up"]["B"]
    rho_xx3 = clean["run3"]["rho_xx"];   rho_xy3 = clean["run3"]["rho_xy"]
    rho_xx4 = clean["run4"]["up"]["rho_xx"]; rho_xy4 = clean["run4"]["up"]["rho_xy"]

    ns_best_l  = float(results["combined"]["ns_best"])
    rho_xy_lin3 = B3 / (ns_best_l * e)   # classical Hall line
    rho_xy_lin4 = B4 / (ns_best_l * e)

    # Build combined two-panel figure
    figc, (axc1, axc2) = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)

    axc1.plot(B3, rho_xy3, linewidth=1.0, label="Run 3")
    axc1.plot(B4, rho_xy4, linewidth=1.0, label="Run 4 (up)")
    axc1.plot(B3, rho_xy_lin3, linestyle="--", linewidth=0.9, label="Linear Hall trend")
    axc1.axvspan(Bmin, Bmax, alpha=0.15)
    axc1.set_ylabel(r"$\rho_{xy}$ (Ω)")
    axc1.legend(loc="best")

    axc2.plot(B3, rho_xx3, linewidth=1.0, label="Run 3")
    axc2.plot(B4, rho_xx4, linewidth=1.0, label="Run 4 (up)")
    axc2.axvspan(Bmin, Bmax, alpha=0.15)
    axc2.set_xlabel(r"$B$ (T)")
    axc2.set_ylabel(r"$\rho_{xx}$ (Ω/sq)")

    savefig_both(figc, outdir, "Fig10_ResistivityTensor_rho_xx_rho_xy")


def fig11_conductivity_tensor_masked(clean: dict, results: dict, outdir: str,
                                     B_mask: float = 0.25) -> None:
    """
    Fig11: conductivity tensor components σxy(B) and σxx(B).
    The low-field region |B| < B_mask is masked because the tensor inversion
    becomes numerically ill-conditioned as ρxy → 0 near B = 0.
    Outside this region, minima in σxx coincide with plateaux in σxy,
    consistent with the localisation picture of the integer QHE.
    """
    require_keys(results["plateau_stats"], ["Bmin_common","Bmax_common"], "results['plateau_stats']")
    require_keys(clean["run3"],            ["B","sigma_xx","sigma_xy"],   "clean['run3']")
    require_keys(clean["run4"]["up"],      ["B","sigma_xx","sigma_xy"],   "clean['run4']['up']")

    Bmin = results["plateau_stats"]["Bmin_common"]
    Bmax = results["plateau_stats"]["Bmax_common"]

    B3      = clean["run3"]["B"]
    sig_xx3 = clean["run3"]["sigma_xx"] / e2_over_h   # in units of e²/h
    sig_xy3 = clean["run3"]["sigma_xy"] / e2_over_h

    B4      = clean["run4"]["up"]["B"]
    sig_xx4 = clean["run4"]["up"]["sigma_xx"] / e2_over_h
    sig_xy4 = clean["run4"]["up"]["sigma_xy"] / e2_over_h

    # Apply low-field mask — exclude ill-conditioned inversion near B = 0
    m3 = (np.abs(B3) >= B_mask) & np.isfinite(sig_xx3) & np.isfinite(sig_xy3)
    m4 = (np.abs(B4) >= B_mask) & np.isfinite(sig_xx4) & np.isfinite(sig_xy4)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)

    ax1.plot(B3[m3], sig_xy3[m3], linewidth=1.0, label="Run 3")
    ax1.plot(B4[m4], sig_xy4[m4], linewidth=1.0, label="Run 4 (up)")
    ax1.axvspan(Bmin, Bmax, alpha=0.15)
    ax1.set_ylabel(r"$\sigma_{xy}$ ($e^2/h$)")
    ax1.set_title(f"Conductivity tensor (|B| < {B_mask} T masked)")
    ax1.legend(loc="best")

    ax2.plot(B3[m3], sig_xx3[m3], linewidth=1.0, label="Run 3")
    ax2.plot(B4[m4], sig_xx4[m4], linewidth=1.0, label="Run 4 (up)")
    ax2.axvspan(Bmin, Bmax, alpha=0.15)
    ax2.set_xlabel(r"$B$ (T)")
    ax2.set_ylabel(r"$\sigma_{xx}$ ($e^2/h$)")

    savefig_both(fig, outdir, "Fig11_ConductivityTensor_sigma_xx_sigma_xy_Masked")


def fig12_mobility_diagnostics(t_arr: np.ndarray, Rxx_arr: np.ndarray,
                                outdir: str) -> None:
    """
    Fig12: two-panel mobility diagnostic from Run 1 (zero-field hold, 3 K).
    Left: Rxx vs time showing stability of the zero-field measurement.
    Right: histogram of Rxx confirming an approximately Gaussian distribution.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.6))

    ax1.plot(t_arr, Rxx_arr, linewidth=1.0)
    ax1.set_xlabel("t (s)")
    ax1.set_ylabel(r"$R_{xx}$ (Ω)")
    ax1.set_title("(a) Low-field stability")

    ax2.hist(Rxx_arr[np.isfinite(Rxx_arr)], bins=30)
    ax2.set_xlabel(r"$R_{xx}$ (Ω)")
    ax2.set_ylabel("Counts")
    ax2.set_title("(b) Distribution")

    savefig_both(fig, outdir, "Fig12_Mobility_LowField_Diagnostics")


def fig13_vxx_overlay_T3_T4_T5K(B3_pos, Vxx3_pos, B6_pos, Vxx6_pos,
                                  B7_pos, Vxx7_pos, outdir: str) -> None:
    """
    Fig13: overlay of Vxx(B) at T = 3 K, 4 K, 5 K in the common field window.
    Demonstrates monotonic decrease in SdH oscillation amplitude with temperature,
    consistent with Lifshitz–Kosevich thermal broadening of Landau levels.
    """
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(B3_pos, Vxx3_pos, linewidth=1.0, label="3 K")
    ax.plot(B6_pos, Vxx6_pos, linewidth=1.0, label="4 K")
    ax.plot(B7_pos, Vxx7_pos, linewidth=1.0, label="5 K")
    ax.set_xlabel(r"$B$ (T)")
    ax.set_ylabel(r"$V_{xx}$ (V)")
    ax.set_title(r"Temperature dependence of SdH oscillations ($V_{xx}$)")
    ax.legend(loc="best")
    savefig_both(fig, outdir, "Fig13_Vxx_Overlay_T3_T4_T5K")


def fig14_sdh_amplitude_vs_temperature(results: dict, outdir: str) -> None:
    """
    Fig14: SdH oscillation amplitude vs temperature.
    Simple peak-to-peak amplitude extracted at the reference field B_peak_ref.
    Only three temperature points available; analysis is qualitative.
    """
    require_keys(results["frozen"],      ["T_run3","T_run6","T_run7"], "results['frozen']")
    require_keys(results["temperature"], ["A3","A4","A5"],             "results['temperature']")

    T = np.array([results["frozen"]["T_run3"],
                  results["frozen"]["T_run6"],
                  results["frozen"]["T_run7"]], dtype=float)
    A = np.array([results["temperature"]["A3"],
                  results["temperature"]["A4"],
                  results["temperature"]["A5"]], dtype=float)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(T, A, marker="o", linewidth=1.0)
    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"SdH amplitude $A$ (V)")
    ax.set_title("SdH amplitude vs temperature (qualitative)")
    savefig_both(fig, outdir, "Fig14_SdH_Amplitude_vs_Temperature")


def fig15_landau_fan(invB_peaks: np.ndarray, a_fit: float, b_fit: float,
                     outdir: str) -> None:
    """
    Fig15: Landau fan diagram — SdH extremum index n vs 1/B_peak.
    The slope of the linear fit gives ns via slope = h/(e ns).
    The intercept has no direct physical meaning here.
    """
    n   = np.arange(invB_peaks.size)
    fit = a_fit * n + b_fit

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(n, invB_peaks, marker="o", linewidth=0.0, label="Extrema")
    ax.plot(n, fit,        linewidth=1.0,              label="Linear fit")
    ax.set_xlabel("Extremum index $n$")
    ax.set_ylabel(r"$1/B$ (T$^{-1}$)")
    ax.set_title("Landau fan diagram from SdH extrema")
    ax.legend(loc="best")

    savefig_both(fig, outdir, "Fig15_LandauFan_SdH_Index_vs_InvB")


def fig16_rxx_vs_second_derivative(B_sdh: np.ndarray, Rxx_sdh: np.ndarray,
                                    d2Rxx_dB2: np.ndarray, peaks_idx: np.ndarray,
                                    outdir: str, show_vlines: bool = True) -> None:
    """
    Fig16: two-panel SdH peak detection diagnostic.
    Top: raw Rxx(B) in the SdH window with vertical markers at detected peaks.
    Bottom: SIGNED d²Rxx/dB² with peak markers.

    Note: peak detection used |d²Rxx/dB²| but the signed version is plotted
    here for visual clarity of the oscillation structure (improvement noted
    in the lab diary on 29 December 2025 compared to earlier scripts).
    """
    B_sdh      = np.asarray(B_sdh,      dtype=float)
    Rxx_sdh    = np.asarray(Rxx_sdh,    dtype=float)
    d2Rxx_dB2  = np.asarray(d2Rxx_dB2,  dtype=float)
    peaks_idx  = np.asarray(peaks_idx,  dtype=int)

    if B_sdh.size != Rxx_sdh.size or B_sdh.size != d2Rxx_dB2.size:
        raise ValueError("B_sdh, Rxx_sdh, and d2Rxx_dB2 must have the same length.")
    if peaks_idx.size > 0:
        peaks_idx = peaks_idx[(peaks_idx >= 0) & (peaks_idx < B_sdh.size)]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True)

    ax1.plot(B_sdh, Rxx_sdh, linewidth=1.0)
    ax1.set_ylabel(r"$R_{xx}$ (Ω)")
    ax1.set_title("SdH peak detection using second derivative")

    if show_vlines and peaks_idx.size > 0:
        for Bp in B_sdh[peaks_idx]:
            ax1.axvline(Bp, linewidth=0.8, alpha=0.25)

    ax2.plot(B_sdh, d2Rxx_dB2, linewidth=1.0)
    if peaks_idx.size > 0:
        ax2.plot(B_sdh[peaks_idx], d2Rxx_dB2[peaks_idx], marker="o", linewidth=0.0)

    if show_vlines and peaks_idx.size > 0:
        for Bp in B_sdh[peaks_idx]:
            ax2.axvline(Bp, linewidth=0.8, alpha=0.25)

    ax2.set_xlabel(r"$B$ (T)")
    ax2.set_ylabel(r"$d^2R_{xx}/dB^2$ (arb.)")

    savefig_both(fig, outdir, "Fig16_Rxx_Raw_vs_SecondDerivative_PeakDetection")


# ============================================================
# 6) APPENDIX FIGURES
# ============================================================

def figA1_sdh_window(B_sdh, Rxx_sdh, outdir: str) -> None:
    """
    FigA1: appendix figure showing the full SdH oscillation window from Run 3.
    Provides context for the peak detection in Fig16.
    """
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(B_sdh, Rxx_sdh, linewidth=1.0)
    ax.set_xlabel(r"$B$ (T)")
    ax.set_ylabel(r"$R_{xx}$ (Ω)")
    ax.set_title("Appendix: SdH oscillation window (Run 3)")
    savefig_both(fig, outdir, "FigA1_SdH_window")


# ============================================================
# 7) EXPORT ALL FIGURES
# ============================================================

def export_report_figures(main_dir: str = "report_figures",
                          appendix_dir: str = "appendix_figures",
                          export_appendix: bool = True) -> None:
    """
    Call all figure functions in sequence, writing PNG and PDF outputs
    to main_dir (report figures) and appendix_dir (appendix figures).

    Directories are created automatically if they do not exist.
    All figure functions consume the pre-computed `clean` and `results`
    dicts populated by the pipeline above; nothing is recomputed here.
    """
    ensure_dir(main_dir)
    ensure_dir(appendix_dir)

    # ---- Main report figures ----
    fig05_overview_Rxy_Rxx_T3K(clean, results, main_dir)
    fig06_carrier_density_comparison(results, main_dir)
    fig07_nu2_rxy_zoom(clean, results, main_dir)
    fig08_nu2_rxx_suppression_zoom(clean, results, main_dir)
    fig09_nu2_updown_deviation(clean, results, main_dir)
    fig10_resistivity_tensor(clean, results, main_dir)
    fig11_conductivity_tensor_masked(clean, results, main_dir,
                                     B_mask=frozen["B0_lowfield"])
    fig12_mobility_diagnostics(t1c, Rxx1, main_dir)
    fig13_vxx_overlay_T3_T4_T5K(B3_pos, Vxx3_pos,
                                  B6_pos, Vxx6_pos,
                                  B7_pos, Vxx7_pos, main_dir)
    fig14_sdh_amplitude_vs_temperature(results, main_dir)

    # Landau fan and peak detection figures use diagnostic arrays from 3D
    fig15_landau_fan(invB_peaks,
                     results["run3"]["a_fit"],
                     results["run3"]["b_fit"],
                     main_dir)
    fig16_rxx_vs_second_derivative(
        results["run3_diag"]["B_sdh"],
        results["run3_diag"]["Rxx_sdh"],
        results["run3_diag"]["d2Rxx_dB2"],
        results["run3_diag"]["peaks_idx"],
        main_dir,
        show_vlines=True
    )

    # ---- Appendix figures ----
    if export_appendix:
        figA1_sdh_window(B3_sdh, Rxx3_sdh, appendix_dir)

    print("\nAll figures exported.")
    print(f"  Main report:   {main_dir}/")
    if export_appendix:
        print(f"  Appendix-only: {appendix_dir}/")


# ---- Run the exporter ----
export_report_figures(
    main_dir       = "report_figures",
    appendix_dir   = "appendix_figures",
    export_appendix = True
)
