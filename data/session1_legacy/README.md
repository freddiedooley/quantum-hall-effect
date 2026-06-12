# Session 1 Data — Legacy Records
**Date:** 6 November 2025  
**Sample:** NU1783 — GaAs/AlGaAs two-dimensional electron gas, Hall bar geometry  
**Cryostat:** Helium insert cryostat, DEWAR 01  

---

## CSV Files

Session 1 was the first cryogenic measurement session of the project. Nine
field sweeps were performed, reaching a base temperature of approximately
1.9–2 K — lower than was achievable in Session 2 due to better cryostat
performance on this occasion.

However, all quantitative results in the final report and in both analysis scripts
are derived exclusively from Session 2 data (20 November 2025).

---

## Eight PNG Screenshots

Eight PNG screenshots of the live acquisition display were captured during
Session 1 and constitute the only surviving record of those measurements.
They are preserved here as a historical record of the session.

These images show the real-time matplotlib plots produced by
`qhe_data_acquisition.py` during each sweep — Vxx and Vxy plotted against
both elapsed time and the field-proportional voltage V_B.

| File | Run | Field range | Temperature | Notes |
|------|-----|-------------|-------------|-------|
| `session1_run2_bipolar.png` | Run 2 | +2 → −2 T | ~1.9 K | SdH oscillations in Vxx visible; Vxy linear Hall slope changes sign at B = 0 as expected |
| `session1_run3_bipolar.png` | Run 3 | −2 → +2 T | ~2.04 K | Mirror image of Run 2; excellent SdH oscillations in both field polarities |
| `session1_run4_partial.png` | Run 4 | +2 → 0 → +2 T | ~2 K | Partial sweep |
| `session1_run5_overrange.png` | Run 5 | 0 → +4 T | ~2 K | Lock-in overrange visible as vertical spikes in Vxy; sensitivity adjusted after this run |
| `session1_run6_vibration.png` | Run 6 | +4 → −4 T | ~2 K | Sharp step discontinuity in both Vxx and Vxy mid-sweep caused by a physical vibration event (footstep near the cryostat); vivid illustration of the vibration sensitivity of cryogenic lock-in measurements |
| `session1_run7_partial.png` | Run 7 | −4 → 0 T | ~2 K | Recovery sweep after vibration event |
| `session1_run8_overrange.png` | Run 8 | 0 → +5 T | ~2 K | Vxy lock-in overloaded at 10 V range near high field; range changed to 20 V mid-run |
| `session1_run9_plateau.png` | Run 9 | +5 → 0 T | ~2 K | Clearest plateau dataset from Session 1; plateau structure near V_B ≈ 2.5 V (≈ 3.3 T) clearly visible with Vxx minimum in the same field range; used for the first circular current estimate in the early analysis |

> **Note on filenames:** The PNG filenames above are descriptive labels
> assigned for this repository. The original filenames from the live
> acquisition display may differ.

---

## Scientific Value of the Session 1 Images

Although no quantitative analysis was performed on Session 1 data, these
images provided important early evidence during the project:

**Run 9** was used in the first analysis session (19 November 2025) to
make an initial circular estimate of the excitation current by assuming
a carrier density and identifying the ν = 2 plateau visually. This
motivated the Session 2 analysis strategy of deriving ns independently
before anchoring the current to the quantised plateau resistance.

**Run 6** (the vibration event) provided a clear practical demonstration
of why cryogenic lock-in measurements require mechanical isolation and
why sudden disturbances during a field sweep can corrupt an entire run.

**Runs 2 and 3** (bipolar sweeps at ~2 K) show SdH oscillations with
greater amplitude than the Session 2 runs at 3 K, consistent with the
lower base temperature achieved in Session 1 before the cryostat leak
degraded thermal performance. This supports the temperature dependence
analysis in the final report.

---

## The Helium Leak

During Session 1, a leak between the sample space and the vacuum chamber
was diagnosed at approximately 4 K. The sample temperature suddenly rose
to 60–70 K and the vacuum gauge tick rate increased rapidly, indicating
helium entering the vacuum space and increasing the heat load on the 4 K
region. A secondary pumping station was quickly connected and the
temperature restabilised.

This leak persisted into Session 2 and prevented stable operation below
approximately 3 K for all subsequent measurements. It is the primary
experimental limitation discussed in Section 5.5 of the report.

---

## Session 1 Run Log (from lab diary, 6 November 2025)

| Run | Field range | Notes |
|-----|-------------|-------|
| Run 1 | 0 → +2 T | First sweep; live plot backend issue; CSV saved but lost |
| Run 2 | +2 → −2 T | Bipolar sweep; T ≈ 1.9 K; oscillations visible |
| Run 3 | −2 → +2 T | T ≈ 2.04 K; excellent SdH oscillations both polarities |
| Run 4 | +2 → 0 → +2 T | Partial sweep |
| Run 5 | 0 → +4 T | Vxy overranged (lock-in saturation); sensitivity adjusted |
| Run 6 | +4 → −4 T | Vibration event mid-sweep; step discontinuity visible |
| Run 7 | −4 → 0 T | Recovery after vibration |
| Run 8 | 0 → +5 T | Vxy overranged at 10 V; range changed to 20 V mid-run |
| Run 9 | +5 → 0 T | Clearest plateau from Session 1; used for initial analysis |

Supervisor sign-off on Session 1 results: *CJ Mellor, 10 November 2025
(discussed during Monday meeting).*
