"""
qhe_data_acquisition.py
=======================
Real-time data acquisition script for the Quantum Hall Effect (QHE) experiment.

Instruments:
    - 2x Stanford Research Systems SR830 Lock-in Amplifiers (via GPIB)
        GPIB0::8  →  Lower SR830, measures longitudinal voltage Vxx
        GPIB0::9  →  Upper SR830, measures Hall voltage Vxy
    - Keithley 2100 Digital Multimeter (via USB)
        Reads a field-proportional voltage VB from the magnet's current monitor output

Experiment overview:
    The superconducting magnet ramp rate, start and stop field are set directly on
    the magnet controller. This script runs continuously from magnet ramp start,
    logging Vxx and Vxy (X, Y, R, theta) alongside VB at a fixed sampling interval.
    Acquisition is stopped manually with Ctrl+C, which triggers a clean shutdown and
    flushes the CSV before closing all instrument connections.

Output:
    CSV file written to ./data/QHE_mergedDATA_<timestamp>.csv
    Live matplotlib plot of Vxx and Vxy vs time and vs VB field voltage.

Usage:
    python qhe_data_acquisition.py

Dependencies:
    pyvisa, numpy, matplotlib
"""

import pyvisa as visa
import numpy as np
import matplotlib.pyplot as plt
import time, csv, datetime, os


# ---------------------------------------------------------------------------
# Connect to VISA resource manager and list available resources
# ---------------------------------------------------------------------------

rm = visa.ResourceManager()
print(rm.list_resources())


# ---------------------------------------------------------------------------
# Open instrument connections
# ---------------------------------------------------------------------------

# Lower SR830 lock-in amplifier — measures longitudinal voltage Vxx
LI_lower = rm.open_resource('GPIB0::8::INSTR', read_termination='\r', send_end=True)
LI_lower.timeout = 5000  # ms

# Upper SR830 lock-in amplifier — measures Hall voltage Vxy
LI_upper = rm.open_resource('GPIB0::9::INSTR', read_termination='\r', send_end=True)
LI_upper.timeout = 5000  # ms

# Keithley 2100 multimeter — reads field-proportional voltage VB via USB
mm1 = rm.open_resource('USB0::0x05E6::0x2100::1174217::INSTR')
mm1.timeout = 5000  # ms


# ---------------------------------------------------------------------------
# Route SR830s to GPIB interface and confirm all instrument IDs
# ---------------------------------------------------------------------------

# OUTX1 switches SR830 output interface to GPIB (required after power cycle)
LI_lower.write('OUTX1')
LI_upper.write('OUTX1')

print("Lower SR830 ID :", LI_lower.query('*IDN?').strip())
print("Upper SR830 ID :", LI_upper.query('*IDN?').strip())
print("Keithley 2100  :", mm1.query('*IDN?').strip())


# ---------------------------------------------------------------------------
# Configure SR830 lock-in amplifiers
# ---------------------------------------------------------------------------

# Reset both to a known state before configuration
LI_lower.write('REST')
LI_upper.write('REST')

# --- Lower SR830 (Vxx) ---
LI_lower.write('FREQ70')       # Set reference frequency to 70 Hz
print(LI_lower.query('FREQ?')) # Confirm frequency setting
LI_lower.write('OFLT9')        # Time constant: code 9 (see SR830 manual for mapping)
print(LI_lower.query('OFLT?')) # Confirm time constant
#LI_lower.write('SENS20')      # Sensitivity: commented out, set manually on front panel

# --- Upper SR830 (Vxy) ---
LI_upper.write('FREQ70')       # Match reference frequency to lower unit
print(LI_upper.query('FREQ?')) # Confirm frequency setting
LI_upper.write('OFLT9')        # Match time constant to lower unit
print(LI_upper.query('OFLT?')) # Confirm time constant
#LI_upper.write('SENS20')      # Sensitivity: commented out, set manually on front panel


# ---------------------------------------------------------------------------
# Configure Keithley 2100 multimeter
# ---------------------------------------------------------------------------

mm1.write('SENS:FUNC "VOLT:DC"')         # DC voltage measurement mode
mm1.write('SENS:VOLT:DC:RANG:AUTO ON')   # Autorange (switch to fixed range if desired)
mm1.write('SENS:VOLT:DC:NPLC 10')        # Integration time = 10 PLC → strong 50 Hz rejection (~0.2 s/reading)
mm1.write('SYST:AZER ONCE')              # Perform autozero once now for accuracy
mm1.write('SYST:AZER:STAT OFF')          # Disable autozero during loop for speed

print("Keithley configured: DCV Autorange, NPLC=10, AZ once")


# ---------------------------------------------------------------------------
# Acquisition parameters
# ---------------------------------------------------------------------------

INTERVAL_S  = 1.0    # Target time between samples (seconds)
NUM_SAMPLES = 2500   # Max samples before loop exits (~42 min at 1 s/sample)
                     # In practice, Ctrl+C is used to stop acquisition manually

print(f"Loop interval = {INTERVAL_S:.3f} s, max rows = {NUM_SAMPLES}")

t0    = time.time()  # Acquisition start time (epoch)
count = 0            # Running sample counter


# ---------------------------------------------------------------------------
# Live plot setup — updates every iteration during acquisition
# ---------------------------------------------------------------------------

plt.ion()  # Enable interactive mode for live updating

# Buffers for live plot data
t_list   = []
vxx_list = []
vxy_list = []
vb_list  = []

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8))

# Top subplot: Vxx and Vxy vs elapsed time
ax1.set_xlim(0, NUM_SAMPLES * INTERVAL_S)
line_vxx_time, = ax1.plot([], [], 'r.-', label='Vxx(t)')
line_vxy_time, = ax1.plot([], [], 'b.-', label='Vxy(t)')
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Voltage (V)")
ax1.set_title("Live: Vxx & Vxy vs Time")
ax1.legend()

# Bottom subplot: Vxx and Vxy vs field proxy voltage VB
line_vxx_vb, = ax2.plot([], [], 'r.-', label='Vxx(VB)')
line_vxy_vb, = ax2.plot([], [], 'b.-', label='Vxy(VB)')
ax2.set_xlabel("VB field (proportional voltage)")
ax2.set_ylabel("Voltage (V)")
ax2.set_title("Live: Vxx & Vxy vs VB")
ax2.legend()

plt.tight_layout()
plt.pause(0.01)


# ---------------------------------------------------------------------------
# Create output CSV file
# ---------------------------------------------------------------------------

os.makedirs('data', exist_ok=True)
stamp    = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
filename = f"QHE_mergedDATA_{stamp}.csv"
filepath = os.path.join('data', filename)

f      = open(filepath, 'w', newline='')
writer = csv.writer(f)

# Column headers — all voltages in volts, theta in degrees, time in seconds
writer.writerow(['t_s', 'Vxx_X_V', 'Vxx_Y_V', 'Vxx_R_V',
                 'Vxx_theta_deg', 'Vxy_X_V', 'Vxy_Y_V',
                 'Vxy_R_V', 'Vxy_theta_deg', 'VB_V'])

print(f"\nLogging to: {filepath}\n")


# ---------------------------------------------------------------------------
# Main acquisition loop
# SNAP queries both lock-ins simultaneously for X, Y, R, theta, then reads VB.
# Loop is paced to INTERVAL_S; stop at any time with Ctrl+C.
# ---------------------------------------------------------------------------

try:
    for n in range(NUM_SAMPLES):
        t_iter_start = time.time()

        # --- Lower SR830 (Vxx): simultaneous read of X, Y, R, theta ---
        sL = LI_lower.query('SNAP?1,2,3,4').strip()  # returns "X,Y,R,theta"
        xL_str, yL_str, rL_str, thL_str = sL.split(',')
        vxx_x  = float(xL_str)   # In-phase component
        vxx_y  = float(yL_str)   # Quadrature component
        vxx_r  = float(rL_str)   # Magnitude
        vxx_th = float(thL_str)  # Phase angle (degrees)

        # --- Upper SR830 (Vxy): simultaneous read of X, Y, R, theta ---
        sU = LI_upper.query('SNAP?1,2,3,4').strip()  # returns "X,Y,R,theta"
        xU_str, yU_str, rU_str, thU_str = sU.split(',')
        vxy_x  = float(xU_str)
        vxy_y  = float(yU_str)
        vxy_r  = float(rU_str)
        vxy_th = float(thU_str)  # degrees

        # --- Keithley 2100: read field-proportional voltage VB ---
        try:
            vb = float(mm1.query('READ?'))
        except Exception as e:
            print("[Keithley read error]", e)
            vb = float('nan')  # Preserve row alignment on read failure

        # Elapsed time since acquisition start
        t = time.time() - t0
        count += 1

        # Print to console
        print(f"{count:4d}  t={t:7.3f}s  "
              f"Vxx[X]={vxx_x:.6e}  Vxx[Y]={vxx_y:.2e}  Vxx[R]={vxx_r:.6e}  Vxx[θ]={vxx_th:.2f}°   "
              f"Vxy[X]={vxy_x:.6e}  Vxy[Y]={vxy_y:.2e}  Vxy[R]={vxy_r:.6e}  Vxy[θ]={vxy_th:.2f}°   "
              f"VB={vb:.6e} V")

        # Write row to CSV
        writer.writerow([f"{t:.6f}",
                         f"{vxx_x:.9e}", f"{vxx_y:.9e}", f"{vxx_r:.9e}", f"{vxx_th:.6f}",
                         f"{vxy_x:.9e}", f"{vxy_y:.9e}", f"{vxy_r:.9e}", f"{vxy_th:.6f}",
                         f"{vb:.9e}"])

        # --- Update live plots ---
        t_list.append(t)
        vxx_list.append(vxx_x)
        vxy_list.append(vxy_x)
        vb_list.append(vb)

        # Time-domain plot
        line_vxx_time.set_xdata(t_list)
        line_vxx_time.set_ydata(vxx_list)
        line_vxy_time.set_xdata(t_list)
        line_vxy_time.set_ydata(vxy_list)
        ax1.relim()
        ax1.autoscale_view()

        # Field-domain plot
        line_vxx_vb.set_xdata(vb_list)
        line_vxx_vb.set_ydata(vxx_list)
        line_vxy_vb.set_xdata(vb_list)
        line_vxy_vb.set_ydata(vxy_list)
        ax2.relim()
        ax2.autoscale_view()

        fig.canvas.draw()
        fig.canvas.flush_events()

        # Pace the loop to INTERVAL_S, accounting for acquisition time
        elapsed  = time.time() - t_iter_start
        to_sleep = INTERVAL_S - elapsed
        if to_sleep > 0:
            time.sleep(to_sleep)

    print("\nSaved:", filepath)

except KeyboardInterrupt:
    print("\n[Interrupted] Saved so far:", filepath)

# ---------------------------------------------------------------------------
# Clean up — flush CSV, close all instruments and plot
# ---------------------------------------------------------------------------

finally:
    try:
        f.flush()
        f.close()
    except:
        pass

    LI_lower.clear(); LI_lower.close()
    LI_upper.clear(); LI_upper.close()
    mm1.clear();      mm1.close()
    rm.close()

    plt.ioff()
    plt.show()


# ---------------------------------------------------------------------------
# Post-run plotting (optional) — uncomment to replay data from saved CSV
# ---------------------------------------------------------------------------

"""
data = np.genfromtxt(filepath, delimiter=",", skip_header=1)

t   = data[:, 0]   # Elapsed time (s)
vxx = data[:, 1]   # Vxx X-component (longitudinal)
vxy = data[:, 5]   # Vxy X-component (Hall)
vb  = data[:, 9]   # Field-proportional voltage VB

plt.figure()
plt.plot(t, vxx, 'r.-', label="Vxx (longitudinal)")
plt.plot(t, vxy, 'b.-', label="Vxy (Hall)")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.title("Raw Hall and Longitudinal Signals vs Time")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.plot(vb, vxx, 'r.', label="Vxx vs VB")
plt.plot(vb, vxy, 'b.', label="Vxy vs VB")
plt.xlabel("VB (field-proportional voltage)")
plt.ylabel("Voltage (V)")
plt.title("Hall and Longitudinal Voltages vs Field Sweep")
plt.legend()
plt.grid()
plt.show()
"""
