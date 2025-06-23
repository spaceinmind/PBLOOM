import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from scipy import integrate, optimize

# --- Parameters ---
mlog = np.linspace(-6, 4, 10000)
masses = 10**mlog * 2.0e+30  # Lens mass in kg
Dt_values = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]  # Time delays in seconds
z = [0.1, 0.5]  # z[0] not used in this version; z[1] = source redshift
Rf = 5  # Flux ratio threshold
SAVEFIG = False

# --- Constants ---
c_ms = 2.99792458e8  # m/s
c_kms = 2.99792458e5  # km/s
G = 6.67408e-11  # m^3 kg^-1 s^-2
H0 = Planck15.H0.value  # km/s/Mpc
Om0 = Planck15.Om0

# --- Core Functions ---
def ymin_func(y, zL, ML, Dt):
    term = y / 2 * np.sqrt(y**2 + 4)
    term += np.log10((np.sqrt(y**2 + 4) + y) / (np.sqrt(y**2 + 4) - y))
    return (4 * G * ML / c_ms**3) * (1 + zL) * term - Dt

def y_min(zL, ML, Dt):
    try:
        return optimize.newton(ymin_func, 1.0, args=(zL, ML, Dt))
    except (RuntimeError, OverflowError):
        return np.nan

def tau0(zL, zS, ML, Dt, fDM=1.0, O_c=0.24):
    DL = Planck15.angular_diameter_distance(zL).value  # Mpc
    DS = Planck15.angular_diameter_distance(zS).value  # Mpc
    DLS = Planck15.angular_diameter_distance_z1z2(zL, zS).value  # Mpc
    H_z = H0 * np.sqrt(Om0 * (1 + zL)**3 + (1 - Om0))  # km/s/Mpc
    ymax = np.sqrt((1 + Rf) / np.sqrt(Rf) - 2)
    y_min_val = y_min(zL, ML, Dt)
    if np.isnan(y_min_val) or y_min_val > ymax:
        return 0.0
    prefactor = (3 / 2) * fDM * O_c * H0**2 * DL * DLS / (c_kms * H_z * DS)
    return prefactor * (1 + zL)**2 * (ymax**2 - y_min_val**2)

def integrated_tau(zS, ML, Dt):
    # Integrate over zL ∈ [0, zS], with zS fixed (matching original script)
    return integrate.quad(tau0, 0, zS, args=(zS, ML, Dt))[0]

# --- Plotting ---
plt.figure(figsize=(19, 9))

for Dt in Dt_values:
    yvals = []
    for ML in masses:
        tau_val = integrated_tau(z[1], ML, Dt)
        yvals.append(tau_val)
    yvals = np.array(yvals)
    plt.plot(masses / 2.0e+30, yvals, linewidth=3,
             label=rf'$\bar{{\Delta t}}={Dt}$ s')

# --- Plot Appearance ---
plt.xscale('log')
# plt.yscale('log')  # Optional, depending on the preferred scale

plt.xlim(1e-8, 1e3)
plt.ylim(0, 0.012)

plt.title(f'z$_{{cut}}$ = {z[1]}', fontsize=30)
plt.xlabel(r'Lensing mass, M$_L$ ($M_{\odot}$)', fontsize=40)
plt.ylabel(r'Integrated optical depth, $\bar{\tau}$', fontsize=40)

plt.legend(loc='upper left', fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

if SAVEFIG:
    plt.savefig(f'tau_vs_mass_zcut_{z[1]}.png', bbox_inches='tight')
else:
    plt.show()
