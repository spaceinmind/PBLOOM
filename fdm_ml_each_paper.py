import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from astropy.cosmology import Planck15
import matplotlib.patheffects as patheffects

# --- Constants ---
c_kms = 2.99792458e5           # km/s
c_ms = 2.99792458e8            # m/s
G = 6.67408e-11                # m^3/kg/s^2
H0 = Planck15.H0.value         # km/s/Mpc
Om0 = Planck15.Om0
Rf = 5                         # flux ratio threshold

SAVEFIG = False  # Change to True if you want to save the figure

# --- Data Loading ---
def load_constraints():
    base = 'dmconstrain/'
    return {
        'EGB': pd.read_csv(base + 'EGB.csv'),
        'EROS': pd.read_csv(base + 'EROS.csv'),
        'FIRAS': pd.read_csv(base + 'FIRAS.csv', engine='python', skipfooter=10),
        'MACHO': pd.read_csv(base + 'MACHO.csv'),
        'NScapture': pd.read_csv(base + 'NScapture.csv', skiprows=range(1, 10)),
        'WB': pd.read_csv(base + 'WB.csv', engine='python', skipfooter=10),
        'WMAP3': pd.read_csv(base + 'WMAP3.csv', engine='python', skipfooter=20),
        'femtolensing': pd.read_csv(base + 'femtolensing.csv'),
        'HSC': pd.read_csv(base + 'HSC.csv')
    }

# --- Core Functions ---
def ymin_func(y, zL, ML, Dt):
    return (4 * G * ML / c_ms**3) * (1 + zL) * (
        y / 2 * np.sqrt(y**2 + 4) + np.log10((np.sqrt(y**2 + 4) + y) / (np.sqrt(y**2 + 4) - y))
    ) - Dt

def y_min(zL, ML, Dt):
    try:
        return optimize.newton(ymin_func, 1.0, args=(zL, ML, Dt))
    except RuntimeError:
        return np.nan

def tau0(zL, zS, ML, Dt, fDM=1.0, Oc=0.24):
    H_z = H0 * np.sqrt(Om0 * (1 + zL)**3 + (1 - Om0))  # Flat ΛCDM
    DL = Planck15.angular_diameter_distance(zL).value
    DS = Planck15.angular_diameter_distance(zS).value
    DLS = Planck15.angular_diameter_distance_z1z2(zL, zS).value
    ymax = np.sqrt((1 + Rf) / np.sqrt(Rf) - 2)
    y_min_val = y_min(zL, ML, Dt)
    if np.isnan(y_min_val):
        return 0.0
    return (3 / 2) * fDM * Oc * H0**2 * DL * DLS / (c_kms * H_z * DS) * (1 + zL)**2 * (ymax**2 - y_min_val**2)

def total_tau(z_list, frb_weights, ML, Dt):
    return sum(
        tau0(z_list[i], z_list[i + 1], ML, Dt) * frb_weights[i]
        for i in range(len(z_list) - 1)
    )

# --- Plotting Functions ---
def plot_tau_constraints(Dt_list, mvals, zvals, frb_weights):
    for Dt in Dt_list:
        yaxis = []
        for ML in mvals:
            tau = total_tau(zvals, frb_weights, ML, Dt)
            val = 1.0 / tau if tau > 0 else np.nan
            yaxis.append(val)
        
        yaxis = np.array(yaxis)
        for i, y in enumerate(yaxis):
            if y > 0 and not np.isnan(y):
                x = np.concatenate([np.full(100, mvals[i]), mvals[i:]]) / 2e30
                y_fill = np.concatenate([np.linspace(y, 1, 100), yaxis[i:]])
                plt.fill_between(x, y_fill, 1, linewidth=2,
                                 label=rf'$\bar{{\Delta t}}={Dt}$s',
                                 alpha=0.9, edgecolor='black')
                break

def plot_existing_constraints(data):
    constraints_main = ['EROS', 'WB', 'MACHO']
    constraints_secondary = ['WMAP3', 'FIRAS']
    color_main = ['black', 'dimgray', 'gray']
    color_secondary = ['darkgrey', 'lightgrey']

    for label, color in zip(constraints_secondary, color_secondary):
        df = data[label]
        mid_idx = int(len(df) / 2 - 12)
        plt.text(df.iloc[mid_idx, 0], df.iloc[mid_idx, 1], label,
                 color=color, fontsize=20,
                 path_effects=[patheffects.withStroke(linewidth=1, foreground='black')])
        plt.fill_between(df.iloc[:, 0], df.iloc[:, 1], 1, linestyle='--',
                         color=color, alpha=0.2, hatch='///', edgecolor='black')

    for label, color in zip(constraints_main, color_main):
        df = data[label]
        mid_idx = int(len(df) / 2)
        plt.text(df.iloc[mid_idx, 0], df.iloc[mid_idx, 1], label,
                 alpha=0.4, color=color, fontsize=20,
                 path_effects=[patheffects.withStroke(linewidth=1, foreground='black')])
        plt.fill_between(df.iloc[:, 0], df.iloc[:, 1], 1, linestyle='--',
                         color=color, alpha=0.2, hatch='///', edgecolor='black')

    hsc = data['HSC']
    plt.text(hsc.iloc[20, 0], hsc.iloc[20, 1], 'HSC', color='whitesmoke',
             fontsize=20, alpha=0.4,
             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')])
    plt.fill_between(hsc.iloc[:, 0], hsc.iloc[:, 1], 1, linestyle='--',
                     color='whitesmoke', alpha=0.2, hatch='///', edgecolor='black')

# --- Main Execution ---
if __name__ == '__main__':
    constraints_data = load_constraints()
    z_values = [0.0682595, 0.5679634, 1.0676673, 1.5673712, 2.0670751]
    frb_weights = [1071.9, 419.4, 139.8, 69.9]  # number of FRBs in each redshift bin

    mlog = np.linspace(-7, 5, 1000)
    mvals = 10**mlog * 2e30  # lensing mass in kg
    Dt_values = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]

    plt.figure(figsize=(19, 9))
    plot_tau_constraints(Dt_values, mvals, z_values, frb_weights)
    plot_existing_constraints(constraints_data)

    # Add text
    plt.text(0.005, 0.014, 'BURSTT-2048', color='red', alpha=0.9, fontsize=50,
             path_effects=[patheffects.withStroke(linewidth=1.5, foreground='black')])

    # Plot settings
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-7, 1e5)
    plt.ylim(1e-2, 1)
    plt.xlabel(r'Lensing mass, M$_L$ ($M_{\odot}$)', fontsize=40)
    plt.ylabel(r'Fraction of PBH in dark matter, f$_{DM}$', fontsize=35)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.legend(loc='lower left', fontsize=20)

    if SAVEFIG:
        plt.savefig('tau_constraints_plot.png', bbox_inches='tight')
    else:
        plt.show()
