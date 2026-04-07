#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# --- Input file ---
input_file = "W_ATLAS_dNNLO_xMuR_1.00_xMuF_1.00.histo"

# --- Load data ---
data = np.loadtxt(input_file)

# --- Unique observables ---
observables = np.unique(data[:, 0])

# --- Plot each observable ---
for obs in observables:
    obs_data = data[data[:, 0] == obs]
    x = obs_data[:, 2]
    y = obs_data[:, 3]
    yerr = obs_data[:, 4]

    # HEP-style plotting
    plt.errorbar(
        x, y, yerr=yerr,
        fmt='o', markersize=0.5,  # small dots
        elinewidth=1, capsize=2, capthick=1,  # thin error bars with caps
        label=f'Observable {int(obs)}'
    )

    plt.xlim(150, 800)  # <-- set the x-axis range here

    plt.xlabel('x')
    plt.ylabel('Value')
    plt.title(f'Observable {int(obs)}')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"observable_{int(obs)}.pdf")
    plt.clf()

print("Plots saved as observable_1.pdf, observable_2.pdf, ...")



