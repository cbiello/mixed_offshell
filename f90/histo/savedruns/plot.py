import numpy as np
import matplotlib.pyplot as plt

filename = "W_ATLAS_lo_ns_na_dynscale_xMuR_1.00_xMuF_1.00.histo"

data = np.loadtxt(filename)

x = data[:, 2]

# Central value
y_central = data[:, 3]

# Lower and upper variations (band)
y_lower = data[:, 5]  # lower
y_upper = data[:, 7]  # upper

plt.figure(figsize=(8, 6))

# Plot central curve with smaller dots
plt.plot(x, y_central, 'o-', color='blue', label='Central', markersize=4)

# Fill between lower and upper for the band
plt.fill_between(x, y_lower, y_upper, color='blue', alpha=0.3, label='Uncertainty band')

plt.ylabel("dσ / dMₜ")
plt.xlabel("Mₜ [GeV]")  # Updated y-axis label
plt.title("Observable 2: massT")
plt.legend()
plt.grid(True)

# Limit x-axis from 100 to 1000
plt.xlim(200, 700)

# Limit y-axis from 100 to 1000
plt.ylim(-1, 200)

plt.tight_layout()

# Save figures
plt.savefig("massT_obs2.pdf")
plt.savefig("massT_obs2.png", dpi=300)

plt.show()
