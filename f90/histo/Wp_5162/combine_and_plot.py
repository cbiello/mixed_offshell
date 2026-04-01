import matplotlib.pyplot as plt
from collections import defaultdict


# --- Parameters ---
scale = 1.00
channel = "ns"
section = "rr"

print(f">>> Processing scale = {scale}")
print(f"  -> Channel = {channel}")
print(f"     * Section = {section}")

# --- List of histogram files ---
hist_files = sorted([f for f in os.listdir(".") if f.endswith(".histo")])
print(" Averaging files")
for i, f in enumerate(hist_files):
    print(f"file {i+1}: {f}")

# --- Read and average histograms ---
hists = [read_histogram(f) for f in hist_files]
avg_hist = hists[0]
for h in hists[1:]:
    avg_hist = add_histograms(avg_hist, h)
write_histogram("W_ATLAS_avg_nnlo_rr_ns_xMuR_1.00_xMuF_1.00.histo", avg_hist)

# --- Sum histograms (if needed) ---
summed_hist = avg_hist  # no rescaling
write_histogram("W_ATLAS_summed_nnlo_ns_xMuR_1.00_xMuF_1.00.histo", summed_hist)

# --- Write final NNLO histogram ---
final_hist = summed_hist
write_histogram("W_ATLAS_dNNLO_xMuR_1.00_xMuF_1.00.histo", final_hist)

print(">>> Done creating histograms")

# --- Plotting ---
print(">>> Creating plots")

data = defaultdict(list)
for line in final_hist.splitlines():
    if line.strip() == "" or line.startswith("#"):
        continue
    cols = line.split()
    obs = int(cols[0])       # Observable ID
    x = float(cols[2])       # x-axis
    y = float(cols[3])       # y-value
    err = float(cols[4])     # y-error
    data[obs].append((x, y, err))

for obs, points in data.items():
    points.sort(key=lambda p: p[0])
    x_vals, y_vals, y_errs = zip(*points)
    
    plt.figure()
    plt.errorbar(x_vals, y_vals, yerr=y_errs, fmt='o', capsize=3)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Observable {obs}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"W_ATLAS_dNNLO_obs{obs}.png")
    plt.close()
    print(f"Saved plot for observable {obs}: W_ATLAS_dNNLO_obs{obs}.png")

print(">>> All plots created")






