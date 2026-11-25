import re

input_file = "input.histo"

# Output files
outputs = {
    "CEN": "output_CEN.top",
    "MIN": "output_MIN.top",
    "MAX": "output_MAX.top"
}

# Storage for each observable
observables = {
    "rate": [],
    "massT": [],
    "pT-lep": [],
    "pT-nu": [],
    "mass": [],
    "y-lep": []
}

current = None

with open(input_file) as f:
    for line in f:
        stripped = line.strip()

        # Detect observable section
        if stripped.startswith("# ----------Obs="):
            name = stripped.split("=")[1].split()[0]
            current = name
            continue
        
        # Skip comments and empty lines
        if stripped.startswith("#") or not stripped:
            continue
        
        if current is None:
            continue  

        parts = line.split()
        if len(parts) < 9:
            continue

        if current == "rate": scale = 1000.0
        if current == "massT": scale = 10000.0
        if current == "pT-lep": scale = 5000.0
        if current == "pT-nu": scale = 5000.0
        if current == "mass": scale = 10000.0
        if current == "y-lep": scale = 200.0
        
        # Read values for CEN, MIN, MAX
        cen_val = float(parts[3]) / scale
        cen_err = float(parts[4]) / scale
        min_val = float(parts[5]) / scale
        min_err = float(parts[6]) / scale
        max_val = float(parts[7]) / scale
        max_err = float(parts[8]) / scale

        observables[current].append({
            "CEN": [float(parts[2]), cen_val, cen_err],
            "MIN": [float(parts[2]), min_val, min_err],
            "MAX": [float(parts[2]), max_val, max_err]
        })

def build_bins(values, key):
    """Convert (L, val, err) to (L, R, val, err)."""
    bins = []
    for i in range(len(values)-1):
        L = values[i][key][0]
        R = values[i+1][key][0]
        V = values[i][key][1]
        E = values[i][key][2]
        bins.append([L, R, V, E])
    # last bin: extrapolate width
    last_L = values[-1][key][0]
    width = values[-1][key][0] - values[-2][key][0]
    bins.append([last_L, last_L + width, values[-1][key][1], values[-1][key][2]])
    return bins

# Loop over CEN, MIN, MAX outputs
for mode, outfile in outputs.items():
    with open(outfile, "w") as out:

        # ----- RATE -----
        out.write("# rate index      0\n")
        rate_bins = build_bins(observables["rate"], mode)
        for L, R, V, E in rate_bins:
            out.write(f" {L:0.8E} {R:0.8E} {V:0.5E} {E:0.3E}\n")

        out.write("\n")

        # ----- MASS T -----
        out.write("# massT index      1\n")
        massT_bins = build_bins(observables["massT"], mode)
        for L, R, V, E in massT_bins:
            out.write(f" {L:0.8E} {R:0.8E} {V:0.5E} {E:0.3E}\n")
            
        out.write("\n")

        # ----- PT lepton -----
        out.write("# pT-lep index      2\n")
        pTlep_bins = build_bins(observables["pT-lep"], mode)
        for L, R, V, E in pTlep_bins:
            out.write(f" {L:0.8E} {R:0.8E} {V:0.5E} {E:0.3E}\n")
            
        out.write("\n")

        # ----- PT neutrino -----
        out.write("# pT-nu index      3\n")
        pTnu_bins = build_bins(observables["pT-nu"], mode)
        for L, R, V, E in pTnu_bins:
            out.write(f" {L:0.8E} {R:0.8E} {V:0.5E} {E:0.3E}\n")

        out.write("\n")

        
        # ----- PT neutrino -----
        out.write("# mass index      4\n")
        mass_bins = build_bins(observables["mass"], mode)
        for L, R, V, E in mass_bins:
            out.write(f" {L:0.8E} {R:0.8E} {V:0.5E} {E:0.3E}\n")

        out.write("\n")

            
        # ----- PT neutrino -----
        out.write("# y-lep index      5\n")
        ylep_bins = build_bins(observables["y-lep"], mode)
        for L, R, V, E in ylep_bins:
            out.write(f" {L:0.8E} {R:0.8E} {V:0.5E} {E:0.3E}\n")

            
            
    print(f"Written {mode} output to {outfile}")
