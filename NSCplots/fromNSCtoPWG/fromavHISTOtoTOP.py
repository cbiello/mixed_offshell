import re

input_file = "input.dat"

# Mapping of integer index to observable name
obs_keys = ["rate", "massT", "pT-lep", "pT-nu", "mass", "y-lep"]
idx_to_name = {i + 1: name for i, name in enumerate(obs_keys)}

# Scaling factors
scales = {
    "rate": 1000.0,
    "massT": 10000.0,
    "pT-lep": 5000.0,
    "pT-nu": 5000.0,
    "mass": 10000.0,
    "y-lep": 200.0
}

outputs = {
    "CEN": "output_CEN.top",
    "MIN": "output_MIN.top",
    "MAX": "output_MAX.top"
}

observables = {name: [] for name in obs_keys}

# --- PARSING SECTION ---
with open(input_file) as f:
    for line in f:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        parts = stripped.split()
        
        # Check for minimum required columns (ObsID, BinID, X, Val, Err)
        if len(parts) < 5:
            continue

        try:
            obs_id = int(parts[0])
            current = idx_to_name.get(obs_id)
        except ValueError:
            continue

        if current:
            scale = scales[current]
            
            # Basic values (always present in 5-column or 9-column format)
            x_val   = float(parts[2])
            cen_val = float(parts[3]) / scale
            cen_err = float(parts[4]) / scale

            # If the file has 9 columns, use the specific MIN/MAX data
            if len(parts) >= 9:
                min_val = float(parts[5]) / scale
                min_err = float(parts[6]) / scale
                max_val = float(parts[7]) / scale
                max_err = float(parts[8]) / scale
            else:
                # Fallback: if MIN/MAX aren't in the file, use CEN values
                min_val, min_err = cen_val, cen_err
                max_val, max_err = cen_val, cen_err

            observables[current].append({
                "CEN": [x_val, cen_val, cen_err],
                "MIN": [x_val, min_val, min_err],
                "MAX": [x_val, max_val, max_err]
            })

# --- HELPER FUNCTION ---
def build_bins(values, key):
    if not values: return []
    bins = []
    for i in range(len(values) - 1):
        L, V, E = values[i][key]
        R = values[i+1][key][0]
        bins.append([L, R, V, E])
    
    # Last bin extrapolation
    if len(values) > 1:
        last_L, V, E = values[-1][key]
        width = last_L - values[-2][key][0]
        bins.append([last_L, last_L + width, V, E])
    else:
        last_L, V, E = values[0][key]
        bins.append([last_L, last_L + 1.0, V, E])
    return bins

# --- WRITING SECTION ---
for mode, outfile in outputs.items():
    with open(outfile, "w") as out:
        for i, name in enumerate(obs_keys):
            # Check if we have data for this observable
            data_list = observables[name]
            
            # Write header
            out.write(f"# {name} index      {i}\n")
            
            if data_list:
                current_bins = build_bins(data_list, mode)
                for L, R, V, E in current_bins:
                    out.write(f" {L:0.8E} {R:0.8E} {V:0.5E} {E:0.3E}\n")
            
            out.write("\n")
    print(f"Done! Created {outfile}")
