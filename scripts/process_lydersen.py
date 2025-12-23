import json
import os
import numpy as np
from scipy.interpolate import PchipInterpolator, UnivariateSpline

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    input_path = os.path.join(project_root, 'data', 'lydersen_extracted_vals_raw.json')
    output_path = os.path.join(project_root, 'data', 'lydersen_normal.json')

    print(f"Reading from {input_path}")
    try:
        with open(input_path, 'r') as f:
            raw_data = json.load(f)
    except FileNotFoundError:
        print(f"Error: Could not find {input_path}")
        return

    processed_data = {}

    # Normalization Rules
    # Tr sets
    start_zero_tr = ["Tr_0.3", "Tr_0.4", "Tr_0.5", "Tr_0.6"]
    end_ten_tr = ["Tr_0.3", "Tr_0.4", "Tr_0.5", "Tr_0.6", "Tr_0.7", "Tr_0.8", "Tr_0.9", "Tr_1.0"]
    end_four_tr = ["Tr_0.97", "Tr_0.99"]
    
    # Process each dataset
    for name, points in raw_data.items():
        if not points:
            continue
            
        pts = np.array(points)
        pr = pts[:, 0]
        rho = pts[:, 1]

        # 1. Normalization (Fixing Endpoints)
        
        # Start at Pr = 0
        if name in start_zero_tr:
            # If the first point is not 0, insert 0.
            if pr[0] > 0.001:
                # Prepend 0. Use the same density as the first point (incompressible assumption)
                pr = np.insert(pr, 0, 0.0)
                rho = np.insert(rho, 0, rho[0])
            else:
                pr[0] = 0.0

        # End at Pr = 10
        if name in end_ten_tr:
            if pr[-1] < 9.99:
                slope = (rho[-1] - rho[-2]) / (pr[-1] - pr[-2])
                new_rho = rho[-1] + slope * (10.0 - pr[-1])
                pr = np.append(pr, 10.0)
                rho = np.append(rho, new_rho)
            else:
                pr[-1] = 10.0

        # End at Pr = 4
        if name in end_four_tr:
            if pr[-1] < 3.99:
                slope = (rho[-1] - rho[-2]) / (pr[-1] - pr[-2])
                new_rho = rho[-1] + slope * (4.0 - pr[-1])
                pr = np.append(pr, 4.0)
                rho = np.append(rho, new_rho)
            else:
                pr[-1] = 4.0

        # Saturation Curve
        if name == "sat":
            # Start at Pr = 0
            if pr[0] > 0.001:
                pr = np.insert(pr, 0, 0.0)
                rho = np.insert(rho, 0, rho[0]) # Extrapolate/Duplicate
            else:
                pr[0] = 0.0
            
            # End at Pr = 1, rho_r = 1
            # Force the last point to be (1, 1)
            # If the last point is close to 1, replace it. Else append.
            if pr[-1] < 0.99:
                pr = np.append(pr, 1.0)
                rho = np.append(rho, 1.0)
            else:
                pr[-1] = 1.0
                rho[-1] = 1.0

        # 2. Smoothing with PCHIP
        # We create a PCHIP interpolator
        # Ensure strictly increasing Pr for PCHIP
        # Remove duplicates if any
        unique_indices = np.unique(pr, return_index=True)[1]
        pr_clean = pr[unique_indices]
        rho_clean = rho[unique_indices]
        
        # Sort again just in case unique messed up order (it shouldn't if sorted)
        sorted_indices = np.argsort(pr_clean)
        pr_clean = pr_clean[sorted_indices]
        rho_clean = rho_clean[sorted_indices]

        try:
            if name == "sat":
                interpolator = PchipInterpolator(pr_clean, rho_clean)
            else:
                # Use UnivariateSpline for smoothing to remove jaggedness
                # Weight endpoints heavily to preserve the fixed boundary conditions
                weights = np.ones(len(pr_clean))
                weights[0] = 1e6
                weights[-1] = 1e6
                # s=0.005 provides mild smoothing. 
                interpolator = UnivariateSpline(pr_clean, rho_clean, w=weights, s=0.005)
        except Exception as e:
            print(f"Failed to interpolate {name}: {e}")
            processed_data[name] = points # Fallback
            continue

        # 3. Resampling
        # Generate a clean set of points.
        # What resolution?
        # The original data has ~80 points for range 0-10.
        # Let's use a fixed step size, e.g., 0.1 for the main curves.
        # For sat curve (0-1), maybe 0.01?
        
        if name == "sat":
            x_new = np.linspace(0, 1, 101) # 0.01 step
        elif name in start_zero_tr:
            # These start at 0 and go to 10
            x_new = np.linspace(0, 10, 201) # 0.05 step
        elif name in end_ten_tr:
            # These (0.7, 0.8, 0.9, 1.0) end at 10 but start at their data start (intersection with sat)
            start_pr = pr_clean[0]
            num_points = int(np.ceil((10.0 - start_pr) / 0.05)) + 1
            x_new = np.linspace(start_pr, 10.0, num_points)
        elif name in end_four_tr:
            # These (0.97, 0.99) end at 4.0
            start_pr = pr_clean[0]
            num_points = int(np.ceil((4.0 - start_pr) / 0.05)) + 1
            x_new = np.linspace(start_pr, 4.0, num_points)
        else:
            # For 0.95, just smooth the existing range
            x_new = np.linspace(pr_clean[0], pr_clean[-1], 100)

        y_new = interpolator(x_new)

        x_new = np.round(x_new, 3)
        y_new = np.round(y_new, 4)

        # Determine the new key
        if name == "sat":
            key = "-1"
        elif name.startswith("Tr_"):
            key = name.replace("Tr_", "")
        else:
            key = name

        processed_points = []
        for px, py in zip(x_new, y_new):
            processed_points.append({"p_r": px, "rho_r": py})
            
        processed_data[key] = processed_points
        print(f"Processed {name} -> {key}: {len(processed_points)} points")

    print(f"Writing to {output_path}")
    with open(output_path, 'w') as f:
        json.dump(processed_data, f, indent=4)
    print("Done.")

if __name__ == "__main__":
    main()
