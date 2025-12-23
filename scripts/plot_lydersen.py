import json
import os
import matplotlib.pyplot as plt
import numpy as np

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    input_path = os.path.join(project_root, 'data', 'lydersen_normal.json')
    output_path = os.path.join(project_root, 'images', 'lydersen_plot.png')

    print(f"Reading from {input_path}")
    try:
        with open(input_path, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: Could not find {input_path}")
        return

    plt.figure(figsize=(10, 8))
    
    # Define a color map or list of colors
    colors = plt.cm.viridis(np.linspace(0, 1, len(data)))

    for i, (name, points) in enumerate(data.items()):
        if not points:
            continue
        
        pr = [p['p_r'] for p in points]
        rho = [p['rho_r'] for p in points]
        
        # Plot
        if name == "-1":
            plt.plot(pr, rho, label="Saturation", color='black', linewidth=2, linestyle='--')
        else:
            plt.plot(pr, rho, label=fr"$T_r = {name}$")

    plt.xlabel(r'Reduced Pressure ($P_r$)')
    plt.ylabel(r'Reduced Density ($\rho_r$)')
    plt.title('Lydersen Chart (Reconstructed)')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlim(0, 10)
    plt.ylim(0, 3.5)

    print(f"Saving plot to {output_path}")
    plt.savefig(output_path, dpi=300)
    print("Done.")

if __name__ == "__main__":
    main()
