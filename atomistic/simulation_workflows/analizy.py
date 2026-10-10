import os
import re
import matplotlib.pyplot as plt
import numpy as np
import glob

def extract_vasp_stress(outcar_path):
    """
    Extract the Z-direction stress (in GPa) from VASP OUTCAR.
    Reads the 'in kB' line below 'FORCE on cell =-STRESS'.
    """
    if not os.path.exists(outcar_path):
        return None
        
    with open(outcar_path, "r") as f:
        lines = f.readlines()
        
    # We want the last occurrence since VASP runs multiple ionic steps
    stress_z_kb = None
    for i in range(len(lines)):
        if "in kB" in lines[i]:
            parts = lines[i].split()
            if len(parts) >= 6 and parts[0] == "in" and parts[1] == "kB":
                # parts[2]=XX, parts[3]=YY, parts[4]=ZZ
                try:
                    stress_z_kb = float(parts[4])
                except ValueError:
                    pass
                    
    if stress_z_kb is not None:
        # VASP "in kB" prints -STRESS (i.e. FORCE on cell = -STRESS).
        # Convention: positive kB value = compressive, negative = tensile.
        # Negate to convert to tensile-positive engineering convention,
        # then divide by 10 to convert kBar -> GPa.
        return -stress_z_kb / 10.0
    return None

def extract_lammps_stress(stress_file_path):
    """
    Extract the Z-direction stress (in GPa) from LAMMPS output file.
    """
    if not os.path.exists(stress_file_path):
        return None
        
    with open(stress_file_path, "r") as f:
        for line in f:
            if "ZZ_STRESS_GPA" in line:
                try:
                    return float(line.split()[1])
                except (ValueError, IndexError):
                    return None
    return None

def analyze_directory(work_dir, orientation):
    print(f"Analyzing {work_dir} (Orientation {orientation})...")
    
    # Detect strains from directory structure
    strain_dirs = glob.glob(os.path.join(work_dir, "strain_*"))
    strains = []
    for d in strain_dirs:
        try:
            name = os.path.basename(d)
            s = int(name.replace("strain_", ""))
            strains.append(s)
        except ValueError:
            continue
    
    if not strains:
        print(f"No strain folders found in {work_dir}")
        return

    strains.sort()
    
    vasp_stresses = []
    eam_stresses = []
    sus2_stresses = []
    
    for s in strains:
        folder = os.path.join(work_dir, f"strain_{s:02d}")
        
        # VASP
        outcar = os.path.join(folder, "OUTCAR")
        v_stress = extract_vasp_stress(outcar)
        vasp_stresses.append(v_stress if v_stress is not None else np.nan)
        
        # EAM
        eam_file = os.path.join(folder, "stress_eam.txt")
        e_stress = extract_lammps_stress(eam_file)
        eam_stresses.append(e_stress if e_stress is not None else np.nan)
        
        # sus2
        sus2_file = os.path.join(folder, "stress_sus2.txt")
        s_stress = extract_lammps_stress(sus2_file)
        sus2_stresses.append(s_stress if s_stress is not None else np.nan)
        
    # Plotting
    plt.figure(figsize=(10, 6))
    
    # Filter out NaNs for cleaner lines, or keep them to show gaps
    # Matplotlib handles NaNs by breaking lines, which is good.
    
    plt.plot(strains, vasp_stresses, 'o-', color='black', label="VASP (DFT Ground Truth)", linewidth=2, markersize=8)
    plt.plot(strains, eam_stresses, 's--', color='blue', label="LAMMPS (EAM)", linewidth=2, markersize=8)
    plt.plot(strains, sus2_stresses, '^-.', color='red', label="LAMMPS (sus2 ML)", linewidth=2, markersize=8)
    
    plt.xlabel("Tensile Strain (%)", fontsize=14)
    plt.ylabel("Tensile Stress (GPa)", fontsize=14)
    plt.title(f"Cu [{orientation}] Uniaxial Tensile Stress-Strain Curve", fontsize=16)
    
    # Adjust X-ticks based on data range
    if strains:
        max_strain = max(strains)
        step = 10 if max_strain > 20 else 2
        plt.xticks(np.arange(0, max_strain + step, step), fontsize=12)
    
    plt.yticks(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=12)
    
    plt.tight_layout()
    output_img = f"stress_strain_cu_{orientation}.png"
    plt.savefig(output_img, dpi=300)
    print(f"Plot saved to {output_img}")

    # Output CSV table for reference
    csv_file = f"stress_strain_data_{orientation}.csv"
    with open(csv_file, "w") as f:
        f.write("Strain(%),VASP_Stress(GPa),EAM_Stress(GPa),sus2_Stress(GPa)\n")
        for strain, v, e, s in zip(strains, vasp_stresses, eam_stresses, sus2_stresses):
            # Handle NaN for string formatting
            v_str = f"{v:.4f}" if not np.isnan(v) else "NaN"
            e_str = f"{e:.4f}" if not np.isnan(e) else "NaN"
            s_str = f"{s:.4f}" if not np.isnan(s) else "NaN"
            f.write(f"{strain},{v_str},{e_str},{s_str}\n")
    print(f"Data saved to {csv_file}")

def main():
    # Check if we are inside a work directory
    current_dir = os.getcwd()
    dir_name = os.path.basename(current_dir)
    
    if dir_name.startswith("work_"):
        # Single directory mode
        orientation = dir_name.replace("work_", "")
        analyze_directory(current_dir, orientation)
    else:
        # Root directory mode: search for work_* folders
        work_dirs = glob.glob("work_*")
        if not work_dirs:
            # Fallback: Check if we have strain_xx folders in current dir (legacy/flat structure)
            if glob.glob("strain_*"):
                 analyze_directory(current_dir, "unknown")
            else:
                print("No work directories (work_*) or strain folders found.")
        else:
            for wd in work_dirs:
                orientation = os.path.basename(wd).replace("work_", "")
                analyze_directory(wd, orientation)

if __name__ == "__main__":
    main()
