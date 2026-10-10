import os
import glob
import matplotlib.pyplot as plt
import csv

def get_area(poscar_path):
    """Calculate the cross-sectional area (xy) from POSCAR."""
    with open(poscar_path, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    scale = float(lines[1])
    lx = float(lines[2].split()[0]) * scale
    ly = float(lines[3].split()[1]) * scale
    return lx * ly

def get_vasp_energy(outcar_path):
    """Extract final total energy from VASP OUTCAR."""
    if not os.path.exists(outcar_path):
        return None
    energy = None
    with open(outcar_path, 'r', encoding='utf-8') as f:
        for line in f:
            if "free  energy   TOTEN" in line:
                energy = float(line.split()[-2])
    return energy

def extract_vasp_data(base_dir):
    calc_dirs = sorted(glob.glob(os.path.join(base_dir, "calc_*")))
    data = []
    
    if not calc_dirs:
        print(f"No calculation directories found in {base_dir}")
        return data

    poscar_files = sorted(glob.glob("stru/*.vasp"))
    if poscar_files:
        area = get_area(poscar_files[0])
    else:
        area = get_area(os.path.join("vasp", "calc_0.0", "POSCAR"))
    
    for cdir in calc_dirs:
        disp_str = os.path.basename(cdir).replace("calc_", "")
        disp = float(disp_str)
        out_file = os.path.join(cdir, "vasp_single.log")
        # Read from log instead if OUTCAR not present
        if not os.path.exists(os.path.join(cdir, "OUTCAR")):
            pass
        e = get_vasp_energy(os.path.join(cdir, "OUTCAR"))
        if e is not None:
            data.append((disp, e))
            
    if not data:
        return []
        
    data.sort(key=lambda x: x[0])
    e0 = data[0][1]
    gsfe_data = []
    for disp, e in data:
        # Divide by 1*area because we use vacuum (1 fault plane)
        gsfe = (e - e0) / area * 16021.7662
        gsfe_data.append((disp, gsfe))
        
    return gsfe_data

def extract_lammps_data(out_path):
    """Extract energies from LAMMPS output where multiple DISP prints exist."""
    if not os.path.exists(out_path):
        print(f"No output file {out_path} found.")
        return []
        
    data = []
    with open(out_path, 'r', encoding='utf-8') as f:
        for line in f:
            if "DISP:" in line and "SFE:" in line:
                # Format: DISP: 0.1 SFE: 12.34 mJ/m^2
                parts = line.split()
                try:
                    disp = float(parts[1])
                    sfe = float(parts[3])
                    data.append((disp, sfe))
                except:
                    pass
    data.sort(key=lambda x: x[0])
    return data

def main():
    print("Extracting VASP data...")
    vasp_data = extract_vasp_data("vasp")
    
    print("Extracting LAMMPS EAM data...")
    eam_data = extract_lammps_data("lammps_eam/lammps.out")
    
    print("Extracting LAMMPS SUS2 data...")
    sus2_data = extract_lammps_data("lammps_sus2/lammps.out")

    # Combine data to save to CSV
    all_disps = set()
    if vasp_data:
        all_disps.update([d[0] for d in vasp_data])
    if eam_data:
        all_disps.update([d[0] for d in eam_data])
    if sus2_data:
        all_disps.update([d[0] for d in sus2_data])
        
    all_disps = sorted(list(all_disps))
    
    vasp_dict = dict(vasp_data) if vasp_data else {}
    eam_dict = dict(eam_data) if eam_data else {}
    sus2_dict = dict(sus2_data) if sus2_data else {}

    csv_file = "gsfe_data.csv"
    with open(csv_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["Displacement_A", "VASP_mJ_m2", "LAMMPS_EAM_mJ_m2", "LAMMPS_SUS2_mJ_m2"])
        for disp in all_disps:
            row = [f"{disp:.2f}"]
            row.append(f"{vasp_dict[disp]:.4f}" if disp in vasp_dict else "")
            row.append(f"{eam_dict[disp]:.4f}" if disp in eam_dict else "")
            row.append(f"{sus2_dict[disp]:.4f}" if disp in sus2_dict else "")
            writer.writerow(row)
    print(f"Data saved to {csv_file}")

    plt.figure(figsize=(8, 6))
    
    if vasp_data:
        x, y = zip(*vasp_data)
        plt.plot(x, y, 'o-', label='VASP (DFT)', linewidth=2, markersize=6)
        
    if eam_data:
        x, y = zip(*eam_data)
        plt.plot(x, y, 's--', label='LAMMPS (EAM)', linewidth=2, markersize=6)
        
    if sus2_data:
        x, y = zip(*sus2_data)
        plt.plot(x, y, '^-.', label='LAMMPS (MLIP SUS2)', linewidth=2, markersize=6)

    plt.xlabel("Displacement (Å)", fontsize=14)
    plt.ylabel("Generalized Stacking Fault Energy (mJ/m$^2$)", fontsize=14)
    plt.title("Generalized Stacking Fault Energy (GSFE) Curve of Cu", fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tick_params(labelsize=12)
    plt.tight_layout()

    out_fig = "gsfe_curve.png"
    plt.savefig(out_fig, dpi=300)
    print(f"Plot saved to {out_fig}")
    # plt.show() # Commented out for automated testing

if __name__ == "__main__":
    main()
