import os
import argparse
import numpy as np
import networkx as nx
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import neighbor_list

def merge_close_atoms(atoms, cutoff):
    """
    Merge atoms that are closer than cutoff.
    New position is the average of the cluster.
    """
    import time
    t0 = time.time()
    print(f"DEBUG: Entering merge_close_atoms with cutoff={cutoff}, atoms={len(atoms)}", flush=True)
    if cutoff <= 0:
        return atoms

    # Use neighbor_list to find all pairs within cutoff
    try:
        print("DEBUG: Calling neighbor_list...", flush=True)
        # neighbor_list can be slow for large systems with large cutoff
        # optimizing: use smaller skin or just ensure cutoff is reasonable
        i_indices, j_indices = neighbor_list('ij', atoms, cutoff)
        print(f"DEBUG: neighbor_list returned {len(i_indices)} pairs. Time: {time.time()-t0:.2f}s", flush=True)
    except Exception as e:
        print(f"  Warning: neighbor_list failed in merge_close_atoms: {e}", flush=True)
        return atoms
    
    # If no close atoms, return original
    if len(i_indices) == 0:
        return atoms
        
    # Build a graph
    print("DEBUG: Building graph...", flush=True)
    g = nx.Graph()
    g.add_nodes_from(range(len(atoms)))
    g.add_edges_from(zip(i_indices, j_indices))
    
    # Find connected components
    print("DEBUG: Finding connected components...", flush=True)
    components = list(nx.connected_components(g))
    
    if len(components) == len(atoms):
        return atoms
        
    print(f"  -> Merging close atoms (tol={cutoff}): reducing {len(atoms)} to {len(components)} atoms.", flush=True)
    
    new_positions = []
    new_numbers = []
    
    # For each component, calculate average position
    print("DEBUG: Calculating new positions...", flush=True)
    for comp in components:
        indices = list(comp)
        ref_idx = indices[0]
        
        if len(indices) == 1:
            new_positions.append(atoms.positions[ref_idx])
            new_numbers.append(atoms.numbers[ref_idx])
            continue
            
        vectors = atoms.get_distances(ref_idx, indices, mic=True, vector=True)
        avg_vector = np.mean(vectors, axis=0)
        new_pos = atoms.positions[ref_idx] + avg_vector
        new_positions.append(new_pos)
        new_numbers.append(atoms.numbers[ref_idx])

    new_atoms = Atoms(numbers=new_numbers, 
                      positions=new_positions, 
                      cell=atoms.cell, 
                      pbc=atoms.pbc)
                      
    print(f"DEBUG: merge_close_atoms done. Time: {time.time()-t0:.2f}s", flush=True)
    return new_atoms

def fast_stack(atoms1, atoms2, axis=2, distance=0.0):
    """
    A faster, simplified stack function avoiding ASE's overhead.
    Assumes cells are compatible perpendicular to the stacking direction.
    """
    print(f"DEBUG: Entering fast_stack axis={axis} distance={distance}", flush=True)
    
    # 1. Base cell from atoms1
    new_cell = atoms1.cell.copy()
    
    # 2. Vector to shift atoms2
    # We take the lattice vector of atoms1 along the axis
    v1 = atoms1.cell[axis]
    v2 = atoms2.cell[axis] # We'll add this to the cell
    
    # Direction for distance
    norm_v1 = np.linalg.norm(v1)
    if norm_v1 > 1e-8:
        dir_v1 = v1 / norm_v1
    else:
        dir_v1 = np.zeros(3)
        dir_v1[axis] = 1.0 # fallback
        
    # Shift vector = v1 + distance * direction
    shift_vec = v1 + dir_v1 * distance
    
    # 3. New cell along axis = shift_vec + v2 (approximately, mimicking stack)
    # Actually, standard stack logic: new_cell[axis] = v1 + v2 + distance_gap
    # But strictly speaking, the shift for atoms2 is v1 + gap.
    # And the new cell vector is v1 + gap + v2.
    
    # Let's verify how ASE does it: 
    # atom2.positions += atom1.cell[axis] + distance * normal
    # new_cell[axis] = atom1.cell[axis] + atom2.cell[axis] + distance * normal
    
    offset = shift_vec
    
    # Update cell
    # Note: This assumes v2 is the vector for atoms2's thickness
    # If atoms2.cell is not aligned or defined differently, this might be simplistic
    # But for VASP files usually OK.
    
    # We add the distance contribution to the cell vector as well
    extra_gap = dir_v1 * distance
    new_cell[axis] = v1 + v2 + extra_gap
    
    # 4. Create new atoms
    pos1 = atoms1.positions
    pos2 = atoms2.positions + offset
    
    new_pos = np.vstack([pos1, pos2])
    new_numbers = np.concatenate([atoms1.numbers, atoms2.numbers])
    
    new_atoms = Atoms(numbers=new_numbers, positions=new_pos, cell=new_cell, pbc=atoms1.pbc)
    
    print("DEBUG: fast_stack done", flush=True)
    return new_atoms

def batch_merge(dir1, dir2, out_dir, axis=1, merge_tol=0.0, distance=None, vacuum=None):
    import sys
    print("DEBUG: Function batch_merge started", flush=True)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"Created output directory: {out_dir}", flush=True)

    files1 = [f for f in os.listdir(dir1) if f.endswith('.vasp')]
    files2 = [f for f in os.listdir(dir2) if f.endswith('.vasp')]
    print(f"DEBUG: Found {len(files1)} files in {dir1} and {len(files2)} files in {dir2}", flush=True)

    map1 = {}
    for f in files1:
        if '_' in f:
            suffix = f.split('_', 1)[1]
            map1[suffix] = f
        else:
            print(f"Warning: File {f} in {dir1} does not contain '_', skipping.", flush=True)

    map2 = {}
    for f in files2:
        if '_' in f:
            suffix = f.split('_', 1)[1]
            map2[suffix] = f
        else:
            print(f"Warning: File {f} in {dir2} does not contain '_', skipping.", flush=True)

    count = 0
    for suffix, f1 in map1.items():
        if suffix in map2:
            f2 = map2[suffix]
            path1 = os.path.join(dir1, f1)
            path2 = os.path.join(dir2, f2)
            
            print(f"Merging {f1} + {f2} along axis {axis}...", flush=True)
            
            try:
                print(f"DEBUG: Reading {path1}...", flush=True)
                atoms1 = read(path1)
                print(f"DEBUG: Reading {path2}...", flush=True)
                atoms2 = read(path2)
                
                # 使用 fast_stack 替代 ase.build.stack
                dist_val = distance if distance is not None else 0.0
                print(f"DEBUG: Stacking with fast_stack...", flush=True)
                merged_atoms = fast_stack(atoms1, atoms2, axis=axis, distance=dist_val)
                print(f"DEBUG: Stacked. Atoms count: {len(merged_atoms)}", flush=True)
                
                # 检查最小原子间距
                try:
                    print("DEBUG: Checking distances...", flush=True)
                    # Use smaller cutoff for check to be faster
                    dists = neighbor_list('d', merged_atoms, 3.0) 
                    if len(dists) > 0:
                        min_dist = np.min(dists)
                        print(f"  -> Min atomic distance after stack: {min_dist:.3f} A", flush=True)
                    else:
                        print(f"  -> Min atomic distance after stack: > 3.0 A", flush=True)
                except Exception as e:
                    print(f"DEBUG: Distance check failed: {e}", flush=True)

                # 如果指定了阈值，进行原子合并
                if merge_tol > 0:
                    print(f"DEBUG: Merging close atoms with tol={merge_tol}...", flush=True)
                    merged_atoms = merge_close_atoms(merged_atoms, merge_tol)
                    print(f"DEBUG: Merge done. New atoms count: {len(merged_atoms)}", flush=True)
                
                # 如果指定了真空层，进行添加
                if vacuum is not None:
                    print(f"DEBUG: Adding vacuum {vacuum}...", flush=True)
                    merged_atoms.center(vacuum=vacuum, axis=axis)
                
                out_name = f"merged_{suffix}"
                out_path = os.path.join(out_dir, out_name)
                print(f"DEBUG: Writing to {out_path}...", flush=True)
                write(out_path, merged_atoms)
                print(f"  -> Saved to {out_path}", flush=True)
                count += 1
            except Exception as e:
                print(f"  Error merging {f1} and {f2}: {e}", flush=True)
                import traceback
                traceback.print_exc()
        else:
            print(f"No match found for {f1} (suffix: {suffix}) in {dir2}", flush=True)

    print(f"\nDone. Merged {count} pairs.", flush=True)

if __name__ == "__main__":
    # ================= 参数修改区域 =================
    
    dir1 = "stru1"
    dir2 = "stru2"
    out_dir = "out"
    
    # 合并轴 (0=x, 1=y, 2=z)
    axis = 1 
    
    # 两层之间的距离
    distance = 0.5
    
    # 原子合并阈值
    merge_tol = 2.5
    
    # 真空层
    vacuum = None
    
    # ===============================================

    print(f"Configuration: stru1={dir1}, stru2={dir2}, out={out_dir}")
    print(f"               axis={axis}, distance={distance}, merge_tol={merge_tol}, vacuum={vacuum}")
    
    batch_merge(dir1, dir2, out_dir, axis=axis, merge_tol=merge_tol, distance=distance, vacuum=vacuum)
