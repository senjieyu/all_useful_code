import numpy as np
import networkx as nx
from ase import Atoms
from ase.neighborlist import neighbor_list

def merge_close_atoms(atoms, cutoff):
    print(f"Testing merge with cutoff {cutoff}")
    i_indices, j_indices = neighbor_list('ij', atoms, cutoff)
    print(f"Found {len(i_indices)} pairs")
    
    if len(i_indices) == 0:
        return atoms
        
    g = nx.Graph()
    g.add_nodes_from(range(len(atoms)))
    g.add_edges_from(zip(i_indices, j_indices))
    
    components = list(nx.connected_components(g))
    print(f"Components: {len(components)}")
    
    return atoms

# Create dummy atoms
atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.5]], cell=[10, 10, 10], pbc=True)
merge_close_atoms(atoms, 0.6)
