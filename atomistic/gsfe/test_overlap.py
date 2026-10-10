import numpy as np
latparam = 3.615
# lattice vectors for FCC oriented as requested
a1 = latparam * np.array([1, 1, 2])/np.sqrt(6)
a2 = latparam * np.array([-1, 1, 0])/np.sqrt(2)
a3 = latparam * np.array([-1, -1, 1])/np.sqrt(3)

print("Vectors:")
print(a1)
print(a2)
print(a3)
