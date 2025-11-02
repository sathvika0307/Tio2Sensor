from ase import Atoms
from ase.build import bulk
from ase.io import write
import numpy as np

# Create bulk ZnO wurtzite structure
# Wurtzite unit cell parameters
a = 3.2495
c = 5.2069
lattice = [[a, 0, 0], [-a/2, a*np.sqrt(3)/2, 0], [0, 0, c]]
# Fractional positions for wurtzite ZnO
scaled_positions = [
    [1/3, 2/3, 0],      # Zn
    [2/3, 1/3, 0.5],    # Zn
    [1/3, 2/3, 0.382],  # O
    [2/3, 1/3, 0.882],  # O
]
symbols = ['Zn', 'Zn', 'O', 'O']
bulk_wurtzite = Atoms(symbols=symbols, scaled_positions=scaled_positions, cell=lattice, pbc=True)

# Function to create spherical nanoparticle by carving from supercell
def create_nanoparticle(bulk, radius, center=None):
    supercell = bulk * (10, 10, 10)  # Large supercell
    if center is None:
        center = supercell.cell.diagonal() / 2
    positions = supercell.positions
    distances = np.linalg.norm(positions - center, axis=1)
    mask = distances <= radius
    nanoparticle = Atoms(symbols=supercell.symbols[mask],
                        positions=supercell.positions[mask],
                        cell=supercell.cell,
                        pbc=False)
    # Remove incomplete molecules or dangling bonds if needed
    return nanoparticle

# Create nanoparticle with radius ~1 nm (diameter 2 nm)
radius = 10  # in Angstroms
np_particle = create_nanoparticle(bulk_wurtzite, radius)

# Save to file
write('zno_nanoparticle.xyz', np_particle)
print(f"Created ZnO nanoparticle with {len(np_particle)} atoms")
