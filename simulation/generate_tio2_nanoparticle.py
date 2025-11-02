from ase import Atoms
from ase.build import bulk
from ase.io import write
import numpy as np

# Create bulk TiO2 anatase structure
# ASE doesn't have anatase in bulk, so create manually
# Anatase unit cell parameters
a = 3.785
c = 9.514
lattice = [[a, 0, 0], [0, a, 0], [0, 0, c]]
# Fractional positions
scaled_positions = [
    [0, 0, 0],  # Ti
    [0.5, 0.5, 0.25],  # Ti
    [0.5, 0, 0.625],  # O
    [0, 0.5, 0.625],  # O
    [0.5, 0.5, 0.875],  # O
    [0, 0, 0.125],  # O
]
symbols = ['Ti', 'Ti', 'O', 'O', 'O', 'O']
bulk_anatase = Atoms(symbols=symbols, scaled_positions=scaled_positions, cell=lattice, pbc=True)

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
np_particle = create_nanoparticle(bulk_anatase, radius)

# Save to file
write('tio2_nanoparticle.xyz', np_particle)
print(f"Created TiO2 nanoparticle with {len(np_particle)} atoms")
