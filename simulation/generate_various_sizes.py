from ase import Atoms
from ase.build import bulk
from ase.io import write
import numpy as np

# Anatase bulk
lattice = [[3.785, 0, 0], [0, 3.785, 0], [0, 0, 9.514]]
scaled_positions = [
    [0, 0, 0], [0.5, 0.5, 0.25],
    [0.5, 0, 0.625], [0, 0.5, 0.625],
    [0.5, 0.5, 0.875], [0, 0, 0.125],
]
symbols = ['Ti', 'Ti', 'O', 'O', 'O', 'O']
bulk_anatase = Atoms(symbols=symbols, scaled_positions=scaled_positions, cell=lattice, pbc=True)

def create_nanoparticle(bulk, radius, center=None):
    supercell = bulk * (10, 10, 10)
    if center is None:
        center = supercell.cell.diagonal() / 2
    positions = supercell.positions
    distances = np.linalg.norm(positions - center, axis=1)
    mask = distances <= radius
    nanoparticle = Atoms(symbols=supercell.symbols[mask],
                        positions=supercell.positions[mask],
                        cell=supercell.cell,
                        pbc=False)
    return nanoparticle

# Different sizes
sizes = [5, 10, 15]  # Angstrom radii
for size in sizes:
    np_part = create_nanoparticle(bulk_anatase, size)
    write(f'tio2_size_{size}A.xyz', np_part)
    print(f"Created nanoparticle with radius {size} A, {len(np_part)} atoms")
