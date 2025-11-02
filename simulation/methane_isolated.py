from ase import Atoms
from ase.io import write
import numpy as np

# Create methane molecule
methane = Atoms('CH4', positions=[[0, 0, 0], [1.09, 0, 0], [-0.363, 1.02, 0], [-0.363, -0.51, 0.89], [-0.363, -0.51, -0.89]])

# Set a cell
methane.set_cell([10, 10, 10])
methane.set_pbc(False)

# Write to LAMMPS data
write('methane_isolated.lammps', methane, format='lammps-data')
print("Created methane isolated data file")
