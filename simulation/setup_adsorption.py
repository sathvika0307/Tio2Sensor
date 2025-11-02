from ase import Atoms
from ase.io import read, write
import numpy as np

# Load the TiO2 nanoparticle
np_particle = read('tio2_nanoparticle.xyz')

# Create methane molecule
methane = Atoms('CH4', positions=[[0, 0, 0], [1.09, 0, 0], [-0.363, 1.02, 0], [-0.363, -0.51, 0.89], [-0.363, -0.51, -0.89]])

# Place methane near the surface, say 5 Angstroms above the center
center = np.mean(np_particle.positions, axis=0)
surface_z = np.max(np_particle.positions[:, 2]) + 5
methane.positions += [center[0], center[1], surface_z]

# Combine structures
adsorption_system = np_particle + methane

# Write to LAMMPS data file
write('adsorption_system.lammps', adsorption_system, format='lammps-data')

print("Created adsorption system with TiO2 nanoparticle and methane molecule")
