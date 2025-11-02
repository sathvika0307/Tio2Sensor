from ase import Atoms
from ase.io import read, write
import subprocess
import os
import numpy as np

# Load doped nanoparticle
doped_np = read('tio2_co_doped.xyz')

# Create methane
methane = Atoms('CH4', positions=[[0, 0, 0], [1.09, 0, 0], [-0.363, 1.02, 0], [-0.363, -0.51, 0.89], [-0.363, -0.51, -0.89]])
from ase import Atoms

# Place methane
center = np.mean(doped_np.positions, axis=0)
surface_z = np.max(doped_np.positions[:, 2]) + 5
methane.positions += [center[0], center[1], surface_z]

# Combine
system = doped_np + methane

# Write LAMMPS data
write('doped_adsorption.lammps', system, format='lammps-data')

# Create input file similar to before, but with masses for Co
lammps_input = """
# LAMMPS input for doped TiO2 methane adsorption

units metal
atom_style atomic
boundary p p p

read_data doped_adsorption.lammps

# Masses
mass 1 12.0  # C
mass 2 1.0   # H
mass 3 16.0  # O
mass 4 48.0  # Ti
mass 5 58.93 # Co

pair_style lj/cut 10.0
pair_coeff * * 0.1 3.0

minimize 1e-4 1e-6 1000 1000

timestep 0.001
fix 1 all nvt temp 100 100 0.1
run 1000

compute pe all pe
thermo_style custom step temp pe
thermo 100
"""

with open('in_doped.lammps', 'w') as f:
    f.write(lammps_input)

# Run LAMMPS
subprocess.run(['/opt/homebrew/Cellar/lammps/20250722-update1/bin/lmp_serial', '-in', 'in_doped.lammps'], cwd=os.getcwd())

print("Simulation for doped system completed")
