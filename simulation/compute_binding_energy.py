import subprocess
import os

# Run sim for doped nanoparticle alone
lammps_input_np = """
units metal
atom_style atomic
boundary p p p

read_data doped_adsorption.lammps

# Remove methane atoms (assume last 5 are methane)
group methane type 1 2
delete_atoms group methane

mass 3 16.0  # O
mass 4 48.0  # Ti
mass 5 58.93 # Co

pair_style lj/cut 10.0
pair_coeff * * 0.1 3.0

minimize 1e-4 1e-6 1000 1000

compute pe all pe
thermo_style custom pe
thermo 1
"""

with open('in_np.lammps', 'w') as f:
    f.write(lammps_input_np)

subprocess.run(['/opt/homebrew/Cellar/lammps/20250722-update1/bin/lmp_serial', '-in', 'in_np.lammps'], cwd=os.getcwd())

# Extract energy from log
with open('log.lammps', 'r') as f:
    lines = f.readlines()
    for line in reversed(lines):
        if 'TotEng' in line:
            e_np = float(line.split()[-1])
            break

# From doped adsorption log
with open('log.lammps', 'r') as f:  # Assuming log.lammps is from doped sim
    lines = f.readlines()
    for line in reversed(lines):
        if 'TotEng' in line:
            e_system = float(line.split()[-1])
            break

# Binding energy (approximate, since methane energy ~0)
binding_energy = e_system - e_np
print(f"Binding energy: {binding_energy} eV")
