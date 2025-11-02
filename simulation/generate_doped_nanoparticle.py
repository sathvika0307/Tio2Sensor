from ase import Atoms
from ase.io import write
import numpy as np

# Load pristine nanoparticle
np_particle = Atoms()  # Placeholder, load from file
# Since file exists, read it
from ase.io import read
np_particle = read('tio2_nanoparticle.xyz')

# Function to dope: replace some Ti atoms with Co
def dope_nanoparticle(atoms, dopant='Co', fraction=0.05):
    symbols = atoms.symbols
    ti_indices = [i for i, sym in enumerate(symbols) if sym == 'Ti']
    num_dope = int(len(ti_indices) * fraction)
    dope_indices = np.random.choice(ti_indices, num_dope, replace=False)
    new_symbols = list(symbols)
    for i in dope_indices:
        new_symbols[i] = dopant
    doped_atoms = Atoms(symbols=new_symbols, positions=atoms.positions, cell=atoms.cell, pbc=False)
    return doped_atoms

# Create 5% Co-doped nanoparticle
doped_np = dope_nanoparticle(np_particle, dopant='Co', fraction=0.05)

# Save
write('tio2_co_doped.xyz', doped_np)
print(f"Created Co-doped TiO2 nanoparticle with {len(doped_np)} atoms, {sum(1 for s in doped_np.symbols if s=='Co')} Co atoms")
