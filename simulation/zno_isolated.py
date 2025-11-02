from ase.io import read, write

# Load ZnO nanoparticle
zno = read('zno_nanoparticle.xyz')

# Write to LAMMPS data
write('zno_isolated.lammps', zno, format='lammps-data')
print("Created ZnO isolated data file")
