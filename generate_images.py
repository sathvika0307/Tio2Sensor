import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase.io import read

# Directory with XYZ files
xyz_dir = 'simulation'

# Function to generate 3D image for XYZ file
def generate_3d_image(xyz_file, output_png):
    atoms = read(os.path.join(xyz_dir, xyz_file))
    positions = atoms.positions
    symbols = atoms.symbols

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Color map for elements
    colors = {'Ti': 'blue', 'O': 'red', 'C': 'black', 'H': 'gray', 'Co': 'green'}

    for pos, sym in zip(positions, symbols):
        color = colors.get(sym, 'purple')
        ax.scatter(pos[0], pos[1], pos[2], color=color, s=50)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'3D Structure: {xyz_file.replace(".xyz", "")}')

    plt.savefig(output_png)
    plt.close()

# List of XYZ files
xyz_files = ['tio2_nanoparticle.xyz', 'tio2_co_doped.xyz', 'tio2_size_5A.xyz', 'tio2_size_10A.xyz', 'tio2_size_15A.xyz']

# Generate images
for xyz_file in xyz_files:
    if os.path.exists(os.path.join(xyz_dir, xyz_file)):
        output_png = xyz_file.replace('.xyz', '.png')
        generate_3d_image(xyz_file, output_png)
        print(f"Generated {output_png}")

# For LAMMPS simulations, since no dump files, note that
print("Note: LAMMPS simulations did not produce dump files for visualization. Use XYZ files for structures.")
